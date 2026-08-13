#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'USAGE'
Usage:
  bash workflow/run_postprocess.sh [--from-module auto|03|04|05] \
    SAMPLE_ID CONTROL_ID SGRNA [config.env]

Consumes completed Phase A outputs for SAMPLE_ID and CONTROL_ID. It does not
rerun FASTQ processing.  The default "auto" mode resumes at the earliest
incomplete or stale module.
USAGE
}

FROM_MODULE=auto
if [[ "${1:-}" == "--from-module" ]]; then
    [[ $# -ge 2 ]] || { usage >&2; exit 2; }
    FROM_MODULE=$2
    shift 2
fi

if [[ $# -lt 3 || $# -gt 4 ]]; then
    usage >&2
    exit 2
fi

sample=$1
control=$2
sgrna=$3

if [[ -n "${INOSEQ_PROJECT_DIR:-}" ]]; then
    PROJECT_DIR=$(cd "${INOSEQ_PROJECT_DIR}" && pwd)
else
    SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
    PROJECT_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
fi
CONFIG_FILE=${4:-"${PROJECT_DIR}/config/inoseq.env"}
[[ -f "$CONFIG_FILE" ]] || { echo "[ERROR] Config not found: $CONFIG_FILE" >&2; exit 2; }
CONFIG_FILE=$(cd "$(dirname "$CONFIG_FILE")" && pwd)/$(basename "$CONFIG_FILE")

# shellcheck disable=SC1090
source "$CONFIG_FILE"
: "${PYTHON_BIN:=python}"
: "${OUTPUT_DIR:=output}"
: "${REFERENCE_FASTA:?REFERENCE_FASTA is not set in $CONFIG_FILE}"
: "${BACKGROUND_WINDOW:=15}"
: "${BACKGROUND_FOLD_CHANGE:=1.5}"
: "${BACKGROUND_PVALUE:=0.05}"
: "${CANDIDATE_MERGE_DISTANCE:=30}"
: "${CANDIDATE_MIN_LENGTH:=30}"
: "${CANDIDATE_MIN_READS:=3}"
: "${OFFTARGET_SEARCH_WINDOW:=25}"
: "${OFFTARGET_MAX_SCORE:=8}"
: "${SPACER_NEIGHBORHOOD:=100}"

if [[ "$OUTPUT_DIR" != /* ]]; then
    OUTPUT_DIR="${PROJECT_DIR}/${OUTPUT_DIR}"
fi
export PYTHONPATH="${PROJECT_DIR}${PYTHONPATH:+:${PYTHONPATH}}"
INOSEQ_VERSION=$(tr -d '[:space:]' < "${PROJECT_DIR}/VERSION")

sample_alignment="${OUTPUT_DIR}/${sample}/alignment"
control_alignment="${OUTPUT_DIR}/${control}/alignment"
sample_end="${sample_alignment}/${sample}.end"
sample_bam="${sample_alignment}/${sample}_end.bam"
control_bam="${control_alignment}/${control}_end.bam"

for path in "$sample_end" "$sample_bam" "${sample_bam}.bai" "$control_bam" "${control_bam}.bai" "$REFERENCE_FASTA"; do
    [[ -r "$path" ]] || { echo "[ERROR] Required input not readable: $path" >&2; exit 2; }
done

post_dir="${OUTPUT_DIR}/${sample}/postprocess"
before_dir="${post_dir}/filt-before"
after_dir="${post_dir}/filt-after"
candidate_dir="${post_dir}/final_merge_distance"
classification_dir="${post_dir}/offtarget"
summary_dir="${post_dir}/summary"
mkdir -p "$before_dir" "$after_dir" "$candidate_dir" "$classification_dir" "$summary_dir"

log_step() { printf '[%s] [%s vs %s] %s\n' "$(date '+%F %T')" "$sample" "$control" "$*" >&2; }
run_module() { "$PYTHON_BIN" "$PROJECT_DIR/workflow/modules/$1" "${@:2}"; }

merged_tsv="${before_dir}/${sample}_merge.tsv"
bed_file="${before_dir}/${sample}.bed"
exp_coverage="${before_dir}/${sample}_coverage.txt"
ctrl_coverage="${before_dir}/${sample}_ctr_${control}_coverage.txt"
background_table="${before_dir}/${sample}.txt"
filtered_table="${after_dir}/${sample}_filted.txt"
query_names="${after_dir}/${sample}_query_name.txt"
filtered_bam="${after_dir}/${sample}_filted.bam"
candidate_table="${candidate_dir}/${sample}_final_merged_distance.txt"
candidate_bam="${candidate_dir}/${sample}_filted.bam"
alignment_table="${classification_dir}/${sample}_align.txt"
dependency_table="${classification_dir}/${sample}_dependent_mark.txt"
spacer_table="${classification_dir}/${sample}_dependent_out_of_spacer.txt"
target_table="${classification_dir}/${sample}_dependent_target.txt"
STATE_TOOL="${PROJECT_DIR}/workflow/utils/stage_state.py"
state_args=(
    --project-dir "$PROJECT_DIR" --config "$CONFIG_FILE"
    --sample "$sample" --control "$control" --sgrna "$sgrna"
)

case "$FROM_MODULE" in
    auto)
        if "$PYTHON_BIN" "$STATE_TOOL" check --stage module05 "${state_args[@]}" --quiet; then
            log_step "Phase B already current; no computation required"
            exit 0
        fi
        if "$PYTHON_BIN" "$STATE_TOOL" check --stage module04 "${state_args[@]}" --quiet; then
            FROM_MODULE=05
        elif "$PYTHON_BIN" "$STATE_TOOL" check --stage module03 "${state_args[@]}" --quiet; then
            FROM_MODULE=04
        else
            FROM_MODULE=03
        fi
        ;;
    03|04|05) ;;
    *)
        echo "[ERROR] --from-module must be auto, 03, 04 or 05: $FROM_MODULE" >&2
        exit 2
        ;;
esac

if [[ "$FROM_MODULE" == "04" ]] \
    && ! "$PYTHON_BIN" "$STATE_TOOL" check --stage module03 "${state_args[@]}" --quiet; then
    echo "[ERROR] Cannot start at module 04: module 03 is incomplete or stale for $sample" >&2
    exit 2
fi
if [[ "$FROM_MODULE" == "05" ]] \
    && ! "$PYTHON_BIN" "$STATE_TOOL" check --stage module04 "${state_args[@]}" --quiet; then
    echo "[ERROR] Cannot start at module 05: module 04 is incomplete or stale for $sample" >&2
    exit 2
fi

"$PYTHON_BIN" "$STATE_TOOL" invalidate \
    --stage "module${FROM_MODULE}" "${state_args[@]}"

if (( 10#$FROM_MODULE <= 3 )); then
    log_step "Phase B 1/10 | module 03a: aggregate signature reads into cleavage-site windows"
    run_module 03_aggregate_sites.py \
        --input "$sample_end" --output-tsv "$merged_tsv" --output-bed "$bed_file" \
        --window "$BACKGROUND_WINDOW"

    log_step "Phase B 2/10 | module 03b: count experimental and matched-control overlaps"
    run_module 03_count_coverage.py --bed "$bed_file" --bam "$sample_bam" --output "$exp_coverage"
    run_module 03_count_coverage.py --bed "$bed_file" --bam "$control_bam" --output "$ctrl_coverage"

    log_step "Phase B 3/10 | module 03c: calculate fold change, p value, and BH-FDR"
    run_module 03_compare_background.py \
        --experiment "$exp_coverage" --control "$ctrl_coverage" --output "$background_table"

    log_step "Phase B 4/10 | module 03d: filter background: fold_change >= ${BACKGROUND_FOLD_CHANGE}; p_value < ${BACKGROUND_PVALUE}"
    run_module 03_filter_background.py \
        --input "$background_table" --input-bam "$sample_bam" \
        --output-table "$filtered_table" --output-query-names "$query_names" \
        --output-bam "$filtered_bam" --min-fold-change "$BACKGROUND_FOLD_CHANGE" \
        --p-value-threshold "$BACKGROUND_PVALUE"
    "$PYTHON_BIN" "$STATE_TOOL" record --stage module03 "${state_args[@]}"
else
    log_step "Phase B | module 03: reused current result"
fi

if (( 10#$FROM_MODULE <= 4 )); then
    log_step "Phase B 5/10 | module 04: merge and standardize candidate intervals"
    run_module 04_candidate_intervals.py \
        --input "$filtered_table" --input-bam "$filtered_bam" \
        --output "$candidate_table" --output-bam "$candidate_bam" \
        --distance "$CANDIDATE_MERGE_DISTANCE" --min-length "$CANDIDATE_MIN_LENGTH" \
        --min-reads "$CANDIDATE_MIN_READS"
    "$PYTHON_BIN" "$STATE_TOOL" record --stage module04 "${state_args[@]}"
else
    log_step "Phase B | module 04: reused current result"
fi

log_step "Phase B 6/10 | module 05a: align sgRNA against candidate windows"
run_module 05_align_sgrna.py \
    --input "$candidate_table" --ref "$REFERENCE_FASTA" --sgrna "$sgrna" \
    --output "$alignment_table" --window "$OFFTARGET_SEARCH_WINDOW" \
    --max-score "$OFFTARGET_MAX_SCORE"

log_step "Phase B 7/10 | module 05b: assign dependency class"
run_module 05_mark_dependency.py --input "$alignment_table" --output "$dependency_table"

log_step "Phase B 8/10 | module 05c: annotate protospacer and neighboring intervals"
run_module 05_annotate_spacer.py \
    --input "$dependency_table" --output "$spacer_table" \
    --neighborhood "$SPACER_NEIGHBORHOOD"

log_step "Phase B 9/10 | module 05d: classify target and non-target strand reads"
run_module 05_classify_strands.py \
    --input "$spacer_table" --bam "$candidate_bam" --output "$target_table" \
    --window "$OFFTARGET_SEARCH_WINDOW"

log_step "Phase B 10/10 | module 05e: write sample-level summaries"
run_module 05_summarize.py \
    --input "$target_table" --sample-id "$sample" \
    --basic-tsv "${summary_dir}/${sample}_offtarget_summary.tsv" \
    --detail-tsv "${summary_dir}/${sample}_strand_summary.tsv" \
    --excel "${summary_dir}/dependent_target_analysis.xlsx"

{
    printf 'parameter\tvalue\n'
    printf 'inoseq_version\t%s\n' "$INOSEQ_VERSION"
    printf 'sample_id\t%s\n' "$sample"
    printf 'control_id\t%s\n' "$control"
    printf 'sgrna\t%s\n' "$sgrna"
    printf 'background_window\t%s\n' "$BACKGROUND_WINDOW"
    printf 'background_fold_change\t%s\n' "$BACKGROUND_FOLD_CHANGE"
    printf 'background_pvalue\t%s\n' "$BACKGROUND_PVALUE"
    printf 'candidate_merge_distance\t%s\n' "$CANDIDATE_MERGE_DISTANCE"
    printf 'candidate_min_length\t%s\n' "$CANDIDATE_MIN_LENGTH"
    printf 'candidate_min_reads\t%s\n' "$CANDIDATE_MIN_READS"
    printf 'offtarget_search_window\t%s\n' "$OFFTARGET_SEARCH_WINDOW"
    printf 'offtarget_max_score\t%s\n' "$OFFTARGET_MAX_SCORE"
    printf 'spacer_neighborhood\t%s\n' "$SPACER_NEIGHBORHOOD"
} > "${post_dir}/run_parameters.tsv"

"$PYTHON_BIN" "$STATE_TOOL" record --stage module05 "${state_args[@]}"
log_step "Phase B complete (modules 03-05); full paired analysis complete"

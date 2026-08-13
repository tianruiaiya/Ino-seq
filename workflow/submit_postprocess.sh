#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'USAGE'
Usage:
  bash workflow/submit_postprocess.sh [pairs.tsv] [config.env]

Maintenance interface for Phase B (modules 03-05). For a standard new run use:
  ./inoseq submit

Environment variables:
  DRY_RUN=1                       Print commands without submitting.
  SKIP_VALIDATION=1               Skip pre-submission validation.
  SUBMIT_OFFTARGET_QC_AFTER=1     Submit cohort summary after all pairs succeed.
USAGE
}

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
PAIR_SHEET=${1:-"${PROJECT_DIR}/config/pairs.tsv"}
CONFIG_FILE=${2:-"${PROJECT_DIR}/config/inoseq.env"}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then usage; exit 0; fi
[[ -f "$PAIR_SHEET" ]] || { echo "[ERROR] Pair sheet not found: $PAIR_SHEET" >&2; exit 2; }
[[ -f "$CONFIG_FILE" ]] || { echo "[ERROR] Config not found: $CONFIG_FILE" >&2; exit 2; }
PAIR_SHEET=$(cd "$(dirname "$PAIR_SHEET")" && pwd)/$(basename "$PAIR_SHEET")
CONFIG_FILE=$(cd "$(dirname "$CONFIG_FILE")" && pwd)/$(basename "$CONFIG_FILE")

if [[ "${SKIP_VALIDATION:-0}" != "1" ]]; then
    "$PROJECT_DIR/workflow/utils/validate_postprocess.sh" "$PAIR_SHEET" "$CONFIG_FILE"
fi

# shellcheck disable=SC1090
source "$CONFIG_FILE"
: "${POSTPROCESS_SLURM_CPUS:=8}"
: "${POSTPROCESS_SLURM_MEM:=32G}"
: "${POSTPROCESS_SLURM_EXCLUDE:=${SLURM_EXCLUDE:-}}"
: "${POSTPROCESS_QC_SLURM_CPUS:=2}"
: "${POSTPROCESS_QC_SLURM_MEM:=8G}"
: "${PYTHON_BIN:=python}"
: "${OUTPUT_DIR:=output}"
: "${LOG_DIR:=logs}"
if [[ "$OUTPUT_DIR" != /* ]]; then OUTPUT_DIR="${PROJECT_DIR}/${OUTPUT_DIR}"; fi
if [[ "$LOG_DIR" != /* ]]; then LOG_DIR="${PROJECT_DIR}/${LOG_DIR}"; fi
mkdir -p "$LOG_DIR"

job_ids=()
submitted=0
while IFS=$'\t' read -r sample_id control_id sgrna extra; do
    [[ -z "${sample_id// }" || "$sample_id" == \#* || "$sample_id" == "sample_id" ]] && continue
    cmd=(
        sbatch --parsable
        --job-name="inoseq_post_${sample_id}"
        --nodes=1 --ntasks=1
        --cpus-per-task="$POSTPROCESS_SLURM_CPUS" --mem="$POSTPROCESS_SLURM_MEM"
        --output="${LOG_DIR}/%x_%j.out" --error="${LOG_DIR}/%x_%j.err"
        --chdir="$PROJECT_DIR" --export="ALL,INOSEQ_PROJECT_DIR=${PROJECT_DIR}"
    )
    [[ -n "$POSTPROCESS_SLURM_EXCLUDE" ]] && cmd+=(--exclude="$POSTPROCESS_SLURM_EXCLUDE")
    cmd+=("${PROJECT_DIR}/workflow/run_postprocess.sh" "$sample_id" "$control_id" "$sgrna" "$CONFIG_FILE")

    if [[ "${DRY_RUN:-0}" == "1" ]]; then
        printf '[DRY-RUN] '; printf '%q ' "${cmd[@]}"; printf '\n'
    else
        job_id=$("${cmd[@]}")
        job_id=${job_id%%;*}
        job_ids+=("$job_id")
        echo "[SUBMITTED] ${sample_id} vs ${control_id}: job ${job_id}"
    fi
    ((submitted+=1))
done < "$PAIR_SHEET"
echo "[INFO] ${submitted} pair job(s) processed from ${PAIR_SHEET}"

if [[ "${DRY_RUN:-0}" != "1" && ${#job_ids[@]} -gt 0 ]]; then
    printf '%s\n' "${job_ids[@]}" > "${LOG_DIR}/inoseq_last_postprocess_job_ids.txt"
fi

if [[ "${SUBMIT_OFFTARGET_QC_AFTER:-1}" == "1" ]]; then
    if [[ "${DRY_RUN:-0}" == "1" ]]; then
        echo "[DRY-RUN] Cohort off-target summary will use afterok dependency"
    elif (( ${#job_ids[@]} > 0 )); then
        dependency=$(IFS=:; echo "${job_ids[*]}")
        qc_dir="${OUTPUT_DIR}/QC"
        mkdir -p "$qc_dir"
        qc_cmd=(
            "$PYTHON_BIN" "$PROJECT_DIR/workflow/qc/collect_offtarget_stats.py"
            --pairs "$PAIR_SHEET" --output-dir "$OUTPUT_DIR" --tsv-dir "$qc_dir"
            --excel "${qc_dir}/dependent_target_analysis.xlsx" --strict
        )
        printf -v qc_wrap '%q ' "${qc_cmd[@]}"
        qc_job=$(sbatch --parsable --job-name=inoseq_offtarget_qc --nodes=1 --ntasks=1 \
            --cpus-per-task="$POSTPROCESS_QC_SLURM_CPUS" --mem="$POSTPROCESS_QC_SLURM_MEM" \
            --dependency="afterok:${dependency}" --output="${LOG_DIR}/%x_%j.out" \
            --error="${LOG_DIR}/%x_%j.err" --chdir="$PROJECT_DIR" \
            --export="ALL,INOSEQ_PROJECT_DIR=${PROJECT_DIR}" --wrap "$qc_wrap")
        echo "[SUBMITTED] cohort off-target summary: job ${qc_job%%;*}"
    fi
fi

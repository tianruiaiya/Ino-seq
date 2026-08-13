#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'USAGE'
Usage:
  bash workflow/submit_full_workflow.sh [--from-stage STAGE] \
    [samples.tsv] [pairs.tsv] [config.env]

Public interface:
  ./inoseq submit [--from-stage STAGE] [samples.tsv] [pairs.tsv] [config.env]

Stages:
  auto       Resume at the earliest incomplete/stale module (default)
  module01   Recompute modules 01-05, cohort summaries and finalization
  module02   Recompute modules 02-05, cohort summaries and finalization
  module03   Recompute modules 03-05, off-target summary and finalization
  module04   Recompute modules 04-05, off-target summary and finalization
  module05   Recompute module 05, off-target summary and finalization
  qc         Recompute cohort summaries and finalization only
  finalize   Re-run final verification only

Aliases: phase-a=module01, phase-b=module03.

Reuse requires a completion marker, all required outputs, and a current stage
fingerprint.  A requested later stage is rejected when its prerequisites are
incomplete or stale.

Environment variables:
  DRY_RUN=1          Print the dynamic job graph without calling sbatch.
  SKIP_VALIDATION=1  Skip full pre-submission validation (development only).
USAGE
}

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
FROM_STAGE=${RESUME_FROM_STAGE:-auto}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --from-stage)
            [[ $# -ge 2 ]] || { echo "[ERROR] --from-stage requires a value" >&2; exit 2; }
            FROM_STAGE=$2
            shift 2
            ;;
        --force)
            FROM_STAGE=module01
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            break
            ;;
        -*)
            echo "[ERROR] Unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
        *) break ;;
    esac
done

case "$FROM_STAGE" in
    auto) ;;
    phase-a|module01|01) FROM_STAGE=module01 ;;
    module02|02) FROM_STAGE=module02 ;;
    phase-b|module03|03) FROM_STAGE=module03 ;;
    module04|04) FROM_STAGE=module04 ;;
    module05|05) FROM_STAGE=module05 ;;
    qc|finalize) ;;
    *)
        echo "[ERROR] Unsupported --from-stage value: $FROM_STAGE" >&2
        usage >&2
        exit 2
        ;;
esac

[[ $# -le 3 ]] || { usage >&2; exit 2; }
SAMPLE_SHEET=${1:-"${PROJECT_DIR}/config/samples.tsv"}
PAIR_SHEET=${2:-"${PROJECT_DIR}/config/pairs.tsv"}
CONFIG_FILE=${3:-"${PROJECT_DIR}/config/inoseq.env"}

for path in "$SAMPLE_SHEET" "$PAIR_SHEET" "$CONFIG_FILE"; do
    [[ -f "$path" ]] || { echo "[ERROR] Required configuration file not found: $path" >&2; exit 2; }
done
SAMPLE_SHEET=$(cd "$(dirname "$SAMPLE_SHEET")" && pwd)/$(basename "$SAMPLE_SHEET")
PAIR_SHEET=$(cd "$(dirname "$PAIR_SHEET")" && pwd)/$(basename "$PAIR_SHEET")
CONFIG_FILE=$(cd "$(dirname "$CONFIG_FILE")" && pwd)/$(basename "$CONFIG_FILE")

if [[ "${SKIP_VALIDATION:-0}" != "1" ]]; then
    "${PROJECT_DIR}/workflow/utils/validate_full_workflow.sh" \
        "$SAMPLE_SHEET" "$PAIR_SHEET" "$CONFIG_FILE"
fi

# shellcheck disable=SC1090
source "$CONFIG_FILE"
: "${PYTHON_BIN:=python}"
: "${OUTPUT_DIR:=output}"
: "${LOG_DIR:=logs}"
: "${SLURM_CPUS:=40}"
: "${SLURM_MEM:=200G}"
: "${SLURM_EXCLUDE:=}"
: "${QC_SLURM_CPUS:=8}"
: "${QC_SLURM_MEM:=16G}"
: "${POSTPROCESS_SLURM_CPUS:=8}"
: "${POSTPROCESS_SLURM_MEM:=32G}"
: "${POSTPROCESS_SLURM_EXCLUDE:=${SLURM_EXCLUDE}}"
: "${POSTPROCESS_QC_SLURM_CPUS:=2}"
: "${POSTPROCESS_QC_SLURM_MEM:=8G}"
: "${FINALIZE_SLURM_CPUS:=1}"
: "${FINALIZE_SLURM_MEM:=2G}"
: "${SUBMIT_QC_AFTER:=1}"
: "${SUBMIT_OFFTARGET_QC_AFTER:=1}"

for toggle_name in SUBMIT_QC_AFTER SUBMIT_OFFTARGET_QC_AFTER; do
    toggle_value=${!toggle_name}
    [[ "$toggle_value" == "0" || "$toggle_value" == "1" ]] \
        || { echo "[ERROR] ${toggle_name} must be 0 or 1: ${toggle_value}" >&2; exit 2; }
done

if [[ "$OUTPUT_DIR" != /* ]]; then OUTPUT_DIR="${PROJECT_DIR}/${OUTPUT_DIR}"; fi
if [[ "$LOG_DIR" != /* ]]; then LOG_DIR="${PROJECT_DIR}/${LOG_DIR}"; fi
mkdir -p "$LOG_DIR"

STATE_TOOL="${PROJECT_DIR}/workflow/utils/stage_state.py"
cohort_state_args=(
    --project-dir "$PROJECT_DIR" --config "$CONFIG_FILE"
    --samples-sheet "$SAMPLE_SHEET" --pairs-sheet "$PAIR_SHEET"
)

declare -a sample_order=() pair_samples=() pair_controls=() pair_sgrnas=()
declare -A sample_r1=() sample_r2=() sample_start=() pair_start=()

while IFS=$'\t' read -r sample_id read1 read2 extra; do
    [[ "$sample_id" == "sample_id" || -z "${sample_id// }" || "$sample_id" == \#* ]] && continue
    [[ -z "${sample_r1[$sample_id]:-}" ]] || { echo "[ERROR] Duplicate sample_id: $sample_id" >&2; exit 2; }
    sample_order+=("$sample_id")
    sample_r1["$sample_id"]=$read1
    sample_r2["$sample_id"]=$read2
done < "$SAMPLE_SHEET"

while IFS=$'\t' read -r sample_id control_id sgrna extra; do
    [[ "$sample_id" == "sample_id" || -z "${sample_id// }" || "$sample_id" == \#* ]] && continue
    [[ -n "${sample_r1[$sample_id]:-}" ]] \
        || { echo "[ERROR] Pair experiment is not in sample sheet: $sample_id" >&2; exit 2; }
    [[ -n "${sample_r1[$control_id]:-}" ]] \
        || { echo "[ERROR] Pair control is not in sample sheet: $control_id" >&2; exit 2; }
    pair_samples+=("$sample_id")
    pair_controls+=("$control_id")
    pair_sgrnas+=("$sgrna")
done < "$PAIR_SHEET"

(( ${#sample_order[@]} > 0 )) || { echo "[ERROR] Sample sheet contains no data rows" >&2; exit 2; }
(( ${#pair_samples[@]} > 0 )) || { echo "[ERROR] Pair sheet contains no data rows" >&2; exit 2; }

sample_state_args() {
    local target=$1
    SAMPLE_STATE_ARGS=(
        --project-dir "$PROJECT_DIR" --config "$CONFIG_FILE"
        --sample "$target" --read1 "${sample_r1[$target]}" --read2 "${sample_r2[$target]}"
    )
}

pair_state_args() {
    local target=$1 control=$2 guide=$3
    PAIR_STATE_ARGS=(
        --project-dir "$PROJECT_DIR" --config "$CONFIG_FILE"
        --sample "$target" --control "$control" --sgrna "$guide"
    )
}

state_valid() {
    local stage=$1
    shift
    "$PYTHON_BIN" "$STATE_TOOL" check --stage "$stage" "$@" --quiet >/dev/null 2>&1
}

prerequisite_error() {
    echo "[ERROR] Cannot start from $1: prerequisite $2 is incomplete or stale for $3" >&2
    echo "[HINT] Use --from-stage $2, or use automatic resume without --from-stage." >&2
    exit 2
}

scheduled_sample_count=0
for sample_id in "${sample_order[@]}"; do
    sample_state_args "$sample_id"
    case "$FROM_STAGE" in
        auto)
            if state_valid module02 "${SAMPLE_STATE_ARGS[@]}"; then
                sample_start["$sample_id"]=complete
            elif state_valid module01 "${SAMPLE_STATE_ARGS[@]}"; then
                sample_start["$sample_id"]=02
                ((scheduled_sample_count+=1))
            else
                sample_start["$sample_id"]=01
                ((scheduled_sample_count+=1))
            fi
            ;;
        module01)
            sample_start["$sample_id"]=01
            ((scheduled_sample_count+=1))
            ;;
        module02)
            state_valid module01 "${SAMPLE_STATE_ARGS[@]}" \
                || prerequisite_error module02 module01 "$sample_id"
            sample_start["$sample_id"]=02
            ((scheduled_sample_count+=1))
            ;;
        *)
            state_valid module02 "${SAMPLE_STATE_ARGS[@]}" \
                || prerequisite_error "$FROM_STAGE" module01 "$sample_id"
            sample_start["$sample_id"]=complete
            ;;
    esac
done

scheduled_pair_count=0
for index in "${!pair_samples[@]}"; do
    sample_id=${pair_samples[$index]}
    control_id=${pair_controls[$index]}
    sgrna=${pair_sgrnas[$index]}
    pair_state_args "$sample_id" "$control_id" "$sgrna"
    key=$index
    case "$FROM_STAGE" in
        auto)
            if [[ "${sample_start[$sample_id]}" != complete \
                || "${sample_start[$control_id]}" != complete ]]; then
                pair_start[$key]=03
            elif state_valid module05 "${PAIR_STATE_ARGS[@]}"; then
                pair_start[$key]=complete
            elif state_valid module04 "${PAIR_STATE_ARGS[@]}"; then
                pair_start[$key]=05
            elif state_valid module03 "${PAIR_STATE_ARGS[@]}"; then
                pair_start[$key]=04
            else
                pair_start[$key]=03
            fi
            ;;
        module01|module02|module03) pair_start[$key]=03 ;;
        module04)
            state_valid module03 "${PAIR_STATE_ARGS[@]}" \
                || prerequisite_error module04 module03 "$sample_id"
            pair_start[$key]=04
            ;;
        module05)
            state_valid module04 "${PAIR_STATE_ARGS[@]}" \
                || prerequisite_error module05 module04 "$sample_id"
            pair_start[$key]=05
            ;;
        qc|finalize)
            state_valid module05 "${PAIR_STATE_ARGS[@]}" \
                || prerequisite_error "$FROM_STAGE" module03 "$sample_id"
            pair_start[$key]=complete
            ;;
    esac
    [[ "${pair_start[$key]}" == complete ]] || ((scheduled_pair_count+=1))
done

phase_a_qc_valid=0
phase_b_qc_valid=0
state_valid phase-a-qc "${cohort_state_args[@]}" && phase_a_qc_valid=1
state_valid phase-b-qc "${cohort_state_args[@]}" && phase_b_qc_valid=1

need_phase_a_qc=0
need_phase_b_qc=0
if [[ "$SUBMIT_QC_AFTER" == "1" ]]; then
    case "$FROM_STAGE" in
        qc|module01|module02) need_phase_a_qc=1 ;;
        finalize)
            (( phase_a_qc_valid == 1 )) \
                || prerequisite_error finalize phase-a-qc cohort
            ;;
        *)
            (( scheduled_sample_count > 0 || phase_a_qc_valid == 0 )) && need_phase_a_qc=1
            ;;
    esac
fi
if [[ "$SUBMIT_OFFTARGET_QC_AFTER" == "1" ]]; then
    case "$FROM_STAGE" in
        qc|module01|module02|module03|module04|module05) need_phase_b_qc=1 ;;
        finalize)
            (( phase_b_qc_valid == 1 )) \
                || prerequisite_error finalize phase-b-qc cohort
            ;;
        auto)
            (( scheduled_pair_count > 0 || phase_b_qc_valid == 0 )) && need_phase_b_qc=1
            ;;
    esac
fi

final_valid=0
state_valid finalize "${cohort_state_args[@]}" && final_valid=1
need_finalize=0
if [[ "$FROM_STAGE" != auto \
    || $scheduled_sample_count -gt 0 \
    || $scheduled_pair_count -gt 0 \
    || $need_phase_a_qc -eq 1 \
    || $need_phase_b_qc -eq 1 \
    || $final_valid -eq 0 ]]; then
    need_finalize=1
fi

echo "[RESUME] requested stage: $FROM_STAGE"
for sample_id in "${sample_order[@]}"; do
    if [[ "${sample_start[$sample_id]}" == complete ]]; then
        echo "[REUSE] Phase A ${sample_id}: modules 01-02 are current"
    else
        echo "[RUN]   Phase A ${sample_id}: start at module ${sample_start[$sample_id]}"
    fi
done
for index in "${!pair_samples[@]}"; do
    if [[ "${pair_start[$index]}" == complete ]]; then
        echo "[REUSE] Phase B ${pair_samples[$index]} vs ${pair_controls[$index]}: modules 03-05 are current"
    else
        echo "[RUN]   Phase B ${pair_samples[$index]} vs ${pair_controls[$index]}: start at module ${pair_start[$index]}"
    fi
done

declare -a graph_rows=() step_job_ids=() post_job_ids=()
declare -A sample_job=()

for sample_id in "${sample_order[@]}"; do
    if [[ "${sample_start[$sample_id]}" == complete ]]; then
        graph_rows+=("SKIP\tphase_a\t${sample_id}\tcurrent")
    fi
done
for index in "${!pair_samples[@]}"; do
    if [[ "${pair_start[$index]}" == complete ]]; then
        graph_rows+=("SKIP\tphase_b\t${pair_samples[$index]}\tcurrent")
    fi
done

if [[ "$need_finalize" == "0" ]]; then
    graph_rows+=("SKIP\tfinalize\tworkflow\tcurrent")
    echo "[RESULT] All requested stages are already complete and current; no Slurm jobs submitted."
    if [[ "${DRY_RUN:-0}" != "1" ]]; then
        rm -f "${LOG_DIR}/inoseq_last_full_workflow_job_id.txt"
        {
            printf 'job_id\tstage\ttarget\tdependency\n'
            printf '%b\n' "${graph_rows[@]}"
        } > "${LOG_DIR}/inoseq_last_full_workflow_jobs.tsv"
    fi
    exit 0
fi

if [[ "${DRY_RUN:-0}" != "1" ]]; then
    for sample_id in "${sample_order[@]}"; do
        [[ "${sample_start[$sample_id]}" == complete ]] && continue
        sample_state_args "$sample_id"
        "$PYTHON_BIN" "$STATE_TOOL" invalidate \
            --stage "module${sample_start[$sample_id]}" "${SAMPLE_STATE_ARGS[@]}" >/dev/null
    done
    for index in "${!pair_samples[@]}"; do
        [[ "${pair_start[$index]}" == complete ]] && continue
        pair_state_args "${pair_samples[$index]}" "${pair_controls[$index]}" "${pair_sgrnas[$index]}"
        "$PYTHON_BIN" "$STATE_TOOL" invalidate \
            --stage "module${pair_start[$index]}" "${PAIR_STATE_ARGS[@]}" >/dev/null
    done
    (( need_phase_a_qc == 0 )) || "$PYTHON_BIN" "$STATE_TOOL" invalidate \
        --stage phase-a-qc "${cohort_state_args[@]}" >/dev/null
    (( need_phase_b_qc == 0 )) || "$PYTHON_BIN" "$STATE_TOOL" invalidate \
        --stage phase-b-qc "${cohort_state_args[@]}" >/dev/null
    "$PYTHON_BIN" "$STATE_TOOL" invalidate \
        --stage finalize "${cohort_state_args[@]}" >/dev/null
fi

dry_counter=0
SUBMITTED_JOB_ID=
submit_job() {
    local label=$1
    shift
    if [[ "${DRY_RUN:-0}" == "1" ]]; then
        ((++dry_counter))
        SUBMITTED_JOB_ID="DRY${dry_counter}"
        printf '[DRY-RUN] %s -> %s\n' "$label" "$SUBMITTED_JOB_ID"
        printf '          '
        printf '%q ' "$@"
        printf '\n'
    else
        local response
        response=$("$@")
        SUBMITTED_JOB_ID=${response%%;*}
        [[ "$SUBMITTED_JOB_ID" =~ ^[0-9]+([_.][0-9]+)?$ ]] \
            || { echo "[ERROR] Unexpected sbatch response for ${label}: ${response}" >&2; exit 1; }
        printf '[SUBMITTED] %s: job %s\n' "$label" "$SUBMITTED_JOB_ID"
    fi
}

join_colon() {
    local IFS=:
    printf '%s' "$*"
}

append_dependency() {
    local -n command_ref=$1
    shift
    if (( $# > 0 )); then
        local dependency
        dependency=$(join_colon "$@")
        command_ref+=(--dependency="afterok:${dependency}")
    fi
}

for sample_id in "${sample_order[@]}"; do
    [[ "${sample_start[$sample_id]}" == complete ]] && continue
    cmd=(
        sbatch --parsable --job-name="inoseq_${sample_id}"
        --nodes=1 --ntasks=1 --cpus-per-task="$SLURM_CPUS" --mem="$SLURM_MEM"
        --output="${LOG_DIR}/%x_%j.out" --error="${LOG_DIR}/%x_%j.err"
        --chdir="$PROJECT_DIR" --export="ALL,INOSEQ_PROJECT_DIR=${PROJECT_DIR}"
    )
    [[ -n "$SLURM_EXCLUDE" ]] && cmd+=(--exclude="$SLURM_EXCLUDE")
    cmd+=(
        "${PROJECT_DIR}/workflow/run_inoseq.sh" --from-module "${sample_start[$sample_id]}"
        "$sample_id" "${sample_r1[$sample_id]}" "${sample_r2[$sample_id]}" "$CONFIG_FILE"
    )
    submit_job "Phase A module ${sample_start[$sample_id]}+ ${sample_id}" "${cmd[@]}"
    sample_job["$sample_id"]=$SUBMITTED_JOB_ID
    step_job_ids+=("$SUBMITTED_JOB_ID")
    graph_rows+=("${SUBMITTED_JOB_ID}\tphase_a\t${sample_id}\t-")
done

step_qc_job=
if (( need_phase_a_qc )); then
    cmd=(
        sbatch --parsable --job-name=inoseq_qc --nodes=1 --ntasks=1
        --cpus-per-task="$QC_SLURM_CPUS" --mem="$QC_SLURM_MEM"
        --output="${LOG_DIR}/%x_%j.out" --error="${LOG_DIR}/%x_%j.err"
        --chdir="$PROJECT_DIR" --export="ALL,INOSEQ_PROJECT_DIR=${PROJECT_DIR}"
    )
    append_dependency cmd "${step_job_ids[@]}"
    [[ -n "$SLURM_EXCLUDE" ]] && cmd+=(--exclude="$SLURM_EXCLUDE")
    cmd+=("${PROJECT_DIR}/workflow/utils/run_cohort_stage.sh" phase-a-qc \
        "$SAMPLE_SHEET" "$PAIR_SHEET" "$CONFIG_FILE")
    submit_job "Phase A cohort QC" "${cmd[@]}"
    step_qc_job=$SUBMITTED_JOB_ID
    dependency=current
    (( ${#step_job_ids[@]} == 0 )) || dependency="afterok:$(join_colon "${step_job_ids[@]}")"
    graph_rows+=("${SUBMITTED_JOB_ID}\tphase_a_qc\tcohort\t${dependency}")
else
    graph_rows+=("SKIP\tphase_a_qc\tcohort\tcurrent")
fi

for index in "${!pair_samples[@]}"; do
    [[ "${pair_start[$index]}" == complete ]] && continue
    sample_id=${pair_samples[$index]}
    control_id=${pair_controls[$index]}
    sgrna=${pair_sgrnas[$index]}
    declare -a pair_deps=()
    [[ -z "${sample_job[$sample_id]:-}" ]] || pair_deps+=("${sample_job[$sample_id]}")
    [[ -z "${sample_job[$control_id]:-}" ]] || pair_deps+=("${sample_job[$control_id]}")
    cmd=(
        sbatch --parsable --job-name="inoseq_post_${sample_id}" --nodes=1 --ntasks=1
        --cpus-per-task="$POSTPROCESS_SLURM_CPUS" --mem="$POSTPROCESS_SLURM_MEM"
        --output="${LOG_DIR}/%x_%j.out" --error="${LOG_DIR}/%x_%j.err"
        --chdir="$PROJECT_DIR" --export="ALL,INOSEQ_PROJECT_DIR=${PROJECT_DIR}"
    )
    append_dependency cmd "${pair_deps[@]}"
    [[ -n "$POSTPROCESS_SLURM_EXCLUDE" ]] && cmd+=(--exclude="$POSTPROCESS_SLURM_EXCLUDE")
    cmd+=(
        "${PROJECT_DIR}/workflow/run_postprocess.sh" --from-module "${pair_start[$index]}"
        "$sample_id" "$control_id" "$sgrna" "$CONFIG_FILE"
    )
    submit_job "Phase B module ${pair_start[$index]}+ ${sample_id} vs ${control_id}" "${cmd[@]}"
    post_job_ids+=("$SUBMITTED_JOB_ID")
    dependency=current
    (( ${#pair_deps[@]} == 0 )) || dependency="afterok:$(join_colon "${pair_deps[@]}")"
    graph_rows+=("${SUBMITTED_JOB_ID}\tphase_b\t${sample_id}\t${dependency}")
done

offtarget_qc_job=
if (( need_phase_b_qc )); then
    cmd=(
        sbatch --parsable --job-name=inoseq_offtarget_qc --nodes=1 --ntasks=1
        --cpus-per-task="$POSTPROCESS_QC_SLURM_CPUS" --mem="$POSTPROCESS_QC_SLURM_MEM"
        --output="${LOG_DIR}/%x_%j.out" --error="${LOG_DIR}/%x_%j.err"
        --chdir="$PROJECT_DIR" --export="ALL,INOSEQ_PROJECT_DIR=${PROJECT_DIR}"
    )
    append_dependency cmd "${post_job_ids[@]}"
    [[ -n "$POSTPROCESS_SLURM_EXCLUDE" ]] && cmd+=(--exclude="$POSTPROCESS_SLURM_EXCLUDE")
    cmd+=("${PROJECT_DIR}/workflow/utils/run_cohort_stage.sh" phase-b-qc \
        "$SAMPLE_SHEET" "$PAIR_SHEET" "$CONFIG_FILE")
    submit_job "off-target cohort summary" "${cmd[@]}"
    offtarget_qc_job=$SUBMITTED_JOB_ID
    dependency=current
    (( ${#post_job_ids[@]} == 0 )) || dependency="afterok:$(join_colon "${post_job_ids[@]}")"
    graph_rows+=("${SUBMITTED_JOB_ID}\tphase_b_qc\tcohort\t${dependency}")
else
    graph_rows+=("SKIP\tphase_b_qc\tcohort\tcurrent")
fi

declare -a final_dependencies=()
if [[ -n "$step_qc_job" ]]; then
    final_dependencies+=("$step_qc_job")
elif [[ "$SUBMIT_QC_AFTER" == "0" ]]; then
    final_dependencies+=("${step_job_ids[@]}")
fi
if [[ -n "$offtarget_qc_job" ]]; then
    final_dependencies+=("$offtarget_qc_job")
elif [[ "$SUBMIT_OFFTARGET_QC_AFTER" == "0" ]]; then
    final_dependencies+=("${post_job_ids[@]}")
fi

cmd=(
    sbatch --parsable --job-name=inoseq_finalize --nodes=1 --ntasks=1
    --cpus-per-task="$FINALIZE_SLURM_CPUS" --mem="$FINALIZE_SLURM_MEM"
    --output="${LOG_DIR}/%x_%j.out" --error="${LOG_DIR}/%x_%j.err"
    --chdir="$PROJECT_DIR" --export="ALL,INOSEQ_PROJECT_DIR=${PROJECT_DIR}"
)
append_dependency cmd "${final_dependencies[@]}"
[[ -n "$SLURM_EXCLUDE" ]] && cmd+=(--exclude="$SLURM_EXCLUDE")
cmd+=(
    "${PROJECT_DIR}/workflow/utils/finalize_full_workflow.sh"
    "$SAMPLE_SHEET" "$PAIR_SHEET" "$CONFIG_FILE"
    "$SUBMIT_QC_AFTER" "$SUBMIT_OFFTARGET_QC_AFTER"
)
submit_job "Phase C full-workflow finalizer" "${cmd[@]}"
final_job=$SUBMITTED_JOB_ID
dependency=current
(( ${#final_dependencies[@]} == 0 )) || dependency="afterok:$(join_colon "${final_dependencies[@]}")"
graph_rows+=("${final_job}\tfinalize\tworkflow\t${dependency}")

if [[ "${DRY_RUN:-0}" != "1" ]]; then
    if (( ${#step_job_ids[@]} > 0 )); then
        printf '%s\n' "${step_job_ids[@]}" > "${LOG_DIR}/inoseq_last_sample_job_ids.txt"
    else
        rm -f "${LOG_DIR}/inoseq_last_sample_job_ids.txt"
    fi
    if (( ${#post_job_ids[@]} > 0 )); then
        printf '%s\n' "${post_job_ids[@]}" > "${LOG_DIR}/inoseq_last_postprocess_job_ids.txt"
    else
        rm -f "${LOG_DIR}/inoseq_last_postprocess_job_ids.txt"
    fi
    [[ -z "$step_qc_job" ]] \
        && rm -f "${LOG_DIR}/inoseq_last_qc_job_id.txt" \
        || printf '%s\n' "$step_qc_job" > "${LOG_DIR}/inoseq_last_qc_job_id.txt"
    [[ -z "$offtarget_qc_job" ]] \
        && rm -f "${LOG_DIR}/inoseq_last_offtarget_qc_job_id.txt" \
        || printf '%s\n' "$offtarget_qc_job" > "${LOG_DIR}/inoseq_last_offtarget_qc_job_id.txt"
    printf '%s\n' "$final_job" > "${LOG_DIR}/inoseq_last_full_workflow_job_id.txt"
    {
        printf 'job_id\tstage\ttarget\tdependency\n'
        printf '%b\n' "${graph_rows[@]}"
    } > "${LOG_DIR}/inoseq_last_full_workflow_jobs.tsv"
fi

echo "[INFO] Submitted graph: ${scheduled_sample_count} Phase A job(s), ${scheduled_pair_count} Phase B job(s), final job ${final_job}"
echo "[INFO] Success marker: ${OUTPUT_DIR}/INOSEQ_WORKFLOW_COMPLETE"

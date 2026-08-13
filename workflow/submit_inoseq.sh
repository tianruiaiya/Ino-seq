#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'USAGE'
Usage:
  bash workflow/submit_inoseq.sh [samples.tsv] [config.env]

Maintenance interface for Phase A (modules 01-02). For a standard new run use:
  ./inoseq submit

Defaults:
  samples.tsv = config/samples.tsv
  config.env  = config/inoseq.env

Environment variables:
  DRY_RUN=1          Print sbatch commands without submitting jobs.
  SKIP_VALIDATION=1  Skip pre-submission validation.
  SUBMIT_QC_AFTER=1  Submit one QC aggregation job after all sample jobs finish.
USAGE
}

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)

SAMPLE_SHEET=${1:-"${PROJECT_DIR}/config/samples.tsv"}
CONFIG_FILE=${2:-"${PROJECT_DIR}/config/inoseq.env"}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

[[ -f "$CONFIG_FILE" ]] || { echo "[ERROR] Configuration file not found: $CONFIG_FILE" >&2; exit 2; }
[[ -f "$SAMPLE_SHEET" ]] || { echo "[ERROR] Sample sheet not found: $SAMPLE_SHEET" >&2; exit 2; }

CONFIG_FILE=$(cd "$(dirname "$CONFIG_FILE")" && pwd)/$(basename "$CONFIG_FILE")
SAMPLE_SHEET=$(cd "$(dirname "$SAMPLE_SHEET")" && pwd)/$(basename "$SAMPLE_SHEET")

if [[ "${SKIP_VALIDATION:-0}" != "1" ]]; then
    "${PROJECT_DIR}/workflow/utils/validate_config.sh" "$SAMPLE_SHEET" "$CONFIG_FILE"
fi

# shellcheck disable=SC1090
source "$CONFIG_FILE"

: "${SLURM_CPUS:=40}"
: "${SLURM_MEM:=200G}"
: "${SLURM_EXCLUDE:=}"
: "${QC_SLURM_CPUS:=8}"
: "${QC_SLURM_MEM:=16G}"
: "${PYTHON_BIN:=python}"
: "${OUTPUT_DIR:=output}"
: "${LOG_DIR:=logs}"

if [[ "$OUTPUT_DIR" != /* ]]; then
    OUTPUT_DIR="${PROJECT_DIR}/${OUTPUT_DIR}"
fi
if [[ "$LOG_DIR" != /* ]]; then
    LOG_DIR="${PROJECT_DIR}/${LOG_DIR}"
fi

mkdir -p "$LOG_DIR"

submitted=0
job_ids=()

while IFS=$'\t' read -r sample_id read1 read2 extra; do
    [[ -z "${sample_id// }" || "$sample_id" == \#* || "$sample_id" == "sample_id" ]] && continue

    sbatch_args=(
        --job-name="inoseq_${sample_id}"
        --nodes=1
        --ntasks=1
        --cpus-per-task="${SLURM_CPUS}"
        --mem="${SLURM_MEM}"
        --output="${LOG_DIR}/%x_%j.out"
        --error="${LOG_DIR}/%x_%j.err"
        --chdir="${PROJECT_DIR}"
        --export="ALL,INOSEQ_PROJECT_DIR=${PROJECT_DIR}"
    )
    [[ -n "${SLURM_EXCLUDE}" ]] && sbatch_args+=(--exclude="${SLURM_EXCLUDE}")

    cmd=(
        sbatch --parsable
        "${sbatch_args[@]}"
        "${PROJECT_DIR}/workflow/run_inoseq.sh"
        "$sample_id" "$read1" "$read2" "$CONFIG_FILE"
    )

    if [[ "${DRY_RUN:-0}" == "1" ]]; then
        printf '[DRY-RUN] '
        printf '%q ' "${cmd[@]}"
        printf '\n'
    else
        job_id=$("${cmd[@]}")
        job_id=${job_id%%;*}
        job_ids+=("$job_id")
        echo "[SUBMITTED] ${sample_id}: job ${job_id}"
    fi
    ((submitted+=1))
done < "$SAMPLE_SHEET"

echo "[INFO] ${submitted} sample job(s) processed from ${SAMPLE_SHEET}"

if [[ "${DRY_RUN:-0}" != "1" && ${#job_ids[@]} -gt 0 ]]; then
    printf '%s\n' "${job_ids[@]}" > "${LOG_DIR}/inoseq_last_sample_job_ids.txt"
fi

if [[ "${SUBMIT_QC_AFTER:-0}" == "1" ]]; then
    if [[ "${DRY_RUN:-0}" == "1" ]]; then
        echo "[DRY-RUN] QC aggregation job will use afterok dependency on all sample jobs"
    elif (( ${#job_ids[@]} > 0 )); then
        dep=$(IFS=:; echo "${job_ids[*]}")
        QC_DIR="${OUTPUT_DIR}/QC"
        mkdir -p "$QC_DIR"

        qc_args=(
            "$PYTHON_BIN"
            "${PROJECT_DIR}/workflow/qc/collect_qc.py"
            --samples "$SAMPLE_SHEET"
            --output-dir "$OUTPUT_DIR"
            --config "$CONFIG_FILE"
            --tsv "${QC_DIR}/inoseq_qc_summary.tsv"
            --csv "${QC_DIR}/inoseq_qc_summary.csv"
            --strict
        )
        printf -v qc_wrap '%q ' "${qc_args[@]}"

        qc_sbatch=(
            sbatch --parsable
            --job-name="inoseq_qc"
            --nodes=1
            --ntasks=1
            --cpus-per-task="${QC_SLURM_CPUS}"
            --mem="${QC_SLURM_MEM}"
            --dependency="afterok:${dep}"
            --output="${LOG_DIR}/%x_%j.out"
            --error="${LOG_DIR}/%x_%j.err"
            --chdir="${PROJECT_DIR}"
            --export="ALL,INOSEQ_PROJECT_DIR=${PROJECT_DIR}"
            --wrap "$qc_wrap"
        )
        [[ -n "${SLURM_EXCLUDE}" ]] && qc_sbatch+=(--exclude="${SLURM_EXCLUDE}")

        qc_job=$("${qc_sbatch[@]}")
        qc_job=${qc_job%%;*}
        echo "[SUBMITTED] aggregate QC: job ${qc_job} (afterok: ${dep})"
        printf '%s\n' "$qc_job" > "${LOG_DIR}/inoseq_last_qc_job_id.txt"
    fi
fi

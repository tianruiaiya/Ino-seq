#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'USAGE'
Usage:
  workflow/utils/finalize_full_workflow.sh \
    SAMPLES_TSV PAIRS_TSV CONFIG_ENV [EXPECT_STEP01_QC] [EXPECT_OFFTARGET_QC]

This command is submitted as the final afterok-dependent Slurm job. It verifies
all stage markers and cohort outputs, then writes the workflow-level status and
completion marker.
USAGE
}

if [[ $# -lt 3 || $# -gt 5 ]]; then
    usage >&2
    exit 2
fi

SAMPLE_SHEET=$1
PAIR_SHEET=$2
CONFIG_FILE=$3
EXPECT_STEP01_QC=${4:-1}
EXPECT_OFFTARGET_QC=${5:-1}

if [[ -n "${INOSEQ_PROJECT_DIR:-}" ]]; then
    PROJECT_DIR=$(cd "$INOSEQ_PROJECT_DIR" && pwd)
else
    SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
    PROJECT_DIR=$(cd "${SCRIPT_DIR}/../.." && pwd)
fi

for path in "$SAMPLE_SHEET" "$PAIR_SHEET" "$CONFIG_FILE"; do
    [[ -r "$path" ]] || { echo "[ERROR] Required file not readable: $path" >&2; exit 2; }
done

# shellcheck disable=SC1090
source "$CONFIG_FILE"
: "${OUTPUT_DIR:=output}"
: "${PYTHON_BIN:=python}"
if [[ "$OUTPUT_DIR" != /* ]]; then OUTPUT_DIR="${PROJECT_DIR}/${OUTPUT_DIR}"; fi
INOSEQ_VERSION=$(tr -d '[:space:]' < "${PROJECT_DIR}/VERSION")

failed=0
sample_count=0
pair_count=0
check_file() {
    [[ -f "$1" ]] || { echo "[ERROR] Missing required completion output: $1" >&2; failed=1; }
}

while IFS=$'\t' read -r sample_id read1 read2 extra; do
    [[ "$sample_id" == "sample_id" || -z "${sample_id// }" || "$sample_id" == \#* ]] && continue
    ((sample_count+=1))
    check_file "${OUTPUT_DIR}/${sample_id}/INOSEQ_COMPLETE"
done < "$SAMPLE_SHEET"

while IFS=$'\t' read -r sample_id control_id sgrna extra; do
    [[ "$sample_id" == "sample_id" || -z "${sample_id// }" || "$sample_id" == \#* ]] && continue
    ((pair_count+=1))
    check_file "${OUTPUT_DIR}/${sample_id}/postprocess/INOSEQ_POSTPROCESS_COMPLETE"
    check_file "${OUTPUT_DIR}/${sample_id}/INOSEQ_FULL_COMPLETE"
done < "$PAIR_SHEET"

if [[ "$EXPECT_STEP01_QC" == "1" ]]; then
    check_file "${OUTPUT_DIR}/QC/inoseq_qc_summary.tsv"
    check_file "${OUTPUT_DIR}/QC/inoseq_qc_summary.csv"
fi
if [[ "$EXPECT_OFFTARGET_QC" == "1" ]]; then
    check_file "${OUTPUT_DIR}/QC/inoseq_offtarget_summary.tsv"
    check_file "${OUTPUT_DIR}/QC/inoseq_strand_summary.tsv"
    check_file "${OUTPUT_DIR}/QC/dependent_target_analysis.xlsx"
fi

if (( failed )); then
    echo "[RESULT] Full-workflow finalization failed" >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}/QC"
status_file="${OUTPUT_DIR}/QC/full_workflow_status.tsv"
{
    printf 'field\tvalue\n'
    printf 'status\tCOMPLETED\n'
    printf 'inoseq_version\t%s\n' "$INOSEQ_VERSION"
    printf 'completed_at_utc\t%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    printf 'sample_count\t%s\n' "$sample_count"
    printf 'pair_count\t%s\n' "$pair_count"
    printf 'step01_qc_expected\t%s\n' "$EXPECT_STEP01_QC"
    printf 'offtarget_qc_expected\t%s\n' "$EXPECT_OFFTARGET_QC"
    printf 'background_fold_change\t%s\n' "${BACKGROUND_FOLD_CHANGE:-1.5}"
    printf 'background_pvalue\t%s\n' "${BACKGROUND_PVALUE:-0.05}"
} > "$status_file"

touch "${OUTPUT_DIR}/INOSEQ_WORKFLOW_COMPLETE"
"$PYTHON_BIN" "${PROJECT_DIR}/workflow/utils/stage_state.py" record \
    --stage finalize --project-dir "$PROJECT_DIR" --config "$CONFIG_FILE" \
    --samples-sheet "$SAMPLE_SHEET" --pairs-sheet "$PAIR_SHEET"
echo "[RESULT] Ino-seq full workflow completed: ${OUTPUT_DIR}"

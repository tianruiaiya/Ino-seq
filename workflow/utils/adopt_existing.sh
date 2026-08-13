#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'USAGE'
Usage:
  ./inoseq adopt-existing [samples.tsv] [pairs.tsv] [config.env]

Create resumable state records for a complete pre-resume Ino-seq output tree.
No analytical module is executed.  Use this only after confirming that the
existing outputs were produced from the supplied inputs and current config.
USAGE
}

if [[ $# -ne 3 ]]; then
    usage >&2
    exit 2
fi

SAMPLE_SHEET=$1
PAIR_SHEET=$2
CONFIG_FILE=$3
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd "${SCRIPT_DIR}/../.." && pwd)

for path in "$SAMPLE_SHEET" "$PAIR_SHEET" "$CONFIG_FILE"; do
    [[ -r "$path" ]] || { echo "[ERROR] Required file not readable: $path" >&2; exit 2; }
done
SAMPLE_SHEET=$(cd "$(dirname "$SAMPLE_SHEET")" && pwd)/$(basename "$SAMPLE_SHEET")
PAIR_SHEET=$(cd "$(dirname "$PAIR_SHEET")" && pwd)/$(basename "$PAIR_SHEET")
CONFIG_FILE=$(cd "$(dirname "$CONFIG_FILE")" && pwd)/$(basename "$CONFIG_FILE")

# shellcheck disable=SC1090
source "$CONFIG_FILE"
: "${PYTHON_BIN:=python}"
: "${OUTPUT_DIR:=output}"
: "${SUBMIT_QC_AFTER:=1}"
: "${SUBMIT_OFFTARGET_QC_AFTER:=1}"
if [[ "$OUTPUT_DIR" != /* ]]; then OUTPUT_DIR="${PROJECT_DIR}/${OUTPUT_DIR}"; fi

while IFS=$'\t' read -r sample_id read1 read2 extra; do
    [[ "$sample_id" == "sample_id" || -z "${sample_id// }" || "$sample_id" == \#* ]] && continue
    [[ -f "${OUTPUT_DIR}/${sample_id}/INOSEQ_COMPLETE" ]] \
        || { echo "[ERROR] Existing Phase A completion marker is missing: $sample_id" >&2; exit 1; }
done < "$SAMPLE_SHEET"
while IFS=$'\t' read -r sample_id control_id sgrna extra; do
    [[ "$sample_id" == "sample_id" || -z "${sample_id// }" || "$sample_id" == \#* ]] && continue
    [[ -f "${OUTPUT_DIR}/${sample_id}/postprocess/INOSEQ_POSTPROCESS_COMPLETE" \
        && -f "${OUTPUT_DIR}/${sample_id}/INOSEQ_FULL_COMPLETE" ]] \
        || { echo "[ERROR] Existing Phase B completion markers are missing: $sample_id" >&2; exit 1; }
done < "$PAIR_SHEET"

STATE_TOOL="${PROJECT_DIR}/workflow/utils/stage_state.py"
while IFS=$'\t' read -r sample_id read1 read2 extra; do
    [[ "$sample_id" == "sample_id" || -z "${sample_id// }" || "$sample_id" == \#* ]] && continue
    args=(--project-dir "$PROJECT_DIR" --config "$CONFIG_FILE" \
        --sample "$sample_id" --read1 "$read1" --read2 "$read2")
    "$PYTHON_BIN" "$STATE_TOOL" record --stage module01 "${args[@]}"
    "$PYTHON_BIN" "$STATE_TOOL" record --stage module02 "${args[@]}"
done < "$SAMPLE_SHEET"

while IFS=$'\t' read -r sample_id control_id sgrna extra; do
    [[ "$sample_id" == "sample_id" || -z "${sample_id// }" || "$sample_id" == \#* ]] && continue
    args=(--project-dir "$PROJECT_DIR" --config "$CONFIG_FILE" \
        --sample "$sample_id" --control "$control_id" --sgrna "$sgrna")
    for stage in module03 module04 module05; do
        "$PYTHON_BIN" "$STATE_TOOL" record --stage "$stage" "${args[@]}"
    done
done < "$PAIR_SHEET"

cohort_args=(--project-dir "$PROJECT_DIR" --config "$CONFIG_FILE" \
    --samples-sheet "$SAMPLE_SHEET" --pairs-sheet "$PAIR_SHEET")
if [[ "$SUBMIT_QC_AFTER" == "1" ]]; then
    "$PYTHON_BIN" "$STATE_TOOL" record --stage phase-a-qc "${cohort_args[@]}"
fi
if [[ "$SUBMIT_OFFTARGET_QC_AFTER" == "1" ]]; then
    "$PYTHON_BIN" "$STATE_TOOL" record --stage phase-b-qc "${cohort_args[@]}"
fi
if [[ -f "${OUTPUT_DIR}/INOSEQ_WORKFLOW_COMPLETE" ]]; then
    "$PYTHON_BIN" "$STATE_TOOL" record --stage finalize "${cohort_args[@]}"
fi

echo "[RESULT] Existing outputs adopted into the resumable state contract."
echo "[NEXT] Run ./inoseq status and ./inoseq plan to review reuse decisions."

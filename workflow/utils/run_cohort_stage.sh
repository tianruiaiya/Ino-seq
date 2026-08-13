#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'USAGE'
Usage:
  workflow/utils/run_cohort_stage.sh \
    phase-a-qc|phase-b-qc SAMPLES_TSV PAIRS_TSV CONFIG_ENV

Internal resumable runner for cohort aggregation stages.  The stage state is
recorded only after the strict collector exits successfully.
USAGE
}

if [[ $# -ne 4 ]]; then
    usage >&2
    exit 2
fi

STAGE=$1
SAMPLE_SHEET=$2
PAIR_SHEET=$3
CONFIG_FILE=$4

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
: "${PYTHON_BIN:=python}"
: "${OUTPUT_DIR:=output}"
if [[ "$OUTPUT_DIR" != /* ]]; then OUTPUT_DIR="${PROJECT_DIR}/${OUTPUT_DIR}"; fi

STATE_TOOL="${PROJECT_DIR}/workflow/utils/stage_state.py"
state_args=(
    --project-dir "$PROJECT_DIR" --config "$CONFIG_FILE"
    --samples-sheet "$SAMPLE_SHEET" --pairs-sheet "$PAIR_SHEET"
)

"$PYTHON_BIN" "$STATE_TOOL" invalidate --stage "$STAGE" "${state_args[@]}"
mkdir -p "${OUTPUT_DIR}/QC"

case "$STAGE" in
    phase-a-qc)
        "$PYTHON_BIN" "${PROJECT_DIR}/workflow/qc/collect_qc.py" \
            --samples "$SAMPLE_SHEET" --output-dir "$OUTPUT_DIR" --config "$CONFIG_FILE" \
            --tsv "${OUTPUT_DIR}/QC/inoseq_qc_summary.tsv" \
            --csv "${OUTPUT_DIR}/QC/inoseq_qc_summary.csv" --strict
        ;;
    phase-b-qc)
        "$PYTHON_BIN" "${PROJECT_DIR}/workflow/qc/collect_offtarget_stats.py" \
            --pairs "$PAIR_SHEET" --output-dir "$OUTPUT_DIR" \
            --tsv-dir "${OUTPUT_DIR}/QC" \
            --excel "${OUTPUT_DIR}/QC/dependent_target_analysis.xlsx" --strict
        ;;
    *)
        echo "[ERROR] Unsupported cohort stage: $STAGE" >&2
        usage >&2
        exit 2
        ;;
esac

"$PYTHON_BIN" "$STATE_TOOL" record --stage "$STAGE" "${state_args[@]}"
echo "[RESULT] Cohort stage completed: $STAGE"

#!/usr/bin/env bash
set -Eeuo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: report_stage_state.sh SAMPLES_TSV PAIRS_TSV CONFIG_ENV" >&2
    exit 2
fi

SAMPLE_SHEET=$1
PAIR_SHEET=$2
CONFIG_FILE=$3
if [[ -n "${INOSEQ_PROJECT_DIR:-}" ]]; then
    PROJECT_DIR=$(cd "$INOSEQ_PROJECT_DIR" && pwd)
else
    SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
    PROJECT_DIR=$(cd "${SCRIPT_DIR}/../.." && pwd)
fi

# shellcheck disable=SC1090
source "$CONFIG_FILE"
: "${PYTHON_BIN:=python}"
STATE_TOOL="${PROJECT_DIR}/workflow/utils/stage_state.py"

status_label() {
    local stage=$1
    shift
    local code=0
    "$PYTHON_BIN" "$STATE_TOOL" check --stage "$stage" "$@" --quiet \
        >/dev/null 2>&1 || code=$?
    case "$code" in
        0) printf 'CURRENT' ;;
        10) printf 'INCOMPLETE' ;;
        11) printf 'STALE' ;;
        *) printf 'ERROR(%s)' "$code" ;;
    esac
}

printf 'scope\ttarget\tstage\tstate\n'
while IFS=$'\t' read -r sample_id read1 read2 extra; do
    [[ "$sample_id" == "sample_id" || -z "${sample_id// }" || "$sample_id" == \#* ]] && continue
    args=(--project-dir "$PROJECT_DIR" --config "$CONFIG_FILE" \
        --sample "$sample_id" --read1 "$read1" --read2 "$read2")
    for stage in module01 module02; do
        printf 'sample\t%s\t%s\t%s\n' "$sample_id" "$stage" \
            "$(status_label "$stage" "${args[@]}")"
    done
done < "$SAMPLE_SHEET"

while IFS=$'\t' read -r sample_id control_id sgrna extra; do
    [[ "$sample_id" == "sample_id" || -z "${sample_id// }" || "$sample_id" == \#* ]] && continue
    args=(--project-dir "$PROJECT_DIR" --config "$CONFIG_FILE" \
        --sample "$sample_id" --control "$control_id" --sgrna "$sgrna")
    for stage in module03 module04 module05; do
        printf 'pair\t%s_vs_%s\t%s\t%s\n' "$sample_id" "$control_id" "$stage" \
            "$(status_label "$stage" "${args[@]}")"
    done
done < "$PAIR_SHEET"

cohort_args=(--project-dir "$PROJECT_DIR" --config "$CONFIG_FILE" \
    --samples-sheet "$SAMPLE_SHEET" --pairs-sheet "$PAIR_SHEET")
for stage in phase-a-qc phase-b-qc finalize; do
    printf 'cohort\tworkflow\t%s\t%s\n' "$stage" \
        "$(status_label "$stage" "${cohort_args[@]}")"
done

#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'USAGE'
Usage:
  workflow/utils/validate_postprocess.sh [--skip-step01-inputs] [pairs.tsv] [config.env]

Defaults:
  pairs.tsv  = config/pairs.tsv
  config.env = config/inoseq.env

Options:
  --skip-step01-inputs  Validate configuration and pair definitions before
                        Phase A has produced .end/BAM files.
USAGE
}

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd "${SCRIPT_DIR}/../.." && pwd)
CHECK_STEP01_INPUTS=1
if [[ "${1:-}" == "--skip-step01-inputs" ]]; then
    CHECK_STEP01_INPUTS=0
    shift
fi
PAIR_SHEET=${1:-"${PROJECT_DIR}/config/pairs.tsv"}
CONFIG_FILE=${2:-"${PROJECT_DIR}/config/inoseq.env"}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then usage; exit 0; fi

failed=0
ok()   { printf '[OK]   %s\n' "$*"; }
miss() { printf '[FAIL] %s\n' "$*" >&2; failed=1; }

[[ -f "$PAIR_SHEET" ]] || { echo "[FAIL] Pair sheet not found: $PAIR_SHEET" >&2; exit 2; }
[[ -f "$CONFIG_FILE" ]] || { echo "[FAIL] Config not found: $CONFIG_FILE" >&2; exit 2; }
# shellcheck disable=SC1090
source "$CONFIG_FILE"
: "${PYTHON_BIN:=python}"
: "${OUTPUT_DIR:=output}"
: "${BACKGROUND_PVALUE:=0.05}"
: "${BACKGROUND_FOLD_CHANGE:=1.5}"

if [[ "$OUTPUT_DIR" != /* ]]; then OUTPUT_DIR="${PROJECT_DIR}/${OUTPUT_DIR}"; fi
command -v "$PYTHON_BIN" >/dev/null 2>&1 && ok "Python: $(command -v "$PYTHON_BIN")" || miss "Python not found: $PYTHON_BIN"
[[ -r "${REFERENCE_FASTA:-}" ]] && ok "reference FASTA: $REFERENCE_FASTA" || miss "REFERENCE_FASTA is unreadable: ${REFERENCE_FASTA:-<unset>}"
[[ -s "${REFERENCE_FASTA:-}.fai" ]] && ok "reference FASTA index" || miss "missing FASTA index: ${REFERENCE_FASTA:-<unset>}.fai"

module_check=$(mktemp)
if PYTHONPATH="$PROJECT_DIR${PYTHONPATH:+:$PYTHONPATH}" "$PYTHON_BIN" - <<'PY' >"$module_check" 2>&1
import Bio
import openpyxl
import pandas
import pyfaidx
import pysam
import regex
import scipy
print("imports ready")
PY
then
    ok "Python postprocessing dependencies"
else
    miss "Python dependency check failed: $(tr '\n' ' ' < "$module_check")"
fi
rm -f "$module_check"

if "$PYTHON_BIN" - "$BACKGROUND_PVALUE" "$BACKGROUND_FOLD_CHANGE" <<'PY'
import sys
p = float(sys.argv[1])
fc = float(sys.argv[2])
assert 0 <= p <= 1
assert fc >= 0
PY
then
    ok "background thresholds: fold_change >= ${BACKGROUND_FOLD_CHANGE}; p_value < ${BACKGROUND_PVALUE}"
else
    miss "invalid background threshold configuration"
fi

header=$(head -n 1 "$PAIR_SHEET" || true)
[[ "$header" == $'sample_id\tcontrol_id\tsgrna' ]] \
    && ok "pair sheet header" \
    || miss "pair sheet header must be: sample_id<TAB>control_id<TAB>sgrna"

pair_count=0
declare -A seen=()
while IFS=$'\t' read -r sample_id control_id sgrna extra; do
    [[ "$sample_id" == "sample_id" || -z "${sample_id// }" || "$sample_id" == \#* ]] && continue
    ((pair_count+=1))
    [[ -z "${extra:-}" ]] || miss "pair '$sample_id': more than 3 columns"
    [[ "$sample_id" =~ ^[A-Za-z0-9._-]+$ ]] || miss "unsupported sample_id: $sample_id"
    [[ "$control_id" =~ ^[A-Za-z0-9._-]+$ ]] || miss "unsupported control_id: $control_id"
    [[ "$sample_id" != "$control_id" ]] || miss "sample and control must differ: $sample_id"
    [[ "$sgrna" =~ ^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]{20,30}$ ]] || miss "invalid sgRNA for '$sample_id': $sgrna"
    [[ -z "${seen[$sample_id]:-}" ]] || miss "duplicate sample_id in pair sheet: $sample_id"
    seen[$sample_id]=1

    if (( CHECK_STEP01_INPUTS )); then
        sample_prefix="${OUTPUT_DIR}/${sample_id}/alignment/${sample_id}"
        control_prefix="${OUTPUT_DIR}/${control_id}/alignment/${control_id}"
        for path in "${sample_prefix}.end" "${sample_prefix}_end.bam" "${sample_prefix}_end.bam.bai" \
                    "${control_prefix}_end.bam" "${control_prefix}_end.bam.bai"; do
            [[ -r "$path" ]] || miss "required Phase A output not readable: $path"
        done
    fi
done < "$PAIR_SHEET"

(( pair_count > 0 )) && ok "pair rows: $pair_count" || miss "pair sheet contains no data rows"
if (( failed )); then echo "[RESULT] Postprocessing validation failed" >&2; exit 1; fi
if (( CHECK_STEP01_INPUTS )); then
    echo "[RESULT] Ino-seq modules 03-05 inputs are ready"
else
    echo "[RESULT] Ino-seq modules 03-05 configuration is ready"
fi

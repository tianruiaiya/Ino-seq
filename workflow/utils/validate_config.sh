#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'USAGE'
Usage:
  workflow/utils/validate_config.sh [samples.tsv] [config.env]

Defaults:
  samples.tsv = config/samples.tsv
  config.env  = config/inoseq.env
USAGE
}

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd "${SCRIPT_DIR}/../.." && pwd)
SAMPLE_SHEET=${1:-"${PROJECT_DIR}/config/samples.tsv"}
CONFIG_FILE=${2:-"${PROJECT_DIR}/config/inoseq.env"}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

failed=0
ok()   { printf '[OK]   %s\n' "$*"; }
warn() { printf '[WARN] %s\n' "$*"; }
miss() { printf '[FAIL] %s\n' "$*" >&2; failed=1; }

[[ -f "$CONFIG_FILE" ]] || { echo "[FAIL] Configuration file not found: $CONFIG_FILE" >&2; exit 2; }
[[ -f "$SAMPLE_SHEET" ]] || { echo "[FAIL] Sample sheet not found: $SAMPLE_SHEET" >&2; exit 2; }
# shellcheck disable=SC1090
source "$CONFIG_FILE"

for exe in fastp cutadapt bwa samtools "${PYTHON_BIN:-python}" sbatch; do
    if command -v "$exe" >/dev/null 2>&1; then
        ok "command: $exe -> $(command -v "$exe")"
    else
        miss "required command not found: $exe"
    fi
done

if [[ -n "${FGBIO_JAR:-}" ]]; then
    command -v java >/dev/null 2>&1 && ok "command: java -> $(command -v java)" || miss "java not found"
    [[ -r "$FGBIO_JAR" ]] && ok "fgbio JAR: $FGBIO_JAR" || miss "fgbio JAR not readable: $FGBIO_JAR"
else
    command -v fgbio >/dev/null 2>&1 && ok "command: fgbio -> $(command -v fgbio)" || miss "fgbio not found in PATH and FGBIO_JAR is empty"
fi


if [[ -n "${READS_TO_PROCESS:-}" ]]; then
    [[ "$READS_TO_PROCESS" =~ ^[0-9]+$ ]] && ok "READS_TO_PROCESS: $READS_TO_PROCESS" || miss "READS_TO_PROCESS must be a non-negative integer"
fi
if [[ -n "${QC_UNIQUE_MAPQ:-}" ]]; then
    [[ "$QC_UNIQUE_MAPQ" =~ ^[0-9]+$ ]] && ok "QC_UNIQUE_MAPQ: $QC_UNIQUE_MAPQ" || miss "QC_UNIQUE_MAPQ must be a non-negative integer"
fi

pysam_check=$(mktemp)
if "${PYTHON_BIN:-python}" - <<'PY' >"$pysam_check" 2>&1
import pysam
print(pysam.__version__)
PY
then
    ok "pysam: $(cat "$pysam_check")"
else
    miss "pysam import failed: $(cat "$pysam_check")"
fi
rm -f "$pysam_check"

if [[ -n "${REFERENCE_FASTA:-}" && -r "${REFERENCE_FASTA:-}" ]]; then
    ok "reference FASTA: $REFERENCE_FASTA"
else
    miss "REFERENCE_FASTA is missing or unreadable: ${REFERENCE_FASTA:-<unset>}"
fi

if [[ -n "${REFERENCE_FASTA:-}" ]]; then
    for ext in amb ann bwt pac sa; do
        [[ -s "${REFERENCE_FASTA}.${ext}" ]] || miss "missing BWA index: ${REFERENCE_FASTA}.${ext}"
    done
    [[ -s "${REFERENCE_FASTA}.fai" ]] && ok "reference FASTA index: ${REFERENCE_FASTA}.fai" || miss "missing FASTA index: ${REFERENCE_FASTA}.fai"
fi

if [[ -n "${REFERENCE_FASTA:-}" ]]; then
    case "$REFERENCE_FASTA" in
        *.fasta) dict_file="${REFERENCE_FASTA%.fasta}.dict" ;;
        *.fa)    dict_file="${REFERENCE_FASTA%.fa}.dict" ;;
        *.fas)   dict_file="${REFERENCE_FASTA%.fas}.dict" ;;
        *)       dict_file="${REFERENCE_FASTA}.dict" ;;
    esac
    [[ -s "$dict_file" ]] && ok "reference sequence dictionary: $dict_file" || miss "missing sequence dictionary: $dict_file"
fi

header=$(head -n 1 "$SAMPLE_SHEET" || true)
if [[ "$header" == $'sample_id\tread1\tread2' ]]; then
    ok "sample sheet header"
else
    miss "sample sheet header must be exactly: sample_id<TAB>read1<TAB>read2"
fi

sample_count=0
declare -A seen=()
while IFS=$'\t' read -r sample_id read1 read2 extra; do
    [[ "$sample_id" == "sample_id" ]] && continue
    [[ -z "${sample_id// }" ]] && continue
    [[ "$sample_id" == \#* ]] && continue
    ((sample_count+=1))

    [[ "$sample_id" =~ ^[A-Za-z0-9._-]+$ ]] || miss "sample_id contains unsupported characters: $sample_id"
    [[ -z "${extra:-}" ]] || miss "sample '$sample_id': more than 3 columns"
    [[ -n "${read1:-}" && -n "${read2:-}" ]] || { miss "sample '$sample_id': missing R1/R2 path"; continue; }

    if [[ -n "${seen[$sample_id]:-}" ]]; then
        miss "duplicate sample_id: $sample_id"
    fi
    seen[$sample_id]=1

    [[ -r "$read1" ]] || miss "sample '$sample_id': R1 not readable: $read1"
    [[ -r "$read2" ]] || miss "sample '$sample_id': R2 not readable: $read2"
done < "$SAMPLE_SHEET"

(( sample_count > 0 )) && ok "sample rows: $sample_count" || miss "sample sheet contains no data rows"

if (( failed )); then
    echo "[RESULT] Validation failed" >&2
    exit 1
fi

echo "[RESULT] Ino-seq environment and inputs are ready"

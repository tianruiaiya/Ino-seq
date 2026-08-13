#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'USAGE'
Usage:
  workflow/utils/prepare_reference.sh [config.env]

Default:
  config.env = config/inoseq.env

This command creates, when missing:
  - BWA index files (.amb/.ann/.bwt/.pac/.sa)
  - FASTA index (.fai)
  - sequence dictionary (.dict)
USAGE
}

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd "${SCRIPT_DIR}/../.." && pwd)
CONFIG_FILE=${1:-"${PROJECT_DIR}/config/inoseq.env"}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

[[ -f "$CONFIG_FILE" ]] || { echo "[ERROR] Configuration file not found: $CONFIG_FILE" >&2; exit 2; }
# shellcheck disable=SC1090
source "$CONFIG_FILE"

: "${REFERENCE_FASTA:?REFERENCE_FASTA is not set in $CONFIG_FILE}"
[[ -r "$REFERENCE_FASTA" ]] || { echo "[ERROR] Reference FASTA not readable: $REFERENCE_FASTA" >&2; exit 2; }
command -v bwa >/dev/null 2>&1 || { echo "[ERROR] bwa not found in PATH" >&2; exit 127; }
command -v samtools >/dev/null 2>&1 || { echo "[ERROR] samtools not found in PATH" >&2; exit 127; }

ref_dir=$(cd "$(dirname "$REFERENCE_FASTA")" && pwd)
[[ -w "$ref_dir" ]] || { echo "[ERROR] Reference directory is not writable: $ref_dir" >&2; exit 2; }

printf '[INFO] Reference: %s\n' "$REFERENCE_FASTA"

bwa_missing=0
for ext in amb ann bwt pac sa; do
    [[ -s "${REFERENCE_FASTA}.${ext}" ]] || bwa_missing=1
done
if (( bwa_missing )); then
    echo "[INFO] Building BWA index..."
    bwa index "$REFERENCE_FASTA"
else
    echo "[OK]   BWA index already present"
fi

if [[ ! -s "${REFERENCE_FASTA}.fai" ]]; then
    echo "[INFO] Building FASTA index (.fai)..."
    samtools faidx "$REFERENCE_FASTA"
else
    echo "[OK]   FASTA index already present"
fi

case "$REFERENCE_FASTA" in
    *.fasta) dict_file="${REFERENCE_FASTA%.fasta}.dict" ;;
    *.fa)    dict_file="${REFERENCE_FASTA%.fa}.dict" ;;
    *.fas)   dict_file="${REFERENCE_FASTA%.fas}.dict" ;;
    *)       dict_file="${REFERENCE_FASTA}.dict" ;;
esac

if [[ ! -s "$dict_file" ]]; then
    echo "[INFO] Building sequence dictionary: $dict_file"
    samtools dict -o "$dict_file" "$REFERENCE_FASTA"
else
    echo "[OK]   Sequence dictionary already present: $dict_file"
fi

echo "[INFO] Reference preparation complete"

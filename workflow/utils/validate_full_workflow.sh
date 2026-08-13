#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'USAGE'
Usage:
  workflow/utils/validate_full_workflow.sh [samples.tsv] [pairs.tsv] [config.env]

Defaults:
  samples.tsv = config/samples.tsv
  pairs.tsv   = config/pairs.tsv
  config.env  = config/inoseq.env

This validates raw FASTQ inputs, the runtime environment, pair definitions,
and sample/control membership before the automatic Slurm dependency graph is
submitted. Phase A output files are intentionally not required at this stage.
USAGE
}

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd "${SCRIPT_DIR}/../.." && pwd)
SAMPLE_SHEET=${1:-"${PROJECT_DIR}/config/samples.tsv"}
PAIR_SHEET=${2:-"${PROJECT_DIR}/config/pairs.tsv"}
CONFIG_FILE=${3:-"${PROJECT_DIR}/config/inoseq.env"}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

"${PROJECT_DIR}/workflow/utils/validate_config.sh" "$SAMPLE_SHEET" "$CONFIG_FILE"
"${PROJECT_DIR}/workflow/utils/validate_postprocess.sh" \
    --skip-step01-inputs "$PAIR_SHEET" "$CONFIG_FILE"

failed=0
ok()   { printf '[OK]   %s\n' "$*"; }
miss() { printf '[FAIL] %s\n' "$*" >&2; failed=1; }

declare -A declared_samples=()
while IFS=$'\t' read -r sample_id read1 read2 extra; do
    [[ "$sample_id" == "sample_id" || -z "${sample_id// }" || "$sample_id" == \#* ]] && continue
    declared_samples["$sample_id"]=1
done < "$SAMPLE_SHEET"

pair_count=0
while IFS=$'\t' read -r sample_id control_id sgrna extra; do
    [[ "$sample_id" == "sample_id" || -z "${sample_id// }" || "$sample_id" == \#* ]] && continue
    ((pair_count+=1))
    [[ -n "${declared_samples[$sample_id]:-}" ]] \
        || miss "pair experiment is absent from sample sheet: $sample_id"
    [[ -n "${declared_samples[$control_id]:-}" ]] \
        || miss "pair control is absent from sample sheet: $control_id"
done < "$PAIR_SHEET"

if (( failed )); then
    echo "[RESULT] Full-workflow validation failed" >&2
    exit 1
fi

ok "all ${pair_count} pair(s) reference samples with raw FASTQ entries"
echo "[RESULT] Ino-seq full workflow is ready for automatic submission"

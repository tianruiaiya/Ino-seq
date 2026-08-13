#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'USAGE'
Usage:
  bash workflow/run_inoseq.sh [--from-module auto|01|02] \
    SAMPLE_ID R1.fastq.gz R2.fastq.gz [config.env]

Internal direct runner for Phase A (modules 01-02).

With "auto" (default), a current module 01 result is reused when module 02 is
the first incomplete/stale stage.  An already current module 02 result causes
the runner to exit without recomputation.

Default config:
  config/inoseq.env
USAGE
}

FROM_MODULE=auto
if [[ "${1:-}" == "--from-module" ]]; then
    [[ $# -ge 2 ]] || { usage >&2; exit 2; }
    FROM_MODULE=$2
    shift 2
fi

if [[ $# -lt 3 || $# -gt 4 ]]; then
    usage >&2
    exit 2
fi

sample=$1
read1=$2
read2=$3

if [[ -n "${INOSEQ_PROJECT_DIR:-}" ]]; then
    PROJECT_DIR=$(cd "${INOSEQ_PROJECT_DIR}" && pwd)
else
    SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
    PROJECT_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
fi

CONFIG_FILE=${4:-"${PROJECT_DIR}/config/inoseq.env"}
[[ -f "$CONFIG_FILE" ]] || { echo "[ERROR] Config not found: $CONFIG_FILE" >&2; exit 2; }
CONFIG_FILE=$(cd "$(dirname "$CONFIG_FILE")" && pwd)/$(basename "$CONFIG_FILE")

export INOSEQ_PROJECT_DIR="${PROJECT_DIR}"

# shellcheck disable=SC1090
source "$CONFIG_FILE"
: "${PYTHON_BIN:=python}"
: "${OUTPUT_DIR:=output}"
: "${MUTATION_LOCATION_QUAL_THRESHOLD:=30}"

if [[ "$OUTPUT_DIR" != /* ]]; then
    OUTPUT_DIR="${PROJECT_DIR}/${OUTPUT_DIR}"
fi

prefix="${OUTPUT_DIR}/${sample}/alignment/${sample}"
STATE_TOOL="${PROJECT_DIR}/workflow/utils/stage_state.py"
state_args=(
    --project-dir "$PROJECT_DIR" --config "$CONFIG_FILE"
    --sample "$sample" --read1 "$read1" --read2 "$read2"
)

# BEGIN INOSEQ MODULE01 AUTO-REPAIR
repair_module01_if_possible() {
    local state_file="${OUTPUT_DIR}/${sample}/.inoseq/module01.json"
    local bam="${prefix}.umi_dedup.bam"
    local bai="${bam}.bai"
    local index_threads="${INDEX_THREADS:-8}"

    [[ ! -e "$state_file" ]] || return 1
    [[ -s "$bam" ]] || return 1

    if ! samtools quickcheck -v "$bam" >/dev/null 2>&1; then
        printf '[%s] [%s] Module 01 recovery skipped: existing BAM is incomplete\n' \
            "$(date '+%F %T')" "$sample" >&2
        return 1
    fi

    if [[ ! -s "$bai" ]] || ! samtools idxstats "$bam" >/dev/null 2>&1; then
        printf '[%s] [%s] Module 01 recovery: rebuilding canonical BAM index\n' \
            "$(date '+%F %T')" "$sample" >&2
        rm -f "$bai"
        samtools index -@ "$index_threads" -o "$bai" "$bam" || return 1
    fi

    samtools idxstats "$bam" >/dev/null 2>&1 || return 1

    if "$PYTHON_BIN" "$STATE_TOOL" record --stage module01 "${state_args[@]}"; then
        printf '[%s] [%s] Module 01 existing outputs repaired and adopted\n' \
            "$(date '+%F %T')" "$sample" >&2
        return 0
    fi

    printf '[%s] [%s] Module 01 recovery could not satisfy the full state contract; recomputation required\n' \
        "$(date '+%F %T')" "$sample" >&2
    return 1
}
# END INOSEQ MODULE01 AUTO-REPAIR

case "$FROM_MODULE" in
    auto)
        if "$PYTHON_BIN" "$STATE_TOOL" check --stage module02 "${state_args[@]}" --quiet; then
            printf '[%s] [%s] Phase A already current; no computation required\n' \
                "$(date '+%F %T')" "$sample" >&2
            exit 0
        fi
if "$PYTHON_BIN" "$STATE_TOOL" check --stage module01 "${state_args[@]}" --quiet; then
    FROM_MODULE=02
elif repair_module01_if_possible \
     && "$PYTHON_BIN" "$STATE_TOOL" check --stage module01 "${state_args[@]}" --quiet; then
    FROM_MODULE=02
else
    FROM_MODULE=01
fi
        ;;
    01|02) ;;
    *)
        echo "[ERROR] --from-module must be auto, 01 or 02: $FROM_MODULE" >&2
        exit 2
        ;;
esac

if [[ "$FROM_MODULE" == "02" ]] \
    && ! "$PYTHON_BIN" "$STATE_TOOL" check --stage module01 "${state_args[@]}" --quiet; then
    echo "[ERROR] Cannot start at module 02: module 01 is incomplete or stale for $sample" >&2
    exit 2
fi

"$PYTHON_BIN" "$STATE_TOOL" invalidate \
    --stage "module${FROM_MODULE}" "${state_args[@]}"

if [[ "$FROM_MODULE" == "01" ]]; then
    printf '[%s] [%s] Phase A | module 01: UMI processing and consensus\n' "$(date '+%F %T')" "$sample" >&2
    bash "${PROJECT_DIR}/workflow/modules/01_umi_consensus.sh" \
        "$sample" "$read1" "$read2" "$CONFIG_FILE"
    "$PYTHON_BIN" "$STATE_TOOL" record --stage module01 "${state_args[@]}"
else
    printf '[%s] [%s] Phase A | module 01: reused current result\n' "$(date '+%F %T')" "$sample" >&2
fi

printf '[%s] [%s] Phase A | module 02: ABE signature-read identification\n' "$(date '+%F %T')" "$sample" >&2
"$PYTHON_BIN" "${PROJECT_DIR}/workflow/modules/02_signature_reads.py" \
    "$prefix" \
    "${prefix}.umi_dedup.bam" \
    "$REFERENCE_FASTA" \
    "$MUTATION_LOCATION_QUAL_THRESHOLD"

"$PYTHON_BIN" "$STATE_TOOL" record --stage module02 "${state_args[@]}"
printf '[%s] [%s] Phase A complete (modules 01-02)\n' "$(date '+%F %T')" "$sample" >&2

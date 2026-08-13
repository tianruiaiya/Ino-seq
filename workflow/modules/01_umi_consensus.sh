#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'USAGE'
Usage:
  workflow/modules/01_umi_consensus.sh SAMPLE_ID R1.fastq.gz R2.fastq.gz [config.env]
USAGE
}

if [[ $# -lt 3 || $# -gt 4 ]]; then
    usage >&2
    exit 2
fi

sample=$1
read1=$2
read2=$3
# Resolve the repository root explicitly for Slurm jobs.  sbatch may execute a
# temporary copy of this script under /var/spool/slurmd, so BASH_SOURCE alone
# is not a reliable way to locate repository files on compute nodes.
if [[ -n "${INOSEQ_PROJECT_DIR:-}" ]]; then
    PROJECT_DIR=$(cd "${INOSEQ_PROJECT_DIR}" && pwd)
    SCRIPT_DIR="${PROJECT_DIR}/workflow/modules"
else
    SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
    PROJECT_DIR=$(cd "${SCRIPT_DIR}/../.." && pwd)
fi
config_file=${4:-"${PROJECT_DIR}/config/inoseq.env"}

[[ -f "$config_file" ]] || { echo "[ERROR] Config not found: $config_file" >&2; exit 2; }
# shellcheck disable=SC1090
source "$config_file"

# Runtime defaults. Analysis parameters correspond to the Step 01 specification.
: "${REFERENCE_FASTA:?REFERENCE_FASTA is not set in $config_file}"
: "${FGBIO_JAR:=}"
: "${PYTHON_BIN:=python}"
: "${OUTPUT_DIR:=output}"
if [[ "$OUTPUT_DIR" != /* ]]; then
    OUTPUT_DIR="${PROJECT_DIR}/${OUTPUT_DIR}"
fi
: "${FASTP_THREADS:=40}"
: "${CUTADAPT_THREADS:=10}"
: "${BWA_THREADS:=40}"
: "${SAMTOOLS_SORT_THREADS:=4}"
: "${SAMTOOLS_SORT_MEM:=5G}"
: "${CONSENSUS_THREADS:=40}"
: "${FINAL_SORT_THREADS:=40}"
: "${STATS_THREADS:=40}"
: "${INDEX_THREADS:=8}"
: "${GROUP_JAVA_XMX:=180g}"
: "${CONSENSUS_JAVA_XMX:=180g}"
: "${ZIPPER_JAVA_XMX:=180g}"
: "${FASTP_QUALIFIED_QUALITY_PHRED:=20}"
: "${FASTP_UNQUALIFIED_PERCENT_LIMIT:=10}"
: "${FASTP_LENGTH_REQUIRED:=50}"
: "${READS_TO_PROCESS:=0}"
: "${QC_UNIQUE_MAPQ:=30}"
: "${UMI_R2_TRIM:=12}"
: "${UMI_ADAPTER:=TGTAGAGCACGCGTGG}"
: "${UMI_STRATEGY:=Adjacency}"
: "${UMI_EDITS:=1}"
: "${CONSENSUS_MIN_READS:=1}"
: "${CONSENSUS_MIN_INPUT_BASE_QUALITY:=20}"
: "${FILTER_MIN_READS:=1}"
: "${FILTER_MIN_BASE_QUALITY:=20}"
: "${FILTER_MAX_BASE_ERROR_RATE:=0.2}"

for exe in fastp cutadapt bwa samtools "$PYTHON_BIN"; do
    command -v "$exe" >/dev/null 2>&1 || { echo "[ERROR] Required command not found: $exe" >&2; exit 127; }
done
[[ -r "$read1" ]] || { echo "[ERROR] R1 not readable: $read1" >&2; exit 2; }
[[ -r "$read2" ]] || { echo "[ERROR] R2 not readable: $read2" >&2; exit 2; }
[[ -r "$REFERENCE_FASTA" ]] || { echo "[ERROR] Reference FASTA not readable: $REFERENCE_FASTA" >&2; exit 2; }
[[ "$READS_TO_PROCESS" =~ ^[0-9]+$ ]] || { echo "[ERROR] READS_TO_PROCESS must be a non-negative integer: $READS_TO_PROCESS" >&2; exit 2; }
[[ "$QC_UNIQUE_MAPQ" =~ ^[0-9]+$ ]] || { echo "[ERROR] QC_UNIQUE_MAPQ must be a non-negative integer: $QC_UNIQUE_MAPQ" >&2; exit 2; }

fgbio_call() {
    local heap=$1
    shift
    if [[ -n "${FGBIO_JAR:-}" ]]; then
        command -v java >/dev/null 2>&1 || { echo "[ERROR] java not found" >&2; return 127; }
        [[ -r "$FGBIO_JAR" ]] || { echo "[ERROR] fgbio JAR not readable: $FGBIO_JAR" >&2; return 2; }
        java -Xmx"${heap}" -jar "$FGBIO_JAR" "$@"
    else
        command -v fgbio >/dev/null 2>&1 || { echo "[ERROR] fgbio not found in PATH" >&2; return 127; }
        env -u JAVA_TOOL_OPTIONS -u _JAVA_OPTIONS -u JDK_JAVA_OPTIONS fgbio -Xms512m -Xmx"${heap}" "$@"
    fi
}

sample_dir="${OUTPUT_DIR}/${sample}"
qc_dir="${sample_dir}/QC"
aln_dir="${sample_dir}/alignment"
prefix="${aln_dir}/${sample}"
mkdir -p "$qc_dir" "$aln_dir"

log_step() { printf '[%s] [%s] %s\n' "$(date '+%F %T')" "$sample" "$*" >&2; }

log_step "Step 1/7: read QC with fastp"
fastp_args=(
    --qualified_quality_phred "$FASTP_QUALIFIED_QUALITY_PHRED"
    --unqualified_percent_limit "$FASTP_UNQUALIFIED_PERCENT_LIMIT"
    --length_required "$FASTP_LENGTH_REQUIRED"
    -w "$FASTP_THREADS"
    --in1 "$read1" --in2 "$read2"
    --out1 "${qc_dir}/${sample}_1.fq.gz"
    --out2 "${qc_dir}/${sample}_2.fq.gz"
    --json "${qc_dir}/${sample}_fastp.json"
    --html "${qc_dir}/${sample}_fastp.html"
)
if (( READS_TO_PROCESS > 0 )); then
    log_step "fastp input cap enabled: ${READS_TO_PROCESS} read pairs"
    fastp_args+=(--reads_to_process "$READS_TO_PROCESS")
else
    log_step "fastp input cap disabled: processing all available read pairs"
fi
fastp "${fastp_args[@]}"

log_step "Step 2/7: UMI extraction with cutadapt"
cutadapt --quiet --discard-untrimmed \
    -j "$CUTADAPT_THREADS" \
    -U "$UMI_R2_TRIM" \
    -G "$UMI_ADAPTER" \
    --rename '{header} RX:Z:{r2.cut_prefix}' \
    -o "${qc_dir}/R1.fq.gz" \
    -p "${qc_dir}/R2.fq.gz" \
    --json "${qc_dir}/${sample}_cutadapt.json" \
    "${qc_dir}/${sample}_1.fq.gz" "${qc_dir}/${sample}_2.fq.gz"

log_step "Step 3/7: initial alignment"
bwa mem -t "$BWA_THREADS" -C -M -Y "$REFERENCE_FASTA" "${qc_dir}/R1.fq.gz" "${qc_dir}/R2.fq.gz" \
    | samtools sort -@ "$SAMTOOLS_SORT_THREADS" -m "$SAMTOOLS_SORT_MEM" -o "${prefix}.mapped.bam"
samtools index "${prefix}.mapped.bam"
# Initial-alignment QC metrics are auxiliary summaries and do not alter the BAM.
samtools flagstat --output-fmt json -@ "$STATS_THREADS" "${prefix}.mapped.bam" > "${prefix}.mapped.flagstat.json"
samtools view -@ "$STATS_THREADS" -c -q "$QC_UNIQUE_MAPQ" -F 2308 "${prefix}.mapped.bam" > "${prefix}.mapped.primary_mapq${QC_UNIQUE_MAPQ}.count.txt"

log_step "Step 4/7: UMI-family grouping"
fgbio_call "$GROUP_JAVA_XMX" --compression 1 --async-io GroupReadsByUmi \
    --input "${prefix}.mapped.bam" \
    --strategy "$UMI_STRATEGY" \
    --edits "$UMI_EDITS" \
    --output "${prefix}.grouped.bam" \
    --family-size-histogram "${prefix}.tag-family-sizes.txt"

log_step "Step 5/7: molecular-consensus construction and filtering"
fgbio_call "$CONSENSUS_JAVA_XMX" --compression 0 CallMolecularConsensusReads \
    --input "${prefix}.grouped.bam" \
    --output /dev/stdout \
    --read-name-prefix "$sample" \
    --min-reads "$CONSENSUS_MIN_READS" \
    --min-input-base-quality "$CONSENSUS_MIN_INPUT_BASE_QUALITY" \
    --threads "$CONSENSUS_THREADS" \
    | fgbio_call "$CONSENSUS_JAVA_XMX" --compression 1 FilterConsensusReads \
        --input /dev/stdin \
        --output "${prefix}.cons.unmapped.bam" \
        --ref "$REFERENCE_FASTA" \
        --min-reads "$FILTER_MIN_READS" \
        --min-base-quality "$FILTER_MIN_BASE_QUALITY" \
        --max-base-error-rate "$FILTER_MAX_BASE_ERROR_RATE"

log_step "Step 6/7: consensus remapping"
samtools fastq "${prefix}.cons.unmapped.bam" \
    | bwa mem -t "$BWA_THREADS" -p -K 150000000 -Y "$REFERENCE_FASTA" - \
    | fgbio_call "$ZIPPER_JAVA_XMX" --compression 0 --async-io ZipperBams \
        --unmapped "${prefix}.cons.unmapped.bam" \
        --ref "$REFERENCE_FASTA" \
        --tags-to-reverse Consensus \
        --tags-to-revcomp Consensus \
    | samtools sort --threads "$FINAL_SORT_THREADS" -o "${prefix}.umi_dedup.bam"

# BEGIN INOSEQ MODULE01 FINAL-BAM CONTRACT REPAIR
final_bam="${prefix}.umi_dedup.bam"
final_bai="${final_bam}.bai"

log_step "Validating final UMI-deduplicated BAM"
[[ -s "$final_bam" ]] || {
    echo "[ERROR] Missing or empty final BAM: $final_bam" >&2
    exit 2
}

samtools quickcheck -v "$final_bam" || {
    echo "[ERROR] Final BAM failed samtools quickcheck: $final_bam" >&2
    exit 2
}

if [[ ! -s "$final_bai" ]] || ! samtools idxstats "$final_bam" >/dev/null 2>&1; then
    log_step "Building canonical BAM index: ${final_bai}"
    rm -f "$final_bai"
    samtools index -@ "$INDEX_THREADS" -o "$final_bai" "$final_bam"
fi

[[ -s "$final_bai" ]] || {
    echo "[ERROR] BAM index missing after indexing: $final_bai" >&2
    exit 2
}

samtools idxstats "$final_bam" >/dev/null || {
    echo "[ERROR] BAM/index validation failed: $final_bam / $final_bai" >&2
    exit 2
}
# END INOSEQ MODULE01 FINAL-BAM CONTRACT REPAIR

log_step "Step 7/7: alignment statistics"
samtools flagstat --output-fmt json -@ "$STATS_THREADS" "${prefix}.umi_dedup.bam" > "${prefix}.flagstat.json"
samtools idxstat "${prefix}.umi_dedup.bam" > "${prefix}.idxstat.txt"

log_step "UMI consensus stage complete"

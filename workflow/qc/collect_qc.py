#!/usr/bin/env python3
"""Collect Ino-seq QC metrics into one tabular summary.

The QC table summarizes fastp sequencing metrics, cutadapt UMI-extraction
retention, initial alignment, fgbio UMI-family distributions, final consensus
alignment, and ABE-signature read yield.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
from collections.abc import Iterable
from pathlib import Path

NA = "NA"

CORE_COLUMNS = [
    "name",
    "raw_bases",
    "raw_reads",
    "raw_reads_Q20_percent",
    "raw_reads_Q30_percent",
    "raw_reads_gc_content",
    "clean_bases",
    "clean_reads",
    "clean_reads_Q20_percent",
    "clean_reads_Q30_percent",
    "clean_reads_gc_content",
    "reads_with_barcode",
    "mapped_reads",
    "mapped_percent",
]

EXTENDED_COLUMNS = [
    "raw_read_pairs",
    "clean_read_pairs",
    "qc_pass_percent",
    "umi_extraction_retention_percent",
    "initial_mapped_reads",
    "initial_mapped_percent",
    "initial_primary_mapped_reads",
    "initial_primary_mapped_percent",
    "high_confidence_primary_mapped_reads",
    "high_confidence_primary_mapped_percent",
    "observed_umi_molecules",
    "umi_family_templates",
    "mean_umi_family_size",
    "singleton_umi_fraction_percent",
    "umi_family_compression_rate_percent",
    "umi_consensus_reads",
    "umi_consensus_read_pairs",
    "primary_mapped_reads",
    "primary_mapped_percent",
    "properly_paired_reads",
    "properly_paired_percent",
    "abe_signature_reads",
    "abe_signature_reads_per_million_consensus_pairs",
    "reads_to_process",
    "input_cap_enabled",
    "qc_status",
]

ALL_COLUMNS = CORE_COLUMNS + EXTENDED_COLUMNS


def load_json(path: Path) -> dict | None:
    if not path.is_file():
        return None
    with path.open() as handle:
        return json.load(handle)


def fmt_pct(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return NA
    return f"{value:.2f}"


def fmt_float(value: float | None, digits: int = 4) -> str:
    if value is None or not math.isfinite(value):
        return NA
    return f"{value:.{digits}f}"


def fmt_bases(value: int | None) -> str:
    if value is None:
        return NA
    return f"{value / 1e9:.2f}G"


def get_nested(d: dict | None, *keys):
    cur = d
    for key in keys:
        if not isinstance(cur, dict) or key not in cur:
            return None
        cur = cur[key]
    return cur


def parse_env(path: Path | None) -> dict[str, str]:
    values: dict[str, str] = {}
    if path is None or not path.is_file():
        return values
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip().strip('"').strip("'")
    return values


def parse_family_histogram(path: Path) -> dict | None:
    if not path.is_file():
        return None
    total_families = 0
    total_templates = 0
    singleton_count = 0
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"family_size", "count"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            return None
        for row in reader:
            size = int(row["family_size"])
            count = int(row["count"])
            total_families += count
            total_templates += size * count
            if size == 1:
                singleton_count += count
    if total_families == 0 or total_templates == 0:
        return None
    return {
        "observed_umi_molecules": total_families,
        "umi_family_templates": total_templates,
        "mean_umi_family_size": total_templates / total_families,
        "singleton_umi_fraction_percent": 100.0 * singleton_count / total_families,
        "umi_family_compression_rate_percent": 100.0 * (1.0 - total_families / total_templates),
    }


def parse_int_file(path: Path) -> int | None:
    if not path.is_file():
        return None
    text = path.read_text().strip().split()
    if not text:
        return None
    try:
        return int(text[0])
    except ValueError:
        return None


def count_data_rows(path: Path) -> int | None:
    if not path.is_file():
        return None
    n = 0
    with path.open() as handle:
        for i, line in enumerate(handle):
            if i == 0:
                continue
            if line.strip():
                n += 1
    return n


def sample_ids_from_sheet(path: Path) -> list[str]:
    samples: list[str] = []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["sample_id", "read1", "read2"]:
            raise ValueError("sample sheet header must be: sample_id<TAB>read1<TAB>read2")
        for row in reader:
            sample = (row.get("sample_id") or "").strip()
            if sample and not sample.startswith("#"):
                samples.append(sample)
    return samples



def ensure_initial_mapping_qc(sample: str, output_dir: Path, env: dict[str, str]) -> None:
    """Generate auxiliary initial-alignment QC summaries when BAM exists.

    This reproduces the behavior previously implemented by the shell QC wrapper.
    It does not modify analytical BAM files.
    """
    aln_dir = output_dir / sample / "alignment"
    mapped_bam = aln_dir / f"{sample}.mapped.bam"
    mapped_flagstat = aln_dir / f"{sample}.mapped.flagstat.json"
    mapq = env.get("QC_UNIQUE_MAPQ", "30") or "30"
    highconf = aln_dir / f"{sample}.mapped.primary_mapq{mapq}.count.txt"
    threads = env.get("STATS_THREADS", "40") or "40"

    if not mapped_bam.is_file():
        return

    samtools = shutil.which("samtools")
    if samtools is None:
        return

    if not mapped_flagstat.is_file():
        with mapped_flagstat.open("w") as handle:
            subprocess.run(
                [samtools, "flagstat", "--output-fmt", "json", "-@", str(threads), str(mapped_bam)],
                check=True,
                stdout=handle,
            )

    if not highconf.is_file():
        result = subprocess.run(
            [
                samtools, "view", "-@", str(threads), "-c",
                "-q", str(mapq), "-F", "2308", str(mapped_bam)
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        highconf.write_text(result.stdout.strip() + "\n")

def collect_sample(sample: str, output_dir: Path, env: dict[str, str]) -> tuple[dict[str, object], list[str]]:
    ensure_initial_mapping_qc(sample, output_dir, env)
    sample_dir = output_dir / sample
    qc_dir = sample_dir / "QC"
    aln_dir = sample_dir / "alignment"

    fastp_path = qc_dir / f"{sample}_fastp.json"
    cutadapt_path = qc_dir / f"{sample}_cutadapt.json"
    family_path = aln_dir / f"{sample}.tag-family-sizes.txt"
    initial_flagstat_path = aln_dir / f"{sample}.mapped.flagstat.json"
    highconf_path = aln_dir / f"{sample}.mapped.primary_mapq{env.get('QC_UNIQUE_MAPQ', '30')}.count.txt"
    final_flagstat_path = aln_dir / f"{sample}.flagstat.json"
    end_path = aln_dir / f"{sample}.end"

    missing: list[str] = []
    for p in [fastp_path, cutadapt_path, family_path, final_flagstat_path]:
        if not p.is_file():
            missing.append(str(p))

    fastp = load_json(fastp_path)
    cutadapt = load_json(cutadapt_path)
    initial_flagstat = load_json(initial_flagstat_path)
    final_flagstat = load_json(final_flagstat_path)
    family = parse_family_histogram(family_path)

    before = get_nested(fastp, "summary", "before_filtering") or {}
    after = get_nested(fastp, "summary", "after_filtering") or {}

    raw_bases = before.get("total_bases")
    raw_reads = before.get("total_reads")
    clean_bases = after.get("total_bases")
    clean_reads = after.get("total_reads")

    cutadapt_output_pairs = get_nested(cutadapt, "read_counts", "output")
    reads_with_barcode = 2 * cutadapt_output_pairs if isinstance(cutadapt_output_pairs, int) else None

    final_qc = get_nested(final_flagstat, "QC-passed reads") or {}
    initial_qc = get_nested(initial_flagstat, "QC-passed reads") or {}

    mapped_reads = final_qc.get("mapped")
    mapped_percent = final_qc.get("mapped %")
    primary_reads = final_qc.get("primary")
    consensus_pairs = final_qc.get("read1")

    highconf = parse_int_file(highconf_path)
    initial_primary = initial_qc.get("primary")

    abe_reads = count_data_rows(end_path)

    reads_to_process_text = env.get("READS_TO_PROCESS", "0") or "0"
    try:
        reads_to_process = int(reads_to_process_text)
    except ValueError:
        reads_to_process = 0

    row: dict[str, object] = {
        "name": sample,
        "raw_bases": fmt_bases(raw_bases),
        "raw_reads": raw_reads if raw_reads is not None else NA,
        "raw_reads_Q20_percent": fmt_pct(100 * before["q20_rate"]) if "q20_rate" in before else NA,
        "raw_reads_Q30_percent": fmt_pct(100 * before["q30_rate"]) if "q30_rate" in before else NA,
        "raw_reads_gc_content": fmt_pct(100 * before["gc_content"]) if "gc_content" in before else NA,
        "clean_bases": fmt_bases(clean_bases),
        "clean_reads": clean_reads if clean_reads is not None else NA,
        "clean_reads_Q20_percent": fmt_pct(100 * after["q20_rate"]) if "q20_rate" in after else NA,
        "clean_reads_Q30_percent": fmt_pct(100 * after["q30_rate"]) if "q30_rate" in after else NA,
        "clean_reads_gc_content": fmt_pct(100 * after["gc_content"]) if "gc_content" in after else NA,
        "reads_with_barcode": reads_with_barcode if reads_with_barcode is not None else NA,
        "mapped_reads": mapped_reads if mapped_reads is not None else NA,
        "mapped_percent": mapped_percent if mapped_percent is not None else NA,
        "raw_read_pairs": raw_reads // 2 if isinstance(raw_reads, int) else NA,
        "clean_read_pairs": clean_reads // 2 if isinstance(clean_reads, int) else NA,
        "qc_pass_percent": fmt_pct(100 * clean_reads / raw_reads) if raw_reads and clean_reads is not None else NA,
        "umi_extraction_retention_percent": fmt_pct(100 * reads_with_barcode / clean_reads) if clean_reads and reads_with_barcode is not None else NA,
        "initial_mapped_reads": initial_qc.get("mapped", NA),
        "initial_mapped_percent": initial_qc.get("mapped %", NA),
        "initial_primary_mapped_reads": initial_qc.get("primary mapped", NA),
        "initial_primary_mapped_percent": initial_qc.get("primary mapped %", NA),
        "high_confidence_primary_mapped_reads": highconf if highconf is not None else NA,
        "high_confidence_primary_mapped_percent": fmt_pct(100 * highconf / initial_primary) if highconf is not None and initial_primary else NA,
        "observed_umi_molecules": family["observed_umi_molecules"] if family else NA,
        "umi_family_templates": family["umi_family_templates"] if family else NA,
        "mean_umi_family_size": fmt_float(family["mean_umi_family_size"]) if family else NA,
        "singleton_umi_fraction_percent": fmt_pct(family["singleton_umi_fraction_percent"]) if family else NA,
        "umi_family_compression_rate_percent": fmt_pct(family["umi_family_compression_rate_percent"]) if family else NA,
        "umi_consensus_reads": primary_reads if primary_reads is not None else NA,
        "umi_consensus_read_pairs": consensus_pairs if consensus_pairs is not None else NA,
        "primary_mapped_reads": final_qc.get("primary mapped", NA),
        "primary_mapped_percent": final_qc.get("primary mapped %", NA),
        "properly_paired_reads": final_qc.get("properly paired", NA),
        "properly_paired_percent": final_qc.get("properly paired %", NA),
        "abe_signature_reads": abe_reads if abe_reads is not None else NA,
        "abe_signature_reads_per_million_consensus_pairs": fmt_float(1e6 * abe_reads / consensus_pairs, 2) if abe_reads is not None and consensus_pairs else NA,
        "reads_to_process": reads_to_process,
        "input_cap_enabled": "yes" if reads_to_process > 0 else "no",
        "qc_status": "PASS" if not missing else "PARTIAL",
    }
    return row, missing


def write_table(rows: Iterable[dict[str, object]], path: Path, delimiter: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ALL_COLUMNS, delimiter=delimiter, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Collect Ino-seq QC metrics")
    parser.add_argument("--samples", required=True, type=Path, help="sample sheet TSV")
    parser.add_argument("--output-dir", required=True, type=Path, help="Ino-seq output directory")
    parser.add_argument("--config", type=Path, help="Ino-seq env config")
    parser.add_argument("--tsv", required=True, type=Path, help="output TSV")
    parser.add_argument("--csv", type=Path, help="optional output CSV")
    parser.add_argument("--strict", action="store_true", help="exit non-zero if required QC files are missing")
    args = parser.parse_args()

    env = parse_env(args.config)
    samples = sample_ids_from_sheet(args.samples)
    rows: list[dict[str, object]] = []
    missing_any = False
    for sample in samples:
        row, missing = collect_sample(sample, args.output_dir, env)
        rows.append(row)
        if missing:
            missing_any = True
            print(f"[WARN] {sample}: missing {len(missing)} required QC file(s)")
            for item in missing:
                print(f"       {item}")

    write_table(rows, args.tsv, "\t")
    if args.csv:
        write_table(rows, args.csv, ",")

    print(f"[OK] QC TSV: {args.tsv}")
    if args.csv:
        print(f"[OK] QC CSV: {args.csv}")
    return 1 if args.strict and missing_any else 0


if __name__ == "__main__":
    raise SystemExit(main())

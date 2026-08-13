#!/usr/bin/env python3
"""Classify reads by target/non-target strand at dependent candidates."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pysam

from workflow.lib.io_utils import ensure_parent, require_columns

REQUIRED_COLUMNS = [
    "chromosome",
    "start",
    "end",
    "max_position",
    "Merge Intervals",
    "dependent_type",
    "Site_SubstitutionsOnly.Strand",
    "Site_GapsAllowed.Strand",
]


def guide_strand(row: pd.Series) -> str | None:
    standard = row["Site_SubstitutionsOnly.Strand"]
    gapped = row["Site_GapsAllowed.Strand"]
    if not pd.isna(standard) and str(standard).strip():
        return str(standard).strip()
    if not pd.isna(gapped) and str(gapped).strip():
        return str(gapped).strip()
    return None


def reads_in_region(
    bam: pysam.AlignmentFile, chromosome: str, start: int, end: int
) -> dict[str, str]:
    # A dictionary intentionally retains one strand per query name, matching
    # the historical implementation when duplicate query names occur.
    reads: dict[str, str] = {}
    for read in bam.fetch(chromosome, start, end):
        reads[read.query_name] = "-" if read.is_reverse else "+"
    return reads


def classify_reads(reads: dict[str, str], sgrna_strand: str) -> tuple[list[str], list[str]]:
    non_target: list[str] = []
    target: list[str] = []
    for read_name, bam_strand in reads.items():
        cleavage_strand = "+" if bam_strand == "-" else "-"
        label = f"{read_name}({cleavage_strand})"
        if sgrna_strand == cleavage_strand:
            non_target.append(label)
        else:
            target.append(label)
    return non_target, target


def classify_strands(
    frame: pd.DataFrame, bam_file: str | Path, window: int = 25
) -> pd.DataFrame:
    columns: dict[str, list[str]] = {
        "NonTarget_reads": [""] * len(frame),
        "NonTarget_reads_ratio": [""] * len(frame),
        "Target_reads": [""] * len(frame),
        "Target_reads_ratio": [""] * len(frame),
    }

    with pysam.AlignmentFile(str(bam_file), "rb") as bam:
        for index, row in frame.iterrows():
            if row["dependent_type"] != "dependent":
                continue
            strand = guide_strand(row)
            if strand is None:
                continue

            if row["Merge Intervals"] == "no":
                center = (int(row["start"]) + int(row["end"])) // 2
            elif row["Merge Intervals"] == "yes":
                center = int(row["max_position"])
            else:
                continue

            reads = reads_in_region(bam, str(row["chromosome"]), center - window, center + window)
            total = len(reads)
            if total == 0:
                columns["NonTarget_reads_ratio"][index] = "0% (0)"
                columns["Target_reads_ratio"][index] = "0% (0)"
                continue

            non_target, target = classify_reads(reads, strand)
            columns["NonTarget_reads"][index] = ",".join(non_target)
            columns["NonTarget_reads_ratio"][index] = (
                f"{len(non_target) / total * 100:.1f}% ({len(non_target)})"
            )
            columns["Target_reads"][index] = ",".join(target)
            columns["Target_reads_ratio"][index] = (
                f"{len(target) / total * 100:.1f}% ({len(target)})"
            )

    result = frame.copy()
    for offset, name in enumerate(
        ["NonTarget_reads", "NonTarget_reads_ratio", "Target_reads", "Target_reads_ratio"]
    ):
        result.insert(12 + offset, name, columns[name])
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--window", type=int, default=25)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.window < 0:
        raise ValueError("--window must be non-negative")
    frame = pd.read_csv(args.input, sep="\t")
    require_columns(frame, REQUIRED_COLUMNS, args.input)
    result = classify_strands(frame, args.bam, args.window)
    output = ensure_parent(args.output)
    result.to_csv(output, sep="\t", index=False)
    dependent = int((result["dependent_type"] == "dependent").sum())
    print(f"[05 strands] {dependent} dependent candidate(s) -> {output}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Append BAM overlap counts to the seven-column Ino-seq candidate BED.

Each overlapping alignment record is counted, matching the historical
``bedtools intersect -a <BED> -b <BAM> -c`` use for valid coordinate-sorted
BAM inputs.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pysam

from workflow.lib.io_utils import ensure_parent


def count_coverage(bed_file: str | Path, bam_file: str | Path, output_file: str | Path) -> int:
    output_path = ensure_parent(output_file)
    row_count = 0

    with pysam.AlignmentFile(str(bam_file), "rb") as bam, open(
        bed_file, newline=""
    ) as source, open(output_path, "w", newline="") as destination:
        reader = csv.reader(source, delimiter="\t")
        writer = csv.writer(destination, delimiter="\t", lineterminator="\n")

        for line_number, row in enumerate(reader, start=1):
            if not row:
                continue
            if len(row) != 7:
                raise ValueError(
                    f"{bed_file}:{line_number}: expected 7 columns, observed {len(row)}"
                )
            chromosome, start_text, end_text = row[:3]
            start = int(start_text)
            end = int(end_text)
            count = sum(1 for _ in bam.fetch(chromosome, start, end))
            writer.writerow([*row, count])
            row_count += 1

    return row_count


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bed", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = count_coverage(args.bed, args.bam, args.output)
    print(f"[03 coverage] {rows} interval(s) -> {args.output}")


if __name__ == "__main__":
    main()

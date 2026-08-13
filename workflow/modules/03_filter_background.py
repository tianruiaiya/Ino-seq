#!/usr/bin/env python3
"""Filter background-comparison sites and retain their BAM alignments.

The historical enrichment rule is preserved except for the user-authorized
unified p-value threshold, now configured as ``p_value < 0.05`` by default.
No FDR or minimum-read filter is applied.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pysam

from workflow.lib.io_utils import ensure_parent, require_columns, split_query_names

REQUIRED_COLUMNS = ("fold_change", "p_value", "query_name")


def filter_table(
    input_file: str | Path,
    min_fold_change: float = 1.5,
    p_value_threshold: float = 0.05,
) -> tuple[pd.DataFrame, list[str]]:
    frame = pd.read_csv(input_file, sep="\t")
    require_columns(frame, REQUIRED_COLUMNS, input_file)
    selected = frame[
        (frame["fold_change"] >= min_fold_change)
        & (frame["p_value"] < p_value_threshold)
    ].copy()
    return selected, split_query_names(selected["query_name"])


def filter_bam(input_bam: str | Path, output_bam: str | Path, query_names: list[str]) -> int:
    output_path = ensure_parent(output_bam)
    selected_names = set(query_names)
    count = 0
    with pysam.AlignmentFile(str(input_bam), "rb") as source, pysam.AlignmentFile(
        str(output_path), "wb", template=source
    ) as destination:
        for read in source.fetch(until_eof=True):
            if read.query_name in selected_names:
                destination.write(read)
                count += 1
    pysam.index(str(output_path))
    return count


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Background comparison TSV")
    parser.add_argument("--input-bam", required=True)
    parser.add_argument("--output-table", required=True)
    parser.add_argument("--output-query-names", required=True)
    parser.add_argument("--output-bam", required=True)
    parser.add_argument("--min-fold-change", type=float, default=1.5)
    parser.add_argument("--p-value-threshold", type=float, default=0.05)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.min_fold_change < 0:
        raise ValueError("--min-fold-change must be non-negative")
    if not 0 <= args.p_value_threshold <= 1:
        raise ValueError("--p-value-threshold must be between 0 and 1")

    selected, query_names = filter_table(
        args.input, args.min_fold_change, args.p_value_threshold
    )
    output_table = ensure_parent(args.output_table)
    output_names = ensure_parent(args.output_query_names)

    # Retain the historical index column in *_filted.txt for output compatibility.
    selected.to_csv(output_table, sep="\t", index=True)
    pd.Series(query_names, dtype="object").to_csv(
        output_names, sep="\t", index=False, header=False
    )
    bam_reads = filter_bam(args.input_bam, args.output_bam, query_names)
    print(
        f"[03 filter] {len(selected)} site(s), {len(query_names)} query name(s), "
        f"{bam_reads} BAM alignment(s) -> {output_table}"
    )


if __name__ == "__main__":
    main()

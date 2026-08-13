#!/usr/bin/env python3
"""Aggregate Step01 signature reads into cleavage-site windows.

This is a structured implementation of the historical ``get_result.py``.
The grouping keys, first-observed strand rule, and +/-15 bp interval definition
are intentionally retained.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from workflow.lib.io_utils import ensure_parent, require_columns

REQUIRED_COLUMNS = ("query_name", "ref_name", "location", "cleavage_direction")


def aggregate_sites(input_file: str | Path, window: int = 15) -> pd.DataFrame:
    frame = pd.read_csv(input_file, sep="\t")
    require_columns(frame, REQUIRED_COLUMNS, input_file)

    grouped = (
        frame.groupby(["ref_name", "location"], sort=True)
        .agg(
            query_name=("query_name", lambda values: ",".join(map(str, values))),
            **{
                "reads count": ("location", "count"),
                "strand": ("cleavage_direction", "first"),
            },
        )
        .reset_index()
    )
    grouped["start"] = (grouped["location"] - window).clip(lower=0)
    grouped["end"] = grouped["location"] + window
    return grouped[
        ["ref_name", "start", "end", "location", "reads count", "strand", "query_name"]
    ]


def write_outputs(frame: pd.DataFrame, output_tsv: str | Path, output_bed: str | Path) -> None:
    output_tsv = ensure_parent(output_tsv)
    output_bed = ensure_parent(output_bed)
    frame.to_csv(output_tsv, sep="\t", index=False)
    frame.to_csv(output_bed, sep="\t", index=False, header=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Step01 <sample>.end file")
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--output-bed", required=True)
    parser.add_argument("--window", type=int, default=15)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.window < 0:
        raise ValueError("--window must be non-negative")
    result = aggregate_sites(args.input, args.window)
    write_outputs(result, args.output_tsv, args.output_bed)
    print(f"[03 aggregate] {len(result)} cleavage site(s) -> {args.output_tsv}")


if __name__ == "__main__":
    main()

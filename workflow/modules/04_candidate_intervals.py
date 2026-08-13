#!/usr/bin/env python3
"""Merge nearby filtered cleavage sites into candidate intervals.

This module retains the historical rules: sort by chromosome/location, merge
successive positions no more than 30 bp apart, center short intervals to 30 bp,
retain intervals with total cleavage reads >=3, and use the first position in a
tie for maximum cleavage-read support.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import pysam

from workflow.lib.io_utils import ensure_parent, require_columns, split_query_names

REQUIRED_COLUMNS = ("chromosome", "location", "clevage_reads", "query_name")
OUTPUT_COLUMNS = [
    "chromosome",
    "start",
    "end",
    "total_reads",
    "max_reads",
    "max_position",
    "Merge Intervals",
    "query_name",
]


@dataclass
class CandidateInterval:
    chromosome: str
    start: int
    end: int
    total_reads: int
    max_reads: int
    max_position: int
    merged: bool
    query_name: str

    def as_row(self) -> list[object]:
        return [
            self.chromosome,
            self.start,
            self.end,
            self.total_reads,
            self.max_reads,
            self.max_position,
            "yes" if self.merged else "no",
            self.query_name,
        ]


def _finalize_group(rows: list[pd.Series]) -> CandidateInterval:
    read_counts = [int(row["clevage_reads"]) for row in rows]
    max_reads = max(read_counts)
    max_index = read_counts.index(max_reads)
    return CandidateInterval(
        chromosome=str(rows[0]["chromosome"]),
        start=int(rows[0]["location"]),
        end=int(rows[-1]["location"]),
        total_reads=sum(read_counts),
        max_reads=max_reads,
        max_position=int(rows[max_index]["location"]),
        merged=len(rows) > 1,
        query_name=",".join(str(row["query_name"]) for row in rows),
    )


def merge_sites(frame: pd.DataFrame, distance: int = 30) -> list[CandidateInterval]:
    if frame.empty:
        return []
    ordered = frame.sort_values(["chromosome", "location"], kind="mergesort").reset_index(
        drop=True
    )
    groups: list[list[pd.Series]] = []
    current = [ordered.iloc[0]]
    for row_index in range(1, len(ordered)):
        row = ordered.iloc[row_index]
        previous = current[-1]
        if (
            str(row["chromosome"]) == str(previous["chromosome"])
            and int(row["location"]) - int(previous["location"]) <= distance
        ):
            current.append(row)
        else:
            groups.append(current)
            current = [row]
    groups.append(current)
    return [_finalize_group(group) for group in groups]


def extend_interval(interval: CandidateInterval, min_length: int = 30) -> CandidateInterval:
    current_length = interval.end - interval.start + 1
    if current_length >= min_length:
        return interval

    center = (interval.start + interval.end) // 2
    new_start = max(1, center - min_length // 2)
    new_end = center + min_length // 2 - 1
    if min_length % 2 == 1:
        new_end += 1
    return CandidateInterval(
        chromosome=interval.chromosome,
        start=new_start,
        end=new_end,
        total_reads=interval.total_reads,
        max_reads=interval.max_reads,
        max_position=interval.max_position,
        merged=interval.merged,
        query_name=interval.query_name,
    )


def build_candidates(
    input_file: str | Path,
    distance: int = 30,
    min_length: int = 30,
    min_reads: int = 3,
) -> pd.DataFrame:
    frame = pd.read_csv(input_file, sep="\t")
    require_columns(frame, REQUIRED_COLUMNS, input_file)
    intervals = [extend_interval(item, min_length) for item in merge_sites(frame, distance)]
    rows = [item.as_row() for item in intervals if item.total_reads >= min_reads]
    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS)


def filter_bam(input_bam: str | Path, output_bam: str | Path, candidates: pd.DataFrame) -> int:
    query_names = split_query_names(candidates.get("query_name", pd.Series(dtype="object")))
    selected = set(query_names)
    output_path = ensure_parent(output_bam)
    written = 0
    with pysam.AlignmentFile(str(input_bam), "rb") as source, pysam.AlignmentFile(
        str(output_path), "wb", template=source
    ) as destination:
        for read in source.fetch(until_eof=True):
            if read.query_name in selected:
                destination.write(read)
                written += 1
    pysam.index(str(output_path))
    return written


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--input-bam", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--output-bam", required=True)
    parser.add_argument("--distance", type=int, default=30)
    parser.add_argument("--min-length", type=int, default=30)
    parser.add_argument("--min-reads", type=int, default=3)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if min(args.distance, args.min_length, args.min_reads) < 0:
        raise ValueError("distance, min-length, and min-reads must be non-negative")
    candidates = build_candidates(
        args.input, args.distance, args.min_length, args.min_reads
    )
    output = ensure_parent(args.output)
    candidates.to_csv(output, sep="\t", index=False)
    alignments = filter_bam(args.input_bam, args.output_bam, candidates)
    print(
        f"[04 candidates] {len(candidates)} interval(s), {alignments} BAM alignment(s) "
        f"-> {output}"
    )


if __name__ == "__main__":
    main()

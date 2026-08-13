#!/usr/bin/env python3
"""Annotate neighboring intervals and cleavage position versus protospacer."""

from __future__ import annotations

import argparse

import pandas as pd

from workflow.lib.io_utils import ensure_parent, require_columns

REQUIRED_COLUMNS = [
    "chromosome",
    "start",
    "end",
    "max_position",
    "dependent_type",
    "Site_SubstitutionsOnly.Start",
    "Site_SubstitutionsOnly.End",
    "Site_GapsAllowed.Start",
    "Site_GapsAllowed.End",
]


def _valid_interval(start: object, end: object) -> bool:
    return (
        not pd.isna(start)
        and not pd.isna(end)
        and str(start).strip() != ""
        and str(end).strip() != ""
    )


def annotate_spacer(frame: pd.DataFrame, neighborhood: int = 100) -> pd.DataFrame:
    overlap_intervals: list[str] = []
    overlap_stats: list[str] = []
    extend_starts: list[object] = []
    extend_ends: list[object] = []
    position_status: list[str] = []
    boundary_distance: list[object] = []

    for index, row in frame.iterrows():
        if row["dependent_type"] != "dependent":
            overlap_intervals.append("")
            overlap_stats.append("")
            extend_starts.append("")
            extend_ends.append("")
            position_status.append("")
            boundary_distance.append("")
            continue

        extended_start = row["start"] - neighborhood
        extended_end = row["end"] + neighborhood
        extend_starts.append(extended_start)
        extend_ends.append(extended_end)

        neighbors: list[str] = []
        dependent_count = 0
        independent_count = 0
        for other_index, other in frame.iterrows():
            if other_index == index or other["chromosome"] != row["chromosome"]:
                continue
            if other["start"] >= extended_start and other["end"] <= extended_end:
                other_type = other["dependent_type"]
                neighbors.append(
                    f"{other['chromosome']}:{other['start']}-{other['end']}({other_type})"
                )
                if other_type == "dependent":
                    dependent_count += 1
                elif other_type == "independent":
                    independent_count += 1

        overlap_intervals.append(";".join(neighbors) if neighbors else "none")
        overlap_stats.append(
            f"dependent({dependent_count});independent({independent_count})"
        )

        if _valid_interval(
            row["Site_SubstitutionsOnly.Start"], row["Site_SubstitutionsOnly.End"]
        ):
            interval_start = float(row["Site_SubstitutionsOnly.Start"])
            interval_end = float(row["Site_SubstitutionsOnly.End"])
        elif _valid_interval(row["Site_GapsAllowed.Start"], row["Site_GapsAllowed.End"]):
            interval_start = float(row["Site_GapsAllowed.Start"])
            interval_end = float(row["Site_GapsAllowed.End"])
        else:
            position_status.append("")
            boundary_distance.append("")
            continue

        max_position = float(row["max_position"])
        if interval_start <= max_position <= interval_end:
            position_status.append("no")
            boundary_distance.append(0)
        else:
            position_status.append("out of protospacer")
            boundary_distance.append(
                min(abs(max_position - interval_start), abs(max_position - interval_end))
            )

    result = frame.copy()
    result.insert(8, "overlap_Intervals", overlap_intervals)
    result.insert(9, "overlap_Intervals_stat", overlap_stats)
    result.insert(10, "extend_start", extend_starts)
    result.insert(11, "extend_end", extend_ends)
    result["max_position_in_protospacer"] = position_status
    result["distance_to_protospacer_boundary"] = boundary_distance
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--neighborhood", type=int, default=100)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.neighborhood < 0:
        raise ValueError("--neighborhood must be non-negative")
    frame = pd.read_csv(args.input, sep="\t")
    require_columns(frame, REQUIRED_COLUMNS, args.input)
    result = annotate_spacer(frame, args.neighborhood)
    output = ensure_parent(args.output)
    result.to_csv(output, sep="\t", index=False)
    print(f"[05 spacer] {len(result)} candidate(s) -> {output}")


if __name__ == "__main__":
    main()

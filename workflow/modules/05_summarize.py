#!/usr/bin/env python3
"""Summarize final Ino-seq target/off-target classifications."""

from __future__ import annotations

import argparse
import re
from collections import defaultdict

import pandas as pd

from workflow.lib.io_utils import ensure_parent, require_columns

REQUIRED_COLUMNS = [
    "dependent_type",
    "total_reads",
    "NonTarget_reads",
    "Target_reads",
    "overlap_Intervals_stat",
    "max_position_in_protospacer",
]


def parse_overlap_stats(value: object) -> dict[str, int]:
    if pd.isna(value) or not str(value).strip() or str(value) == "none":
        return {}
    counts: defaultdict[str, int] = defaultdict(int)
    for part in str(value).split(";"):
        match = re.fullmatch(r"(\w+)\((\d+)\)", part.strip())
        if match:
            counts[match.group(1)] += int(match.group(2))
    return dict(counts)


def _not_empty(series: pd.Series) -> pd.Series:
    return series.notna() & series.astype(str).ne("")


def summarize(frame: pd.DataFrame, sample_id: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    require_columns(frame, REQUIRED_COLUMNS, "final classification table")
    on_target = frame[frame["dependent_type"] == "onTarget"]
    dependent = frame[frame["dependent_type"] == "dependent"]
    independent = frame[frame["dependent_type"] == "independent"]

    basic = pd.DataFrame(
        [
            {
                "样本名": sample_id,
                "OnTarget_number": len(on_target),
                "OnTarget_reads_number": on_target["total_reads"].sum() if len(on_target) else 0,
                "Dependent_offTarget_number": len(dependent),
                "Independent_offTarget_number": len(independent),
                "Total_offTarget_number": len(dependent) + len(independent),
            }
        ]
    )

    non_target_present = _not_empty(dependent["NonTarget_reads"])
    target_present = _not_empty(dependent["Target_reads"])
    overlap_totals: defaultdict[str, int] = defaultdict(int)
    for value in dependent["overlap_Intervals_stat"]:
        for label, count in parse_overlap_stats(value).items():
            overlap_totals[label] += count
    overlap_text = ";".join(
        f"{label}({count})" for label, count in overlap_totals.items() if count > 0
    ) or "none"

    detail = pd.DataFrame(
        [
            {
                "样本名": sample_id,
                "Dependent_offTarget_number": len(dependent),
                "Only NonTarget strand": int((non_target_present & ~target_present).sum()),
                "Only Target strand": int((~non_target_present & target_present).sum()),
                "NonTarget strand & Target strand": int(
                    (non_target_present & target_present).sum()
                ),
                "out of protospacer_max_position": int(
                    (dependent["max_position_in_protospacer"] == "out of protospacer").sum()
                ),
                "Out of protospacer": overlap_text,
            }
        ]
    )
    return basic, detail


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--basic-tsv", required=True)
    parser.add_argument("--detail-tsv", required=True)
    parser.add_argument("--excel", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    frame = pd.read_csv(args.input, sep="\t", low_memory=False)
    basic, detail = summarize(frame, args.sample_id)
    basic_path = ensure_parent(args.basic_tsv)
    detail_path = ensure_parent(args.detail_tsv)
    excel_path = ensure_parent(args.excel)
    basic.to_csv(basic_path, sep="\t", index=False)
    detail.to_csv(detail_path, sep="\t", index=False)
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        basic.to_excel(writer, sheet_name="基本统计", index=False)
        detail.to_excel(writer, sheet_name="详细分析", index=False)
    print(f"[05 summary] sample={args.sample_id} -> {excel_path}")


if __name__ == "__main__":
    main()

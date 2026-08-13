#!/usr/bin/env python3
"""Collect per-sample Ino-seq off-target summaries into cohort outputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd


def sample_ids(pair_sheet: str | Path) -> list[str]:
    values: list[str] = []
    with open(pair_sheet, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["sample_id", "control_id", "sgrna"]:
            raise ValueError(
                "pair sheet header must be exactly: sample_id<TAB>control_id<TAB>sgrna"
            )
        for row in reader:
            if row["sample_id"] and not row["sample_id"].startswith("#"):
                values.append(row["sample_id"])
    return values


def collect(pair_sheet: str | Path, output_dir: str | Path, strict: bool = False):
    basic_frames: list[pd.DataFrame] = []
    detail_frames: list[pd.DataFrame] = []
    missing: list[str] = []
    root = Path(output_dir)
    for sample_id in sample_ids(pair_sheet):
        summary_dir = root / sample_id / "postprocess" / "summary"
        basic = summary_dir / f"{sample_id}_offtarget_summary.tsv"
        detail = summary_dir / f"{sample_id}_strand_summary.tsv"
        if not basic.is_file() or not detail.is_file():
            missing.append(sample_id)
            continue
        basic_frames.append(pd.read_csv(basic, sep="\t"))
        detail_frames.append(pd.read_csv(detail, sep="\t"))
    if strict and missing:
        raise FileNotFoundError(f"missing per-sample summaries: {', '.join(missing)}")
    return (
        pd.concat(basic_frames, ignore_index=True) if basic_frames else pd.DataFrame(),
        pd.concat(detail_frames, ignore_index=True) if detail_frames else pd.DataFrame(),
        missing,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--tsv-dir", required=True)
    parser.add_argument("--excel", required=True)
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    basic, detail, missing = collect(args.pairs, args.output_dir, args.strict)
    tsv_dir = Path(args.tsv_dir)
    tsv_dir.mkdir(parents=True, exist_ok=True)
    excel = Path(args.excel)
    excel.parent.mkdir(parents=True, exist_ok=True)
    basic.to_csv(tsv_dir / "inoseq_offtarget_summary.tsv", sep="\t", index=False)
    detail.to_csv(tsv_dir / "inoseq_strand_summary.tsv", sep="\t", index=False)
    with pd.ExcelWriter(excel, engine="openpyxl") as writer:
        basic.to_excel(writer, sheet_name="基本统计", index=False)
        detail.to_excel(writer, sheet_name="详细分析", index=False)
    print(
        f"[off-target QC] {len(basic)} sample(s) collected; "
        f"missing={len(missing)} -> {excel}"
    )


if __name__ == "__main__":
    main()

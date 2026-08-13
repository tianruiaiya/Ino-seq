#!/usr/bin/env python3
"""Assign on-target, sgRNA-dependent, or independent candidate labels."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from workflow.lib.io_utils import ensure_parent, require_columns

STANDARD_MISMATCHES = "Site_SubstitutionsOnly.NumSubstitutions"
GAPPED_SUBSTITUTIONS = "Site_GapsAllowed.Substitutions"


def _has_value(value: object) -> bool:
    return not pd.isna(value) and str(value).strip() != ""


def classify_row(row: pd.Series) -> str:
    standard = row[STANDARD_MISMATCHES]
    gapped = row[GAPPED_SUBSTITUTIONS]
    if not pd.isna(standard) and float(standard) == 0:
        return "onTarget"
    if (not pd.isna(standard) and float(standard) > 0) or _has_value(gapped):
        return "dependent"
    return "independent"


def mark_dependency(input_file: str | Path) -> pd.DataFrame:
    frame = pd.read_csv(input_file, sep="\t")
    require_columns(frame, [STANDARD_MISMATCHES, GAPPED_SUBSTITUTIONS], input_file)
    labels = frame.apply(classify_row, axis=1)
    frame.insert(7, "dependent_type", labels)
    return frame


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = mark_dependency(args.input)
    output = ensure_parent(args.output)
    result.to_csv(output, sep="\t", index=False)
    counts = result["dependent_type"].value_counts().to_dict()
    print(f"[05 dependency] {len(result)} candidate(s), labels={counts} -> {output}")


if __name__ == "__main__":
    main()

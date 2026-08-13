#!/usr/bin/env python3
"""Compare experimental and control overlap counts.

The statistical equations are a direct port of
``nonormalized_all_calculate_FDR.R``: unnormalized fold change, a 0.5
pseudocount only for zero-control fold changes, the same Poisson-tail cases,
and Benjamini-Hochberg FDR adjustment. The FDR is reported but is not used by
the downstream historical filter.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import poisson

from workflow.lib.io_utils import ensure_parent

COVERAGE_COLUMNS = [
    "chromosome",
    "start",
    "end",
    "location",
    "clevage_reads",
    "strand",
    "query_name",
    "reads",
]

OUTPUT_COLUMNS = [
    "chromosome",
    "start",
    "end",
    "exp_reads",
    "ctrl_reads",
    "original_ctrl_reads",
    "fold_change",
    "p_value",
    "FDR",
    "location",
    "clevage_reads",
    "strand",
    "query_name",
]


def read_coverage(path: str | Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", header=None)
    if frame.shape[1] != 8:
        raise ValueError(f"{path}: expected 8 coverage columns, observed {frame.shape[1]}")
    frame.columns = COVERAGE_COLUMNS
    return frame


def poisson_p_value(exp_reads: int, ctrl_reads: int) -> float:
    """Return the historical unnormalized Poisson-tail comparison p-value."""
    if exp_reads == 0 and ctrl_reads == 0:
        return 1.0
    if ctrl_reads == 0:
        return float(poisson.sf(exp_reads - 1, mu=0.1))
    if exp_reads == 0:
        return float(2 * poisson.cdf(0, mu=ctrl_reads))

    expected_each = (exp_reads + ctrl_reads) / 2
    if exp_reads > expected_each:
        p_value = 2 * poisson.sf(exp_reads - 1, mu=expected_each)
    else:
        p_value = 2 * poisson.cdf(exp_reads, mu=expected_each)
    return float(min(p_value, 1.0))


def benjamini_hochberg(values: pd.Series) -> np.ndarray:
    """Benjamini-Hochberg adjustment equivalent to R ``p.adjust(..., 'BH')``."""
    p_values = values.to_numpy(dtype=float)
    adjusted = np.full(p_values.shape, np.nan, dtype=float)
    valid_indices = np.flatnonzero(~np.isnan(p_values))
    if valid_indices.size == 0:
        return adjusted

    valid = p_values[valid_indices]
    order = np.argsort(valid, kind="mergesort")
    ranked = valid[order]
    count = len(ranked)
    scaled = ranked * count / np.arange(1, count + 1)
    scaled = np.minimum.accumulate(scaled[::-1])[::-1]
    scaled = np.minimum(scaled, 1.0)
    reverse_order = np.empty_like(order)
    reverse_order[order] = np.arange(count)
    adjusted[valid_indices] = scaled[reverse_order]
    return adjusted


def compare_background(exp_file: str | Path, ctrl_file: str | Path) -> pd.DataFrame:
    exp = read_coverage(exp_file)
    ctrl = read_coverage(ctrl_file)
    if len(exp) != len(ctrl):
        raise ValueError(
            f"coverage row-count mismatch: experiment={len(exp)}, control={len(ctrl)}"
        )

    keys = ["chromosome", "start", "end"]
    mismatch = (exp[keys].reset_index(drop=True) != ctrl[keys].reset_index(drop=True)).any(axis=1)
    if mismatch.any():
        first = int(np.flatnonzero(mismatch.to_numpy())[0]) + 1
        raise ValueError(f"experiment/control interval mismatch at data row {first}")

    result = pd.DataFrame(
        {
            "chromosome": exp["chromosome"],
            "start": exp["start"],
            "end": exp["end"],
            "exp_reads": exp["reads"].astype(int),
            "ctrl_reads": ctrl["reads"].astype(int),
            "location": exp["location"],
            "strand": exp["strand"],
            "query_name": exp["query_name"],
            "clevage_reads": exp["clevage_reads"],
        }
    )
    result["original_ctrl_reads"] = result["ctrl_reads"]
    adjusted_control = result["ctrl_reads"].where(result["ctrl_reads"] != 0, 0.5)
    result["fold_change"] = result["exp_reads"] / adjusted_control
    result["p_value"] = [
        poisson_p_value(int(exp_count), int(ctrl_count))
        for exp_count, ctrl_count in zip(result["exp_reads"], result["ctrl_reads"])
    ]
    result["FDR"] = benjamini_hochberg(result["p_value"])
    result["fold_change"] = result["fold_change"].round(6)

    return result[OUTPUT_COLUMNS].sort_values(
        ["FDR", "fold_change"], ascending=[True, False], kind="mergesort"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment", required=True)
    parser.add_argument("--control", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = compare_background(args.experiment, args.control)
    output = ensure_parent(args.output)
    result.to_csv(output, sep="\t", index=False, float_format="%.6g")
    print(f"[03 compare] {len(result)} site(s) -> {output}")


if __name__ == "__main__":
    main()

from __future__ import annotations

import math

import pandas as pd


def test_aggregate_sites_retains_historical_columns(tmp_path, load_module):
    module = load_module("03_aggregate_sites.py")
    source = tmp_path / "sample.end"
    pd.DataFrame(
        [
            ["q1", "chr1", 1, 50, 100, "+"],
            ["q2", "chr1", 2, 51, 100, "+"],
            ["q3", "chr2", 3, 52, 5, "-"],
        ],
        columns=[
            "query_name",
            "ref_name",
            "ref_start",
            "ref_end",
            "location",
            "cleavage_direction",
        ],
    ).to_csv(source, sep="\t", index=False)

    result = module.aggregate_sites(source, window=15)
    assert result.columns.tolist() == [
        "ref_name",
        "start",
        "end",
        "location",
        "reads count",
        "strand",
        "query_name",
    ]
    assert result.iloc[0].to_dict() == {
        "ref_name": "chr1",
        "start": 85,
        "end": 115,
        "location": 100,
        "reads count": 2,
        "strand": "+",
        "query_name": "q1,q2",
    }
    assert result.iloc[1]["start"] == 0


def test_poisson_cases_and_bh_adjustment(load_module):
    module = load_module("03_compare_background.py")
    assert module.poisson_p_value(0, 0) == 1.0
    assert math.isclose(module.poisson_p_value(1, 0), 1 - math.exp(-0.1), rel_tol=1e-12)
    assert math.isclose(module.poisson_p_value(0, 1), 2 * math.exp(-1), rel_tol=1e-12)
    adjusted = module.benjamini_hochberg(pd.Series([0.01, 0.04, 0.03]))
    assert adjusted.tolist() == [0.03, 0.04, 0.04]


def test_background_filter_uses_strict_p_less_than_005(tmp_path, load_module):
    module = load_module("03_filter_background.py")
    source = tmp_path / "background.tsv"
    pd.DataFrame(
        {
            "fold_change": [2.0, 2.0, 1.49, 3.0],
            "p_value": [0.049, 0.05, 0.001, 0.051],
            "query_name": ["pass", "equal", "low_fc", "high_p"],
        }
    ).to_csv(source, sep="\t", index=False)
    selected, names = module.filter_table(source)
    assert selected["query_name"].tolist() == ["pass"]
    assert names == ["pass"]

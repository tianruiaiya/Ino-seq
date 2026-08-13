from __future__ import annotations

import pandas as pd


def test_candidate_merge_extension_and_minimum_reads(tmp_path, load_module):
    module = load_module("04_candidate_intervals.py")
    source = tmp_path / "filtered.tsv"
    pd.DataFrame(
        {
            "chromosome": ["chr1", "chr1", "chr1"],
            "location": [100, 120, 500],
            "clevage_reads": [2, 1, 2],
            "query_name": ["q1,q2", "q3", "q4,q5"],
        }
    ).to_csv(source, sep="\t", index=False)
    result = module.build_candidates(source)
    assert len(result) == 1
    row = result.iloc[0]
    assert (row["start"], row["end"]) == (95, 124)
    assert row["total_reads"] == 3
    assert row["max_reads"] == 2
    assert row["max_position"] == 100
    assert row["Merge Intervals"] == "yes"
    assert row["query_name"] == "q1,q2,q3"

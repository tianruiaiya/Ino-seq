from __future__ import annotations

import pandas as pd


def test_dependency_labels_use_legacy_fields(tmp_path, load_module):
    module = load_module("05_mark_dependency.py")
    source = tmp_path / "aligned.tsv"
    base = {
        "chromosome": ["chr1"] * 4,
        "start": [1] * 4,
        "end": [30] * 4,
        "total_reads": [3] * 4,
        "max_reads": [3] * 4,
        "max_position": [15] * 4,
        "Merge Intervals": ["no"] * 4,
        "query_name": ["a", "b", "c", "d"],
        "Site_SubstitutionsOnly.NumSubstitutions": [0, 2, None, None],
        "Site_GapsAllowed.Substitutions": [None, None, None, 1],
    }
    pd.DataFrame(base).to_csv(source, sep="\t", index=False)
    result = module.mark_dependency(source)
    assert result.columns[7] == "dependent_type"
    assert result["dependent_type"].tolist() == [
        "onTarget",
        "dependent",
        "independent",
        "dependent",
    ]


def test_summary_counts_named_columns(load_module):
    module = load_module("05_summarize.py")
    frame = pd.DataFrame(
        {
            "dependent_type": ["onTarget", "dependent", "independent"],
            "total_reads": [5, 4, 3],
            "NonTarget_reads": ["", "n1(+)", ""],
            "Target_reads": ["", "t1(-)", ""],
            "overlap_Intervals_stat": ["", "dependent(1);independent(2)", ""],
            "max_position_in_protospacer": ["", "out of protospacer", ""],
        }
    )
    basic, detail = module.summarize(frame, "sampleA")
    assert basic.iloc[0]["OnTarget_reads_number"] == 5
    assert basic.iloc[0]["Total_offTarget_number"] == 2
    assert detail.iloc[0]["NonTarget strand & Target strand"] == 1
    assert detail.iloc[0]["out of protospacer_max_position"] == 1
    assert detail.iloc[0]["Out of protospacer"] == "dependent(1);independent(2)"

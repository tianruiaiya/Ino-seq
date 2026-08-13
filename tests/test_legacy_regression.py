from __future__ import annotations

import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
from test_end_to_end import write_bam, write_reference

PROJECT_DIR = Path(__file__).resolve().parents[1]
LEGACY_DIR = Path(os.environ.get("INOSEQ_LEGACY_DIR", "/missing/legacy/scripts"))

pytestmark = pytest.mark.skipif(
    not LEGACY_DIR.is_dir(), reason="set INOSEQ_LEGACY_DIR to run source-regression tests"
)


def run_python(script: Path, *arguments: object) -> None:
    subprocess.run(
        [os.sys.executable, str(script), *map(str, arguments)],
        check=True,
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONPATH": str(PROJECT_DIR)},
    )


def import_script(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def assert_tsv_equal(left: Path, right: Path) -> None:
    pd.testing.assert_frame_equal(
        pd.read_csv(left, sep="\t"),
        pd.read_csv(right, sep="\t"),
        check_dtype=False,
    )


def test_fuzzy_alignment_matches_legacy_for_substitution_and_bulges():
    legacy = import_script(LEGACY_DIR / "identify.py", "legacy_identify")
    refactored = import_script(
        PROJECT_DIR / "workflow/modules/05_align_sgrna.py", "refactored_identify"
    )
    target = "ACGTACGTACGTACGTACGTNGG"
    explicit = target.replace("N", "A")
    windows = [
        "TTTT" + explicit + "TTTT",
        "TTTT" + explicit[:7] + "T" + explicit[8:] + "TTTT",
        "TTTT" + explicit[:9] + "A" + explicit[9:] + "TTTT",
        "TTTT" + explicit[:9] + explicit[10:] + "TTTT",
    ]
    for window in windows:
        assert refactored.align_sequences(target, window, 8) == legacy.alignSequences(
            target, window, 8
        )


def test_refactored_core_matches_legacy_valid_input(tmp_path):
    sample = "sample"
    sgrna = "ACGTACGTACGTACGTACGTNGG"
    reference = tmp_path / "reference.fa"
    write_reference(reference, sgrna.replace("N", "A"))
    bam = tmp_path / "input.bam"
    read_names = [f"q{index}" for index in range(8)]
    write_bam(bam, read_names)

    end_file = tmp_path / f"{sample}.end"
    pd.DataFrame(
        {
            "query_name": read_names,
            "ref_name": ["chr1"] * 8,
            "location": [500] * 8,
            "cleavage_direction": ["+"] * 8,
        }
    ).to_csv(end_file, sep="\t", index=False)

    legacy_aggregate = tmp_path / "legacy_aggregate"
    new_aggregate = tmp_path / "new_aggregate"
    legacy_aggregate.mkdir()
    new_aggregate.mkdir()
    run_python(LEGACY_DIR / "get_result.py", end_file, legacy_aggregate)
    run_python(
        PROJECT_DIR / "workflow/modules/03_aggregate_sites.py",
        "--input",
        end_file,
        "--output-tsv",
        new_aggregate / f"{sample}_merge.tsv",
        "--output-bed",
        new_aggregate / f"{sample}.bed",
    )
    assert_tsv_equal(
        legacy_aggregate / f"{sample}_merge.tsv", new_aggregate / f"{sample}_merge.tsv"
    )
    assert (legacy_aggregate / f"{sample}.bed").read_bytes() == (
        new_aggregate / f"{sample}.bed"
    ).read_bytes()

    filtered = tmp_path / "filtered.tsv"
    pd.DataFrame(
        {
            "chromosome": ["chr1"],
            "location": [500],
            "clevage_reads": [8],
            "query_name": [",".join(read_names)],
        }
    ).to_csv(filtered, sep="\t", index=False)
    legacy_candidates = tmp_path / "legacy_candidates.tsv"
    new_candidates = tmp_path / "new_candidates.tsv"
    legacy_bam = tmp_path / "legacy_candidates.bam"
    new_bam = tmp_path / "new_candidates.bam"
    run_python(
        LEGACY_DIR / "final_merged_distance.py",
        "--input",
        filtered,
        "--output",
        legacy_candidates,
        "--bam",
        bam,
        "--output_bam",
        legacy_bam,
    )
    run_python(
        PROJECT_DIR / "workflow/modules/04_candidate_intervals.py",
        "--input",
        filtered,
        "--output",
        new_candidates,
        "--input-bam",
        bam,
        "--output-bam",
        new_bam,
    )
    assert_tsv_equal(legacy_candidates, new_candidates)

    legacy_align = tmp_path / "legacy_align.tsv"
    new_align = tmp_path / "new_align.tsv"
    run_python(
        LEGACY_DIR / "identify.py",
        "--input",
        legacy_candidates,
        "--ref",
        reference,
        "--output",
        legacy_align,
        "--window",
        25,
        "--max_score",
        8,
        "--sgrna",
        sgrna,
    )
    run_python(
        PROJECT_DIR / "workflow/modules/05_align_sgrna.py",
        "--input",
        new_candidates,
        "--ref",
        reference,
        "--output",
        new_align,
        "--window",
        25,
        "--max-score",
        8,
        "--sgrna",
        sgrna,
    )
    assert_tsv_equal(legacy_align, new_align)

    legacy_dependency = tmp_path / "legacy_dependency.tsv"
    new_dependency = tmp_path / "new_dependency.tsv"
    run_python(LEGACY_DIR / "dependent_mark.py", legacy_align, legacy_dependency)
    run_python(
        PROJECT_DIR / "workflow/modules/05_mark_dependency.py",
        "--input",
        new_align,
        "--output",
        new_dependency,
    )
    assert_tsv_equal(legacy_dependency, new_dependency)

    legacy_spacer = tmp_path / "legacy_spacer.tsv"
    new_spacer = tmp_path / "new_spacer.tsv"
    run_python(
        LEGACY_DIR / "detect_out_of_spacer.py",
        "--input",
        legacy_dependency,
        "--output",
        legacy_spacer,
    )
    run_python(
        PROJECT_DIR / "workflow/modules/05_annotate_spacer.py",
        "--input",
        new_dependency,
        "--output",
        new_spacer,
    )
    assert_tsv_equal(legacy_spacer, new_spacer)

    legacy_target = tmp_path / "legacy_target.tsv"
    new_target = tmp_path / "new_target.tsv"
    run_python(
        LEGACY_DIR / "detect_target.py",
        "--input",
        legacy_spacer,
        "--bam",
        legacy_bam,
        "--output",
        legacy_target,
    )
    run_python(
        PROJECT_DIR / "workflow/modules/05_classify_strands.py",
        "--input",
        new_spacer,
        "--bam",
        new_bam,
        "--output",
        new_target,
    )
    assert_tsv_equal(legacy_target, new_target)

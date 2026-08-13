from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pandas as pd
import pysam

PROJECT_DIR = Path(__file__).resolve().parents[1]


def write_reference(path: Path, sgrna_explicit: str) -> None:
    sequence = list("A" * 1000)
    sequence[488 : 488 + len(sgrna_explicit)] = sgrna_explicit
    path.write_text(">chr1\n" + "".join(sequence) + "\n")
    pysam.faidx(str(path))


def write_bam(path: Path, read_names: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for index, read_name in enumerate(read_names):
            read = pysam.AlignedSegment()
            read.query_name = read_name
            read.query_sequence = "A" * 50
            read.flag = 16 if index % 2 else 0
            read.reference_id = 0
            read.reference_start = 475
            read.mapping_quality = 60
            read.cigar = ((0, 50),)
            read.query_qualities = pysam.qualitystring_to_array("I" * 50)
            bam.write(read)
    pysam.index(str(path))


def test_synthetic_pair_runs_from_step01_outputs(tmp_path):
    sample = "experiment"
    control = "blank"
    sgrna = "ACGTACGTACGTACGTACGTNGG"
    explicit_target = sgrna.replace("N", "A")
    output_dir = tmp_path / "output"
    reference = tmp_path / "reference.fa"
    write_reference(reference, explicit_target)

    sample_alignment = output_dir / sample / "alignment"
    control_alignment = output_dir / control / "alignment"
    read_names = [f"q{index}" for index in range(8)]
    write_bam(sample_alignment / f"{sample}_end.bam", read_names)
    write_bam(control_alignment / f"{control}_end.bam", [])
    sample_alignment.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "query_name": read_names,
            "ref_name": ["chr1"] * 8,
            "ref_start": [476] * 8,
            "ref_end": [525] * 8,
            "location": [500] * 8,
            "cleavage_direction": ["+"] * 8,
            "rbase": ["A"] * 8,
            "qbase": ["G"] * 8,
            "mutations": ["500:A>G"] * 8,
            "flank_seq": ["A" * 41] * 8,
            "Poly": [1] * 8,
        }
    ).to_csv(sample_alignment / f"{sample}.end", sep="\t", index=False)

    config = tmp_path / "inoseq.env"
    config.write_text(
        f"REFERENCE_FASTA={reference}\n"
        f"PYTHON_BIN={os.sys.executable}\n"
        f"OUTPUT_DIR={output_dir}\n"
        "BACKGROUND_PVALUE=0.05\n"
    )
    environment = os.environ.copy()
    environment["INOSEQ_PROJECT_DIR"] = str(PROJECT_DIR)
    completed = subprocess.run(
        [
            "bash",
            str(PROJECT_DIR / "workflow" / "run_postprocess.sh"),
            sample,
            control,
            sgrna,
            str(config),
        ],
        check=True,
        capture_output=True,
        text=True,
        env=environment,
    )
    assert "full paired analysis complete" in completed.stderr

    post = output_dir / sample / "postprocess"
    assert (post / "INOSEQ_POSTPROCESS_COMPLETE").is_file()
    assert (output_dir / sample / "INOSEQ_FULL_COMPLETE").is_file()
    parameters = pd.read_csv(post / "run_parameters.tsv", sep="\t")
    parameter_values = parameters.set_index("parameter")["value"]
    assert parameter_values.loc["inoseq_version"] == "1.0.0"
    assert float(parameter_values.loc["background_pvalue"]) == 0.05
    summary = pd.read_csv(post / "summary" / f"{sample}_offtarget_summary.tsv", sep="\t")
    assert summary.iloc[0]["OnTarget_number"] == 1
    assert summary.iloc[0]["OnTarget_reads_number"] == 8
    final_table = pd.read_csv(post / "offtarget" / f"{sample}_dependent_target.txt", sep="\t")
    assert final_table["dependent_type"].tolist() == ["onTarget"]

    pair_sheet = tmp_path / "pairs.tsv"
    pair_sheet.write_text(
        f"sample_id\tcontrol_id\tsgrna\n{sample}\t{control}\t{sgrna}\n"
    )
    cohort_qc = output_dir / "QC"
    subprocess.run(
        [
            os.sys.executable,
            str(PROJECT_DIR / "workflow/qc/collect_offtarget_stats.py"),
            "--pairs",
            str(pair_sheet),
            "--output-dir",
            str(output_dir),
            "--tsv-dir",
            str(cohort_qc),
            "--excel",
            str(cohort_qc / "dependent_target_analysis.xlsx"),
            "--strict",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert (cohort_qc / "inoseq_offtarget_summary.tsv").is_file()
    assert (cohort_qc / "dependent_target_analysis.xlsx").is_file()

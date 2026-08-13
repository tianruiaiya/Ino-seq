from __future__ import annotations

import os
import subprocess
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]
SUBMIT_SCRIPT = PROJECT_DIR / "workflow" / "submit_full_workflow.sh"
FINALIZE_SCRIPT = PROJECT_DIR / "workflow" / "utils" / "finalize_full_workflow.sh"


def write_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample_id\tread1\tread2\n"
        "blank\t/data/blank_R1.fastq.gz\t/data/blank_R2.fastq.gz\n"
        "experiment\t/data/experiment_R1.fastq.gz\t/data/experiment_R2.fastq.gz\n"
    )
    pairs = tmp_path / "pairs.tsv"
    pairs.write_text(
        "sample_id\tcontrol_id\tsgrna\n"
        "experiment\tblank\tACGTACGTACGTACGTACGTNGG\n"
    )
    config = tmp_path / "inoseq.env"
    config.write_text(
        f"OUTPUT_DIR={tmp_path / 'output'}\n"
        f"LOG_DIR={tmp_path / 'logs'}\n"
        "PYTHON_BIN=python\n"
    )
    return samples, pairs, config


def test_full_workflow_dry_run_has_complete_afterok_graph(tmp_path: Path):
    samples, pairs, config = write_inputs(tmp_path)
    env = os.environ.copy()
    env.update({"DRY_RUN": "1", "SKIP_VALIDATION": "1"})
    completed = subprocess.run(
        ["bash", str(SUBMIT_SCRIPT), str(samples), str(pairs), str(config)],
        check=True,
        capture_output=True,
        text=True,
        env=env,
    )

    output = completed.stdout
    assert "Phase A module 01+ blank -> DRY1" in output
    assert "Phase A module 01+ experiment -> DRY2" in output
    assert "Phase A cohort QC -> DRY3" in output
    assert "--dependency=afterok:DRY1:DRY2" in output
    assert "Phase B module 03+ experiment vs blank -> DRY4" in output
    assert "--dependency=afterok:DRY2:DRY1" in output
    assert "off-target cohort summary -> DRY5" in output
    assert "--dependency=afterok:DRY4" in output
    assert "Phase C full-workflow finalizer -> DRY6" in output
    assert "--dependency=afterok:DRY3:DRY5" in output
    assert "INOSEQ_WORKFLOW_COMPLETE" in output


def test_full_workflow_records_real_sbatch_job_graph(tmp_path: Path):
    samples, pairs, config = write_inputs(tmp_path)
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    counter = tmp_path / "sbatch.counter"
    calls = tmp_path / "sbatch.calls"
    fake_sbatch = fake_bin / "sbatch"
    fake_sbatch.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        f"counter={counter!s}\n"
        f"calls={calls!s}\n"
        "value=1000\n"
        "[[ ! -f \"$counter\" ]] || value=$(cat \"$counter\")\n"
        "value=$((value + 1))\n"
        "printf '%s\\n' \"$value\" > \"$counter\"\n"
        "printf '%q ' \"$@\" >> \"$calls\"\n"
        "printf '\\n' >> \"$calls\"\n"
        "printf '%s\\n' \"$value\"\n"
    )
    fake_sbatch.chmod(0o755)

    output_root = tmp_path / "output"
    output_root.mkdir()
    stale_marker = output_root / "INOSEQ_WORKFLOW_COMPLETE"
    stale_marker.touch()

    env = os.environ.copy()
    env.update(
        {
            "SKIP_VALIDATION": "1",
            "PATH": f"{fake_bin}{os.pathsep}{env['PATH']}",
        }
    )
    completed = subprocess.run(
        ["bash", str(SUBMIT_SCRIPT), str(samples), str(pairs), str(config)],
        check=True,
        capture_output=True,
        text=True,
        env=env,
    )

    assert "final job 1006" in completed.stdout
    assert not stale_marker.exists()
    log_dir = tmp_path / "logs"
    graph = (log_dir / "inoseq_last_full_workflow_jobs.tsv").read_text()
    assert "1001\tphase_a\tblank\t-" in graph
    assert "1002\tphase_a\texperiment\t-" in graph
    assert "1003\tphase_a_qc\tcohort\tafterok:1001:1002" in graph
    assert "1004\tphase_b\texperiment\tafterok:1002:1001" in graph
    assert "1005\tphase_b_qc\tcohort\tafterok:1004" in graph
    assert "1006\tfinalize\tworkflow\tafterok:1003:1005" in graph
    assert (log_dir / "inoseq_last_full_workflow_job_id.txt").read_text().strip() == "1006"


def test_finalizer_requires_and_records_all_enabled_outputs(tmp_path: Path):
    samples, pairs, config = write_inputs(tmp_path)
    output_root = tmp_path / "output"
    for sample in ["blank", "experiment"]:
        (output_root / sample).mkdir(parents=True, exist_ok=True)
        (output_root / sample / "INOSEQ_COMPLETE").touch()
    postprocess = output_root / "experiment" / "postprocess"
    postprocess.mkdir()
    (postprocess / "INOSEQ_POSTPROCESS_COMPLETE").touch()
    (output_root / "experiment" / "INOSEQ_FULL_COMPLETE").touch()
    qc_dir = output_root / "QC"
    qc_dir.mkdir()
    for filename in [
        "inoseq_qc_summary.tsv",
        "inoseq_qc_summary.csv",
        "inoseq_offtarget_summary.tsv",
        "inoseq_strand_summary.tsv",
        "dependent_target_analysis.xlsx",
    ]:
        (qc_dir / filename).touch()

    completed = subprocess.run(
        [
            "bash",
            str(FINALIZE_SCRIPT),
            str(samples),
            str(pairs),
            str(config),
            "1",
            "1",
        ],
        check=True,
        capture_output=True,
        text=True,
        env={**os.environ, "INOSEQ_PROJECT_DIR": str(PROJECT_DIR)},
    )

    assert "full workflow completed" in completed.stdout
    assert (output_root / "INOSEQ_WORKFLOW_COMPLETE").is_file()
    status = (qc_dir / "full_workflow_status.tsv").read_text()
    assert "status\tCOMPLETED" in status
    assert "inoseq_version\t1.0.0" in status
    assert "sample_count\t2" in status
    assert "pair_count\t1" in status
    assert "background_pvalue\t0.05" in status

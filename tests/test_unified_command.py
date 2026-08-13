from __future__ import annotations

import os
import subprocess
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]
CLI = PROJECT_DIR / "inoseq"


def run_cli(*args: object, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(CLI), *map(str, args)],
        check=True,
        capture_output=True,
        text=True,
        env=env,
    )


def test_help_and_version_present_one_public_interface():
    help_result = run_cli("help")
    assert "Official full-workflow path" in help_result.stdout
    assert "./inoseq submit" in help_result.stdout
    assert "Phase A" in help_result.stdout
    assert run_cli("version").stdout.strip() == "Ino-seq v1.0.0"


def test_init_creates_three_configs_without_overwriting(tmp_path: Path):
    result = run_cli("init", tmp_path)
    assert "[CREATE]" in result.stdout
    for name in ["inoseq.env", "samples.tsv", "pairs.tsv"]:
        assert (tmp_path / name).is_file()

    protected = tmp_path / "inoseq.env"
    protected.write_text("DO_NOT_OVERWRITE=1\n")
    second = run_cli("init", tmp_path)
    assert "[KEEP]" in second.stdout
    assert protected.read_text() == "DO_NOT_OVERWRITE=1\n"


def test_plan_routes_to_complete_dependency_graph(tmp_path: Path):
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
    environment = {**os.environ, "SKIP_VALIDATION": "1"}
    result = run_cli("plan", samples, pairs, config, env=environment)
    assert "Phase A module 01+ blank -> DRY1" in result.stdout
    assert "Phase B module 03+ experiment vs blank -> DRY4" in result.stdout
    assert "Phase C full-workflow finalizer -> DRY6" in result.stdout


def test_status_reads_workflow_level_completion(tmp_path: Path):
    output = tmp_path / "output"
    logs = tmp_path / "logs"
    (output / "QC").mkdir(parents=True)
    logs.mkdir()
    (output / "INOSEQ_WORKFLOW_COMPLETE").touch()
    (output / "QC" / "full_workflow_status.tsv").write_text(
        "field\tvalue\nstatus\tCOMPLETED\ninoseq_version\t1.0.0\n"
    )
    (logs / "inoseq_last_full_workflow_job_id.txt").write_text("12345\n")
    config = tmp_path / "inoseq.env"
    config.write_text(f"OUTPUT_DIR={output}\nLOG_DIR={logs}\n")

    result = run_cli("status", config)
    assert "[STATUS] COMPLETED" in result.stdout
    assert "[JOB] Final Slurm job: 12345" in result.stdout

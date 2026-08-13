from __future__ import annotations

import importlib.util
import os
import subprocess
import sys
from pathlib import Path
from types import ModuleType

PROJECT_DIR = Path(__file__).resolve().parents[1]
SUBMIT_SCRIPT = PROJECT_DIR / "workflow" / "submit_full_workflow.sh"
STATE_SCRIPT = PROJECT_DIR / "workflow" / "utils" / "stage_state.py"


def load_state_module() -> ModuleType:
    spec = importlib.util.spec_from_file_location("inoseq_stage_state", STATE_SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


STATE = load_state_module()


def write_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample_id\tread1\tread2\n"
        f"blank\t{tmp_path / 'blank_R1.fastq.gz'}\t{tmp_path / 'blank_R2.fastq.gz'}\n"
        f"experiment\t{tmp_path / 'experiment_R1.fastq.gz'}\t"
        f"{tmp_path / 'experiment_R2.fastq.gz'}\n"
    )
    for name in [
        "blank_R1.fastq.gz",
        "blank_R2.fastq.gz",
        "experiment_R1.fastq.gz",
        "experiment_R2.fastq.gz",
    ]:
        (tmp_path / name).write_bytes(b"FASTQ")
    pairs = tmp_path / "pairs.tsv"
    pairs.write_text(
        "sample_id\tcontrol_id\tsgrna\n"
        "experiment\tblank\tACGTACGTACGTACGTACGTNGG\n"
    )
    reference = tmp_path / "hg38.fa"
    reference.write_text(">chr1\nACGT\n")
    config = tmp_path / "inoseq.env"
    config.write_text(
        f"REFERENCE_FASTA={reference}\n"
        f"OUTPUT_DIR={tmp_path / 'output'}\n"
        f"LOG_DIR={tmp_path / 'logs'}\n"
        "PYTHON_BIN=python\n"
        "SUBMIT_QC_AFTER=1\n"
        "SUBMIT_OFFTARGET_QC_AFTER=1\n"
    )
    return samples, pairs, config


def make_context(
    config: Path,
    *,
    sample: str | None = None,
    control: str | None = None,
    sgrna: str | None = None,
    read1: str | None = None,
    read2: str | None = None,
    samples: Path | None = None,
    pairs: Path | None = None,
):
    env = STATE.load_shell_environment(config)
    return STATE.Context(
        project=PROJECT_DIR,
        config=config,
        env=env,
        output=Path(env["OUTPUT_DIR"]),
        sample=sample,
        control=control,
        sgrna=sgrna,
        read1=read1,
        read2=read2,
        samples_sheet=samples,
        pairs_sheet=pairs,
    )


def create_output_files(stage: str, context) -> None:
    for path in STATE.required_outputs(stage, context):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"x")
    STATE.record(stage, context)


def create_current_run(tmp_path: Path) -> tuple[Path, Path, Path]:
    samples, pairs, config = write_inputs(tmp_path)
    for sample in ["blank", "experiment"]:
        context = make_context(
            config,
            sample=sample,
            read1=str(tmp_path / f"{sample}_R1.fastq.gz"),
            read2=str(tmp_path / f"{sample}_R2.fastq.gz"),
        )
        create_output_files("module01", context)
        create_output_files("module02", context)
    pair_context = make_context(
        config,
        sample="experiment",
        control="blank",
        sgrna="ACGTACGTACGTACGTACGTNGG",
    )
    for stage in ["module03", "module04", "module05"]:
        create_output_files(stage, pair_context)
    cohort_context = make_context(config, samples=samples, pairs=pairs)
    create_output_files("phase-a-qc", cohort_context)
    create_output_files("phase-b-qc", cohort_context)
    create_output_files("finalize", cohort_context)
    return samples, pairs, config


def plan(samples: Path, pairs: Path, config: Path, *options: str) -> str:
    completed = subprocess.run(
        ["bash", str(SUBMIT_SCRIPT), *options, str(samples), str(pairs), str(config)],
        check=True,
        capture_output=True,
        text=True,
        env={**os.environ, "DRY_RUN": "1", "SKIP_VALIDATION": "1"},
    )
    return completed.stdout


def test_auto_resume_submits_nothing_when_every_stage_is_current(tmp_path: Path):
    samples, pairs, config = create_current_run(tmp_path)
    output = plan(samples, pairs, config)
    assert "[REUSE] Phase A blank" in output
    assert "[REUSE] Phase B experiment vs blank" in output
    assert "no Slurm jobs submitted" in output
    assert "[DRY-RUN]" not in output


def test_auto_resume_starts_at_module02_and_rebuilds_dynamic_dependencies(tmp_path: Path):
    samples, pairs, config = create_current_run(tmp_path)
    experiment = make_context(
        config,
        sample="experiment",
        read1=str(tmp_path / "experiment_R1.fastq.gz"),
        read2=str(tmp_path / "experiment_R2.fastq.gz"),
    )
    STATE.invalidate("module02", experiment)

    output = plan(samples, pairs, config)
    assert "[REUSE] Phase A blank" in output
    assert "Phase A experiment: start at module 02" in output
    assert "Phase A module 02+ experiment -> DRY1" in output
    assert "Phase A module 01+ blank" not in output
    assert "--dependency=afterok:DRY1" in output
    assert "Phase B module 03+ experiment vs blank -> DRY3" in output
    assert "Phase C full-workflow finalizer -> DRY5" in output


def test_explicit_module04_preserves_phase_a_and_module03(tmp_path: Path):
    samples, pairs, config = create_current_run(tmp_path)
    output = plan(samples, pairs, config, "--from-stage", "module04")
    assert "[REUSE] Phase A blank" in output
    assert "Phase A module" not in output
    assert "Phase B module 04+ experiment vs blank -> DRY1" in output
    assert "off-target cohort summary -> DRY2" in output
    assert "Phase C full-workflow finalizer -> DRY3" in output


def test_sgrna_change_invalidates_only_module05_and_downstream(tmp_path: Path):
    samples, pairs, config = create_current_run(tmp_path)
    pairs.write_text(
        "sample_id\tcontrol_id\tsgrna\n"
        "experiment\tblank\tTTTTACGTACGTACGTACGTNGG\n"
    )
    output = plan(samples, pairs, config)
    assert "[REUSE] Phase A blank" in output
    assert "Phase B experiment vs blank: start at module 05" in output
    assert "Phase B module 05+ experiment vs blank -> DRY1" in output
    assert "Phase B module 04+" not in output


def test_module03_parameter_change_invalidates_module03_chain(tmp_path: Path):
    samples, pairs, config = create_current_run(tmp_path)
    with config.open("a") as handle:
        handle.write("BACKGROUND_WINDOW=21\n")
    output = plan(samples, pairs, config)
    assert "[REUSE] Phase A blank" in output
    assert "Phase B experiment vs blank: start at module 03" in output
    assert "Phase B module 03+ experiment vs blank -> DRY1" in output

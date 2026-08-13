#!/usr/bin/env python3
"""Validate and record resumable Ino-seq stage state.

The state contract deliberately combines three signals before a stage can be
reused: a completion marker, all required outputs, and a reproducibility
fingerprint.  Fingerprints include the relevant input metadata, configuration
values, upstream fingerprints, workflow code, and Ino-seq version.

Only inexpensive file metadata are used for large FASTQ/reference inputs.  The
workflow source files and small sample/pair sheets are content-hashed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

SCHEMA_VERSION = 1
VALID = 0
INCOMPLETE = 10
STALE = 11

STAGES = (
    "module01",
    "module02",
    "module03",
    "module04",
    "module05",
    "phase-a-qc",
    "phase-b-qc",
    "finalize",
)

CONFIG_DEFAULTS = {
    "FASTP_THREADS": "40",
    "CUTADAPT_THREADS": "10",
    "BWA_THREADS": "40",
    "SAMTOOLS_SORT_THREADS": "4",
    "SAMTOOLS_SORT_MEM": "5G",
    "CONSENSUS_THREADS": "40",
    "FINAL_SORT_THREADS": "40",
    "STATS_THREADS": "40",
    "GROUP_JAVA_XMX": "180g",
    "CONSENSUS_JAVA_XMX": "180g",
    "ZIPPER_JAVA_XMX": "180g",
    "FASTP_QUALIFIED_QUALITY_PHRED": "20",
    "FASTP_UNQUALIFIED_PERCENT_LIMIT": "10",
    "FASTP_LENGTH_REQUIRED": "50",
    "READS_TO_PROCESS": "0",
    "QC_UNIQUE_MAPQ": "30",
    "UMI_R2_TRIM": "12",
    "UMI_ADAPTER": "TGTAGAGCACGCGTGG",
    "UMI_STRATEGY": "Adjacency",
    "UMI_EDITS": "1",
    "CONSENSUS_MIN_READS": "1",
    "CONSENSUS_MIN_INPUT_BASE_QUALITY": "20",
    "FILTER_MIN_READS": "1",
    "FILTER_MIN_BASE_QUALITY": "20",
    "FILTER_MAX_BASE_ERROR_RATE": "0.2",
    "MUTATION_LOCATION_QUAL_THRESHOLD": "30",
    "BACKGROUND_WINDOW": "15",
    "BACKGROUND_FOLD_CHANGE": "1.5",
    "BACKGROUND_PVALUE": "0.05",
    "CANDIDATE_MERGE_DISTANCE": "30",
    "CANDIDATE_MIN_LENGTH": "30",
    "CANDIDATE_MIN_READS": "3",
    "OFFTARGET_SEARCH_WINDOW": "25",
    "OFFTARGET_MAX_SCORE": "8",
    "SPACER_NEIGHBORHOOD": "100",
}

STAGE_CONFIG_KEYS = {
    "module01": (
        "REFERENCE_FASTA",
        "FGBIO_JAR",
        "FASTP_THREADS",
        "CUTADAPT_THREADS",
        "BWA_THREADS",
        "SAMTOOLS_SORT_THREADS",
        "SAMTOOLS_SORT_MEM",
        "CONSENSUS_THREADS",
        "FINAL_SORT_THREADS",
        "STATS_THREADS",
        "GROUP_JAVA_XMX",
        "CONSENSUS_JAVA_XMX",
        "ZIPPER_JAVA_XMX",
        "FASTP_QUALIFIED_QUALITY_PHRED",
        "FASTP_UNQUALIFIED_PERCENT_LIMIT",
        "FASTP_LENGTH_REQUIRED",
        "READS_TO_PROCESS",
        "QC_UNIQUE_MAPQ",
        "UMI_R2_TRIM",
        "UMI_ADAPTER",
        "UMI_STRATEGY",
        "UMI_EDITS",
        "CONSENSUS_MIN_READS",
        "CONSENSUS_MIN_INPUT_BASE_QUALITY",
        "FILTER_MIN_READS",
        "FILTER_MIN_BASE_QUALITY",
        "FILTER_MAX_BASE_ERROR_RATE",
    ),
    "module02": ("REFERENCE_FASTA", "MUTATION_LOCATION_QUAL_THRESHOLD"),
    "module03": ("BACKGROUND_WINDOW", "BACKGROUND_FOLD_CHANGE", "BACKGROUND_PVALUE"),
    "module04": (
        "CANDIDATE_MERGE_DISTANCE",
        "CANDIDATE_MIN_LENGTH",
        "CANDIDATE_MIN_READS",
    ),
    "module05": (
        "REFERENCE_FASTA",
        "OFFTARGET_SEARCH_WINDOW",
        "OFFTARGET_MAX_SCORE",
        "SPACER_NEIGHBORHOOD",
    ),
    "phase-a-qc": ("QC_UNIQUE_MAPQ", "READS_TO_PROCESS"),
    "phase-b-qc": (),
    "finalize": (
        "BACKGROUND_FOLD_CHANGE",
        "BACKGROUND_PVALUE",
        "SUBMIT_QC_AFTER",
        "SUBMIT_OFFTARGET_QC_AFTER",
    ),
}

STAGE_CODE_FILES = {
    "module01": ("workflow/modules/01_umi_consensus.sh",),
    "module02": ("workflow/modules/02_signature_reads.py",),
    "module03": (
        "workflow/modules/03_aggregate_sites.py",
        "workflow/modules/03_count_coverage.py",
        "workflow/modules/03_compare_background.py",
        "workflow/modules/03_filter_background.py",
        "workflow/lib/io_utils.py",
    ),
    "module04": (
        "workflow/modules/04_candidate_intervals.py",
        "workflow/lib/io_utils.py",
    ),
    "module05": (
        "workflow/modules/05_align_sgrna.py",
        "workflow/modules/05_mark_dependency.py",
        "workflow/modules/05_annotate_spacer.py",
        "workflow/modules/05_classify_strands.py",
        "workflow/modules/05_summarize.py",
        "workflow/lib/io_utils.py",
    ),
    "phase-a-qc": ("workflow/qc/collect_qc.py",),
    "phase-b-qc": ("workflow/qc/collect_offtarget_stats.py",),
    "finalize": ("workflow/utils/finalize_full_workflow.sh",),
}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_shell_environment(config: Path) -> dict[str, str]:
    command = 'set -a; source "$1"; env -0'
    completed = subprocess.run(
        ["bash", "-c", command, "stage-state", str(config)],
        check=True,
        capture_output=True,
    )
    result: dict[str, str] = {}
    for item in completed.stdout.split(b"\0"):
        if not item or b"=" not in item:
            continue
        key, value = item.split(b"=", 1)
        result[key.decode()] = value.decode()
    return result


def file_metadata(path_value: str | Path | None) -> dict[str, Any]:
    if not path_value:
        return {"path": "", "exists": False}
    path = Path(path_value).expanduser().resolve()
    if not path.exists():
        return {"path": str(path), "exists": False}
    stat = path.stat()
    return {
        "path": str(path),
        "exists": True,
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def content_descriptor(path: Path, label: str | None = None) -> dict[str, Any]:
    path = path.resolve()
    return {
        "path": label or str(path),
        "exists": path.is_file(),
        "sha256": sha256_file(path) if path.is_file() else None,
    }


def reference_descriptors(reference: str | None, include_alignment_indexes: bool) -> list[dict[str, Any]]:
    if not reference:
        return [file_metadata(reference)]
    fasta = Path(reference).expanduser().resolve()
    paths = [fasta, Path(f"{fasta}.fai")]
    if include_alignment_indexes:
        paths.extend(Path(f"{fasta}.{extension}") for extension in ("amb", "ann", "bwt", "pac", "sa"))
        if fasta.suffix in {".fa", ".fas", ".fasta"}:
            paths.append(fasta.with_suffix(".dict"))
        else:
            paths.append(Path(f"{fasta}.dict"))
    return [file_metadata(path) for path in paths]


def read_tsv(path: Path, expected: list[str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != expected:
            raise ValueError(f"Unexpected header in {path}: {reader.fieldnames}")
        for row in reader:
            first = (row.get(expected[0]) or "").strip()
            if first and not first.startswith("#"):
                rows.append({key: (row.get(key) or "").strip() for key in expected})
    return rows


@dataclass
class Context:
    project: Path
    config: Path
    env: dict[str, str]
    output: Path
    sample: str | None = None
    control: str | None = None
    sgrna: str | None = None
    read1: str | None = None
    read2: str | None = None
    samples_sheet: Path | None = None
    pairs_sheet: Path | None = None

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> Context:
        project = Path(args.project_dir).expanduser().resolve()
        config = Path(args.config).expanduser().resolve()
        env = load_shell_environment(config)
        output = Path(env.get("OUTPUT_DIR", "output")).expanduser()
        if not output.is_absolute():
            output = project / output
        return cls(
            project=project,
            config=config,
            env=env,
            output=output.resolve(),
            sample=args.sample,
            control=args.control,
            sgrna=args.sgrna,
            read1=args.read1,
            read2=args.read2,
            samples_sheet=Path(args.samples_sheet).resolve() if args.samples_sheet else None,
            pairs_sheet=Path(args.pairs_sheet).resolve() if args.pairs_sheet else None,
        )

    @property
    def version(self) -> str:
        return (self.project / "VERSION").read_text().strip()


def require(value: str | Path | None, label: str) -> str | Path:
    if value is None or str(value) == "":
        raise ValueError(f"{label} is required for this stage")
    return value


def manifest_path(stage: str, ctx: Context) -> Path:
    if stage in {"module01", "module02"}:
        sample = str(require(ctx.sample, "--sample"))
        return ctx.output / sample / ".inoseq" / f"{stage}.json"
    if stage in {"module03", "module04", "module05"}:
        sample = str(require(ctx.sample, "--sample"))
        return ctx.output / sample / "postprocess" / ".inoseq" / f"{stage}.json"
    if stage in {"phase-a-qc", "phase-b-qc"}:
        return ctx.output / "QC" / ".inoseq" / f"{stage}.json"
    return ctx.output / ".inoseq" / "finalize.json"


def marker_paths(stage: str, ctx: Context) -> list[Path]:
    sample = ctx.sample or ""
    mapping = {
        "module01": [ctx.output / sample / ".inoseq" / "MODULE01_COMPLETE"],
        "module02": [ctx.output / sample / "INOSEQ_COMPLETE"],
        "module03": [ctx.output / sample / "postprocess" / ".inoseq" / "MODULE03_COMPLETE"],
        "module04": [ctx.output / sample / "postprocess" / ".inoseq" / "MODULE04_COMPLETE"],
        "module05": [
            ctx.output / sample / "postprocess" / "INOSEQ_POSTPROCESS_COMPLETE",
            ctx.output / sample / "INOSEQ_FULL_COMPLETE",
        ],
        "phase-a-qc": [ctx.output / "QC" / ".inoseq" / "PHASE_A_QC_COMPLETE"],
        "phase-b-qc": [ctx.output / "QC" / ".inoseq" / "PHASE_B_QC_COMPLETE"],
        "finalize": [ctx.output / "INOSEQ_WORKFLOW_COMPLETE"],
    }
    return mapping[stage]


def required_outputs(stage: str, ctx: Context) -> list[Path]:
    sample = ctx.sample or ""
    sample_dir = ctx.output / sample
    qc = sample_dir / "QC"
    aln = sample_dir / "alignment"
    prefix = aln / sample
    post = sample_dir / "postprocess"
    control = ctx.control or ""
    mapping = {
        "module01": [
            qc / f"{sample}_fastp.json",
            qc / f"{sample}_fastp.html",
            qc / f"{sample}_cutadapt.json",
            Path(f"{prefix}.mapped.bam"),
            Path(f"{prefix}.mapped.bam.bai"),
            Path(f"{prefix}.mapped.flagstat.json"),
            Path(f"{prefix}.tag-family-sizes.txt"),
            Path(f"{prefix}.cons.unmapped.bam"),
            Path(f"{prefix}.umi_dedup.bam"),
            Path(f"{prefix}.umi_dedup.bam.bai"),
            Path(f"{prefix}.flagstat.json"),
            Path(f"{prefix}.idxstat.txt"),
        ],
        "module02": [
            Path(f"{prefix}.end"),
            Path(f"{prefix}_end.bam"),
            Path(f"{prefix}_end.bam.bai"),
        ],
        "module03": [
            post / "filt-before" / f"{sample}_merge.tsv",
            post / "filt-before" / f"{sample}.bed",
            post / "filt-before" / f"{sample}_coverage.txt",
            post / "filt-before" / f"{sample}_ctr_{control}_coverage.txt",
            post / "filt-before" / f"{sample}.txt",
            post / "filt-after" / f"{sample}_filted.txt",
            post / "filt-after" / f"{sample}_query_name.txt",
            post / "filt-after" / f"{sample}_filted.bam",
            post / "filt-after" / f"{sample}_filted.bam.bai",
        ],
        "module04": [
            post / "final_merge_distance" / f"{sample}_final_merged_distance.txt",
            post / "final_merge_distance" / f"{sample}_filted.bam",
            post / "final_merge_distance" / f"{sample}_filted.bam.bai",
        ],
        "module05": [
            post / "offtarget" / f"{sample}_align.txt",
            post / "offtarget" / f"{sample}_dependent_mark.txt",
            post / "offtarget" / f"{sample}_dependent_out_of_spacer.txt",
            post / "offtarget" / f"{sample}_dependent_target.txt",
            post / "summary" / f"{sample}_offtarget_summary.tsv",
            post / "summary" / f"{sample}_strand_summary.tsv",
            post / "summary" / "dependent_target_analysis.xlsx",
            post / "run_parameters.tsv",
        ],
        "phase-a-qc": [
            ctx.output / "QC" / "inoseq_qc_summary.tsv",
            ctx.output / "QC" / "inoseq_qc_summary.csv",
        ],
        "phase-b-qc": [
            ctx.output / "QC" / "inoseq_offtarget_summary.tsv",
            ctx.output / "QC" / "inoseq_strand_summary.tsv",
            ctx.output / "QC" / "dependent_target_analysis.xlsx",
        ],
        "finalize": [ctx.output / "QC" / "full_workflow_status.tsv"],
    }
    return mapping[stage]


def upstream_fingerprint(stage: str, ctx: Context, sample: str | None = None) -> str:
    prior = {
        "module02": "module01",
        "module04": "module03",
        "module05": "module04",
    }[stage]
    prior_ctx = Context(**{**ctx.__dict__, "sample": sample or ctx.sample})
    path = manifest_path(prior, prior_ctx)
    if not path.is_file():
        return "MISSING"
    try:
        return str(json.loads(path.read_text())["fingerprint"])
    except (KeyError, json.JSONDecodeError):
        return "INVALID"


def manifest_fingerprints(stage: str, ctx: Context) -> list[dict[str, str]]:
    result: list[dict[str, str]] = []
    if stage == "phase-a-qc":
        sheet = Path(require(ctx.samples_sheet, "--samples-sheet"))
        for row in read_tsv(sheet, ["sample_id", "read1", "read2"]):
            sample_ctx = Context(**{**ctx.__dict__, "sample": row["sample_id"]})
            path = manifest_path("module02", sample_ctx)
            value = "MISSING"
            if path.is_file():
                try:
                    value = str(json.loads(path.read_text())["fingerprint"])
                except (KeyError, json.JSONDecodeError):
                    value = "INVALID"
            result.append({"sample": row["sample_id"], "fingerprint": value})
    elif stage == "phase-b-qc":
        sheet = Path(require(ctx.pairs_sheet, "--pairs-sheet"))
        for row in read_tsv(sheet, ["sample_id", "control_id", "sgrna"]):
            sample_ctx = Context(**{**ctx.__dict__, "sample": row["sample_id"]})
            path = manifest_path("module05", sample_ctx)
            value = "MISSING"
            if path.is_file():
                try:
                    value = str(json.loads(path.read_text())["fingerprint"])
                except (KeyError, json.JSONDecodeError):
                    value = "INVALID"
            result.append({"sample": row["sample_id"], "fingerprint": value})
    return result


def fingerprint_basis(stage: str, ctx: Context) -> dict[str, Any]:
    config = {
        key: ctx.env.get(key, CONFIG_DEFAULTS.get(key, ""))
        for key in STAGE_CONFIG_KEYS[stage]
    }
    codes = [
        content_descriptor(ctx.project / item, label=item) for item in STAGE_CODE_FILES[stage]
    ]
    basis: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "inoseq_version": ctx.version,
        "stage": stage,
        "config": config,
        "code": codes,
        "environment": content_descriptor(
            ctx.project / "envs/inoseq.yml", label="envs/inoseq.yml"
        ),
    }
    reference = ctx.env.get("REFERENCE_FASTA")
    if stage in {"module01", "module02", "module05"}:
        basis["reference"] = reference_descriptors(
            reference, include_alignment_indexes=stage == "module01"
        )

    if stage == "module01":
        basis.update(
            {
                "sample": str(require(ctx.sample, "--sample")),
                "read1": file_metadata(require(ctx.read1, "--read1")),
                "read2": file_metadata(require(ctx.read2, "--read2")),
                "fgbio_jar": file_metadata(ctx.env.get("FGBIO_JAR")),
            }
        )
    elif stage == "module02":
        basis.update(
            {
                "sample": str(require(ctx.sample, "--sample")),
                "upstream": upstream_fingerprint(stage, ctx),
            }
        )
    elif stage == "module03":
        sample = str(require(ctx.sample, "--sample"))
        control = str(require(ctx.control, "--control"))
        sample_ctx = Context(**{**ctx.__dict__, "sample": sample})
        control_ctx = Context(**{**ctx.__dict__, "sample": control})
        basis.update(
            {
                "sample": sample,
                "control": control,
                "sample_module02": manifest_value("module02", sample_ctx),
                "control_module02": manifest_value("module02", control_ctx),
            }
        )
    elif stage == "module04":
        basis.update(
            {
                "sample": str(require(ctx.sample, "--sample")),
                "control": str(require(ctx.control, "--control")),
                "upstream": upstream_fingerprint(stage, ctx),
            }
        )
    elif stage == "module05":
        basis.update(
            {
                "sample": str(require(ctx.sample, "--sample")),
                "control": str(require(ctx.control, "--control")),
                "sgrna": str(require(ctx.sgrna, "--sgrna")),
                "upstream": upstream_fingerprint(stage, ctx),
            }
        )
    elif stage in {"phase-a-qc", "phase-b-qc"}:
        sheet = ctx.samples_sheet if stage == "phase-a-qc" else ctx.pairs_sheet
        basis["sheet"] = content_descriptor(
            Path(require(sheet, "cohort sheet")), label=f"{stage}-input-sheet"
        )
        basis["upstream"] = manifest_fingerprints(stage, ctx)
    elif stage == "finalize":
        samples = Path(require(ctx.samples_sheet, "--samples-sheet"))
        pairs = Path(require(ctx.pairs_sheet, "--pairs-sheet"))
        expect_a = ctx.env.get("SUBMIT_QC_AFTER", "1") == "1"
        expect_b = ctx.env.get("SUBMIT_OFFTARGET_QC_AFTER", "1") == "1"
        deps: list[dict[str, str]] = []
        for row in read_tsv(samples, ["sample_id", "read1", "read2"]):
            sample_ctx = Context(**{**ctx.__dict__, "sample": row["sample_id"]})
            deps.append(
                {
                    "stage": "module02",
                    "target": row["sample_id"],
                    "fingerprint": manifest_value("module02", sample_ctx),
                }
            )
        for row in read_tsv(pairs, ["sample_id", "control_id", "sgrna"]):
            pair_ctx = Context(**{**ctx.__dict__, "sample": row["sample_id"]})
            deps.append(
                {
                    "stage": "module05",
                    "target": row["sample_id"],
                    "fingerprint": manifest_value("module05", pair_ctx),
                }
            )
        if expect_a:
            deps.append(
                {
                    "stage": "phase-a-qc",
                    "target": "cohort",
                    "fingerprint": manifest_value("phase-a-qc", ctx),
                }
            )
        if expect_b:
            deps.append(
                {
                    "stage": "phase-b-qc",
                    "target": "cohort",
                    "fingerprint": manifest_value("phase-b-qc", ctx),
                }
            )
        basis.update(
            {
                "samples_sheet": content_descriptor(samples, label="samples-sheet"),
                "pairs_sheet": content_descriptor(pairs, label="pairs-sheet"),
                "upstream": deps,
            }
        )
    return basis


def manifest_value(stage: str, ctx: Context) -> str:
    path = manifest_path(stage, ctx)
    if not path.is_file():
        return "MISSING"
    try:
        return str(json.loads(path.read_text())["fingerprint"])
    except (KeyError, json.JSONDecodeError):
        return "INVALID"


def expected_fingerprint(stage: str, ctx: Context) -> tuple[str, dict[str, Any]]:
    basis = fingerprint_basis(stage, ctx)
    raw = json.dumps(basis, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()
    return sha256_bytes(raw), basis


def inspect(stage: str, ctx: Context) -> tuple[int, str]:
    direct_upstream = {
        "module02": "module01",
        "module04": "module03",
        "module05": "module04",
    }.get(stage)
    if direct_upstream is not None:
        upstream_code, upstream_detail = inspect(direct_upstream, ctx)
        if upstream_code != VALID:
            return STALE, f"upstream {direct_upstream} is not current: {upstream_detail}"
    manifest = manifest_path(stage, ctx)
    markers = marker_paths(stage, ctx)
    outputs = required_outputs(stage, ctx)
    missing = [str(path) for path in [*markers, *outputs] if not path.is_file()]
    if not manifest.is_file() or missing:
        details = []
        if not manifest.is_file():
            details.append(f"manifest={manifest}")
        if missing:
            details.append("missing=" + ",".join(missing))
        return INCOMPLETE, "; ".join(details)
    try:
        recorded = json.loads(manifest.read_text())
    except json.JSONDecodeError:
        return STALE, f"invalid manifest JSON: {manifest}"
    expected, _ = expected_fingerprint(stage, ctx)
    if recorded.get("schema_version") != SCHEMA_VERSION:
        return STALE, f"state schema changed: {manifest}"
    if recorded.get("fingerprint") != expected:
        return STALE, f"fingerprint changed: {manifest}"
    return VALID, f"current: {manifest}"


def record(stage: str, ctx: Context) -> None:
    outputs = required_outputs(stage, ctx)
    missing = [str(path) for path in outputs if not path.is_file()]
    if missing:
        raise FileNotFoundError("Cannot record stage; missing outputs: " + ", ".join(missing))
    fingerprint, basis = expected_fingerprint(stage, ctx)
    manifest = manifest_path(stage, ctx)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "schema_version": SCHEMA_VERSION,
        "stage": stage,
        "fingerprint": fingerprint,
        "inoseq_version": ctx.version,
        "recorded_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "basis": basis,
        "outputs": [str(path) for path in outputs],
    }
    temporary = manifest.with_suffix(manifest.suffix + ".tmp")
    temporary.write_text(json.dumps(data, ensure_ascii=False, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, manifest)
    for marker in marker_paths(stage, ctx):
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.touch()


def invalidate(stage: str, ctx: Context) -> None:
    chains = {
        "module01": ("module01", "module02"),
        "module02": ("module02",),
        "module03": ("module03", "module04", "module05"),
        "module04": ("module04", "module05"),
        "module05": ("module05",),
        "phase-a-qc": ("phase-a-qc",),
        "phase-b-qc": ("phase-b-qc",),
        "finalize": ("finalize",),
    }
    for item in chains[stage]:
        manifest_path(item, ctx).unlink(missing_ok=True)
        for marker in marker_paths(item, ctx):
            marker.unlink(missing_ok=True)
    if stage == "finalize":
        (ctx.output / "QC" / "full_workflow_status.tsv").unlink(missing_ok=True)


def add_context_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--stage", choices=STAGES, required=True)
    parser.add_argument("--project-dir", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--sample")
    parser.add_argument("--control")
    parser.add_argument("--sgrna")
    parser.add_argument("--read1")
    parser.add_argument("--read2")
    parser.add_argument("--samples-sheet")
    parser.add_argument("--pairs-sheet")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ("check", "record", "invalidate"):
        child = subparsers.add_parser(command)
        add_context_arguments(child)
        if command == "check":
            child.add_argument("--quiet", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    ctx = Context.from_args(args)
    try:
        if args.command == "check":
            code, detail = inspect(args.stage, ctx)
            if not args.quiet:
                label = "VALID" if code == VALID else "INCOMPLETE" if code == INCOMPLETE else "STALE"
                print(f"[{label}] {args.stage}: {detail}")
            return code
        if args.command == "record":
            record(args.stage, ctx)
            print(f"[STATE] recorded {args.stage}: {manifest_path(args.stage, ctx)}")
            return 0
        invalidate(args.stage, ctx)
        print(f"[STATE] invalidated {args.stage}")
        return 0
    except (FileNotFoundError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

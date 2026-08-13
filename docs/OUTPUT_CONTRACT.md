# Ino-seq v1.0.0 complete output contract

This document defines the complete output surface from modules 01–05 through
cohort reporting and final workflow verification. Historical basenames are
retained where downstream compatibility requires them.

## Output lifecycle

| Level | Completion evidence | Meaning |
|---|---|---|
| Module 01 | `<sample>/.inoseq/module01.json` plus `MODULE01_COMPLETE` | UMI consensus outputs passed state recording |
| Phase A sample | `<sample>/INOSEQ_COMPLETE` plus `.inoseq/module02.json` | modules 01–02 finished for one experimental or control sample |
| Modules 03–04 | `<sample>/postprocess/.inoseq/module03.json` / `module04.json` | the corresponding paired stage passed state recording |
| Phase B pair | `<sample>/postprocess/INOSEQ_POSTPROCESS_COMPLETE` plus `.inoseq/module05.json` | modules 03–05 finished for one experimental/control pair |
| Full experimental sample | `<sample>/INOSEQ_FULL_COMPLETE` | Phase A and Phase B products exist for that experimental sample |
| Cohort summaries | `QC/.inoseq/phase-a-qc.json` and `phase-b-qc.json` | enabled cohort reports were generated from current upstream states |
| Whole run | `INOSEQ_WORKFLOW_COMPLETE`, `QC/full_workflow_status.tsv` and `.inoseq/finalize.json` | every enabled branch and cohort summary passed final verification |

An intermediate file or marker alone is not a reuse criterion. Automatic resume
requires the marker, every required output, and a matching state fingerprint.

## Resume-state contract

Each `.inoseq/*.json` record contains the stage name, schema, Ino-seq version,
successful recording time, required output paths, and a SHA256 fingerprint of
the stage's reproducibility basis. Depending on the stage, that basis includes:

- FASTQ and reference path/size/mtime metadata;
- relevant configuration values and sgRNA;
- upstream stage fingerprints;
- content hashes of workflow code and small sample/pair sheets.

Large sequencing/reference files are not content-hashed on every submission.
A changed path, size or modification time invalidates their dependent stage.
State JSON and marker files are orchestration/provenance metadata; they do not
change analytical tables, thresholds or biological interpretation.

## Directory map

```text
<OUTPUT_DIR>/
├── <sample>/
│   ├── .inoseq/                     Module 01–02 resume state
│   ├── QC/                           Phase A read-processing files
│   ├── alignment/                    Phase A alignment and signature reads
│   ├── postprocess/                  Phase B paired detection; experiments only
│   │   ├── .inoseq/                  Module 03–05 resume state
│   │   ├── filt-before/
│   │   ├── filt-after/
│   │   ├── final_merge_distance/
│   │   ├── offtarget/
│   │   ├── summary/
│   │   ├── run_parameters.tsv
│   │   └── INOSEQ_POSTPROCESS_COMPLETE
│   ├── INOSEQ_COMPLETE
│   └── INOSEQ_FULL_COMPLETE          experiments only
├── QC/                               Phase C cohort reports and run status
│   └── .inoseq/                      Cohort QC resume state
├── .inoseq/finalize.json             Run-level resume state
└── INOSEQ_WORKFLOW_COMPLETE
```

Control samples normally contain Phase A outputs but no `postprocess/`
directory unless they are also declared as an experimental sample in a pair.

## Phase A — modules 01–02

### Read processing and UMI products

| Filename | Role |
|---|---|
| `<sample>/QC/<sample>_fastp.json` | machine-readable raw/clean read QC |
| `<sample>/QC/<sample>_fastp.html` | fastp visual report |
| `<sample>/QC/<sample>_cutadapt.json` | UMI-adapter extraction statistics |
| `<sample>/alignment/<sample>.mapped.bam` | coordinate-sorted initial alignment |
| `<sample>/alignment/<sample>.mapped.flagstat.json` | initial mapping summary |
| `<sample>/alignment/<sample>.mapped.primary_mapq<q>.count.txt` | high-confidence primary mapped-read count at configured MAPQ |
| `<sample>/alignment/<sample>.tag-family-sizes.txt` | fgbio UMI family-size histogram |
| `<sample>/alignment/<sample>.cons.unmapped.bam` | filtered molecular consensus before remapping |
| `<sample>/alignment/<sample>.umi_dedup.bam` | final remapped UMI-consensus BAM |
| `<sample>/alignment/<sample>.flagstat.json` | final consensus alignment summary |
| `<sample>/alignment/<sample>.idxstat.txt` | final per-contig alignment counts |

Every coordinate-sorted BAM that is used downstream has a BAM index.
Cleaned FASTQ and other implementation intermediates under `<sample>/QC/` are
runtime products; they are not biological result tables.

### ABE signature-read products

| Filename | Role |
|---|---|
| `<sample>/alignment/<sample>.end` | read-level ABE signature sites and cleavage direction |
| `<sample>/alignment/<sample>_end.bam` | alignments passing the signature-read definition |
| `<sample>/alignment/<sample>_end.bam.bai` | index for the signature-read BAM |

The `.end` table is the Phase A input for experimental samples in Phase B. A
matched control contributes `_end.bam` coverage to the paired background
comparison.

## Phase B — module 03 background filtering

All paths below are relative to `<sample>/postprocess/`.

| Filename | Description |
|---|---|
| `filt-before/<sample>_merge.tsv` | aggregated cleavage sites with header |
| `filt-before/<sample>.bed` | the same seven aggregation fields without header |
| `filt-before/<sample>_coverage.txt` | experimental BAM overlap count appended |
| `filt-before/<sample>_ctr_<control>_coverage.txt` | matched-control overlap count |
| `filt-before/<sample>.txt` | fold change, raw P value and BH-FDR |
| `filt-after/<sample>_filted.txt` | sites passing `fold_change >=1.5` and raw `P<0.05` |
| `filt-after/<sample>_query_name.txt` | query names retained by the site filter |
| `filt-after/<sample>_filted.bam` | retained signature-read alignments |

BH-FDR is a reported field and is not a filter in v1.0.0.

## Phase B — module 04 candidate intervals

| Filename | Description |
|---|---|
| `final_merge_distance/<sample>_final_merged_distance.txt` | merged and standardized candidate intervals |
| `final_merge_distance/<sample>_filted.bam` | alignments assigned to retained candidate intervals |

The historical misspelling `filted` is preserved as part of the compatibility
interface.

## Phase B — module 05 classification

| Filename | Description |
|---|---|
| `offtarget/<sample>_align.txt` | sgRNA fuzzy alignment results |
| `offtarget/<sample>_dependent_mark.txt` | `onTarget`, `dependent` or `independent` label |
| `offtarget/<sample>_dependent_out_of_spacer.txt` | neighboring-site and protospacer annotations |
| `offtarget/<sample>_dependent_target.txt` | final strand-aware site-level result table |

Classification order is frozen:

1. `onTarget`: zero-mismatch substitution-only match to the supplied sgRNA;
2. `dependent`: an accepted sgRNA alignment that is not `onTarget`;
3. `independent`: no accepted sgRNA alignment.

The `onTarget` definition is sequence-based. It does not require a separately
provided expected coordinate, so multiple perfect genomic matches can yield
multiple `onTarget` candidates.

## Per-sample summaries and parameters

| Filename | Description |
|---|---|
| `summary/<sample>_offtarget_summary.tsv` | site/read counts by primary dependency class |
| `summary/<sample>_strand_summary.tsv` | strand and protospacer-detail summary |
| `summary/dependent_target_analysis.xlsx` | two-sheet legacy-compatible sample workbook |
| `run_parameters.tsv` | version, IDs, sgRNA and exact frozen parameters used for the pair |

`run_parameters.tsv` is the provenance record for interpreting a sample result.

## Phase C — cohort outputs

All files are under `<OUTPUT_DIR>/QC/`.

| Filename | Description |
|---|---|
| `inoseq_qc_summary.tsv` | complete sample-level sequencing, mapping, UMI and signature-read QC |
| `inoseq_qc_summary.csv` | CSV representation of the same QC table |
| `inoseq_offtarget_summary.tsv` | concatenated sample-level off-target summary |
| `inoseq_strand_summary.tsv` | concatenated strand/protospacer summary |
| `dependent_target_analysis.xlsx` | cohort workbook with basic and detailed sheets |
| `full_workflow_status.tsv` | final status, version, time, sample/pair counts and core thresholds |

`full_workflow_status.tsv` describes the run as a whole. It does not replace
the per-sample `run_parameters.tsv` files.

## Scheduler provenance

The following files are written under `<LOG_DIR>` by the complete submission
path:

| Filename | Description |
|---|---|
| `inoseq_last_full_workflow_job_id.txt` | finalizer job ID |
| `inoseq_last_full_workflow_jobs.tsv` | job ID, stage, target and `afterok` dependency graph |
| `inoseq_last_sample_job_ids.txt` | Phase A sample job IDs |
| `inoseq_last_postprocess_job_ids.txt` | Phase B pair job IDs |
| `inoseq_last_qc_job_id.txt` | Phase A cohort QC job ID when enabled |
| `inoseq_last_offtarget_qc_job_id.txt` | off-target cohort summary job ID when enabled |

The graph may contain `job_id=SKIP` rows with `dependency=current`; these are
verified stages reused by the most recent resume plan. Only newly submitted
jobs receive numeric Slurm IDs.

Slurm `.out` and `.err` files are execution logs, not scientific outputs.

## Coordinate and compatibility notes

- The historical `.end` table retains its 1-based reported `location` field.
- BAM and FASTA access uses 0-based half-open coordinates internally.
- The historical `clevage_reads` spelling and the `*_filted.txt` index column
  are retained for downstream compatibility.
- See [ALGORITHM_CONTRACT.md](ALGORITHM_CONTRACT.md) for frozen equations,
  boundaries and source-script mapping.

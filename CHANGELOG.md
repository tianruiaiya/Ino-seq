# Changelog

## 1.0.0 — 2026-08-10

- Standardized the complete Ino-seq core workflow from paired-end FASTQ to
  sample- and cohort-level off-target summaries.
- Preserved the validated modules 01–02 for UMI consensus and ABE
  signature-read detection.
- Added paired background comparison, candidate-site merging, sgRNA alignment,
  dependency annotation, protospacer annotation and strand classification as
  modules 03–05.
- Frozen background filtering at `fold_change >= 1.5` and raw
  `p_value < 0.05`; BH-FDR remains a reported value and is not a filter.
- Added one-command Slurm submission with explicit `afterok` dependencies from
  raw FASTQ through final workflow verification.
- Added configuration validation, run-parameter records, completion markers,
  deterministic output directories and cohort TSV/XLSX summaries.
- Added unit, synthetic integration, source-regression and Slurm dependency
  tests.
- Added the root-level `./inoseq` command as the single public interface for
  initialization, validation, planning, submission, status and testing.
- Unified documentation around Phase A (modules 01–02), Phase B (modules
  03–05) and Phase C (cohort reporting), and expanded the output contract to
  cover the complete workflow.
- Added module-level resume state with completion/output/fingerprint checks,
  automatic skip of current work, dynamic Slurm dependency rebuilding,
  explicit `--from-stage` boundaries, detailed stage status, and controlled
  adoption of complete pre-resume outputs.

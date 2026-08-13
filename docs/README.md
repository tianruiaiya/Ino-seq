# Ino-seq v1.0.0 documentation

This directory separates operating instructions from versioned scientific
contracts. Start from the repository [README](../README.md); use this page when
you need to choose the authoritative reference for a specific question.

## Reading paths

### Running an analysis

1. Read [QUICKSTART_CN.md](QUICKSTART_CN.md).
2. Use the root-level `./inoseq` command.
3. Consult [OUTPUT_CONTRACT.md](OUTPUT_CONTRACT.md) after completion.

### Reviewing or reproducing the method

1. Read [ALGORITHM_CONTRACT.md](ALGORITHM_CONTRACT.md).
2. Confirm parameters in `config/inoseq.env` and each sample's
   `postprocess/run_parameters.tsv`.
3. Use [OUTPUT_CONTRACT.md](OUTPUT_CONTRACT.md) for table and path definitions.

## Canonical terminology

| Public term | Scope | Historical/internal name retained for compatibility |
|---|---|---|
| Phase A — sample processing | modules 01–02 for every sample | Step01, `run_inoseq.sh`, `INOSEQ_COMPLETE` |
| Phase B — paired detection | modules 03–05 for each experimental/control pair | postprocess, `run_postprocess.sh` |
| Phase C — reporting | cohort QC and final verification | QC jobs, finalizer |

The phases are scheduling units within one complete workflow. They must not be
presented as separate alternative pipelines.

## Sources of truth

| Question | Authoritative location |
|---|---|
| Which command should a user run? | root `./inoseq` help and `QUICKSTART_CN.md` |
| How are completed stages reused? | `QUICKSTART_CN.md` sections 6 and 10, plus `.inoseq/*.json` state records |
| What do configuration variables mean? | annotated `config/inoseq.env.example` |
| What calculation and threshold are frozen? | `ALGORITHM_CONTRACT.md` |
| What does an output file mean? | `OUTPUT_CONTRACT.md` |
| Did the whole run finish? | `QC/full_workflow_status.tsv` plus `INOSEQ_WORKFLOW_COMPLETE` |
| Which exact parameters produced a sample result? | `<sample>/postprocess/run_parameters.tsv` |

Low-level scripts under `workflow/` remain callable for testing and controlled
partial reruns, but they are implementation interfaces. A standard user should
not need to assemble those scripts manually.

# Ino-seq v1.0.0 algorithm compatibility contract

Version `1.0.0` restructures the original Ino-seq research scripts into a
reproducible workflow. Valid-input analytical decisions are preserved, with
the released background-filter threshold frozen at `p_value < 0.05`.

The execution model has three phases within one workflow: Phase A contains
modules 01–02, Phase B contains modules 03–05, and Phase C performs cohort
reporting and final verification. The phase boundary changes scheduling only;
it does not introduce an alternative calculation path.

## Source mapping

| Original script | Refactored module | Frozen behavior |
|---|---|---|
| `get_result.py` | `03_aggregate_sites.py` | group by chromosome/location; first strand; +/-15 bp |
| `ABE_filt_merge.sh` coverage blocks | `03_count_coverage.py` | count BAM alignment records overlapping each experimental interval |
| `nonormalized_all_calculate_FDR.R` | `03_compare_background.py` | raw-count fold change, zero-control pseudocount, Poisson-tail cases, BH-FDR |
| `filt.py` | `03_filter_background.py` | fold change and raw P-value filter; query-name BAM extraction |
| `final_merged_distance.py` | `04_candidate_intervals.py` | <=30 bp merge, 30 bp extension, total cleavage reads >=3 |
| `identify.py` | `05_align_sgrna.py` | +/-25 bp window, both strands, substitutions/bulges, max score 8 |
| `dependent_mark.py` | `05_mark_dependency.py` | onTarget/dependent/independent decision order |
| `detect_out_of_spacer.py` | `05_annotate_spacer.py` | +/-100 bp neighbor scan and protospacer-boundary annotation |
| `detect_target.py` | `05_classify_strands.py` | cleavage-strand inversion and target/non-target read assignment |
| `get_stat.py` | `05_summarize.py` | basic and strand/protospacer summary definitions |

## Background statistics

Let `E` and `C` be experimental and control overlap counts.

- Fold change is `E / C`, except `C=0` uses `0.5` only in the denominator.
- If `E=C=0`, P is 1.
- If `C=0<E`, P is the upper Poisson tail with lambda 0.1.
- If `E=0<C`, P is twice the lower Poisson probability at 0 with lambda C.
- Otherwise, the historical two-sided Poisson-tail approximation around
  `(E+C)/2` is retained and capped at 1.
- BH-FDR is calculated across sites and reported.
- Candidate filtering uses `fold_change >= 1.5` and raw `P < 0.05`.

This method is an unnormalized paired count comparison. It must not be
described as Fisher's exact test in a manuscript or methods section.

## Compatibility safeguards

- Historical output basenames are retained inside per-sample stage folders.
- The historical misspelling `clevage_reads` remains part of the tabular API.
- The `*_filted.txt` index column is retained for backward compatibility.
- The first maximum-support position is retained when cleavage-read counts tie.
- The historical bulge realignment parameter behavior is frozen and regression
  tested; it is not silently reinterpreted in this engineering release.
- Valid inputs produce the same downstream tables in source-regression tests.

## Safer failure behavior

The refactor intentionally fails early for malformed tables, mismatched
experimental/control intervals, unreadable inputs and BAM fetch errors. These
checks do not change valid-input calculations; they prevent technical failures
from being reported as biological zeroes.

## Coordinate convention

The historical `.end` table retains its 1-based reported `location` field,
while BAM/FASTA interval access uses 0-based half-open coordinates internally.
The refactored modules preserve the validated conversion behavior. For a new
reference build or assay design, users should still confirm at least one known
on-target locus against the reference sequence before interpreting novel sites.

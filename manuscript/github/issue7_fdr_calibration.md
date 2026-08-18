Permutation-null calibration of `ema switch diff --strategy fisher` shows it is **decisively anti-conservative** — its q-values cannot support publishable claims as currently defaulted.

## Protocol

TRUE run: the cohort TREG per-celltype h5ad (3,790 cells × 115,450 PAS, `--cluster-key stage`, unified pasbed, GTF), fisher, BH `--fdr 0.05` → 56 tests across 6 stage pairs, **40 q<0.05 hits**.
NULL: 20 runs with `obs['stage']` shuffled across cells (label counts preserved, numpy seed=k), **full pipeline rerun per permutation** including marker selection, so selection effects are captured.

## Result

| statistic | expectation under null | observed |
|---|---|---|
| null p-values < 0.05 | ~5% | **38.6%** (of 2,640) |
| null tests with q < 0.05 | ~0 | **31.9%** (841/2,640) |
| null runs with ≥1 q<0.05 hit | ~5% of runs | **20/20 (100%)** |
| mean hits per null run | ~0 | **42** (range 4–67) — *the TRUE run has 40* |
| min null p | ~1/2,640 | **1e-20** |

A null run reports as many "significant switches" as the real labels do. Reproducible: `PeakATail_wd/results/fdr_calibration/` + `scripts/manuscript_figures/fdr_calibration.py` (figure in `manuscript/figures/fdr_calibration.*`).

## Causes (verified in source)

1. **Read-count pseudoreplication**: FisherStrategy builds 2×2 tables from read counts; per-cell reads are correlated. ema's own docstring flags fisher as anti-conservative — this quantifies *how* anti-conservative.
2. **Marker pre-selection on the tested labels** (top-200/cluster) — classic double dipping; the null captures it because markers were re-selected per permutation.

## Recommendations

- **Expose D4** (`count_mode="cells"`, already landed per HANDOFF) **as a CLI flag** — it is currently unreachable from `ema switch diff` (verified `--help`), so no pipeline can use the fix.
- Re-run this calibration for `nb_pairwise` (in progress on our side) and for fisher+cells-mode once exposed; whichever calibrates becomes the manuscript's test.
- **Affects #67**: the sweep's SWITCH_CELLTYPE fix should not rerun with `--strategy fisher` defaults for publication results — either add nb, or wait for D4 exposure, or treat fisher output as a ranking screen only (say so in Methods).

@TRextabat

# Issue: `switch diff` default marker pre-selection makes every test anti-conservative (label double-dip); also changes the Fisher denominator

**Labels:** bug, statistics, switch-test. **Blocks:** reliable cell-type switch calls for the paper.

## Summary

On a correctly keyed count matrix (post-#92 clip-seeded caller, testis mouse 1, 1,230 cells with
GEX-derived SPC/RS/ES labels, 20 label permutations, full `switch diff` pipeline per run) the
shipped defaults are not FDR-controlled under the label-permutation null:

| test | marker pre-selection | null p<0.05 | null runs (of 20) with ≥1 q<0.05 hit |
|---|---|---|---|
| fisher `--count-mode reads` | top-200 (default) | 20.3% | 20 |
| fisher `--count-mode cells` | top-200 (default) | 13.0% | 19 |
| nb_pairwise | top-200 (default) | 24.7% | 19 |
| fisher `--count-mode cells` | **`--marker-top-n 0`** | **3.0%** | **0** |
| fisher `--count-mode reads` | `--marker-top-n 0` | 9.6% | 20 |
| nb_pairwise | `--marker-top-n 0` | 5.1% | 20 (dispersion-floor tail) |

Expected for a calibrated test: ~5% null p<0.05 and a q<0.05 hit in roughly 1–3 of 20 runs.

## Two causes, both in the marker path

1. **Label double-dip.** `switch diff` selects the top-200 wilcoxon markers per cluster on `.X`
   *using the labels that are then tested*, and tests only those. Restricting the no-marker
   cells-mode Fisher p-values to the PAS that marker mode selected on *permuted* labels gives
   17.4% p<0.05 — selection alone reproduces the inflation.
2. **The Fisher denominator silently changes.** With pre-selection on, the within-gene "gene
   total" is computed from the marker-restricted matrix (`ema/switch_test/strategies/fisher.py`
   ~L176-179: `cm1[pas_in_gene]` of the restricted count matrix), i.e. only that gene's *marker*
   PAS. Joining marker-on vs marker-off runs on identical permutation/PAS/pair, only 629/6,453
   p-values agree (example: `n_reads_gene` 1,746 vs 5,289 for the same PAS).

Separately: reads-mode Fisher is inflated by pseudoreplication (UMIs within a cell are not
independent) regardless of pre-selection, and nb_pairwise's per-PAS dispersion hits the 1e-4
floor in 8.5% of null tests but produces 67% of its false q<0.05 hits.

## Proposed change

- [ ] Default `--marker-top-n 0` (test all PAS that pass a **label-independent** filter, e.g.
      min expressing cells per group), or replace wilcoxon pre-selection with such a filter.
- [ ] Compute the Fisher gene denominator from the **full** count matrix in every mode; the
      restricted matrix must only restrict *which* PAS are tested, never the gene totals.
- [ ] Make `--count-mode cells` the default for fisher (reads mode is pseudoreplicated).
- [ ] nb_pairwise: raise/regularise the dispersion floor or shrink toward a trend; until then
      document that it needs permutation-calibrated q.
- [ ] Docs: the per-pair TSV column list on the CLI page is stale (`statistic` is gone; current
      columns end `odds_ratio delta_proportion log2fc`), and note that `delta_proportion =
      prop(c1) - prop(c2)` while `log2fc = log2(prop(c2)/prop(c1))` — opposite sign conventions.
- [ ] Run manifests should record the tool commit hash (the Stage-1b manifest records none).

## Reproduce

`results/fdr_calibration_v2/` on biolab: `build_input.py` → `run_all_arms.sh` / `run_diag_arms.sh`
(resumable, `DONE.ok` markers); harvest + figure: `scripts/manuscript_figures/fdr_calibration_v2.py`.
Write-up: `manuscript/14_switch_calibration_v2.md`. Figure: `manuscript/figures/fdr_calibration_v2.png`.

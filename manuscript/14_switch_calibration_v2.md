# Switch-test calibration v2 — which `ema switch diff` test can ship as "reliable"

**Date:** 2026-08-21. **Status:** analysis + adversarial verification complete → **SOUND**
(verifier reproduced every headline number from the raw per-permutation TSVs with an
independent parser; two wording/caveat corrections folded in below).
**Serves:** PI directive in [13_reliability_positioning.md](13_reliability_positioning.md) §2 —
"cell type specific switches should be reliable."

Supersedes the v1 calibration (`manuscript/figures/fdr_calibration.png`), which ran on the
bug-0a mis-keyed matrix. Everything here is on a correctly keyed count matrix.

## Input

- Stage-1b acceptance run, testis mouse 1, clip-seeded caller, `07_clustering/default/clusters.h5ad`
  (post-0a fix; `layers['counts']` integral, 32.7 M UMI over the subset; stage-marker genes show the
  expected pattern — Prm1 CP10k SPC 20 / RS 54 / ES 1,064).
- Cells: the 1,294 STARsolo cells with **GEX-derived** stage labels (marker-panel argmax on the
  STARsolo gene matrix — orthogonal to the PAS matrix under test). SPG dropped (64 cells, fragile per
  [11_verifier_corrections.md](11_verifier_corrections.md)) → **1,230 cells** (SPC 370 / RS 504 /
  ES 356) × 76,711 PAS, 3 stage pairs per run. Agreement with PAS-derived labels 98.8%.
- Null: 20 stage-label permutations (seeds 1–20), **identical permutations shared by all arms**; full
  `switch diff` pipeline per run; 6 arms × 21 runs, 0 failures.

## Arms and verdicts

Pre-registered rule (anti-conservative checks only; KS reported, not gated): *calibrated iff*
null p<0.05 ≤ 7% AND null p<0.01 ≤ 1.5% AND null q<0.05 ≤ 5% AND ≤ 25% of null runs carry any
q<0.05 hit.

| arm | marker pre-selection | null p<0.05 | null q<0.05 | null runs with ≥1 q<0.05 hit | mean hits / null run | TRUE hits | verdict |
|---|---|---|---|---|---|---|---|
| A fisher, `--count-mode reads` (legacy) | top-200 wilcoxon (default) | 20.3% | 7.66% | 20/20 | 24.7 | 307 | anti-conservative |
| B fisher, `--count-mode cells` (#86) | top-200 wilcoxon (default) | 13.0% | 0.96% | 19/20 | 3.1 (max 11) | 668 | anti-conservative |
| C nb_pairwise | top-200 wilcoxon (default) | 24.7% | 1.27% | 19/20 | 22.9 (max 141) | 1,754 | anti-conservative |
| A0 fisher reads | **off** | 9.6% | 0.69% | 20/20 | 1,487 | 49,923 | anti-conservative (pseudoreplication) |
| **B0 fisher cells** | **off** | **3.0%** | **0.00%** | **0/20** (0/60 BH families) | **0** | 60,332 | **valid FDR control — conservative** |
| C0 nb_pairwise | **off** | 5.1% | 0.13% | 20/20 | 139 | 43,127 | bulk OK, **tail inflated** |

Per-pair breakdown: B0 is 2.94–3.13% null p<0.05 with 0/20 families hit in *every* stage pair;
B is 12.6–13.7% in all three; no pair drives the result. Stratifying B0's null p-values by a
label-independent expression measure (6 strata, <50 to ≥600 expressing cells) gives 2.6–4.0%
everywhere — the test itself is never anti-conservative at any expression level.

## Mechanism (verified)

1. **Marker pre-selection is a label double-dip.** The shipped default selects the top-200
   wilcoxon markers per cluster *on the same labels* being tested, then tests them. Restricting
   B0's own null p-values to the PAS that B selected on the permuted labels gives 17.4% p<0.05 —
   selection alone reproduces the anti-conservativeness. Switching it off takes cells-mode Fisher
   from 13.0% → 3.0%.
2. **Verifier addition — marker mode also changes the test's denominator.** With pre-selection
   on, the within-gene Fisher "gene total" is computed from the marker-restricted matrix
   (`fisher.py:176-179`), i.e. only the marker PAS of that gene. Joining B and B0 on the same
   permutation/PAS/pair, only 629/6,453 p-values are identical (e.g. one PAS: `n_reads_gene`
   1,746 vs 5,289). So "marker mode" is not just a subset of the same tests — it is a different
   test.
3. **reads-mode Fisher is inflated by the test itself** (pseudoreplication of UMIs within cells):
   9.6% null p<0.05 and ~1,500 false q<0.05 hits per null run even with pre-selection off.
4. **nb_pairwise is tail-inflated by the dispersion floor**: per-PAS plug-in dispersion hits the
   1e-4 floor in 8.5% of null tests but accounts for 67% of null q<0.05 hits (C0); 13% in C.
   Min null p 4e-182.

## Verifier corrections (already applied)

- B0 is **conservative, not exactly calibrated**: an exactly calibrated test with 3 BH families
  per run would show a hit in ~1–3 of 20 null runs; 0/20 plus a 3.0% null p<0.05 rate (24.8% of
  null p exactly 1, discrete Fisher) is conservative. Figure title now reads "valid FDR control
  … (conservative)"; the `conservative()` rule in the figure script was widened accordingly.
- Denominator caveat (point 2 above) added here and in the figure index.
- Provenance of "correctly keyed": supported by run timestamp post-dating the 0a fix (51cbe67,
  ancestor of tool HEAD 6697b0d), the 85,993/85,993 keying check in
  [12_stage2_gate.md](12_stage2_gate.md), and the biological spot check — **not** by a stored
  commit hash (the run manifest records none). Final-run manifests must record the tool commit.

## Permutation-calibrated q (the alternative if pre-selection is kept)

Empirical p = rank within the pooled same-pair null p-values, BH per pair. TRUE hits surviving
q_perm<0.05: A 186/307 (85 with |Δprop| ≥ 0.1), B 640/668 (572 with floor), C 1,733/1,754
(1,588 with |log2FC| ≥ 1). For B0 the permutation null is conservative so q_perm *adds* hits
(66,630 vs 60,332 nominal; 46,230 with the |Δprop| ≥ 0.1 floor). Resolution is 1/(N_null+1)
(~0.00046 for the marker arms).

## What ships (decision)

1. **Default reliable switch test = Fisher, `--count-mode cells`, marker pre-selection OFF
   (`--marker-top-n 0`)**, nominal BH q. It controls FDR (conservatively) in every stratum and
   stage pair tested, and it needs no permutation machinery.
2. Whenever marker pre-selection or `nb_pairwise` is used, **permutation-calibrated q is
   required** — and the paper says so.
3. On top of (1): the replication filter (≥2 samples, same direction) and the effect-size floor
   (|Δprop| or |ΔPDUI| ≥ 0.1) from [13](13_reliability_positioning.md) §2, implemented in
   `scripts/reliability/replication_filter.py`.
4. Tool change requested from Amir (issue draft: [github/issue9_switch_marker_preselection.md](github/issue9_switch_marker_preselection.md)):
   default `--marker-top-n 0`, or a **label-independent** pre-filter (e.g. minimum expressing
   cells), and compute the Fisher gene denominator from the full matrix regardless of mode.

## Caveats that travel with the figure

- Single mouse, single tissue, very large true stage effects: TRUE hit counts are upper bounds on
  detectability, not precision estimates. Null-run hit *counts* are not comparable to TRUE counts
  (permuted labels change the marker set and the ≥10-cells filter); the calibrated quantities are
  the per-test null rates.
- The label-permutation null tests exchangeability under the global null; it does not model
  within-stage heterogeneity or batch. Cross-dataset confirmation (Laughney TREG, per-patient)
  is outstanding and is part of the final run.
- Two-sided KS-vs-uniform rejects for every arm (10^4–10^6 tests, discrete Fisher mass at p=1);
  reported, not gated.
- "Calibrated" is our operational rule above, pre-stated in the figure script, not a standard.

## Files

- Figure: `manuscript/figures/fdr_calibration_v2.{png,pdf}`; script
  `scripts/manuscript_figures/fdr_calibration_v2.py` (harvest → metrics → perm-q → figure → report).
- Tables: `results/figures/manuscript/fdr_calibration_v2{,_stats,_per_pair,_hist,_null_pvalues,_true_permcal}.tsv`.
- Machine-readable verdict: `results/fdr_calibration_v2/report.json`. Runs (resumable, `DONE.ok`
  markers): `results/fdr_calibration_v2/<arm>/{true,null/perm_01..20}/` (7.0 GB; the 20 permuted
  h5ads = 4 GB are deletable).

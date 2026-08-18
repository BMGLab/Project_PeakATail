The atlas benchmark as currently computed (`ema/benchmark/metrics.py`, surfaced as `benchmark_vs_polyasite_v3.json`) reports numbers that will not survive peer review. A full null-control analysis (BioLab: `PeakATail_wd/scripts/manuscript_figures/null_control.py`, adversarially verified) shows:

- The reference used has **18,432,135 entries** (mean 9 bp clusters) — one every 168 bp of genome. At a 100 bp cutoff, **33% of the genome** lies within reach of a strand-matched entry. This is not the curated ~570k-cluster PolyASite atlas.
- `metrics.py` scores the **whole predicted interval** (mean peak width 821 bp, max 10.2 kb) via `predicted.window(reference)` with **no strand matching** — a peak is credited if it merely contains any entry on either strand.
- Consequence: width/chrom/count-matched **shuffled peaks score 0.62–0.75 precision** under the identical procedure. The reported 0.9986 sits on that floor, not on 0.
- The reported "recall" (0.0326 @100 bp) is coverage of a genome-wide reference including unexpressed genes — not sensitivity — and the derived F1 is meaningless.

The signal is real once scored properly — that's the good news:

| test | real PAS | shuffled null |
|---|---|---|
| point-mode (peak 3′ base), strand-matched, 100 bp | **0.994** | 0.35–0.46 |
| terminal-exon fraction of matches | **0.310** (10×) | ~0.03 (background) |
| high-confidence subset (atlas avg TPM ≥ 1) | **0.37** | 0.002 |

## Proposed change (PR to follow)

1. Point-mode scoring: strand-aware 3′-most base of each peak (BED end exclusive: `+` → end−1, `−` → start).
2. Strand-matched matching.
3. Optional shuffled-null baseline (N seeds, width/chrom-preserving, optional gene-body constraint) reported alongside.
4. Rename the recall-like quantity to `reference_coverage` in the new output block; support a restricted-reference BED for cohort-relevant recall. Legacy JSON keys stay unchanged for back-compat.

@TRextabat — flagging before the PR lands so the design isn't a surprise; happy to adjust the output schema to whatever the sweep's downstream expects.

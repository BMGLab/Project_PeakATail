# Positioning decision — reliability first (PI, 2026-08-21)

**Decision (Ebru, PI):** drop the "PAS-based clustering beats/extends gene-expression clustering" claim
(ablation showed it is expression signal re-encoded; kept only as an optional utility). The tool's
claims are: (1) trustworthy PAS detection with **low false positives**, (2) **reliable cell-type-specific
APA switches**, (3) 3'UTR shortening/lengthening, (4) **trustworthy de novo (atlas-novel) PAS**.

This document pre-registers what those words mean operationally, BEFORE the final Stage-2 run exists.
Numbers that informed the choices come from the now-stale Stage-2 runs and are disclosed as such; the
choices are made on product grounds (precision over recall), not on the final result.

## 1. Pre-registered default output (precision-first)

`ema run --peak-strategy clip_seeded --ip-filter --ip-filter-mode filter` and the default output =
**clip-supported tier only (coverage-only tier written to a separate `low_confidence` file), internal-priming
filter ON, ≥2 distinct molecules of clip support**. Rationale: every component independently removes a
measured false-positive class — tier-2 (P ≈ 0.06 vs null 0.02), IP-flagged sites (30% of old calls per
long reads), single-molecule clusters (23.7% of ≥2-read sites were one duplicated molecule).
Expected cost, from stale numbers: recall ≈ 0.19 (vs polyApipe 0.20), precision ≈ 0.5–0.6.
**Gate for the default output (added to, not replacing, the original): P@100 ≥ 0.50 on pbmc_10k_v3
and P ≥ 0.50 on both testis mice; F1 reported but not gated.** The original two-sided gate stays
reported for the ≥1-molecule output so the comparison to the earlier numbers is never hidden.

## 2. Reliable cell-type-specific switches — Stage 3 reliability program

| requirement | measured state | action |
|---|---|---|
| FDR-calibrated test | fisher 38.6% / nb 21.3% null p<0.05 (NOT calibrated) | re-calibrate on correctly keyed matrices with `--count-mode cells` (#86) and nb; if still off, ship **permutation-calibrated** q-values (label-shuffle null, 20+ perms) as the default |
| replication | none required today | **replication filter**: a switch is reported only if called (same direction) in ≥2 independent samples (Laughney 17 samples; testis 2 mice) — the single strongest FP control we have |
| strand/counts correctness | fixed (#82, Stage 0) | re-validate on spermatogenesis SPC>RS>ES with the new caller |
| effect-size floor | none | require |ΔPDUI| ≥ 0.1 (or Δproportion) in addition to q — report both |

## 3. Trusted de novo PAS — pre-registered definition

A PAS is "trusted novel" iff: clip-supported with ≥2 distinct molecules AND not IP-flagged AND
AATAAA/ATTAAA (or one of the 12 canonical hexamers) within −40..−5 nt of the cleavage site AND ≥ 100 bp
from any PolyASite 2.0 / PolyA_DB site. Validation (atlas-independent): fraction of trusted-novel PBMC PAS
with a Kinnex long-read 3' end within 25 bp (we have 104 M poly(A)-verified molecules) — target ≥ 70%,
reported against a shuffled-position null.

## 4. What changes in the manuscript
- Title/abstract lead with reliable APA analysis and trusted PAS calling, not clustering.
- Benchmark figure reports the precision-first default as PeakATail's primary arm, the ≥1-molecule
  output as the sensitivity arm, both tiers, both gates.
- Clustering moves to a supplementary utility section, with the ablation stated honestly.

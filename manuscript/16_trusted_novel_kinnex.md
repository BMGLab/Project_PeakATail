# Trusted de novo PAS — pre-registered definition vs Kinnex long-read truth (verified FIXED, 2026-08-21)

**Outcome: the pre-registered definition does NOT meet its target.** It was set in
[13](13_reliability_positioning.md) §3 before any number existed: trusted-novel = clip-supported ∧ ≥2
molecules ∧ not internal-priming ∧ canonical hexamer in −40..−5 ∧ ≥100 bp from any atlas site; target
≥70% within 25 bp (strand-matched) of a poly(A)-verified Kinnex long-read 3′ end. Input = the final
PBMC precision default (IP arm, code `4efeb125`, 44,413 sites). Verifier reproduced every count
exactly; one README fix (IP-flag column lives in `03_gtf_annotation/default/annotatedpas.bed:10`, not
the 9-column top-level file) and one robustness table were added.

## Funnel (verified)

44,394 on primary contigs → clip-supported 44,394 → ≥2 molecules 44,394 → not-IP 44,394 (the run was
already IP-filtered; re-deriving IP from the genome would flag 389–736, ≤3% change) → hexamer (any
of 12) **34,192** (77.0%) → ≥100 bp from PolyASite 2.0 **6,628 trusted-novel** (14.9% of input) →
AATAAA/ATTAAA **4,037 strong** (9.1%). 3,350 genes. PolyA_DB could not be used (hg19 only).

## Long-read concordance (≤25 bp, same strand; Kinnex x3p truth; Wilson 95% CI)

| set | n | ≥5 UMI | ≥20 | ≥100 | ≥500 | target 0.70 |
|---|---:|---:|---:|---:|---:|---|
| **trusted-novel (pre-registered)** | 6,628 | **0.521** [0.509–0.533] | 0.261 | 0.073 | 0.013 | **not met** |
| trusted-novel, strong hexamer | 4,037 | 0.525 | 0.256 | 0.074 | 0.014 | not met |
| atlas-known complement, hexamer-pass (calibration) | 27,564 | 0.902 | 0.715 | 0.393 | 0.144 | — |
| full precision default | 44,394 | 0.786 | 0.576 | 0.289 | 0.101 | — |
| atlas-novel, hexamer-FAIL (diagnostic) | 5,972 | 0.528 | 0.281 | 0.078 | 0.016 | — |
| sensitivity: ≥1 molecule trusted-novel | 37,438 | 0.309 | 0.103 | 0.017 | 0.003 | — |

Gene-body-shuffled null ≈ 0.007 at ≥5 UMI (10 seeds) → enrichment ~70–250×, empirical p at the 1/11
floor: the sites are far from random, but the target is an absolute fraction. Robustness (verifier):
GEM-X truth 0.554; x3p+GEM-X pooled 0.623 at 25 bp (0.698 at 50 bp, 0.739 at 100 bp) — the
pre-registered 25-bp metric is missed under every truth choice. Donor mismatch cannot explain the gap:
the atlas-known complement scores 0.90–0.94 under the same truths.

## What the data say about *why* (diagnostic, not a re-definition)

- The atlas-novelty stage costs the most concordance (0.828 → 0.521); the hexamer stage adds +0.04 on
  the whole set but **nothing among atlas-novel sites** (hexamer-fail 0.528 vs trusted-novel 0.521).
- 68.7% of trusted-novel sites are **intronic** (atlas-known complement: 68.6% in 3′ UTRs); only 15.6%
  are in annotated 3′ UTRs.
- **24.5% of trusted-novel sites lie within 25 bp of a Kinnex internal-priming decoy terminus** (6.6%
  for the atlas-known reference). The caller's IP rule (−10..+30, ≥6 A or ≥70% A) is looser than
  Kinnex's (+1..+18, ≥12 of 18 A) — the residual false-positive class is internal priming the filter
  does not catch. → Stage-1d follow-up (issue 10 addendum): a `--ip-rule kinnex` option.
- Stratification (exploratory only): 3′-UTR trusted-novel sites 0.78 at ≥5 UMI / 0.56 at ≥20;
  ≥5-molecule sites 0.73–0.76 / 0.53–0.66; 2-molecule sites 0.43 / 0.16.
- Signed distances are asymmetric: Kinnex termini sit a few bp upstream of PeakATail cleavage points
  (1,228 hits at −5..0 vs 210 at 0..5) — consistent with the clip-site convention; no effect at 25 bp.

## What the paper says (reliability-first)

1. PeakATail's atlas-novel calls are strongly enriched for real cleavage sites (70–250× over null) but
   **the pre-registered "trusted de novo PAS" definition did not reach the 70% long-read target
   (52%)**; we report this as a negative result with the funnel figure.
2. The residual error is dominated by intronic sites and by internal priming the current filter misses.
   A stricter definition (3′-UTR-restricted and/or ≥5 molecules and/or a Kinnex-style IP rule) reaches
   the target in exploratory stratification, but **was chosen after seeing the data**; it may be
   pre-registered as v2 and must be validated on data not used here (e.g. after the Stage-1d re-run, on
   a held-out truth set) before the paper calls any de novo site "trusted".
3. Until then, de novo sites are reported with their support (molecules), hexamer tier, feature class
   and IP status, and the paper makes no "trusted novel" claim.

Files: `results/reliability/trusted_novel_final_pbmc/` (REPORT_PROVISIONAL.md §1–9 incl. verifier
robustness table; call_primary/, call_sens_min1/, validate_*/). Tool: `scripts/reliability/trusted_novel_pas.py`.

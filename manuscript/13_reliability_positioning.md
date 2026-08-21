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

## Addendum 2026-08-21 05:30 — Stage-3 Laughney replication: pre-registered design (written BEFORE any cross-sample number exists)

Triggered by the Stage-3 dry-run verifier (PROBLEM verdict; scripts/stage3/, dry run GSM3516665):

1. **PAS-id space = design A.** One cohort-mode `ema run -c` over all 17 GSM datasets on the frozen
   snapshot `tools/pa-polya-run-4efeb125` (clip_seeded + IP filter, mode filter), so every sample's
   switch table shares the tool's own unified PAS ids (`unified/multi_sample_merged.bed`). No post-hoc
   coordinate harmonisation. Per-sample independent runs (design B) are NOT used for replication.
2. **Unit of replication = patient, not GSM.** 17 GSMs = 14 patients (LX675, LX682, LX684 each have a
   tumour + normal GSM). A switch "replicates" iff called (q<0.05, same direction) in ≥ 2 *patients*;
   for a patient with two GSMs, either GSM counts once. Discordant policy = exclude (any opposite-
   direction call in another patient vetoes). Effect floor |Δproportion| ≥ 0.1 reported alongside.
3. **Tested PAS universe = the pre-registered precision default** (tier-1 ∩ IP-pass ∩ ≥2 molecules;
   `TIER_FILTER=tier1_ge2`). Coverage-only tier-2 PAS are never switch-tested in the paper.
4. **Test = Fisher, `--count-mode cells`, `--marker-top-n 0`** (manuscript/14), nominal BH q per pair,
   ≥ 10 label-shuffle nulls per GSM through the identical pipeline; the replication filter is run on the
   nulls to give an empirical false-replication rate.
5. **Label policy (pre-registered; the dry run showed the curated argmax "B_cell" class is mostly
   DC/monocyte in most samples).** A cell enters the switch test only if its curated label is
   *confirmed* by an independent labeller: (a) immune types (B_cell, T_cell, NK, Monocyte, Macrophage,
   Dendritic, Mast, Plasma_IG, Treg) require agreement with the CellTypist call
   (`B2_celltypist/percell_labels.csv`, conf ≥ 0.5; mapping table fixed in
   `scripts/stage3/LABEL_POLICY.md` before use); (b) non-immune types (Epithelial_Tumor, Fibroblast,
   Endothelial, Pericyte), which CellTypist's immune model cannot confirm, require a canonical-marker
   gate on the cell's own GEX (own-type marker score highest among the four AND above the
   cohort-wide 25th percentile of confirmed cells of that type; markers: EPCAM/KRT8/KRT18;
   COL1A2/DCN/LUM; PECAM1/VWF/CLDN5; RGS5/ACTA2/PDGFRB). Unconfirmed cells are dropped, not
   relabelled. Cell types with < 20 confirmed cells in a GSM are dropped for that GSM. Per-GSM
   retained counts per type are reported before any switch statistic is viewed.
6. "Epithelial_Tumor" in Normal GSMs is normal epithelium and is reported as "Epithelial".
7. Nothing from the dry run (GSM3516665, tier filter none, unconfirmed labels) is a result.

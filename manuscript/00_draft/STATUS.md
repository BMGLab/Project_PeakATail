<!--
STATUS.md — state of the PeakATail submission draft after the 2026-09-02 fix pass
(number-trace audit -> simulated referee report -> fixes). MAIN.md and SUPPLEMENTARY.md are the
authoritative texts; the sec_*.md section files predate the audit and fix passes and are stale
(except sec_front.md / sec_declarations.md, which the PI's title/author update touched directly).
-->

# Draft status — 2026-09-02, post fix pass

## 1. What happened in this pass

Inputs: NUMBERS_LEDGER.tsv (audit verdict FIXED: 139/139 claim clusters traced, 14 fixes applied,
2 PI flags) and REVIEW_SIMULATION.md (simulated GB referee: enthusiastic major revision — 12 majors,
14 minors, 2 binding-rule defects, 1 structural slip). During the pass the PI independently set the
title and author list in MAIN.md (reviewer minor 12 thereby resolved; A.A.T. ORCID still a
placeholder).

**Applied in-draft (evidence-supported):**
- Audit PI-flag 1 — Fig S12 over-pointing: the "atlas-shaped priors" Results paragraph no longer
  cites Fig S12 in its heading; the figure is cited only for the byte-for-byte identity it carries.
- Audit PI-flag 2 — abstract R_det: range split by comparator scope ("0.175–0.208 on the three
  benchmarked libraries — below polyApipe's — and 0.110 on the second donor"), denominators named.
- Review defect 1 — bare "precision" for the P@100 metric in the Results Fig 3 paragraph and the
  Fig 3b legend: now "atlas-agreement precision" (plus one tightening in the S12 legend).
- Review defect 2 — 20,000-vs-20,320: standardised to the 25 §8.7-licensed "from 20,320 (calls)
  upward" in abstract, key points, Results one-liner; Fig 2 legend already carried it.
- Review structure slip — the two-truth hard-FP composite (15.10% vs 49.50%, 35.5% residual class)
  now debuts in Results (25 §8.6 slots it there) with call counts; Discussion keeps a one-sentence
  back-reference.
- M2 (partly) — R_det denominator provenance stated exactly (PeakATail `lambda_gradient` run for
  PBMC; pbmc4k's own run as a declared deviation; one fixed denominator shared by every tool) with
  the measured +1.26% construction-sensitivity control and its direction of bias (against PeakATail).
- M3 (partly) — one connecting sentence ties the 58.9 M library-level null tests to the 2,128,711
  patient-level (pair, PAS) hypotheses.
- M4 (partly) — Fig 4 section now states the harness has been run on no competitor's test, so the
  field-level concern remains an untested hypothesis for other tools.
- M5a/b — the R_det 0.34 "hard ceiling" is scoped to PeakATail's candidate generator, not the data;
  the 44.6% Kinnex-corroboration fraction is reframed as context for the denominator, not a bound on
  what is real (donor mismatch + un-deduplicated support named in-sentence).
- M5c — Sierra/scAPAtrap matched-N truncation now stated in Results (below PeakATail on both axes at
  every budget tested, both datasets) and in Methods (their own shipped rankings; scAPAtrap's
  19,410-call cap caveat; scAPAtrap's shipped N not reachable).
- M6a (partly) — the 20-shuffle z-scores now carry their empirical floor (p ≤ 1/21) and the pending
  200-permutation null is named; "0 of 40 null draws" is defined (both mice's 20 shuffles pooled;
  derivation re-verified from the fig5 null-shuffles TSV).
- M6b — the PDUI proximal/distal pair rule (first/last PAS of the gene in transcription order),
  depth-weighted pseudobulk statistic, ≥50-UMI-per-stage guard and label recipe are now in Methods.
- M6c (partly) — the literature-panel negative is scoped "under our PAS-to-gene assignment", with
  issue-#99 mis-assignment named as an unexcluded alternative (Results + Fig 5d legend).
- M7 (partly) — Fig 6's gene-keyed quantities (2,883 genes, 10,415 combinations, 57.8%) are flagged
  as assignment-dependent in Results and the legend; PAS-level statistics explicitly are not.
- M8 (partly) — a [[CITE]] slot for the Kinnex x3p/GEM-X data citation (with accession) was added to
  Methods.
- M9i — the `--seq-len` acceptance rule is defined (reference-span test; excludes spliced alignments;
  13.74% of valid-CB reads, 96% spliced; measured cost of relaxing it).
- M9iv — BH scope stated: one BH family per (library, cell-type pair).
- M9v — two-library patients: counts-once and the veto across a patient's own two libraries, with
  the audited numbers (89 collapsed; 2 vetoed).
- M11 — Discussion ¶3 condensed (anecdote details trimmed; the disclosures remain in Methods and
  Fig S11).
- M12 — a Methods paragraph on the cell-independence assumption and the structural defence
  (patient-level replication; nulls that preserve within-library dependence).
- Minors 1, 2, 4 ("dominates" defined in-sentence: higher on both axes at the same call count),
  5 (92-fold measurement defined), 6 ("0 of 40" defined), 10 (already fixed by audit), 11 (pair
  universe: 59 pairs, 12 single-patient, 47 replicable — in Methods).

**Routed to PI_DECISIONS.md (new work or PI call):** M1+M8 core (panel-wide de-duplicated Kinnex
axis — the review's top request), M3 (more permutations), M4 (competitor calibration run), M5
(depth downsampling), M6a (200 shuffles), M6c (16-gene manual check), M7 (non-overlap recompute),
M2 residual (tool-independent denominator rebuild), M9ii/iii (marker panel list; curated-label
provenance), M9vi/A3 (repo public + Zenodo), M10 (T7, S5/S6), minors 7 (gate CIs), 8 (genotype
check), 9 (Fig 3c panel marking), 13 (key-point overlap), 14 (`--polya-min-umis` default).

## 2. Verification after editing

- Both mandated verbatim claim blocks re-verified word-for-word against 18 §v2 and 20
  ("Recommended claim wording") after all edits.
- Forbidden-phrasing greps re-run, clean: no beats/outperforms/superior (the Background "can beat a
  simple rule on atlas agreement" describes our own model vs our own rule, as audited); no 0.9986,
  no "104 M", no 72.5%, no ~17×; SCAPTURE mouse 2 never "invalid" (the word appears only inside the
  mandated never-say-it caveat); 1.152% only as the corrected/superseded estimate; "none in 10
  label-shuffle nulls (empirical p ≤ 0.091, the 10-permutation floor)" never an FDR; cohort
  phrasing "17 libraries from 14 patients" / "15 libraries from 12 patients" intact; no positive
  trusted-novel claim; no clustering-novelty claim; every remaining bare "precision" is either the
  licensed matched-N formulation, "precision-first"/"precision default" (names), or generic
  non-metric use; no 20,000-call leftovers; nothing cited from S5/S6.

## 3. Word counts (LC_ALL=C tokens, comments excluded; table also at the end of MAIN.md)

| Section | Words |
|---|---|
| Title page block | 52 |
| Abstract + keywords | 293 (abstract alone 276) |
| Key points | 350 |
| Background | 924 |
| Results | 4,415 |
| Discussion | 690 |
| Conclusions | 134 |
| Methods | 3,236 |
| Declarations | 347 |
| Figure legends (Figs 1–6) | 1,963 |
| Supplementary legends section | 566 |
| **TOTAL (MAIN.md)** | **12,970** |

SUPPLEMENTARY.md (Additional file 1) is unchanged except one language tightening in the S12 legend.

## 4. Ledger totals

NUMBERS_LEDGER.tsv now holds **148 rows** (139 from the audit + 9 added by this pass), all
**traced**. Rows updated this pass: 5 (abstract R_det split), 6/16/38 (20,320 harmonisation),
19 (atlas-agreement wording), 23 (92-fold definition), 65 (0-of-40 derivation + 1/21 floor),
85 (two-truth moved to Results). Rows added: 140–148 (two-truth in Results; Sierra/scAPAtrap
matched-N; span rule; denominator provenance/control; pair universe; replication audits; PDUI
definition; Fig 6e assignment caution; null-test link).

## 5. Placeholder and citation inventory

**[[PLACEHOLDER]] (6, all PI-decision/author items):** A.A.T. ORCID (optional); Zenodo DOI
(blocked on repo publication); competing interests; funding; author contributions (CRediT —
suggested split drafted in-comment); acknowledgements.

**[[CITE]] (47 instances, 26 distinct targets, all in MAIN.md):** tool/resource papers with
unambiguous identity (PolyASite 2.0 ×4, polyApipe ×4, Kinnex ×4, Sierra ×3, scAPAtrap ×3,
SCAPTURE ×3, scTail ×3, scUTRquant ×3, CellRanger ×2, STARsolo ×2, bedtools, Signac TF-IDF,
CellTypist, Benjamini–Hochberg ×2 phrasings — the reference manager should merge these two);
literature placeholders (APA reviews/mechanisms ×3, internal-priming artefacts, benchmarks of
single-cell APA tools, spermatogenesis panel sources); dataset citations (GSE104556 ×2 phrasings —
merge; Laughney ×2 phrasings — merge; **new:** the public Kinnex PBMC x3p/GEM-X datasets with
accession — must be filled from the truth-build record, review M8).

## 6. What stands between this draft and submission

1. **A1 (blocking, top reviewer request):** de-duplicate the Kinnex truth, recompute the four
   inheriting numbers, extend long-read concordance to the whole panel at matched N (also completes
   S5 and fills the Kinnex data accession).
2. **A2 (blocking):** Supplementary Table T7 collated from the run registry; S6 built or the S5/S6
   slots dropped with a single renumbering pass.
3. **A3 (blocking):** companion repository public + Zenodo DOI into the placeholder.
4. **PI placeholders:** competing interests, funding, CRediT, acknowledgements, A.A.T. ORCID.
5. **Reference pass:** resolve all 26 [[CITE]] targets (merge the double-phrased BH/GSE104556/
   Laughney entries).
6. **PI decisions B1–B7 and C1–C9** (PI_DECISIONS.md): scheduled runs or declared declines for the
   reviewer's strongly-requested analyses; confirmations of the two applied wording calls (C1, C2).
7. **Figure regeneration items:** Fig 3c v1-era row marking (m9); gate-value CIs if C5 is accepted
   (Fig 2/S4 scripts); any Fig 2 panel change from A1.
8. **Housekeeping:** re-sync or retire the stale sec_*.md section files against MAIN.md before the
   next assembly; word-count comment is current as of this pass.

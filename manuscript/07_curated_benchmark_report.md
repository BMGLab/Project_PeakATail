# Curated-reference benchmark + atlas-independent motif validation — round 2 report

Generated 2026-08-18. Companion to `05_figure_index.md` (round 1). This round executes priority items
2, 3 and 4 of the round-1 "next analyses" list: re-benchmark against a curated reference with point-mode
strand-matched matching and a single reconciled null, plus an orthogonal, atlas-independent sequence
validation of the calls. Both analyses were independently re-verified by an adversarial pass
(verdict: **SOUND** for both); all numbers below are post-verification.

Scripts: `/mnt/ssd1/Projects/PeakATail_wd/scripts/manuscript_figures/benchmark_curated.py`,
`.../motif_validation.py`. Data: `/mnt/ssd1/Projects/PeakATail_wd/results/figures/manuscript/
benchmark_curated.tsv` (405 rows), `.../motif_validation.tsv` (6,455 rows). Figures shipped to
`manuscript/figures/` as PNG (300 dpi) + PDF (fonttype 42). All sorts and bedtools steps under
`LC_ALL=C` (tr_TR locale hazard).

---

## 1. What changed versus the discredited 18.4M benchmark, and why

Round 1 established (figs 4–5 of the index) that the 0.9986 precision headline was a reference-density
artifact: the reference was an **18,432,135-entry unfiltered PolyASite dump** (92% IN/DI/UI entries,
±100 bp windows tiling a third of the genome), matched **strand-agnostically** with **whole-interval
window padding** (median peak width 403–811 bp), on which randomly placed width-matched peaks already
scored 0.62–0.75. Round 2 replaces every one of those inflators:

| Component | Old (figs 4–5) | New (this round) |
|---|---|---|
| Reference | 18.4M-entry unfiltered dump | **Curated PolyASite 2.0 atlas, 569,005 clusters (GRCh38.96)** via representative point sites — ~32× sparser, biologically curated |
| Secondary reference | none | Protein-coding TES, 73,439 unique ends (Ensembl 99, GRCh38.p13) |
| Query geometry | whole peak interval (contains-a-site = perfect hit) | **point mode: peak 3'-most base only** |
| Strand | agnostic (`window`, no `-s`) | **strand-matched** (`bedtools closest -s -d`) |
| Null | two conflicting definitions across figs 4/5 | **one reconciled null**: width-preserving shuffle constrained to merged annotated gene bodies (1.80 Gb), same-chromosome, non-overlapping, 3 seeds × 5 arms |
| Chrom naming / sort | mixed | Ensembl naming throughout; `LC_ALL=C` verified on both references |

Reference provenance (all md5-logged in `data/references/atlases/README.md`):
- **Primary**: `polyasite2.GRCh38.96.rep_sites.bed6` — 569,005 representative point sites (score = TPM),
  built from the raw PolyASite 2.0 download (in the expected 500–600k range, *not* the 18.4M dump).
  Sanity checks passed: 0 representative sites outside their own cluster, 0 strand/chrom mismatches.
- **Secondary**: `tes.protein_coding.GRCh38.99.bed6` — 73,439 unique protein-coding transcript 3' ends
  from 84,077 transcripts (score = n transcripts sharing the end, not expression).
- **On hand but unusable as-is**: PolyA_DB v3.2 (`polyadb3.2.human.hg19.pas_sites.bed6`, 311,594 PAS) —
  **hg19-only** and UCSC `chr` naming; liftOver/CrossMap are not installed on this machine, so no GRCh38
  version exists. Do not intersect with GRCh38 calls without lifting over first.

**Consequence:** the new numbers are *not comparable* to the old 0.9986 — they use a deliberately
stricter matcher than the pipeline's own (round 1's `null_control` reconciles the two matchers).
Precision at 100 bp drops from 0.9986 to 0.36–0.49; what survives is a ~16–23× enrichment over the
matched null at 100 bp (~20× typical). The enrichment fold is cutoff-dependent: 25–42× at 10 bp,
19–32× at 50 bp, 16–23× at 100 bp, but only 10–14× at 200 bp — quote the fold with its cutoff.

---

## 2. Results per arm and reference, versus the roadmap bars

Roadmap publication bars (`06_roadmap.md`): **precision ≥ 0.70, recall ≥ 0.60, F1 ≥ 0.65 at 100 bp**.

### Primary table — curated PolyASite 2.0 representative sites, 100 bp, strand-matched point mode

| Arm | n sites | Precision | Null precision (3 seeds) | Recall (a) full atlas (n=569,005) | Recall (b) detected-gene atlas (n=285,220) | F1 (b) | Pass P/R/F1? |
|---|---|---|---|---|---|---|---|
| lg_annotate (= lg_ip_off) | 22,629 | 0.447 | 0.021–0.024 | 0.033 | 0.066 | 0.115 | **FAIL / FAIL / FAIL** |
| lp_annotate | 21,470 | **0.492** | 0.020–0.022 | 0.036 | 0.072 | 0.125 | **FAIL / FAIL / FAIL** |
| si_annotate | 43,023 | 0.355 | 0.022–0.023 | **0.048** | **0.095** | **0.150** | **FAIL / FAIL / FAIL** |
| lg_ip_filter | 20,605 | 0.457 | 0.020–0.022 | 0.031 | 0.062 | 0.109 | **FAIL / FAIL / FAIL** |
| B1_cohort_full | 16,497 | 0.475 | 0.022–0.023 | 0.025 | 0.050 | 0.091 | **FAIL / FAIL / FAIL** |

F1 against the full atlas (flavor a) is lower still: 0.048–0.084. At the most permissive cutoff
tested (200 bp): best precision 0.557 (lp/B1), best recall-b 0.140 (si), best F1-b 0.210 (si) — still
nowhere near any bar.

### Secondary table — protein-coding TES (n=73,439), 100 bp

| Arm | Precision | Recall | Pass? |
|---|---|---|---|
| lg_annotate | 0.280 | 0.149 | **FAIL** |
| lp_annotate | 0.313 | 0.159 | **FAIL** |
| si_annotate | 0.204 | 0.202 | **FAIL** |
| lg_ip_filter | 0.293 | 0.142 | **FAIL** |
| B1_cohort_full | 0.303 | 0.120 | **FAIL** |

TES maxima anywhere in 10–200 bp: precision 0.348, recall 0.239 (both @200 bp) — never approaches a bar.

### Verdict versus the bars

**0 of 150 roadmap rows pass. Every arm fails all three bars at every cutoff from 10 to 200 bp, on
both recall flavors and on both references.** This is not marginal: the best F1 anywhere (0.210,
si @200 bp, favourable denominator) is less than a third of the 0.65 bar.

Interpretation constraints that must accompany this table:
- **Recall is arithmetically capped by set size.** With 16.5–43k called sites against a 285–569k-site
  reference, per-arm recall ceilings (n_called / 285,220) are 0.058 (B1), 0.072 (lg_ip), 0.075 (lp),
  0.079 (lg), 0.151 (si). The recall bar of 0.60 is unreachable *by construction* for a single-cohort
  3' assay against a pan-tissue atlas; recall-b rises roughly linearly with n called (si best at 0.095).
- Recall flavor (a) additionally includes 397 atlas sites on scaffolds absent from the genome file,
  unreachable by construction; flavor (b) is the favourable denominator and still misses the bar ~6×.
- The matcher is deliberately stricter than the pipeline's (strand-matched, point mode); these numbers
  measure single-base 3'-end placement, not whole-peak overlap.
- The arms now separate on precision (0.355–0.492, spread 0.14 — versus 0.004 on the old dump), with the
  usual precision/recall trade: si calls ~2× the sites, lowest precision, highest recall.
- Set sizes are after dropping 4/4/12/4/3 scaffold peaks; `lg_ip_off` is byte-identical to `lg_annotate`
  (md5 `20cb4550…`) and was scored once — the round-1 duplicate-arm caveat still applies.

---

## 3. Atlas-independent evidence: the calls flank genuine polyadenylation sites

`motif_validation` tests the sequence around each peak's 3'-most base against the same genic-shuffle
null, using no atlas at all (real n=22,629; ip-filtered n=20,605; 3 null seeds).

**The peak 3' end is not the cleavage site — it stops ~90–105 nt short.** Windows that assume
peak end = cleavage contain the signal at near-null rates: AATAAA|ATTAAA in the canonical −40..−5
window is 7.7% (real) vs 5.4% (null), nowhere near the 70–80% literature rate for true PAS. But the
signal is not absent — it is displaced downstream:

- **AATAAA|ATTAAA within +0..+100 nt downstream of the peak end: 37.1% real vs 14.7% null (2.5×);
  any-of-12 canonical hexamers: 65.8% vs 43.7%.**
- AATAAA positional density is enriched >2× null across +29..+94 nt, mode **+75 nt** (bimodal
  sub-peaks +48/+75 — at least two offset populations, cause not investigated).
- Base composition shows the full canonical architecture: A-fraction climbs to **0.426 at +98 nt then
  cliffs to background (~0.30)** — the modal cleavage position, where the poly(A) tail would begin;
  T overtakes A past the modal cleavage position (directionally supported but weak, ~1 pp above null —
  do not lean on "U-rich DSE" language).
- Given canonical 15–30 nt AATAAA-to-cleavage spacing, implied cleavage is **+90..+105 nt past the
  called peak end** — consistent with 10x R2 coverage ending before the poly(A) junction. This offset
  also rationalizes the section-2 geometry: point-mode precision rises 0.30 → 0.45 → 0.52–0.56 across
  50/100/200 bp as the cutoff grows past the offset (interpretation, not separately proven).
- Even re-anchored, 37% is below the 70–80% of curated PAS — per-site offset variability and residual
  false calls are not separable here.
- **The internal-priming filter demonstrably works:** A-rich downstream tracts (≥6 consecutive A or
  ≥70% A at +10..+30) drop from 6.3% (real) to **3.0%** (ip-filtered; 628/20,605 = 3.05% — do not
  write 3.1%) against a 2.2% null; run-A6 0.060→0.028, fracA70 0.024→0.012. The filter removes A-tract
  sites without changing positional geometry (ip tracks real in every other panel).

Together with section 2 this is the evidential core: enrichment ~20× over null on a curated atlas,
plus an atlas-independent canonical PAS sequence architecture at the expected (offset) position, plus
a mechanistically sensible filter effect.

---

## 4. What the manuscript can now honestly say

**Retired permanently:** "precision 0.9986" (and any number against the 18.4M dump), atlas recall/F1
against the full dump, and any claim that the atlas benchmark ranks the five strategies (round 1);
"99.4% of PAS in annotated 3'UTRs" unconditioned (round 1, fig 1).

**Usable claims (verified wording):**

1. *Benchmark reframe:* "Benchmarking targets the curated PolyASite 2.0 atlas (569,005 clusters,
   GRCh38.96) via its representative point sites — a ~32-fold sparser, biologically curated reference
   versus the 18.4M-entry dump that inflated precision to 0.9986 — with annotated protein-coding TES
   (73,439 sites, Ensembl 99) as an independent secondary reference."
2. *Benchmark result:* "Against the curated atlas, strand-matched point-mode precision at 100 bp is
   0.36–0.49 across all pipeline arms — roughly 20-fold above a width-preserving genic-shuffle null
   (0.02) — but recall is at most 0.05 against the full atlas and 0.09 restricted to the 14,851
   detected genes, so under this strict matcher no arm reaches conventional benchmark bars
   (P ≥ 0.70, R ≥ 0.60, F1 ≥ 0.65). The honest claim is strong enrichment over null, not
   benchmark-passing accuracy." (Tie the ~20× fold to the 100 bp cutoff; at 200 bp it is 10–14×.)
3. *Sequence validation:* "An atlas-independent sequence test shows the called peaks flank genuine
   polyadenylation sites — AATAAA/ATTAAA occurs within 100 nt downstream of the peak 3' end in 37% of
   sites versus 15% in gene-body-shuffled nulls, with the full canonical architecture (AATAAA density
   mode +75 nt; A-fraction cresting at 43% at +98 nt then collapsing to background where the poly(A)
   tail would begin; T overtaking A past the modal cleavage position) — but the peak 3' end itself is
   not the cleavage site: it systematically stops ~90–105 nt short, consistent with 10x R2 coverage
   ending before the poly(A) junction."
4. *IP filter:* "The internal-priming filter acts as designed, halving A-rich downstream tracts
   (6.3% → 3.0%; null 2.2%)."
5. *Framing:* describe called peaks as "peaks whose 3' ends flank/precede cleavage sites", not as
   cleavage-site calls; describe recall against pan-tissue references as size-capped (ceilings
   0.06–0.15), not as sensitivity.

**Still not sayable:** that any arm meets a published accuracy benchmark; that 37% downstream-hexamer
incidence matches curated-PAS rates; anything based on PolyA_DB v3.2 (hg19-only, un-lifted);
"U-rich downstream element" as a strong claim; and everything in the round-1 "not established" list,
above all any APA/PDUI biology (SWITCH_CELLTYPE still not rerun).

---

## 5. Open items

1. **SWITCH_CELLTYPE rerun** (round-1 priority 1, unchanged, still the only item that adds biology).
   Fix documented in `scripts/pipeline/rerun_switch_celltype_fix.md`; ~24 short tasks, upstream cached.
2. **The +90–105 nt offset question.** Decide whether to model the systematic 3'-truncation
   (shift anchors / extend peaks 3' by the modal offset, or re-call cleavage from read 3' ends) and
   re-benchmark; the bimodal offset populations (+48/+75) are uninvestigated and may be strategy- or
   class-specific. This is the single change most likely to move precision materially.
3. **Roadmap bars need a decision.** P ≥ 0.70 / R ≥ 0.60 / F1 ≥ 0.65 at 100 bp is arithmetically
   unreachable on recall for this design (ceilings 0.06–0.15); either restate the bars for a
   single-cohort assay (e.g. precision + fold-over-null + motif architecture) in `06_roadmap.md`,
   or the manuscript ships with "fails its own bars".
4. **PolyA_DB v3.2 as a third reference** requires installing liftOver/CrossMap, lifting hg19→GRCh38,
   and harmonizing UCSC↔Ensembl chromosome naming before any intersection.
5. **Reconcile round-1 figs 4/5 with this round:** `benchmark_strategies` (fig 4) still shows
   dump-based 0.996–1.000 curves and its own null; retire it or move to supplement with the curated
   benchmark as the primary figure, and propagate the reconciled null definition everywhere.
6. Minor/latent: `benchmark_curated.py` sh pipelines run without pipefail (a bedtools failure could be
   masked and cached — this run was proven correct by independent recomputation; fix before reuse);
   null shuffle lets 95/22,629 anchors (0.42%) hang just off gene-body ends (negligible); TES `score`
   is transcript-sharing count, not expression; delete `.cache_benchmark_curated/` to force recompute.

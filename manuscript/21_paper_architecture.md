<!-- SUPERSEDED FRAMING NOTE (2026-08-22): this document was written before the matched-call-count
analysis in 25_competitive_position.md and before the second-donor validation in 26. Where it states the
head-to-head position as a single operating point, 25 supersedes it; where it treats one human donor as
the deepest gap, 26 (pbmc4k, gate PASS at P@100 0.8279, cross-donor concordance 84.2% within 100 bp)
partially closes it. The claim ladder and figure scheme below stand. -->

# 21 — Paper architecture: the story, the figures, the claim ladder

**Status: PROVISIONAL (architecture decisions, not new measurements) — 2026-08-21, written for the PI to steer the write-up from.**
**Adversarially verified 2026-08-21, verdict FIXED: twelve corrections applied into the text; see [§11](#11-adversarial-verification-pass-2026-08-21-verdict-fixed).**
This document supersedes the *structure* of `01_outline_and_journals.md` v0.4 (title/abstract/figure plan/journal
ranking). It does **not** supersede `01`'s number policy, its Methods outline, or any verified result document.
Where `01` and this file disagree on structure, this file wins; where they disagree on a number, the cited
verified document wins over both.

**Number policy for this file.** Every number below is either (a) quoted from an adversarially verified
document — `19_final_gate_v2.md` (FIXED), `16_trusted_novel_kinnex.md` §v2 (FIXED), `14_switch_calibration_v2.md`
(SOUND), `18_spermatogenesis_final.md` (FIXED), `20_stage3_replication.md` (FIXED, v1 code), `15_final_gate.md` §3
(SOUND, competitor rows), `10_caller_fix_plan.md`, `09_headtohead_results.md`, `05_figure_index.md`, and the
generated `figures/*.caption.md` sidecars — or (b) computed by me in this pass and labelled **[checked here]**.
Nothing else. Four things were re-derived directly from primary files while writing this (§10).

---

## 0. The decisions, on one screen

| # | Decision | Why |
|---|---|---|
| D1 | The paper is a **calibrated-reliability methods paper**, not a sensitivity paper and not a discovery paper. | We lead every de novo caller on precision on both datasets and trail polyApipe on recall on both; the biology we can name is blocked by a gene-assignment defect (issue #99). Reliability is the only claim the evidence carries end to end. |
| D2 | **Six main figures.** Fig 1 overview (new), Fig 2 accuracy, Fig 3 trade surface, Fig 4 calibration, Fig 5 cohort replication, Fig 6 the negative de novo result. | Genome Biology Method allows 6–8; six is what the ladder supports without padding. |
| D3 | **The spermatogenesis control goes to Supplementary S1**, with its own Results paragraph and a line in the abstract — but that line is added **only after the v2 re-run lands** (gap 2). The §3 draft carries no spermatogenesis number on purpose: `18` is v1 code, exactly like the Stage-3 counts the abstract also omits. | It is a positive control, its literature-panel leg does not reproduce, its across-gene index reverses, and it is the one asset that must be re-run anyway. Its main-text payload is two sentences. |
| D4 | **The negative result gets a main figure (Fig 6).** | It answers the question the title asks ("how many can be trusted?"), it is the paper's most quotable field-level lesson, and a pre-registration whose failures live in the supplement is decoration. |
| D5 | **Retire `cohort_qc`, `parameter_sweep`, `benchmark_strategies` as figures; rebuild `null_control` and `clustering_concordance`; regenerate `kinnex_truth_validation` and `motif_validation` as supplements.** | §4.3 argues each one individually. |
| D6 | **Rename everything once, after the two v2 regenerations land** (`fig01_overview` … `figS11_compute`). | Renaming before the last regeneration means renaming twice and desynchronising captions from scripts. |
| D7 | **Journal: Genome Biology (Method) → Genome Research → Briefings in Bioinformatics.** Nature Communications is off the top list until issue #99 lets us name genes. | §9. |
| D8 | **Do not wait for `--ip-rule kinnex` or a trusted-novel v2 definition.** They are follow-up work, named as such in the Discussion. | Chasing the failed gate turns a clean negative result into an unfinished one. |

---

## 1. The one-sentence claim, and what carries it

> **PeakATail makes single-cell alternative-polyadenylation analysis auditable rather than merely productive:
> seeding poly(A) sites from non-templated poly(A)-tail read evidence and shipping a pre-registered,
> internal-priming-filtered, ≥2-molecule default produces the most atlas-concordant de novo PAS call set among
> the de novo callers benchmarked, on two datasets, corroborated by donor-independent long reads; and pairing
> it with the only one of six differential-test configurations that controls the false-discovery rate and with a
> patient-level replication filter yields cell-type APA switches that survive label-shuffle nulls — at a recall
> at a deliberately conservative point on a precision/recall curve that, at matched call count, lies above every de novo
> competitor's on both axes (see [25_competitive_position.md](25_competitive_position.md); the recall half of that claim
> is robust at every matched N tested, the precision half is not established for 20,320 <= N <= 35,759), and with one of
> the four pre-registered claims reported as a failure.**

This is a paper about **calibrated reliability**, and it must say so in its own abstract. It is *not* a paper
that claims best-in-class sensitivity (we are below polyApipe's recall on both datasets and far below
scAPAtrap's), *not* a paper that claims overall accuracy leadership (the F1_det lead over polyApipe is +0.020
PBMC / +0.013–0.014 mouse, which is not a margin, and catalog-based scUTRquant is above every de novo arm), and
*not* a discovery paper (no named-gene biology can be published until PAS→gene re-assignment, `20` disclosure 1 /
issue #99).

### The four load-bearing results

| # | Result | Verified source |
|---|---|---|
| **L1** | The pre-registered precision-first default reaches **atlas-agreement precision @100 bp of 0.7062 (PBMC 10k v3, n 46,524), 0.7450 (testis mouse 1, n 26,255) and 0.7572 (mouse 2, n 26,526)** against curated PolyASite 2.0 representative sites — gate P ≥ 0.50 **PASS on all three**, 33–55× the 3-seed gene-body-shuffled null (0.0217; 0.0127–0.0150) — at detected-gene recall 0.1754 / 0.2048 / 0.2080. It is the most atlas-concordant **de novo** call set in a six-tool panel on both datasets (next best de novo: SCAPTURE 0.652 PBMC, 0.694/0.672 mouse). | `19_final_gate_v2.md` §1–2 (verifier: FIXED, every value reproduced to six decimals by an independent bedtools pipeline); competitor rows `15_final_gate.md` §3 |
| **L2** | **Atlas-independent corroboration:** 0.7647 [0.7608–0.7685] of the default PBMC sites have a poly(A)-verified PacBio Kinnex long-read 3′ end within 25 bp at ≥5 UMI (35,575 / 46,524), against a gene-body-shuffled null of 0.0072 (106.5×); the atlas-known, hexamer-pass complement scores 0.8941, so the metric is calibrated on sites we already trust. | `16_trusted_novel_kinnex.md` §v2 (FIXED) + `results/reliability/trusted_novel_final_v2_pbmc/REPORT_PROVISIONAL.md` (status VERIFIED) |
| **L3** | **Only one of six shipped/candidate differential-APA test configurations controls FDR.** Fisher, `--count-mode cells`, marker pre-selection off: 3.0% of label-shuffled tests p < 0.05, 0.00% q < 0.05, **0/20 null runs (0/60 BH families)** with any hit — conservative, not exactly calibrated. The shipped defaults are anti-conservative (reads-Fisher 20.3%, cells-Fisher with top-200 marker pre-selection 13.0%, `nb_pairwise` 24.7%), with three named mechanisms: UMI pseudoreplication, marker double-dipping (which also changes the Fisher denominator), and the 1e-4 NB dispersion floor carrying 67% of that arm's false hits. | `14_switch_calibration_v2.md` (SOUND) |
| **L4** | **Switches are only reported when they replicate across independent biological units, and nothing replicates under the null.** Laughney lung adenocarcinoma, 12 patients / 15 libraries, one cohort PAS space, 74,954-site precision-first universe, 59 cell-type pairs: 14,480 of 1,977,134 tested (pair, PAS) hypotheses replicate in ≥2 patients in the same direction with an opposite-direction veto (0.73%), 13,826 of them over \|Δproportion\| ≥ 0.1, covering 5,395 PAS in 2,687 genes across 47 pairs — and **0 replicated features in each of 10 patient-wise label-shuffle nulls**. Independently, on two testis mice: 6,049 / 9,469 / 11,263 same-direction replicated PAS per stage pair, 99.7–99.8% sign agreement, **0 in all 15 null pairings**, per-gene effect ρ = 0.655 (n = 917; null 0.001 ± 0.032). | `20_stage3_replication.md` (FIXED; **v1 code `4efeb125` — the Stage-3 v2 chain on `9dfdefb` is finishing tonight and supersedes these counts**) and `18_spermatogenesis_final.md` (FIXED; also v1 code) |

**The honest fifth result, which is a limitation and a finding at once (and belongs in the Results, not only the
Discussion):** the discriminator between callers is **evidence type, not peak-calling sophistication** — our own
coverage-only caller scored P@100 0.118 on the same PBMC BAM and ranked last among de novo tools; re-seeding the
identical pipeline on poly(A) clips took it to 0.7062 (`09_headtohead_results.md` §"Why PeakATail underperforms",
`19` §1). The price is a hard sensitivity budget: only **1.152%** of CB-bearing reads in that BAM carry a
non-templated poly(A) soft clip and only **28.77%** of detected-gene atlas sites carry any clip read within
100 bp (`10_caller_fix_plan.md` §R1), which is why **no clip-seeded arm can pass R_det ≈ 0.29** and why 0.4 is
out of reach for this evidence type at this chemistry. **Do not stretch that budget to explain the gap to
polyApipe** (Gap 0): polyApipe uses the same clip evidence and reaches 0.199, and our own ≥1-molecule arm
reaches 0.2685 — the further fall to the default's 0.1754 is the pre-registered ≥2-molecule threshold, a
chosen operating point, not the budget.

---

## 2. Title — five candidates, ranked

| Rank | Title | What it promises | Can we deliver it? |
|---|---|---|---|
| **1** | **How much of a single-cell poly(A)-site call set can be trusted? Precision-first calling, calibrated switch testing and pre-registered evaluation in PeakATail** | A reliability audit of the field plus a tool that answers its own question. | **Yes, fully.** Every element is measured: precision vs two truths, the FDR audit of six test configurations, replication, and a pre-registered gate we failed and report. This is the only title where the failed gate is an asset rather than an embarrassment. |
| **2** | **Precision-first poly(A)-site calling and calibrated cell-type APA switch testing from 3′ single-cell RNA-seq** | A method with two components, both validated. | **Yes.** Descriptive, safe, no over-promise. Weakness: it is forgettable and does not signal the benchmark, which is half our value. |
| **3** | **Poly(A)-tail read evidence, not peak shape, determines what a single-cell APA caller can be trusted to find** | A field-level finding with a tool attached. | **Mostly.** The evidence-type result is verified on two datasets and one internal before/after (0.118 → 0.706 on the same BAM), but "determines" is stronger than a six-tool, two-dataset panel supports; "determines" would have to become "dominates" or "separates". Best title for Briefings in Bioinformatics. |
| **4** | **PeakATail: tiered, internal-priming-filtered poly(A)-site detection with replicated cell-type APA switches in single cells** | A tool paper naming its mechanisms. | **Yes**, but it buries the calibration and the pre-registration, which are the parts reviewers cannot get from any competitor. |
| **5** | **Low-false-positive poly(A)-site discovery and replicated 3′UTR switches in single cells with PeakATail** | Compact, journal-title-length; leads on false positives. | **Yes on the words, no on the frame.** "Discovery" invites the de novo-site question we answer with a negative result; the reader arrives expecting novel sites and finds a failed gate. Keep only if a length limit forces it. |

**Retired for good:** anything leading with clustering, "cell-state discovery from PAS usage alone",
"expression-independent", "beyond gene counts", or any form of "most accurate".

---

## 3. Abstract (draft, verified numbers only)

> Alternative polyadenylation remodels 3′UTRs across cell types, and its signal is present in every 3′
> scRNA-seq dataset, yet poly(A)-site (PAS) calls are rarely scored against independent truth and
> differential tests are rarely calibrated. We present PeakATail, which seeds PAS from non-templated
> poly(A) soft-clipped reads, emits clip-supported and coverage-only sites as separate, molecule-counted tiers,
> filters internal priming, and tests cell-type APA switches with calibrated,
> replication-filtered statistics. Under gates registered before the final run, its default output —
> clip-supported, ≥2 distinct molecules, internal-priming-filtered — reaches atlas-agreement precision at 100 bp
> of 0.706 on PBMC and 0.745/0.757 on two testis replicates against curated PolyASite 2.0 sites
> (gene-body-shuffled nulls 0.013–0.022), the highest of the de novo callers benchmarked on both datasets, at
> detected-gene recall 0.175–0.208, below polyApipe's. Independently, 76.5% of default PBMC sites lie within
> 25 bp of a poly(A)-verified Kinnex long-read 3′ end at ≥5 UMI. Of six differential-test configurations, only
> Fisher on cells without marker pre-selection controls the false-discovery rate (3.0% of label-shuffled tests
> p < 0.05), and switches are reported only when replicated across patients, where label-shuffle nulls yield
> none. A pre-registered "trusted novel PAS" definition reached 48.5% long-read support against a 70% target and
> is reported as a negative result. Code: github.com/BMGLab/PeakATail.

**200 words [checked here: `wc -w` on the blockquote body, re-counted after the verifier's edits].** Sources:
`19` §1 (0.706 / 0.745 / 0.757, recall range, nulls), `16` §v2 + the v2 REPORT (0.7647 → "76.5%", 0.4853 →
"48.5%"), `14` (3.0%, six arms), `20` (the replication sentence, wording only — **no Stage-3 count appears in
the abstract on purpose**, because the v2 chain supersedes those counts tonight and an abstract number must not
move after the figures freeze). **Two verifier corrections are already folded in.** (i) The Kinnex sentence now
carries **"at ≥5 UMI"**: the fraction is threshold-dependent (0.7647 at ≥5, 0.5584 at ≥20, 0.2807 at ≥100), so
the threshold is load-bearing and may never be dropped for brevity. (ii) **"MIT-licensed" was removed.** `01`
v0.4 asserts MIT, but gap 4 of §8 records that `Project_PeakATail` has **no `LICENSE` file** and that the
licence is still a PI decision — an abstract must not state a licence the repository does not carry. Restore
the word only after the file exists. **The abstract also carries no spermatogenesis number**, deliberately, for
the same reason it carries no Stage-3 count: `18` is v1 code and gap 2 supersedes it (see D3).
For Nature Communications' 150-word limit, cut the testis replicate values and the null range. For Genome
Biology's structured abstract, split at "We present PeakATail" (Background / Results / Conclusions).

---

## 4. Final figure scheme

### 4.1 Naming convention (adopt now, execute once — see §4.4)

> **[EXECUTED 2026-09-02 — with amendments.]** The rename pass ran on 2026-09-02 under the approved
> definitive map (recorded in `figures/FIGURE_MAP.tsv` and `figures/FIGURES_MANIFEST.md`), which amends
> this section: (a) stems are **unpadded** (`fig1_…`/`figS1_…`) because the built Fig 1 already shipped as
> `fig1_overview` and D6's rename-once rule forbids repadding it; (b) D3 is amended — spermatogenesis is
> **main Fig 5** (`fig5_spermatogenesis`) and the cohort figure is **Fig 6** (`fig6_cohort`); (c) D4 is
> amended — the trusted-novel negative result moves to **`figS7_novelfunnel`** with its full text
> prominence kept (own subsection R6, abstract sentence, Author Summary line); (d) the supplement is
> renumbered S1–S12 (sperm leaves for main; novelfunnel enters as S7; new slots for the second-donor
> validation `figS4_seconddonor` and the four-way version benchmark `figS12_versions`; §4.3's proposed
> `figS10_ipfilter` has no slot in the executed map — the IP-filter measurement lives in Fig 1 panel d and
> the caveat text). No numbers were
> changed in the pass: all renamed figures' audit TSVs verified byte-identical (`cmp`). Env-var switches
> keep their pre-rename names. The tables below are preserved as the historical proposal.


```
main:           fig01_overview   fig02_accuracy   fig03_tradeoff
                fig04_calibration  fig05_cohort   fig06_novelfunnel
supplementary:  figS01_spermatogenesis … figS11_compute
```

One **stem** owns every artefact of a figure, with no exceptions:

| Artefact | Path |
|---|---|
| generating script | `scripts/manuscript_figures/<stem>.py` (each script already exposes a single `NAME = "<stem>"` constant — `cohort_qc.py` is the only one that does not, and it is retired) |
| raster + vector | `manuscript/figures/<stem>.png` (300 dpi) and `<stem>.pdf` (fonttype 42, vector) |
| caption sidecar | `manuscript/figures/<stem>.caption.md` (written by the script, never by hand) |
| every plotted value | `results/figures/manuscript/<stem>*.tsv` |

Rules that make the convention worth having: the number in the stem **is** the number in the manuscript, so a
figure cannot be reordered without renaming its script; the caption sidecar is generated by the script, so a
caption cannot drift from the numbers; and `manuscript/figures/FIGURE_MAP.tsv` (new, 5 columns: `old_stem`,
`new_stem`, `script`, `source_docs`, `status`) is the only place the old names survive, so that every verified
caveat in `05_figure_index.md` can be traced to its new home.

### 4.2 Main figures

---

**Fig 1 — `fig01_overview` — What PeakATail does and what its default output is.**
- **Proves:** nothing on its own; it makes claims C1, C3, C5 legible and it is where the pre-registration is
  stated in the main text.
- **Supports:** the whole ladder (orientation figure).
- **Panels.** **(a)** Data flow: tagged BAM (CellRanger / STARsolo / Alevin-fry) → per-strand read ingestion with
  UMI dedup and `-F 3844` flag filtering → clip-site detection (`min_clip`, `min_purity`) → 25 bp single-linkage
  clustering with the read-weighted mode as cleavage position → **tier 1** (clip-supported, molecule-counted) and
  **tier 2** (coverage summits without clip support, written to a separate `low_confidence` file) → strand-aware
  internal-priming filter (−10..+30, ≥6 A or ≥70% A) → tiered gene-end/3′UTR annotation → PAS × cell matrix →
  `switch diff / length / trend / geneview`. **(b)** Why coverage shape alone is not enough on deep data: clip
  reads sit ~90–105 nt past the coverage-peak end (`07_curated_benchmark_report.md` §3), and the honest
  sensitivity budget drawn as two bars — 1.152% of CB reads carry a poly(A) clip, 28.77% of detected-gene atlas
  sites carry any clip read within 100 bp (`10` §R1). **(c)** The default as a decision tree with real PBMC v2
  counts: tier 1 IP-pass 167,565 → ≥2 distinct molecules **46,524** (the default), with the ≥1-molecule arm
  labelled "sensitivity arm" and tier 2 labelled "not switch-tested in this paper"; the pre-registration stamp
  printed on the panel (`0e27b1a`, 2026-08-21 01:18:59; original arms 02:59:36; v2 re-run 16:14). **(d)** One
  worked gene: a real `ema switch geneview` track from the v2 run, two cell types, proximal/distal PAS, per-site
  clip-molecule counts on the track.
- **Script / who builds it:** panels a–c are a hand-authored SVG (no clip art, no external assets, colourblind
  palette per `figures/README.md`) drafted by me and approved by the PI on wording; panel d is generated by a new
  `scripts/manuscript_figures/fig01_overview.py` that shells out to `ema switch geneview` on the v2 run and
  writes `fig01_overview_paneld*.tsv`; final assembly SVG → PDF/PNG in the same script. **Hard rule: every number
  printed on the schematic cites `19`, `16`, `13` or `10` — no illustrative or rounded-for-beauty numbers.**
- **Status: NOT YET BUILT.** Blocks submission.
- **Caveats that travel with it:** the tier label reflects gene assignment, not 3′-end proximity (`09` §5);
  panel b's clip rate is CellRanger 3.0.0 / 10x v3 / 91 bp R2 and any pipeline that trims poly(A) before
  alignment destroys the evidence entirely (`10` §R2).

---

**Fig 2 — `fig02_accuracy` (was `final_benchmark`) — The default output is a low-false-positive PAS set.** *[executed 2026-09-02 as `fig2_accuracy`]*
- **Proves:** L1. Precision default P@100 0.7062 / 0.7450 / 0.7572, gate PASS on all three; the highest de novo
  precision in the panel on both datasets; recall lower than polyApipe's AT OUR OWN CALL COUNT but higher at polyApipe's own N (25 §2); the tier decomposition that explains
  where the precision comes from (IP arm: tier 2 0.0558, tier 1 ≥1 mol 0.3520, default 0.7062).
- **Supports:** C1, C2, C5, C6.
- **Script:** `scripts/manuscript_figures/final_benchmark.py` → rename to `fig02_accuracy.py`. *[executed: `fig2_accuracy.py`]*
- **Status: DONE on v2** (regenerated 2026-08-21 21:25; `05_figure_index.md` entry written; caption sidecar
  current). Needs only the rename.
- **Caveats that travel with it (verbatim from `figures/final_benchmark.caption.md`, now `figures/fig2_accuracy.caption.md`):** atlas-agreement precision
  is agreement with a curated atlas, not ground truth — atlas-novel true sites count as false positives; the
  recall denominator is the detected-gene atlas and full-atlas recall (0.089 PBMC; 0.098 / 0.100 mice) must
  appear in the caption; single PBMC donor and single CellRanger BAM, and the two mice are one study and one
  chemistry; the non-IP ≥2-molecule file (`pas_tier1_ge2mol_noIP_POSTHOC.bed`) is not the default and is not
  shown; scUTRquant is catalog-based and is drawn but not ranked.

---

**Fig 3 — `fig03_tradeoff` (was `trade_reproducibility`) — There is a trade surface, not a ranking.** *[executed 2026-09-02 as `fig3_tradeoff`]*
- **Proves:** C5, C6, C8. The molecule-threshold sweep (≥1 → ≥10) moves PBMC precision 0.706 → 0.941 while
  R_det falls 0.175 → 0.083, and **only the ≥2 point was pre-registered**; the de novo precision ordering is
  unchanged at a 10 bp matching window (P@10/P@100 = 0.74 for our default and 0.77 for polyApipe against 0.46
  SCAPTURE / 0.33 Sierra / 0.24 scAPAtrap), so the precision lead is not a loose-window artefact; the two testis
  mice agree on 0.776–0.779 of default sites at 100 bp (0.729–0.737 at 25 bp) against a ≤0.010 chance level, and
  the ≥2-molecule threshold buys 12.5–14.3 points of replicate agreement over the ≥1-molecule arm; compute falls
  293.7 GB → 12.53 GB peak RSS on the PBMC BAM.
- **Script:** `trade_reproducibility.py` → `fig03_tradeoff.py`. *[executed: `fig3_tradeoff.py`]*
- **Status: DONE** (2026-08-21 20:29, caption current). Rename only.
- **Caveats:** the ≥1/≥3/≥5/≥10 points are descriptive and were never gated — never present ≥10 as a proposal;
  panel c compares a v2 PeakATail arm against v1-era competitor measurements (no competitor was re-run on v2);
  our replicate agreement is above polyApipe and above our own shipped caller but **below Sierra (0.790/0.834)
  and scAPAtrap (0.884/0.793)**, and scAPAtrap's lead is measured after its own `reducePeaks` depth cleaning;
  PBMC is a single donor, so no human biological replicate exists in this benchmark; peak RSS is the production
  number and wall time must be quoted with the four-arm concurrency disclosed (~28–35 min; 27m43s uncontended).

---

**Fig 4 — `fig04_calibration` (was `fdr_calibration_v2`) — A switch test you can trust, and five you cannot.** *[executed 2026-09-02 as `fig4_calibration`; the script now writes the previously missing caption sidecar]*
- **Proves:** L3 / C3. Six arms × 21 runs, 20 shared label permutations, 0 failures; only Fisher cells-mode with
  marker pre-selection off controls FDR (3.0% null p < 0.05, 0/20 null runs with any q < 0.05 hit); the three
  mechanisms are named and each is separately demonstrated.
- **Script:** `fdr_calibration_v2.py` → `fig04_calibration.py`. *[executed: `fig4_calibration.py`]*
- **Status: DONE** (2026-08-21 02:30). Rename only. **One addition worth its cost:** overlay the final-caller
  re-validation from `18` (null p < 0.05 3.169% / 3.166%, 0 q < 0.05 hits in all 30 null BH families) as a
  fifth marker in panel a, so the calibration is not resting on a single run of a superseded caller.
  **Take that overlay from the v2 spermatogenesis re-run (gap 2), not from the current `18` record:** those two
  values are v1 code (`4efeb125`), and §7 item 4 forbids quoting a v1 spermatogenesis number once v2 is
  verified — and forbids a v1/v2 mix inside one panel, which is exactly what an overlay onto a v2 figure would
  be. If gap 2 slips, ship Fig 4 without the overlay rather than with a v1 marker.
- **Caveats:** single mouse, single tissue, very large true stage effects, so TRUE hit counts are upper bounds on
  detectability and are not comparable to null-run hit counts; B0 is **conservative, not exactly calibrated**;
  marker-on and marker-off are different tests, not subsets, because pre-selection also changes the Fisher gene
  denominator; the label-permutation null tests exchangeability under the global null only; two-sided KS versus
  uniform rejects for every arm (discrete Fisher mass at p = 1) and is reported, not gated; "calibrated" is our
  operational rule, pre-stated in the figure script, not a standard.

---

**Fig 5 — `fig05_cohort` (was `laughney_switches`) — Replicated cell-type APA switches in a tumour cohort.** *[executed 2026-09-02 as **Fig 6** `fig6_cohort` — the approved map's amendment (b): spermatogenesis holds the Fig 5 slot]*
- **Proves:** L4 / C4. The replication funnel against the mean of 10 label-shuffle nulls; per-cell-type-pair
  yields; the effect-size distribution against the pre-registered \|Δproportion\| ≥ 0.1 floor; patient support;
  and the honesty panel (genomic context of the replicated sites: 58.6% own-gene 3′UTR, 96.8% inside a
  same-strand gene body).
- **Script:** `laughney_switches.py` → `fig05_cohort.py`, re-run with `LAUGHNEY_SWITCHES_VERSION=v2_code`. *[executed: `fig6_cohort.py`; v2_code is the default]*
- **Status: NEEDS REGENERATION from the Stage-3 v2 chain.** [checked here] as of 2026-08-21 22:31 the chain at
  `/mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v3/` has finished `replication/primary_noMetBone/pas_K2`
  (22:29) and is running `pas_K3`; `gene_K2` follows. The figure is regenerable within minutes of the chain's
  `replication/ALL.done`, but **every count in it changes**, so the caption, `20`, this file's L4 and the claim
  ladder are all downstream of one verification pass.
- **Caveats:** the null resolves only to empirical p ≤ 0.091 (the 10-permutation floor) — say "none replicated in
  10 label-shuffle nulls", **never** an FDR or a false-discovery estimate; the unit of replication is the
  **patient**, not the GSM; the MetBone exclusion is pre-registered but its stated premise (a 0.015% clip rate)
  was a sampling artefact of the caller reading only the first 200,000 CB reads, and the 13-patient sensitivity
  run differs by 0.4%; **no ranked list of named top-switch genes may be published** until PAS→gene
  re-assignment (issue #99: 33% of the top-30 gene rows name a spanning or readthrough model, and 6.4% of
  3′UTR PAS sit in a different gene's 3′UTR than the one assigned); pairs tested in more patients replicate
  more, so the per-pair panel tracks cohort composition as much as biology; the between-donor same-cell-type
  specificity null of `02_validation_plan.md` is a *different*, still-unrun control and must not be conflated
  with the label-shuffle null.

---

**Fig 6 — `fig06_novelfunnel` (was `trusted_novel_funnel`) — The pre-registered de novo test, and its failure.** *[executed 2026-09-02 as **`figS7_novelfunnel`** — the approved map's amendment (c) demotes the figure to the supplement; D4's prominence argument is honoured in text (R6, abstract, Author Summary)]*
- **Proves:** C11 (the negative result) and, positively, C2's boundary: the funnel shows exactly where long-read
  concordance is lost. 46,524 → hexamer 35,712 (76.7%) → atlas-novel 7,259 (15.6%) → strong hexamer 4,329; the
  per-stage concordance 0.765 → 0.811 → **0.485** → 0.502 shows the loss is at the atlas-novelty stage and that
  the hexamer adds nothing among atlas-novel sites; 69.4% of the trusted-novel set is intronic and 27.0% sits
  within 25 bp of a Kinnex internal-priming decoy terminus against 7.6% for the atlas-known reference.
- **Script:** `trusted_novel_funnel.py` → `fig06_novelfunnel.py`. *[executed: `figS7_novelfunnel.py`]*
- **Status: DONE on v2** (2026-08-21 17:56). Rename only.
- **Caveats:** this is a negative result and **no site may be called "trusted"**; the truth set comes from other
  donors, so it is site-level truth (robustness: GEM-X 0.5154, pooled 0.5816 at 25 bp; the atlas-known
  hexamer-pass complement scores **0.8941 on the same x3p truth**, so donor mismatch does not explain the gap —
  the GEM-X/pooled atlas-known calibration exists only for **v1** (0.9200 / 0.9436) and must never be reported
  as a v2 value); **the pre-registration (`13` §3) defined novelty as ≥100 bp from any PolyASite 2.0 *or
  PolyA_DB* site, and PolyA_DB could not be used (hg19 only), so "atlas-novel" here means PolyASite 2.0 only —
  a deviation that must be printed with the negative result**; the 50 bp
  (0.6688) and 100 bp (0.7143) values — both on the **pooled x3p+GEM-X** truth, not on x3p — are **not** the
  pre-registered window and the 100 bp value must never be
  quoted as meeting the target; panel d's stratification (3′UTR 0.769, ≥5 molecules 0.692–0.728) was chosen
  after seeing the data and is exploratory, not a definition; single PBMC donor; the caller's IP rule is looser
  than Kinnex's and that gap is a named follow-up (`--ip-rule kinnex`, issue #95 addendum), not a fix applied here.

---

### 4.3 Supplementary figures — every legacy asset placed

| # | Stem | Content | Source / script | Decision & status |
|---|---|---|---|---|
| S1 | `figS01_spermatogenesis` | Testis 3′UTR control on the final caller: per-gene monotone shortening 31.5% / 30.2% vs shuffle nulls 17.5% / 21.0% (z 10.5 / 5.8); composition-controlled per-cell distal residual falls at every step (Cliff's δ SPC–ES 0.56 / 0.61); cross-mouse replication 6,049 / 9,469 / 11,263 PAS, 0 in all 15 null pairings, per-gene ρ 0.655; boxed negatives. | `spermatogenesis_final.py`, `18` | **KEEP, REGENERATE on `9dfdefb`** (currently the `4efeb125` record). Demoted from main per D3; it gets its own Results paragraph now, and an abstract line **only once the v2 re-run is verified** — the §3 draft quotes no spermatogenesis number because `18` is v1 code. |
| S2 | `figS02_datasets` | Dataset table made visual: cells, depth, chemistry, aligner, QC, **and the per-BAM poly(A) clip rate** (the number that predicts whether the method works at all). | NEW script; sources = the four v2 `run_config.json` / `run_manifest.json`, the cohort run manifest, `04_datasets.md`, `10` §R2 | **NEW — replaces the retired `cohort_qc`.** Cheap, and it pre-empts the first methods question. |
| S3 | `figS03_peakqc` | Peak-call QC on the v2 arms: width, per-gene PAS counts, tier composition, molecule-count distributions, and the **cleavage-offset** analysis (clip sites ~90–105 nt past coverage-peak ends) that motivated the redesign. | NEW script; `07` §3, v2 `pas_support.tsv` | **NEW.** This is where "72.1% of tier-1 sites are single-molecule (no-IP arm; 72.2% on the IP arm)" gets its distribution instead of a bare percentage — §7 item 8 requires the arm to be named every time. |
| S4 | `figS04_nulldesign` | Why the benchmark is built the way it is: reference density, **point vs interval matching** (the single biggest inflator), match-class composition, and the fact that the older 18.4M-entry dump could not rank five peak-calling strategies (spread 0.004) where the curated point-site reference can (0.355–0.492; `07` §"the arms now separate"). | `null_control.py` + the surviving methodological residue of `benchmark_strategies.py` | **REBUILD as one figure against the curated PolyASite 2.0 point reference.** Two mandatory preconditions. (i) The two scripts currently compute *different* nulls (0.753 vs 0.615 for the same arm) and both cannot ship — reconcile to the `score_tool.py` gene-body-shuffled null used everywhere else. (ii) **The published density figures (23/33/45/61/73% of the genome within 50/100/200/500/1,000 bp; "one entry every 168 bp") are properties of the retired 18,432,135-entry PolyASite 3.0 dump, not of the curated PolyASite 2.0 point reference** (`05` §5, binding caveat). They must be recomputed on the curated reference before they can appear anywhere; carrying them over would put a dump-derived number back in the paper, which §7 forbids. |
| S5 | `figS05_longread_rerank` | Per-tool long-read re-ranking: Spearman ρ between the atlas and Kinnex precision orderings (0.90 at 7 of 8 arms, 0.80 at one) and the genuine / internal-priming / unsupported decomposition per tool — including the finding that Sierra's atlas precision was ~half internal priming. | `kinnex_truth_validation.py` (**quarantined**) → rebuild from `results/reliability/trusted_novel_final_v2_pbmc/` + competitor call sets | **REGENERATE** under the binding corrections of `11_verifier_corrections.md` (never "every stringency"; last only at the ≥5/≥20 support thresholds; and the truth-set support is an **alignment-record** count, *not* a de-duplicated UMI count — de-duplicate before the rebuild or relabel the axis). **Until the rebuild these ρ values and the Sierra decomposition are quarantined and not quotable** (`05` QUARANTINE NOTE). Strengthens R3; does not block. |
| S6 | `figS06_motif` | Sequence architecture: canonical hexamer at −40..−5 of the clip-seeded cleavage point in 76.7% of default sites (35,712 / 46,524); A-fraction profile; hexamer adds +0.046 to long-read concordance overall and **nothing** among atlas-novel sites. | `motif_validation.py` → rebuild re-anchored on cleavage points | **REBUILD, and it is blocked on one missing measurement:** a shuffled-position null for the re-anchored −40..−5 window does not exist. Cheap to compute; without it the panel cannot ship. |
| S7 | `figS07_clustering` | The **dropped** claim, stated honestly: PAS-profile clustering recovers GEX cell types (PBMC AMI 0.708 / ARI 0.502; Laughney median AMI 0.662 / ARI 0.463 in 17/17 samples) **and the ablation that killed it** — collapsing 275,370 sites to 14,891 gene totals gives AMI 0.698, so the site resolution adds nothing. | `clustering_concordance.py` + the ablation | **REBUILD with the ablation as a panel.** A figure of the positive result without the ablation on the same axes must never ship. Offered as an optional utility for label-free exploration, nothing more. |
| S8 | `figS08_calibration_extended` | Calibration detail: per-stage-pair and per-expression-stratum null rates (2.6–4.0% everywhere), permutation-calibrated q, `count_mode` reads vs cells, the NB dispersion-floor diagnostic, runtime of the permutation harness. | `fdr_calibration_v2*.tsv` (already written) | **NEW script over existing verified TSVs.** Low cost. |
| S9 | `figS09_gatehistory` | The pre-registration timeline as a figure: the original two-sided gate, the stale Stage-2 failure, the post-hoc sweep that informed the ≥2-molecule threshold (disclosed, not used as a result), the pre-registered default at 01:18:59, the arms at 02:59:36, the v2 re-run at 16:14, and the outcome of every gate on every arm. | NEW; `12` CORRECTION, `13`, `15`, `19` | **NEW.** This is the figure that makes the pre-registration auditable rather than asserted. Mostly a timeline plus table T3. |
| S10 | `figS10_ipfilter` | Internal priming: the rule (−10..+30, ≥6 A or ≥70% A) vs Kinnex's (+1..+18, ≥12/18 A); the minus-strand window defect and its correction (+ strand 0 changes; − strand 5,599 newly removed / 11,217 newly kept on PBMC); the filter's cost/benefit (+5.3 pp P for −1.2 pp R_det, mouse 1); the residual 27.0% vs 7.6% decoy rate. | NEW; `19` §1, `16` §v2, `08_ipfilter_live_check.md` | **NEW.** Reviewer-bait in the good sense: it shows a bug found by our own verifier, fixed, re-run and re-scored, with the gate unchanged. |
| S11 | `figS11_compute` | Compute scaling and the resolved Stage-1d item: 293.7 GB → 12.53 GB peak RSS and 3:45:53 → 34:37 on the PBMC BAM (v2 wall time measured with all four arms concurrent; 27m43s uncontended — §7 item 5 applies to every panel of this figure); the 17-library cohort run 9:06:26 → 1:06:26 at identical output scale (505,197 unified PAS); competitor runtimes with their retry/skipped-stage caveats. | `trade_reproducibility_compute.tsv` + `19` §3/§5 | **NEW script over existing verified TSVs** (panel d of Fig 3 is the summary; this is the detail). |

**Retired outright** (do not regenerate, do not cite as figures; their surviving findings live in the text or in
S4): `cohort_qc` (its headline is a tautology of the upstream TIER filter and its 16,500-site "cohort PAS set"
is a 17/17 intersection that no analysis uses; replaced by S2); `parameter_sweep` (it sweeps annotation-trim
knobs on the pre-clip-seeded lambda_gradient cohort run and its real finding — Leiden resolution dominates
cluster count — belongs to the dropped clustering claim; the parameter sweep that matters now is the
molecule-threshold sweep in Fig 3a); `benchmark_strategies` (superseded by curated-atlas scoring; its one
durable lesson moves into S4); `benchmark_curated`, `benchmark_headtohead`, `benchmark_tools_running`,
`pbmc_novelty`, `fdr_calibration` (v1, mis-keyed matrix), `spermatogenesis_control` (quarantined v1).

**Supplementary tables** (carried from `01`, unchanged in intent): T1 19-tool feature comparison with an
evidence-type column; T2 per-tool × per-dataset × per-window precision / R_det / full-atlas recall / F1_det with
nulls; T3 all pre-registered gates with values and pass/fail; T4 replicated switches per contrast with per-patient
direction calls (**PAS-id keyed, no gene-name column until issue #99**); T5 the 7,259 atlas-novel sites with
molecule support, hexamer tier, feature class, IP status and Kinnex outcome — **no "trusted" label**; T6 PDUI /
length per gene × cell type; T7 software versions, parameters, frozen commit hashes and manifest digests per run.

### 4.4 The rename pass (execute once, after the two v2 regenerations land)

1. Freeze: Stage-3 v2 verified (Fig 5) and spermatogenesis v2 verified (S1). **Do not rename before this** — a
   rename plus a number change in the same pass makes verification ambiguous.
2. For each figure: `git mv` the script, change its `NAME` constant, its module docstring's output paths, the
   caption header line, and any TSV names embedded in the caption text (each script writes its own caption, so
   these are the only four places a stem appears).
3. Re-run every script from scratch under `export LC_ALL=C` and confirm: the new PNG/PDF/caption/TSV set is
   complete, the old stems appear nowhere under `manuscript/figures/` or `results/figures/manuscript/`, and
   `05_figure_index.md` + `figures/README.md` + `FIGURE_MAP.tsv` all resolve.
4. Grep the manuscript tree for every old stem and fix the references in `01`, `05`, `11`, `14`, `16`, `18`,
   `19`, `20` and this file. A figure referenced by a name that no longer generates is the failure mode this
   convention exists to prevent.

---

## 5. The claim ladder

Ordered strongest → weakest. "Effect size" is the number the sentence carries; "Control" is the measurement
that rules out the obvious alternative explanation. The verbatim caveat sentence for each claim is in §5.1.

| # | Claim, as it will appear | Evidence | Effect size | Control that kills the obvious alternative |
|---|---|---|---|---|
| **C1** | The pre-registered precision-first default is a low-false-positive PAS set on both datasets. | `19` §1 (FIXED) | P@100 **0.7062** / **0.7450** / **0.7572**; gate ≥ 0.50 PASS ×3; R_det 0.1754 / 0.2048 / 0.2080 | *Alternative: reference density.* 3-seed gene-body-shuffled nulls 0.0217 / 0.0127–0.0150, i.e. 33–55× below; strand-matched point matching, not interval matching; every value reproduced to six decimals by an independent `bedtools closest`/`window` pipeline. |
| **C2** | The default set is corroborated by truth that owes nothing to the atlas. | `16` §v2 (FIXED) + v2 REPORT | **0.7647** [0.7608–0.7685] within 25 bp of a Kinnex 3′ end at ≥5 UMI (35,575 / 46,524) | *Alternative: the long-read truth is dense enough to hit anything.* Gene-body-shuffled null 0.0072 (106.5×); the atlas-known hexamer-pass complement scores 0.8941 on the same truth, calibrating the ceiling; GEM-X and pooled truths reproduce the ordering. |
| **C3** | Of six differential-APA test configurations, exactly one controls FDR; the shipped defaults did not. | `14` (SOUND) | 3.0% null p < 0.05 and **0/20** null runs with any q < 0.05 hit, vs 20.3% / 13.0% / 24.7% for the shipped arms | *Alternative: our null is too weak.* 20 label permutations shared identically across all six arms, full pipeline per run, 6 × 21 runs with 0 failures; the three mechanisms are each demonstrated separately (restricting B0's nulls to the permuted-label marker set reproduces the anti-conservativeness at 17.4%). |
| **C4** | Reported switches replicate across independent biological units and do not appear under label-shuffle nulls. | `20` (FIXED, **v1 code — v2 pending**) *[superseded: `20` §v2 verified SOUND 2026-08-22 on the v2_code chain — the record is now **15,942 / 2,128,711** replicate, 15,212 over the effect floor, **0 in each of 10 nulls** (Fig 6 `fig6_cohort`, `fig6_cohort_funnel.tsv`); the row's 14,480 / 1,977,134 / 13,826 are the v1-code values]* + `18` (FIXED) | Laughney: **14,480 / 1,977,134** replicate (0.73%), 13,826 over the effect floor, **0 in each of 10 nulls**. Testis: 6,049 / 9,469 / 11,263 per pair, 99.7–99.8% sign agreement, **0 in all 15 null pairings**, ρ 0.655 | *Alternative: replication is an artefact of shared cells, shared labels or a shared PAS space.* One cohort PAS space so ids are comparable by construction; unit = patient with a two-GSM patient counted once; any opposite-direction patient vetoes; the nulls run through the identical pipeline with patient-wise label shuffles. |
| **C5** | Evidence type, not peak-calling sophistication, separates single-cell PAS callers. | `09`, `19` §1, `15` §3 | Same BAM, same framework: coverage-only 0.118 → clip-seeded default **0.7062**; panel ordering (SCAPTURE 0.652 > polyApipe 0.380 > Sierra 0.256 > scAPAtrap 0.100) tracks evidence type, not algorithm class | *Alternatives tested and excluded:* not a coordinate offset (relaxing to 200 bp lifts 0.118 → 0.166 while polyApipe goes 0.380 → 0.410; `09` §2); not gene-distance bias (all tools ≥95% gene-proximal; `12` §5, verified). **The third alternative — "not a threshold problem" — is PROVISIONAL and may not be printed**: the top-22k/36k/106k UMI-rank triple (0.118–0.120) was computed on a mis-keyed matrix, and even the re-keyed triple is flagged provisional in `09` §"Why PeakATail underperforms" and forbidden by §7 item 11 of this file. Gap 13 re-runs it on v2; until then C5 stands on two excluded alternatives, not three. |
| **C6** | Molecule support is the knob that buys reliability, and it is a trade, not a free lunch. | `19` §4, Fig 3 caption | ≥1 → ≥10 molecules moves P 0.352 → 0.941 as R_det falls 0.269 → 0.083; the ≥2 default discards **72.1%** of tier-1 sites (no-IP arm; 72.2% on the IP arm); replicate agreement 0.636–0.651 → 0.776–0.779 | *Alternative: the threshold is post-hoc tuning.* Only the ≥2 point was pre-registered (`13` §1, commit `0e27b1a` 01:18:59, arms 02:59:36); the sweep points are labelled descriptive on the figure and the post-hoc sweep that informed the threshold is disclosed (`12` CORRECTION). |
| **C7** | The precision lead is not an artefact of a loose matching window. | Fig 3 caption (from verified per-tool score TSVs) | P@10 / P@100 = **0.74** (our default) and 0.77 (polyApipe) vs 0.46 / 0.33 / 0.24 for SCAPTURE / Sierra / scAPAtrap; the de novo ordering is unchanged at 10 bp | *Alternative: 100 bp flatters a caller whose sites are systematically offset.* The ordering is recomputed at 10 / 25 / 50 / 100 bp from the same call sets and the same scorer. |
| **C8** | The method is now cheap enough to be run by anyone. | `19` §3, §5 | PBMC BAM **12.53 GB** peak RSS, ~28–35 min (27m43s uncontended), from 293.7 GB / 3:45:53; the 17-library cohort 9:06:26 → **1:06:26** at identical output scale (505,197 unified PAS) | *Alternative: the speed-up changed the results.* The `--peak-workers` change was verified byte-identical where it promised to be (19/19 files `cmp`-equal on the no-IP arm) and the cohort produced the same 505,197 unified PAS. |
| **C9** | 3′UTR usage shifts progressively proximal along spermatocyte → round → elongating spermatid. | `18` (FIXED, **v1 code — v2 pending**) | 31.5% / 30.2% of depth-guarded genes shorten monotonically vs nulls 17.5% / 21.0% (z 10.5 / 5.8); shortening exceeds lengthening 489 vs 325 and 361 vs 296; per-cell composition-controlled residual falls at every step (Cliff's δ 0.56 / 0.61) | *Alternative: label noise or composition.* 20 label shuffles per mouse; GEX-derived labels orthogonal to the PAS matrix under test; the per-cell index is composition-controlled; both mice independently. |
| **C10** | Single-cell PAS benchmarking must be reported per dataset, and F1 only with its truth set named. | `09`, `05` §8 | Our own shipped caller spans P 0.118 (PBMC) → 0.369 (testis) → 0.447 (Laughney); F1 rank reorders every tool but one across truth sets differing 53-fold in size | *Alternative: one number would do.* The same fixed call sets are scored against two truth sets and two datasets by one scorer with per-dataset denominators. |
| **C11** | **NEGATIVE.** The pre-registered "trusted de novo PAS" definition did not reach its target. | `16` §v2 (FIXED) | **0.4853** [0.474–0.497] vs a target of 0.70; enrichment 65–179× over a positional null; 69.4% intronic; 27.0% at Kinnex IP decoys vs 7.6% | *Alternative: the fixes, the donor, or the truth set caused it.* Re-running the unchanged definition after #96/#97 moved it **down** (0.5211 → 0.4853); missed under GEM-X (0.5154) and pooled (0.5816) truths at the pre-registered window; on the **x3p** truth the atlas-known hexamer-pass complement scores **0.8941** (hexamer-strong 0.9031), so donor mismatch is excluded. **No v2 atlas-known calibration exists under the GEM-X or pooled truths** — that pairing was measured only on v1 (0.9200 / 0.9436) and must not be quoted as a v2 number. |

### 5.1 The caveat sentence that must travel with each claim (verbatim)

- **C1** — "Precision here is atlas-agreement precision at 100 bp against curated PolyASite 2.0 representative
  sites, not ground-truth precision: atlas-novel true sites count as false positives, the recall denominator is
  the atlas restricted to genes detected in the dataset (full-atlas recall 0.089 and 0.098/0.100), and the
  evaluation rests on one PBMC donor and two mice from a single study and chemistry."
- **C2** — "The Kinnex truth set comes from different donors, so this is site-level and not donor-matched
  truth, and long-read support is counted in **alignment records at the terminus**: `11_verifier_corrections.md`
  established that the truth set's '≥5 UMI' threshold is a record count, with ~10.5% of records repeating a
  (CB, UMI) at the same terminus, so the support thresholds are slightly optimistic. **The v2 truth BEDs are
  md5-identical to v1 and were never de-duplicated** (`REPORT_PROVISIONAL.md` §0; `scripts/benchmark_tools/kinnex/build_truth2.sh`
  counts termini rows), so this caveat is outstanding, not fixed." [verifier correction: the previous wording,
  "deduplicated UMIs", asserted the opposite of the source disclosure.]
- **C3** — "The calibrated arm is conservative rather than exactly calibrated, the assessment rests on one mouse
  and one tissue with very large true stage effects, and 'calibrated' is our pre-stated operational rule
  (null p < 0.05 ≤ 7% AND null p < 0.01 ≤ 1.5% AND null q < 0.05 ≤ 5% AND ≤ 25% of null runs with any hit),
  not a community standard." [verifier correction: the p < 0.01 clause of `14`'s pre-registered rule was
  missing from this sentence.]
- **C4** — "None of the tested features replicated in any of ten label-shuffle nulls; because ten permutations
  floor the empirical p at 0.091, this is a null control and not a false-discovery-rate estimate, and no ranked
  list of named genes is reported because PAS-to-gene assignment in overlapping loci is under revision."
- **C5** — "This is a six-tool, two-dataset comparison under one scorer and one atlas-agreement metric; it ranks
  the evidence available to each caller on these data, not the callers' best achievable performance."
- **C6** — "Only the ≥2-molecule operating point was pre-registered; the remaining points on the sweep are
  descriptive, and the threshold itself was informed by a post-hoc sweep on a superseded run, which we disclose."
- **C7** — "Window sensitivity was computed on the same fixed call sets, so it tests the sharpness of each
  caller's coordinates and not its ability to find sites the atlas does not contain."
- **C8** — "Peak resident memory is the production number; the quoted wall times were measured with all four
  arms running concurrently on one machine, with the uncontended single-run time given alongside."
- **C9** — "The gradient is a per-gene, equal-weight statement across spermatocyte → round → elongating
  spermatid: the per-gene medians are not monotone, the across-gene UMI-weighted index reverses at the final
  step because protamine transcripts at ceiling PDUI dominate elongating-spermatid UMIs, the spermatogonia step
  is excluded as fragile, and a 16-gene literature panel does not reproduce on this caller. Monotone
  *lengthening* is also above its own null, so the direction-specific evidence is the shortening-minus-
  lengthening excess, which is strong in mouse 1 (binomial p 9.9e-9) but only **modest in mouse 2**
  (p 0.0125; one of 40 null draws reached the same excess in absolute value)." [verifier addition: `18`
  "Cautions bound to the claim" and `05`'s "mouse-2 direction excess modest" were missing from this sentence.]
- **C10** — "F1 is reported only with its truth set and denominator named, because F1 rank is not stable across
  truth sets of different sizes."
- **C11** — "We report this as a negative result: no site in this paper is called 'trusted novel', and the
  stricter strata that would reach the target (3′UTR-restricted, ≥5 molecules) were chosen after seeing the data
  and must be pre-registered and validated on held-out truth before they mean anything."

---

## 6. The negative-results section

The paper carries **four** negative or self-refuting results. They are not scattered through the Discussion;
they get one Results subsection ("R6 — What did not work"), one main figure (Fig 6), and one paragraph in the
Discussion. This is the paragraph:

> **We pre-registered four claims and their acceptance gates before the final run existed, and we report the
> outcome of all four, including the one that failed.** Three cleared. The fourth did not: our pre-registered
> definition of a *trusted de novo* poly(A) site — clip-supported, ≥2 distinct molecules, internal-priming
> filtered, a canonical hexamer at −40 to −5 nt, and ≥100 bp from any curated atlas site (PolyASite 2.0; the
> pre-registration also named PolyA_DB, which is hg19-only and could not be used) — was required to show
> ≥70% concordance with poly(A)-verified long-read 3′ ends within 25 bp. It reached 48.5% (95% CI 47.4–49.7).
> Re-running the unchanged definition after two caller bug fixes moved the value down rather than up
> (52.1% → 48.5%), so the failure is a property of the definition, not of the bugs; and it is missed under every
> alternative truth set at the pre-registered window, while the atlas-known sites scored on the same x3p truth reach
> 89.4%, which excludes donor mismatch as the explanation. The residual is diagnosable — 69.4% of these sites
> are intronic and 27.0% sit within 25 bp of a long-read internal-priming decoy against 7.6% of atlas-known
> sites — but a diagnosis chosen after seeing the data is a hypothesis, not a result, so we report the funnel
> and make no trusted-novel claim. Three further results ran against our own expectations. **Clustering:**
> cells cluster on poly(A)-site profiles into the expression-derived cell types (adjusted mutual information
> 0.708), but collapsing 275,370 sites to 14,891 per-gene totals recovers them just as well (0.698), so the
> signal is 3′-end expression re-encoded rather than isoform choice, and we withdraw the claim and ship the
> capability as a convenience. **Literature panel:** a 16-gene panel of spermatogenesis 3′UTR shortening does
> not reproduce on the final caller — 4 genes shorten in both mice, 5 lengthen, 5 are discordant, 2 are
> uninformative — so we drop that leg and keep only the per-gene, null-referenced statement it was meant to
> illustrate. **Our own defaults:** the differential-APA settings we shipped do not control the false-discovery
> rate, with 13.0% to 24.7% of label-shuffled tests reaching p < 0.05, a failure we found only because we built
> the calibration harness that our field's tools, ours included, did not have. We regard the disclosure of all
> four outcomes as a result in its own right. A pre-registration that can only be confirmed is decoration; a
> literature in which no single-cell APA caller has ever reported a failed acceptance gate is a literature in
> which acceptance gates are not being set.

**Framing rules for the write-up.** (i) Each negative appears *first* in its own section, before the positive
result it constrains — the reader must never discover a limitation after being sold the claim it limits.
(ii) Every negative names the alternative explanation that was tested and excluded, so it reads as a
measurement and not as an apology. (iii) The word "unfortunately" does not appear in the manuscript.
(iv) The four negatives are listed in the Author Summary / Key Points box, not only in the body.

---

## 7. MUST-NOT-CLAIM (carried forward from `01` §(f), extended)

Carried forward unchanged: accuracy leadership beyond what the gates show; a single cross-dataset accuracy
number; "precision" without "atlas-agreement"; F1 superiority; F1 without its truth set named; the non-IP
≥2-molecule PBMC file as the default; "trusted novel PAS" as a positive claim; recall parity; reproducibility
leadership; clustering novelty or cell identity from PAS usage alone; FDR-controlled discovery from the shipped
Fisher/NB q-values; the retired dump-based precision figure, anything scored against the 18.4M-entry dump, or
atlas snap-rate as precision; stale Stage-2 numbers; quarantined spermatogenesis/Kinnex numbers; monotone
four-stage spermatogenesis shortening; Read1 cleavage precision, sequence-model filtering or spatial support.

**Added here:**

1. **No named top-switch gene list, anywhere, until issue #99 lands** — 33% of the top-30 gene-level rows name a
   spanning or readthrough model, and seven of them have zero assigned PAS cohort-wide. Switch results are
   PAS-id keyed until re-assignment.
2. **Never call the 10-permutation label-shuffle result an FDR, a q-value or a false-positive rate.** The only
   permitted wording is "none replicated in ten label-shuffle nulls (empirical p ≤ 0.091, the ten-permutation
   floor)".
3. **Never repeat the MetBone "no clip evidence" justification.** The exclusion stands because it was
   pre-registered; its stated premise was a clip-rate sampling artefact (issue #99), and the sensitivity run
   including it differs by 0.4%.
4. **Never quote a Stage-3 or spermatogenesis number from v1 code once the v2 chain is verified**, and never mix
   v1 and v2 numbers in one sentence, table or figure.
5. **Never quote wall time without the concurrency disclosure**; peak RSS is the production number.
6. **Never quote the 50 bp (0.6688) or 100 bp (0.7143) trusted-novel values as meeting the 70% target** — the
   pre-registered window is 25 bp, and both of those values are on the **pooled x3p+GEM-X** truth, not on the
   pre-registered x3p truth, so they differ from the metric in two ways at once.
7. **Never quote "104 M poly(A)-verified molecules"** (`13` §3 correction note): the verified truth sets are the
   x3p and GEM-X 3′-end point BEDs at ≥5/20/100/500 UMI.
8. **Never call the Kinnex support threshold a de-duplicated UMI count.** `11` established that the truth
   set's "UMI" support is an **alignment-record** count (~10.5% of records repeat a (CB, UMI) at the same
   terminus), and the correction was never applied: the v2 truth BEDs are md5-identical to v1 and
   `build_truth2.sh` counts termini rows. Every "≥5 UMI" in `16`, `19` §4 and Fig 6 inherits this — write
   "≥5 long-read records at the terminus", or write "≥5 UMI" with the record-count caveat attached. Also never
   quote 72.5% single-molecule (72.1% no-IP arm / 72.2% IP arm; always name the arm).
9. **Never present the ≥5 or ≥10 molecule operating points as a recommendation** — they are the shape of the
   trade surface, and only ≥2 was pre-registered.
10. **Never call the between-donor same-cell-type specificity null "done"** — it is a different, still-unrun
    control from the label-shuffle null (`02_validation_plan.md`).
11. **Never quote the top-N UMI-ranking triple (0.118–0.120)** as a current result until it is re-run on the v2
    caller; the published triple was computed on a mis-keyed matrix and even the re-keyed one is flagged
    provisional in `09`.
12. **Never claim the internal-priming problem is solved** — the corrected filter is strand-aware but still
    looser than the long-read rule, and 27.0% of atlas-novel calls sit at decoys.

---

## 8. What is still missing for submission

**Gap 0 — the limitation this list must not hide: the recall deficit is real, and nothing below closes it.**
The default's detected-gene recall (0.1754 / 0.2048 / 0.2080) is below polyApipe's (0.199 / 0.250 / 0.251) on
**both** datasets — −12% PBMC, −18% mouse (`final_benchmark.caption.md`, now `fig2_accuracy.caption.md`) — and far below scAPAtrap's 0.300 on
PBMC; full-atlas recall is 0.089 / 0.098 / 0.100. **The clip-rate budget is not an answer to this, and must
never be offered as one.** The 1.152% clip rate and the 28.77% ceiling (`10` §R1) bound *every* clip-seeded
caller equally, and polyApipe also requires non-templated poly(A) soft clips (`09` §4 "the discriminator is
evidence type"): it gets **closer** to that shared ceiling than our default does, and our own ≥1-molecule arm
reaches R_det 0.2685 (`19` §2), i.e. nearly the ceiling, at P 0.3520. So the deficit at the default is the
price of the pre-registered ≥2-molecule threshold, which discards 72.1% of tier-1 sites (no-IP arm; 72.2% IP
arm) — **a chosen operating point on the trade surface of Fig 3, not a physical limit of the evidence.** Write
it that way in the Abstract (already done: "below polyApipe's"), in Results R1 *before* the precision claim,
and in the Discussion. Closing it would require recovering true single-molecule sites without their false
positives — tool work that nobody has scheduled and that no gap below delivers. This item cannot be closed
before submission; it can only be stated honestly, and it is listed first so that no one mistakes its absence
from the P0/P1/P2 tables for its absence from the paper.

**P0 — blocks submission.**

| # | Gap | What closes it | Cost | Note |
|---|---|---|---|---|
| 1 | **Stage-3 v2 numbers** (Fig 5, `20`, L4, C4). | The chain at `stage3_laughney_v3/` finishes `replication/ALL.done`; re-run `laughney_switches.py` (now `fig6_cohort.py`) with `LAUGHNEY_SWITCHES_VERSION=v2_code`; write `20` §v2; one adversarial verification pass. | Chain: finishing tonight ([checked here] `pas_K2` done 22:29, `pas_K3` running at 22:31, `gene_K2` to follow). Figure: minutes. **Verification: ~half a day of a verifier's attention — that is the real cost.** | Every count in Fig 5 and C4 changes. Nothing downstream may be frozen before this. |
| 2 | **Spermatogenesis v2 re-run** (S1, `18`, C9). | Re-run both mice on `9dfdefb` from a frozen worktree, redo `switch diff/length` + the replication analysis, regenerate `spermatogenesis_final.py` (now `fig5_spermatogenesis.py`), verify. | 2 caller runs (~15 min each at v2 speeds) + switch/length + figure + verification ≈ 1 day. | The current record is `4efeb125`. Also the chance to re-check whether the literature panel behaves differently — if it does, that changes D3. |
| 3 | **Figure 1 does not exist.** | Draft panels a–c as hand-authored SVG; build panel d from a real `ema switch geneview` track on the v2 run; assemble in `fig01_overview.py`. | 1–2 days including PI review of every printed number. | The only main figure with no asset at all. |
| 4 | **Data and code availability.** | `Project_PeakATail` blockers: placeholder authors in `CITATION.cff`, no `LICENSE` file, benchmark pin `18678ef` living on the unpushed `biolab-manuscript` branch; repo must go public before Zenodo can mint a DOI; Docker image; bioRxiv preprint. | ~1 day of packaging plus PI decisions on authorship/licence. | Genome Biology and NAR both make this mandatory at submission, not at acceptance. |
| 5 | **The re-anchored hexamer null** (S6, and the motif sentence in Results). | Shuffled-position null for the canonical hexamer in −40..−5 of clip-seeded cleavage points. | ~1 hour of compute with `scripts/reliability/trusted_novel_pas.py` machinery. | Without it, "76.7% carry a canonical hexamer" has no baseline and must not be presented as evidence. |
| 6 | **The rename pass and caption/script sync** (§4.4). | Execute after gaps 1 and 2. | Half a day. | Cosmetic in isolation, load-bearing for reproducibility: a figure whose script no longer generates it is an unverifiable figure. |
| 7 | **The null definitions in `benchmark_strategies` and `null_control` disagree** (0.753 vs 0.615 for the same arm). | Reconcile to the `score_tool.py` gene-body-shuffled null and rebuild as S4. | ~half a day. | Two different nulls cannot ship in one paper. |

**P1 — strengthens materially; do these if there is any time.**

| # | Gap | What closes it | Cost | Why it matters |
|---|---|---|---|---|
| 8 | **A second human dataset.** | Run the pre-registered default, unchanged and un-tuned, on a second public 10x BAM (pbmc4k/8k are already listed in `04_datasets.md`) and report P@100 / R_det with the same scorer and nulls. | ~1 h of compute per BAM plus scoring. | **The single highest value-per-hour item in this list.** Our precision claim currently rests on one human donor and one CellRanger BAM; a second dataset with no re-tuning converts "we tuned to one dataset" from an open question into a measured answer, and it is the first thing a Genome Biology reviewer will ask for. |
| 9 | **Competitors have never been run on v2 code paths / one machine snapshot.** | Re-measure runtime and peak RSS for every tool in one sitting, with retries and skipped stages disclosed. | ~1 day of wall time, mostly unattended. | Fig 3d currently mixes a v2 PeakATail point with v1-era competitor runs, and scAPAtrap's 4.07 h is a resumed run. |
| 10 | **Long-read re-ranking (S5) is quarantined.** | Regenerate under `11`'s binding corrections from the v2 truth sets. | ~half a day. | It is the strongest answer to "your precision metric is circular", and it carries the Sierra internal-priming finding, which is the paper's most useful gift to the field. |
| 11 | **Between-donor same-cell-type specificity null.** | Contrast the same cell type between donors in the Laughney normals through the identical switch pipeline; report the false-positive yield. | ~1 day. | This is the control a statistical reviewer will ask for that the label-shuffle null does not provide. |
| 12 | **PAS→gene re-assignment (issue #99).** | Tool fix, then re-derive T4 and the biology paragraph. | Tool-side; not ours to schedule. | Without it there is no named biology, which is exactly why Nature Communications is off the top list (§9). |
| 13 | **The top-N ranking test is stale.** | Re-run the UMI-rank top-N precision test on the v2 caller. | ~1 h. | It is one of the three alternative explanations excluded under C5, and it currently rests on a provisional re-keyed measurement. |

**P2 — reviewer bait; only if the schedule allows, and never at the cost of P0.**

14. **Downsampling / saturation curve** — precision and recall vs read depth on the PBMC BAM; pre-empts "does
    your recall just reflect depth?" (~half a day).
15. **Chemistry and aligner dependence of the clip rate** — measure it across the BAMs we already have
    (CellRanger v3, STARsolo v2) and state plainly that poly(A)-trimming pipelines destroy the evidence (~2 h;
    partially in S2 already).
16. **A `--ip-rule kinnex` arm** as an exploratory row, pre-registered as v2 and validated on the held-out
    GEM-X truth (~1 day). **This must never become a blocker** (D8): the negative result is the result.
17. **Literature sweep of 2025–26 single-cell PAS callers using clip evidence or long-read truth** — positioning
    risk, not novelty risk (~half a day).
18. **A simulation arm** with known ground-truth PAS. Cheap to over-claim from, so only as a supplement, and
    only if it is generated by a model we did not fit to our own caller.

---

## 9. Journal fit — revised ranking

The ranking in `01` was written when Fig 6 was expected to deliver a biological program. It cannot, until
issue #99. The revision below reflects what the results actually turned out to be: a precision-leading,
recall-lagging caller; a strong benchmark; a genuine calibration contribution; and a disclosed failed gate.

| Rank | Journal | Change from `01` | Fit |
|---|---|---|---|
| **1** | **Genome Biology** (Method) | unchanged | Still the best fit and the fit improved: pre-registered gates, gene-body-shuffled nulls, long-read truth, replicate reproducibility and a disclosed negative result are exactly what the Method article type rewards, and Genome Biology published APAeval. 6–8 figures, structured abstract ≤ 350 words, mandatory availability statement. |
| **2** | **Genome Research** (Methods/Resource) | ↑ from 3 | Accuracy-centred work with nulls and long-read truth is squarely in its lineage. Degrades from a Genome Biology submission with no restructuring. |
| **3** | **Briefings in Bioinformatics** | ↑ from 5 | The six-tool, two-dataset, one-scorer benchmark with nulls and a reproducibility metric is now one of the paper's two strongest assets, and BiB's centre of gravity is exactly that. The safest high-visibility home if reviewers weigh the recall cost heavily. |
| 4 | **NAR / NAR Genomics & Bioinformatics** | ↓ from 4 (unchanged rank, weaker case) | The right readership for an atlas-referenced benchmark with the circularity handled explicitly, but NAR methods papers tend to want a resource or webserver flavour we do not have. |
| 5 | **Bioinformatics** (OUP) | unchanged | Viable if descoped to "precision-first caller + calibration"; the length limit forces the benchmark into supplements. |
| 6 | **PLoS Computational Biology** | unchanged | Safety net; its open-science ethos fits the pre-registration and the disclosed failure better than any other venue on this list. |
| — | **Nature Communications** | ↓ **off the top list** | Without a named, orthogonally supported biological program it is a methods/benchmark paper, and NatComms takes those only when the benchmark changes field practice. Reconsider only if issue #99 lands *and* the v2 cohort run yields a replicated per-cell-type program with orthogonal support. |

**What a reviewer at the top three will demand that we do not yet have.**

- **Genome Biology.** (i) A **second human dataset** scored with the unchanged pre-registered default — the
  precision claim currently rests on one donor and one BAM (gap 8). (ii) Some engagement with **APAeval**: either
  run their harness or state precisely why our scorer differs (point vs interval matching, per-dataset
  denominators, gene-body nulls). (iii) A **complete availability package** at submission — Zenodo DOI, Docker
  image, and the manifests as supplement (gap 4). (iv) They will press on **recall**, and the clip-rate
  budget is **not** a sufficient answer — see **Gap 0**. The 1.152% clip rate and the 28.77% ceiling bound
  every clip-seeded caller equally; polyApipe uses the same poly(A)-clip evidence and gets *closer* to that
  ceiling than our default does; and our own ≥1-molecule arm reaches R_det 0.2685. The answer that survives
  review, stated up front in the Results rather than defended in the rebuttal, is: the ceiling bounds what any
  clip-seeded caller can reach, the gap to polyApipe is the price of the pre-registered ≥2-molecule threshold,
  the sensitivity arm sits on the same figure, and we regard the false-positive reduction as worth it — as a
  **choice**, never as a physical constraint.
- **Genome Research.** (i) A **depth/saturation analysis** (gap 14) — they will not accept a recall number
  without knowing how it moves with depth. (ii) **Chemistry dependence** of the evidence (gap 15), because the
  method's premise is a read feature. (iii) Ideally one dataset with **bulk 3′-seq truth** in the same tissue;
  we have none, and the honest answer is that Kinnex long reads substitute for it at the site level.
- **Briefings in Bioinformatics.** (i) **More tools** — they will want SCAPE, scDaPars, MAAPER or scAPA added to
  the panel; we have the harness, and each addition is bounded work. (ii) **More datasets**, same argument as
  gap 8. (iii) A **released benchmark harness**, not just released tool code: `score_tool.py`, the null
  construction and the per-tool environment fixes packaged so a reader can re-rank the panel themselves. That
  last one is cheap for us and is probably worth doing regardless of venue.

---

## 10. Provenance of this document

Verified documents read in full for this pass: `01` (v0.4), `05`, `09`, `11`, `13`, `14`, `15` §1/§3, `16`,
`18`, `19`, `20`, `github/PUBLISHED.md`, and the four generated caption sidecars
(`final_benchmark`, `trade_reproducibility`, `laughney_switches` and, via `05`, `trusted_novel_funnel` and
`spermatogenesis_final`; since 2026-09-02 these stems are `fig2_accuracy`, `fig3_tradeoff`, `fig6_cohort`,
`figS7_novelfunnel` and `fig5_spermatogenesis`).

Four things I re-derived from primary files rather than quoting **[checked here]**:

1. **Call-set sizes.** `wc -l` on the v2 BEDs gives 46,544 (PBMC default), 26,263 / 26,533 (mice) and 167,629
   (PBMC tier-1 IP) — each 7–64 lines above the scored n of `19` §1 (46,524 / 26,255 / 26,526 / 167,565). The
   difference is the non-primary contigs the scorer excludes, exactly as `16` §v2 describes for PBMC
   ("46,544 sites" → "46,524 on primary contigs"). No discrepancy; worth one sentence in Methods so that a
   reader comparing the shipped BED to the reported n is not confused.
2. **The trusted-novel funnel counts.** `call_primary/trusted_novel.bed` has 7,260 lines of which the first is a
   `#chrom` header → **7,259 sites**, and `trusted_novel.strong.bed` → **4,329**; both match `16` §v2 exactly.
3. **The Kinnex concordance of the full default.** `REPORT_PROVISIONAL.md` row `00_input`/t5:
   35,575 / 46,524 = **0.7647** [0.7608–0.7685], null 0.00718, enrichment 106.5× — the source of the abstract's
   "76.5%".
4. **Stage-3 v2 chain state at 22:31 on 2026-08-21:** `replication/primary_noMetBone/pas_K2` complete (22:29),
   `pas_K3` running, `gene_K2` pending; the chain driver (`run_replication_chain.sh`) and the
   `replication_filter.py` job for `pas_K3` were both live. The tree was inspected read-only (`ls`, `ps`) and
   nothing in it was touched.

**Everything in this file that is a number is either from a verified document or from that list. Nothing here is
a new measurement, and every architectural decision above is reversible by the PI — but each one is a decision,
not an option.**

---

## 11. Adversarial verification pass (2026-08-21, verdict FIXED)

An adversarial verifier re-derived every number in this file against its cited source (`19` + the
`final_v2_verify/VERIFIED_v2.md` table, `16` §v2 + `trusted_novel_final_v2_pbmc/REPORT_PROVISIONAL.md` and
`verifier_crosscheck/ADVERSARIAL_VERIFY_v2.txt`, `14`, `18`, `20` + `laughney_switches.caption.md` (now `fig6_cohort.caption.md`), `15` §3,
`13`, `12`, `11`, `10`, `09`, `07`, `05`, `04`, `02`, `01`, the four generated caption sidecars) and, where the
file claimed `[checked here]`, against the primary files. Twelve corrections were applied **into the text
above**; each tightens a caveat and none softens one. Nothing else in the architecture was changed.

1. **Mouse null range.** `0.0129–0.0150` → **`0.0127–0.0150`** (L1 and C1). `VERIFIED_v2.md` §1 gives mouse-1
   seeds 0.0149 / 0.0129 / 0.0139 and mouse-2 seeds 0.0137 / **0.0127** / 0.0150.
2. **Kinnex support unit** (C2 caveat, §7 item 8, S5). The file said support is "counted in deduplicated UMIs".
   `11` established the opposite — the truth set's "UMI" support is an **alignment-record** count, ~10.5% of
   records repeating a (CB, UMI) at the same terminus — and the correction was never applied: the v2 truth
   BEDs are md5-identical to v1 and `scripts/benchmark_tools/kinnex/build_truth2.sh` counts termini rows.
3. **Atlas-known calibration under GEM-X / pooled truths** (C11, Fig 6 caveats, §6). "0.89–0.90 under those
   same truths" does not exist for v2: only the **x3p** value (0.8941 hexamer-pass, 0.9031 hexamer-strong) was
   measured. The GEM-X / pooled atlas-known calibration is a **v1** measurement (0.9200 / 0.9436).
4. **Abstract: the ≥5-UMI threshold restored.** 76.5% is 0.7647 at ≥5 UMI; the same set scores 0.5584 at ≥20
   and 0.2807 at ≥100, so the threshold is load-bearing.
5. **Abstract: "MIT-licensed" removed.** Gap 4 records that the repository has no `LICENSE` file and that the
   licence is an open PI decision; `01` v0.4 asserts MIT, but an abstract may not state a licence the
   repository does not carry. Recount after both edits: **200 words**.
6. **D3 / S1 vs the drafted abstract.** D3 promised the spermatogenesis control "a line in the abstract" and §3
   contained none. Reconciled in favour of §3: `18` is v1 code and gap 2 supersedes it, so the abstract carries
   no spermatogenesis number until the v2 re-run lands — the same rule the file already applies to Stage-3.
7. **C9 caveat completed** with `18`'s own bound cautions: monotone *lengthening* is also above its null, so the
   direction-specific evidence is the shortening-minus-lengthening excess, strong in mouse 1 (p 9.9e-9) but
   only **modest in mouse 2** (p 0.0125; 1 of 40 null draws matched it in absolute value).
8. **C5's third alternative demoted.** "Not a threshold problem (0.118–0.120)" is the exact number §7 item 11
   of this file forbids quoting; C5 now stands on two excluded alternatives and cites `09` §2 and `12` §5 for
   the two that survive.
9. **Fig 4's proposed overlay** must come from the v2 spermatogenesis re-run: `18` is v1 code, and §7 item 4
   forbids both a v1 spermatogenesis number after v2 and a v1/v2 mix inside one panel.
10. **S4's reference-density figures flagged as dump-derived.** 23/33/45/61/73% and "one entry every 168 bp"
    are properties of the retired 18,432,135-entry PolyASite 3.0 dump (`05` §5), not of the curated PolyASite
    2.0 point reference S4 is to be rebuilt against; carrying them over would return a dump-derived number to
    the paper.
11. **Two of this file's own must-not-claim rules were being broken inside it** and are now obeyed: S3 quoted
    72.1% without naming the arm (item 8) and S11 quoted wall times without the concurrency disclosure (item 5).
    Separately, C3's caveat sentence was missing the `null p < 0.01 ≤ 1.5%` clause of `14`'s pre-registered
    calibration rule, and Fig 6's caveat list was missing the pre-registration deviation that PolyA_DB
    (hg19-only) could not be used, so "atlas-novel" means PolyASite 2.0 only (`13` §3 vs `16`); both added.
12. **Gap 0 added, and the Genome Biology recall answer rewritten.** The gap list contained no entry for the
    recall deficit, and §9 offered the 1.152% clip rate as "the answer" to it. That answer does not hold:
    polyApipe uses the same clip evidence (`09` §4) and reaches R_det 0.199 against our default's 0.1754, and
    our own ≥1-molecule arm reaches 0.2685 — so the budget bounds the **ceiling** (≈0.29), while the gap to
    polyApipe is the price of the pre-registered ≥2-molecule threshold, a chosen operating point. The same
    over-reach in §1's "honest fifth result" was tightened.

**Checked against primary files and confirmed, not changed:** the four `[checked here]` items of §10 (BED line
counts 46,544 / 26,263 / 26,533 / 167,629; `trusted_novel.bed` 7,260 lines − 1 header = 7,259 and
`trusted_novel.strong.bed` 4,330 − 1 = 4,329; `REPORT_PROVISIONAL.md` row `00_input`/t5 = 35,575 / 46,524 =
0.7647, null 0.00718, 106.5×); `cohort_qc.py` is indeed the only figure script without a `NAME` constant
(checked across all 19); the figure scheme places all 19 existing stems exactly once (6 main incl. one new,
5 existing in supplements, 9 retired, `benchmark_strategies` retired-as-figure with its lesson moved to S4);
the trusted-novel funnel, decoy (0.2699 / 0.0763), intronic (0.6939), enrichment (65.5–179.1×), stratification
(0.7687 3′UTR; 0.6920 / 0.7281 support bins) and per-stage concordances (0.7647 → 0.8110 → 0.4853 → 0.5024);
`14`'s six arms and mechanisms; `20`'s counts and the 58.6% own-gene 3′UTR figure (which correctly supersedes
`20` disclosure 1's 62.5% via the generated caption); `12` §5's ≥95% gene-proximal; `07`'s 0.355–0.492 and the
90–105 nt offset; `04`'s pbmc4k/8k; `02` line 110 / V11's still-unrun between-donor control.

**One unresolved inconsistency outside this file, for the rename pass to reconcile:** `final_benchmark.caption.md`
panel b says "SCAPTURE is mouse 1 only", but `15` §3 reports SCAPTURE on both mice (0.694 / 0.672) and
`results/benchmark_tools/gse104556/scapture/DONE.mouse2.ok` records that mouse-2 **PAS calling completed**
(24,076 points, P@100 0.672) and only the per-cell `PASquant` step failed. The value quoted in §1 and §4.2 of
this file is correct; the caption line is stale and must be fixed when Fig 2 is renamed.
**[RESOLVED 2026-09-02: fixed at the rename — `fig2_accuracy.caption.md` panel b now records that the mouse-2 site-level run completed (24,076 points, P@100 0.672, only PASquant failed) while the plot still shows mouse 1 only; the on-figure note reads "SCAPTURE: mouse 1 plotted (mouse-2 sites scored 0.672; 15 §3)".]**
 **[CORRECTED 2026-08-21: the genome-wide clip rate is 0.573%, not 1.15%. The 1.15% figure came from the caller's head-sampling QC estimator, which reads only the first 200,000 CB reads of a coordinate-sorted BAM (the head of chr1) and reports 2.26% where the truth is 0.536%; see results/algo_headroom/VERIFY/ and manuscript/23_algorithm_roadmap.md §3 Step 0. The conclusion the figure supported — that clip evidence is a small, highly specific channel — is unchanged and in fact strengthened.]**

# PeakATail manuscript figures — round 1 index

Generated 2026-08-12/13. Scripts: `/mnt/ssd1/Projects/PeakATail_wd/scripts/manuscript_figures/<name>.py`
(each re-runnable end to end). Outputs: `<name>.png` (300 dpi) + `<name>.tsv` (every plotted value).
Every figure below was independently re-verified from source by an adversarial pass; the numbers quoted
here are the post-verification numbers. Where a verdict says a wording still has to change, that wording
is given in **Must travel with it**.

All five figures re-run under `LC_ALL=C` (the machine locale is tr_TR; GNU `sort` misorders BED files
without it and `bedtools closest` then returns silently wrong distances).

---

## 1. `cohort_qc.png` / `.tsv`

**Finding.** PeakATail retains 55,422 cells across 17 lung samples and defines a cohort PAS set of
16,500 sites — which is the strict 17/17-sample intersection, 9.0% of the 183,510 PAS detected in at
least one sample.

**Headline numbers.**
- 17 samples, 55,422 cells retained (1,036–4,923/sample, median 3,418; 4.8-fold spread). 98.1% of
  barcode-filtered cells survive PAS QC (min 87.2%).
- 16,500 cohort PAS, 100% gene-assigned; 8,039 genes; 3,507 (43.6%) multi-PAS (max 32/gene).
- PAS peak width median 811 bp, IQR 544–1,276, max 15,195 bp.
- Called-peak universe 301,860 → 184,295 (61.1%) survive the TIER_1/TIER_2 annotation filter;
  38.9% dropped. Within the cohort set, 16,409 TIER_1 / 91 TIER_2.
- 3'UTR length median 1,092 bp (single-PAS genes) vs 1,856 bp (multi-PAS genes).
- Stage composition: normal 4, primary stage I–II 7, primary stage IV 1, metastasis 5.

**Verdict: FIXED.** All counts recomputed from source (`clustering_stats.json`, `pasbed.bed`,
`annotatedpas.bed`, `annotated_pas_ids.tsv`) with 0 mismatches. The 17/17-intersection identity was
proved, not assumed. Three defects fixed in the figure (undisclosed selection rule; circular TIER panel;
"universe" mislabel).

**Must travel with it.**
- Do not write "99.4% of PAS fall in annotated 3'UTRs" unconditioned — it is a tautology of the upstream
  TIER_1+TIER_2 filter. Correct: *of the PAS retained by the TIER_1/TIER_2 filter, 99.4% are TIER_1; that
  filter discarded 38.9% of the 301,860 called peaks.* TIER_1 is 60.3% of called peaks.
- Do not write that the 3,507 multi-PAS genes are "the gene set on which APA/PDUI analysis can act."
  The switch analysis runs in the unified 115,450-PAS space; `pasbed.bed` covers 16,405/115,450 = 14% of
  it. 3,507 is a floor from a 17/17 subset.
- Per-sample `final_pas` (11,946–36,232) exceeds 16,500 in 13/17 samples — never compare per-sample PAS
  counts to the cohort number.
- `filtered_cb.tsv` is pre-QC, not retained cells (plotted value is `final_cells`).
- Cell yield is confounded with clinical category (metastases 2,589 vs primary I–II 3,694 mean cells), so
  any met-vs-primary comparison is also a power comparison. Stage IV primary is n=1.
- PAS are broad intervals, not cleavage sites; wide peaks probably merge several true sites.
- Multi-PAS status is confounded by 3'UTR length — the APA-testable set is biased to long-3'UTR genes.

---

## 2. `clustering_concordance.png` / `.tsv`

**Finding.** Clustering cells on poly(A)-site profiles alone, with no gene-level expression matrix,
recovers the GEX-derived cell-type partition in all 17 samples at comparable, slightly finer granularity.

**Headline numbers.**
- Median AMI 0.662 (0.367–0.780); median ARI 0.463 (0.200–0.671); AMI > ARI in 17/17.
- Worst GSM3516664-MetBone AMI 0.367 / ARI 0.226; best GSM3516672-StageIB 0.780 / ARI 0.661.
- Against raw GEX Leiden instead of cell types: median AMI 0.657 / ARI 0.456 (essentially identical).
- Granularity: mean 14.1 PAS clusters vs 12.5 GEX cell types (above y=x in 12/17); Spearman rho 0.58, p=0.014.
- Coverage (GEX-labelled / PAS cells) median 0.498, range 0.274–0.869; 9/17 below 0.500.
- No stage effect detectable: Kruskal-Wallis AMI H=3.84 p=0.147, ARI H=3.32 p=0.190.
- Coverage does not drive agreement: rho(coverage, AMI) = −0.34 (p=0.18).

**Verdict: FIXED.** AMI/ARI recomputed from the raw per-cell labels in all 34 source `.h5ad` files (not
from the figure TSV); all nine headline numbers reproduced exactly. No fabrication. Five presentation
defects fixed (feature-space mislabel, truncated x-axis, missing caveats on-figure, 2-dp rounding that
contradicted a panel title, undisclosed cluster-count provenance).

**Must travel with it.**
- **Say "poly(A)-site profiles", not "usage".** The feature space is per-site *counts*, TF-IDF-transformed
  with an LSI embedding — abundance features, not within-gene usage fractions. A PAS count proxies its
  gene's expression, so agreement with GEX cell types is partly expected by construction. This figure
  cannot separate isoform-choice signal from expression-level signal.
- The reference partition is **not ground truth**: marker-signature labels assigned per GEX Leiden cluster.
  Report as agreement with an automated GEX partition.
- Quote ARI 0.46 as the conservative number; AMI exceeds ARI in 17/17 because PAS splits GEX clusters.
- 13–73% (median 50%) of PAS cells carry no GEX label and are entirely unassessed. `matched_cells ==
  gex_cells` in all 17 — the GEX side is always limiting.
- Both cluster counts are Leiden-resolution-dependent; "comparable granularity" is partly a parameter choice.
- Stage: write "no stage effect was detectable at this n" (n=4/7/5, stage IV primary n=1 excluded), not
  "holds uniformly across stage". The descriptive stage ordering is not significant.
- The 17 samples come from a smaller number of donors (normal/tumour/met sets likely paired); the 17
  values are not independent and no patient/batch effect was tested.

---

## 3. `parameter_sweep.png` / `.tsv` (+ `parameter_sweep_per_dataset.tsv`, 221 rows)

**Finding.** PAS and cell recovery are insensitive to the annotation-trim parameters swept here, while
cluster granularity is governed almost entirely by one knob — Leiden resolution.

**Headline numbers.**
- Leiden resolution 0.5 → 2.0: **+112.5%** median clusters (11 → 23; per-dataset +71% to +200%).
- Every annotation knob is ≥8× smaller: `include_extended` off→on +13.2% PAS; `max_gene_distance`
  1,000→10,000 (a 10× window) only +6.4%; `utr_multiplier` 1.5→3.0 +1.1%.
- `n_neighbors` 15→50: −14.3% clusters. tfidf→libsize: +5.9% median clusters.
- Cell recovery effectively invariant: median 0.14%, worst case 3.8% across all 13 branches.
- Median `final_pas`: d1000_ext 26,596 → d10000_ext 28,417; pipeline default (ext=off) 24,136.
- Median PAS retention 23.1% of called peaks at the pipeline default (16.8–28.6%).
- Between-dataset IQR (~10,800 PAS) is ~6× the entire `max_gene_distance` effect (~1,800 PAS).

**Verdict: FIXED.** All 221 branch×dataset rows re-extracted independently with 0 mismatches; every
panel-d contrast audited from `branch_manifest.json` as a clean single-knob span; the claim that
clustering branches reuse a bit-identical matrix confirmed by md5 of `preprocessed.h5ad`. Seven defects
fixed, chiefly a false symmetry claim and median-only bars that hid per-dataset movement.

**Must travel with it.**
- **Orthogonality is one-directional.** Clustering knobs leave PAS/cell counts bit-identical (0/17
  datasets move, md5-confirmed). Annotation knobs are only *median*-neutral on cluster count: 7–9 of 17
  datasets shift by 1–2 clusters (−13% to +13%). Write "orthogonal in one direction by construction,
  median-neutral in the other."
- Scope the robustness claim to "insensitive to the annotation-trim parameters swept here." The filter
  thresholds (`min_read=1500`, `min_cells=3`, `min_pas_per_cell=50`) were held fixed in all 13 branches
  and are plausibly the larger lever — only ~23% of called peaks survive filtering.
- `n_clusters` is a granularity proxy only; no ARI between branches was computed, so nothing here says
  cell *assignments* are stable.
- The tfidf→libsize swap is stable on the median (+5.9%) but ranges −41.7% to +93.3% per dataset — do not
  quote its median alone.
- Confound: the five distance branches all have `include_extended=True`, the pipeline default does not.
  `A2_trim_default` is not the d5000 point of that curve (drawn detached on purpose).

---

## 4. `benchmark_strategies.png` / `.tsv`

**Finding.** Five peak-calling strategies differ ~2-fold in how many PAS they call but are
indistinguishable against the PolyASite reference (precision 0.996–1.000, full spread < 0.004), so this
atlas benchmark cannot rank them — it only shows that no strategy produces off-target calls.

**Headline numbers.**
- PAS called: sierra_iterative 43,035; lambda_gradient 22,633; lambda_poisson 21,474; lg + IP-filter
  20,609. sierra_iterative calls 1.90–2.00× the lambda methods. The internal-priming filter removes
  2,024 peaks (−8.9%).
- Precision @ ±50 bp: lg 0.9974, lg_ip_filter 0.9972, si 0.9964, lp 0.9959; all ≥0.9997 by ±5,000 bp.
- **Chance floor:** the same peaks shuffled genome-wide with widths preserved score 0.671–0.735 @ ±50 bp
  and ~0.928 @ ±5,000 bp. The 0.999 sits on a ~0.73 floor, not on 0.
- Reference density: atlas ±50 bp windows already tile 38.6% of the 3.088 Gb genome (92.7% @ ±5 kb).
- Atlas coverage ("recall") @ ±50 bp: si 0.0448, lg 0.0298, lg_ip_filter 0.0270, lp 0.0172.
- Called-peak width median: lg 640 bp, si 522 bp, lp 403 bp.

**Verdict: FIXED.** Every precision and recall numerator recomputed from source through an independent
pipeline (awk pad → `bedtools intersect -sorted -u`) — exact match on all cells; atlas coverage
re-derived with `bedtools merge`; the chance floor reproduced byte-for-byte on one seed and confirmed on
two fresh seeds (0.739, 0.740). Three defects fixed.

**Must travel with it.**
- **Duplicate arm.** `lg_annotate/pasbed.bed` and `lg_ip_off/pasbed.bed` are byte-identical
  (md5 `20cb4550…`). The "IP filter off" arm is not an independent run — the default *already* has the
  IP filter off. The only real internal-priming contrast is 20,609 vs 22,633.
- The reference is **18,432,135 PolyASite 3.0 3'-end clusters (mean 9 bp wide, 165 Mb total)**, not
  "single-nucleotide sites" — that earlier wording is wrong and must not reach the manuscript.
- Atlas coverage rises with **bases searched**, n × (peak width + 2 × cutoff) (Spearman rho = 1.00 at
  every cutoff), **not** with peak count (rho 0.80; lambda_poisson calls more peaks than lg_ip_filter yet
  covers less atlas). It is not sensitivity and a perfect lung PAS caller could not approach 1.0.
- Never quote precision ≈0.999 as an accuracy claim without the 0.67–0.74 random-placement floor.
- Effective window is wider than the nominal cutoff: `bedtools window` pads whole intervals (median
  403–640 bp), so "±50 bp" is really a ~0.5–0.7 kb search window.
- Matching is **strand-agnostic** (`metrics.py` uses `window(w=cutoff)` with no `-s`).
- **Do not report F1** from these JSONs — it is built on the uninterpretable recall term.
- No replicates: one cohort-level run per arm, no CIs on the arm curves.

---

## 5. `null_control.png` / `.tsv` (502 rows)

**Finding.** The 0.9986 precision headline is largely reference-density artifact, but PeakATail's calls
are still not explained by density: scored at the inferred cleavage site against matched shuffled
controls, real PAS separate from the null by ~0.54 and show a poly(A)-class signature the null cannot fake.

**Headline numbers.**
- Reference density: 23% / 33% / 45% / 61% / 73% of the genome lies within 50/100/200/500/1,000 bp of a
  strand-matched atlas entry (one entry every 168 bp).
- Interval matching (what the pipeline does), @100 bp: real 0.9985 vs null 0.615 (chromosome-shuffle) /
  0.749 (gene-body shuffle). At 1,000 bp the gap collapses to +0.087 — the benchmark is near-meaningless there.
- **3'-end point matching @100 bp: real 0.9940 vs null 0.348 / 0.458 (gap +0.536); @50 bp 0.9564 vs
  0.241 / 0.327.** This is the discriminating measurement.
- High-confidence subsets (point, @100 bp): TE/AL/EX classes (1,549,168 sites) real 0.485 vs null
  0.019/0.035; atlas average TPM ≥ 1 (40,243 sites) real 0.374 vs null 0.0013/0.0021.
- Match-class composition @100 bp — real: TE 0.310, AL 0.115, IN 0.483, UI 0.004. Null (gene-body):
  TE 0.022, IN 0.612, UI 0.128. Atlas background: TE 0.031, IN 0.570, UI 0.133. Real PAS are 10× enriched
  for terminal-exon and 32× depleted for upstream/intergenic; the null tracks background at 0.7–1.1×.
- Corrected recall @100 bp, interval mode: full atlas 0.0239 → within detected gene bodies 0.0444 →
  + TE/AL/EX 0.1779 → + TPM ≥ 1 **0.5534**. In point mode the same four are 0.005 / 0.010 / 0.038 / **0.290**.
- Replicate spread across 5 shuffles < 0.006.

**Verdict: FIXED.** Density, precision, class composition and restricted recall all reproduced exactly by
an independent implementation; an independent gene-body-constrained shuffle with its own RNG gave
0.466–0.472 (point) vs the figure's 0.458, i.e. the plotted null is if anything slightly conservative.
Six defects fixed, including a proven misattribution of strand-matching to the pipeline.

**Must travel with it.**
- **The reference is not the curated PolyASite 3.0 atlas.** 18,432,135 entries vs the published ~570k
  human clusters; column 8 (supporting protocols) is 1 for every row; 92% of entries are IN/DI/UI. It is
  an unfiltered dump. Either re-benchmark against the filtered cluster set or state the filter used —
  *every* precision number in the manuscript inherits this.
- **The pipeline does not strand-match.** `metrics.py` uses `window(w=cutoff)`, no `-s`; the
  strand-agnostic recomputation reproduces the pipeline's JSON exactly at all five cutoffs, the
  strand-matched one does not. This figure's strand-matched null is the better control, but it is *our*
  control, not the pipeline's procedure.
- **Interval vs point matching is the single biggest inflator and the pipeline uses the inflating
  version** — `closest -d` returns 0 whenever a peak merely *contains* an atlas site, so a 10 kb "PAS"
  scores as a perfect hit. Report the 3'-end point number (0.994 @100 bp), or report both.
- If quoting recall 55.3%, say "whole-peak interval matching"; the point-mode equivalent is 29.0%.
  The pipeline's own published full-atlas recall is 0.033, not the 0.024 computed here — do not put the
  two in the same table unexplained. Its F1 of 0.0632 is meaningless.
- Say the null's class composition is "within 0.7–1.1× of the atlas background", **not** "statistically
  indistinguishable" — chi-square 89.2, 5 df, p ≈ 1e-17.
- Call the TPM ≥ 1 subset "high-confidence (atlas average TPM ≥ 1)", not "expressed" — column 5 is
  PolyASite's cross-sample average, not expression in this cohort. The threshold is arbitrary
  (TPM ≥ 0.1 → 320,188 sites; ≥ 10 → 9,445).
- Restricting the reference to detected genes is mildly circular (which is why it alone lifts recall only
  0.024 → 0.044); the class and TPM filters do the real work and are not circular.
- Neither null models read pileup, so both are optimistic nulls for the caller.
- **Unreconciled between figures:** `benchmark_strategies.py` computes its null without `-s` and without
  `-chrom` (0.753 for lg_annotate @100 bp) where `null_control.py` gives 0.615. Both cannot ship with
  different null definitions.

---

## What these figures do and do not establish

**Established.**
1. The pipeline runs end to end on 17 samples and yields a QC-clean, reproducible cohort (fig 1), with
   PAS and cell recovery insensitive to the annotation-trim parameters swept (fig 3).
2. Called peaks are genuinely poly(A)-site-like and not explained by reference density: at the inferred
   cleavage site they beat width-, count-, chromosome- and strand-matched controls by ~0.54 in precision,
   and their PolyASite class composition is 10× enriched for terminal-exon and 32× depleted for
   upstream/intergenic sites — an orthogonal signature that distance cannot fake (fig 5).
3. PAS-space cell clustering carries real cell-type structure (median ARI 0.46, no sample near chance,
   17/17 samples) (fig 2).
4. Cluster granularity is a deliberate parameter choice (Leiden resolution, +112%), not an emergent
   property of the data (fig 3).

**Not established — the atlas-benchmark limitation.** The precision numbers in figs 4 and 5 cannot
support an accuracy claim in the form they are usually quoted. The reference is a low-stringency
18.4M-entry file whose ±100 bp windows tile a third of the genome; the pipeline matches strand-agnostically
and pads whole multi-hundred-bp intervals; randomly placed peaks of the same widths already score
0.62–0.75. Consequently: precision at ≥500 bp cutoffs is uninformative; "recall"/F1 against the full atlas
are uninterpretable; and the atlas benchmark **cannot rank the five peak-calling strategies** — they span
0.9959–0.9998, a spread of 0.004 (fig 4). Only the point-mode, null-referenced, class-composition analysis
of fig 5 carries evidential weight, and even it inherits the wrong reference file.

**Also not established.** (i) That the fig-2 concordance reflects *isoform choice* — the features are
per-site counts carrying gene abundance, so expression-level signal is not excluded. (ii) That agreement
is with true cell identity — the reference partition is automated marker-signature labels. (iii) Anything
about clinical stage — n=4/7/5 with stage IV primary n=1, and samples are donor-correlated. (iv) That
cluster *assignments* (not counts) are parameter-stable. (v) Robustness to the filtering stage, which was
never swept and removes ~77% of called peaks.

**Blocked on the per-celltype switch outputs.** The headline biology — differential poly(A)-site usage
between cell types and between normal/primary/metastasis — is not in this round at all. All 24
`SWITCH_CELLTYPE` tasks failed (exit 2) in the 2026-08-11/12 sweep: `pipeline/main.nf` line 215 builds
`"${h5.parent}/pasbed.bed"` from a staged single-segment path, so `h5.parent` is `null` and the command
requests the literal `null/pasbed.bed`. Diagnosis and the fix (point it at
`${cohort}/unified/multi_sample_merged.bed`, **not** `${cohort}/pasbed.bed`, which would silently drop 86%
of PAS) are in `/mnt/ssd1/Projects/PeakATail_wd/scripts/pipeline/rerun_switch_celltype_fix.md`.
Everything upstream is cached; only ~24 short tasks need to re-execute. Until they do, there is no
APA/PDUI result, no per-celltype switch figure, and fig 1's "APA-testable gene set" number (3,507) is a
placeholder from the wrong PAS space.

---

## Next analyses, in priority order

1. **Rerun the 24 `SWITCH_CELLTYPE` tasks** (one-line fix, everything upstream cached). This is the only
   item that adds new biology rather than new controls, and three claims across figs 1 and 2 stay
   provisional until it lands.
2. **Re-benchmark against the curated, filtered PolyASite 3.0 cluster set** (~570k clusters, or a stated
   TPM/protocol filter) and make point (3'-most base) matching the primary metric, with interval matching
   reported alongside. This retires the largest single caveat in the round.
3. **Reconcile the two null definitions** so figs 4 and 5 ship one control (recommend strand-matched,
   chromosome-preserving, gene-body-constrained), and drop recall/F1 against the full atlas everywhere.
4. **Orthogonal, atlas-independent validation of the calls:** poly(A) signal motif (AATAAA/ATTAAA)
   enrichment in the 40 bp upstream of the 3'-most base versus the matched shuffled controls, plus a
   genomic A-richness / internal-priming scan. This is the cleanest way to support an accuracy claim
   without the reference-density problem, and it also gives the internal-priming filter an evidence base
   beyond "removes 8.9% of peaks".
5. **Re-run fig 2's concordance on within-gene usage fractions** (per-gene-normalised PAS composition)
   instead of TF-IDF'd counts. If AMI/ARI survive, the isoform-choice claim becomes defensible; if they
   collapse, fig 2 must be reframed as an expression-proxy result.
6. **Sweep the filter thresholds** (`min_read`, `min_cells`, `min_pas_per_cell`) — the unswept stage that
   discards ~77% of called peaks and is plausibly a larger lever on `final_pas` than every annotation knob
   in fig 3 combined.
7. **Cluster-assignment stability across parameter branches** (ARI between branches, not just cluster
   counts), which is what "robust clustering" actually means and what fig 3 currently cannot say.
8. **Choose and justify the PAS set used for APA testing.** The 17/17 intersection (16,500) is one
   defensible extreme; the unified space (115,450) is another. A prevalence rule (detected in ≥ k of 17)
   with k chosen on a stated criterion is probably the right answer — and fig 1 needs whichever it is.
9. **Peak-width diagnostics:** test whether the wide peaks (cohort median 811 bp, max 15 kb) merge
   distinct cleavage sites — sub-peak decomposition, or per-peak read-3'-end profiles. Merged sites blunt
   proximal/distal separation and directly limit PDUI sensitivity.
10. **Donor and batch structure:** map the 17 samples to donors, run an integrated-cohort clustering, and
    test patient effects. Every per-sample n=17 test in this round assumes an independence the design
    probably does not have; this also determines whether any stage comparison is worth powering.

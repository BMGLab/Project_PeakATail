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

---

# Round 2 addendum — curated benchmark + motif validation (2026-08-18)

Full report: `07_curated_benchmark_report.md`. Both figures independently re-verified (adversarial
pass); numbers below are post-verification. Executes round-1 next-analyses items 2, 3 and 4.

## 6. `benchmark_curated.png` / `.tsv`

**Finding.** Re-benchmarked against the curated PolyASite 2.0 atlas (569,005 representative sites,
GRCh38.96) with strand-matched point-mode matching and one reconciled genic-shuffle null: precision at
100 bp is 0.36–0.49 across all five arms — ~16–23× the null — but recall is size-capped at 0.03–0.09,
so **every arm fails all three roadmap bars (P ≥ 0.70, R ≥ 0.60, F1 ≥ 0.65) at every cutoff 10–200 bp,
on both recall flavors and on the TES reference** (0 of 150 roadmap rows pass).

**Headline numbers.**
- Precision @100 bp vs curated atlas: lg 0.447, lp 0.492, si 0.355, lg_ip_filter 0.457, B1 0.475;
  max anywhere 0.557 (@200 bp). Genic-shuffle null 0.020–0.024 (5 arms × 3 seeds).
- Fold over null by cutoff: 25–42× @10 bp, 19–32× @50 bp, 16–23× @100 bp, 10–14× @200 bp — quote the
  fold with its cutoff (~20× is a 100 bp number).
- Recall @100 bp: full atlas 0.025–0.048; restricted to strand-matched bodies of the 14,851 detected
  genes (n=285,220) 0.050–0.095. F1 (flavor b) 0.091–0.150 @100 bp; max 0.210 (si @200 bp).
- TES secondary reference (73,439 protein-coding ends): precision @100 bp 0.204–0.313; never near a bar.
- Set sizes after scaffold drop: lg 22,629 / lp 21,470 / si 43,023 / lg_ip 20,605 / B1 16,497.

**Verdict: SOUND.** Every precision/recall/F1/null cell reproduced exactly by independent pipelines
(`bedtools window -sm`, fresh shuffle seed 777 landed inside the null band); point-mode strand
reduction proven formally and empirically (swapped 5'-base control collapses precision 0.447→0.130).

**Must travel with it.**
- **Not comparable to the retired 0.9986** — deliberately stricter matcher (strand-matched, point-mode,
  curated ~32×-sparser reference); `null_control` (fig 5) reconciles the matchers. The 0.9986 was a
  reference-density artifact and must not reach the manuscript.
- **Recall is arithmetically capped by set size**: per-arm ceilings (n_called/285,220) 0.058–0.151. The
  0.60 recall bar is unreachable by construction for a single-cohort assay vs a pan-tissue atlas;
  recall-b rises roughly linearly with n called. Flavor (a) also counts 397 unreachable scaffold sites.
- The honest claim is strong enrichment over null, not benchmark-passing accuracy.
- `lg_ip_off` is byte-identical to `lg_annotate` (md5 `20cb4550…`) — scored once, still not an
  independent arm.
- Latent (not triggered): script pipelines lack pipefail and results are cached under
  `.cache_benchmark_curated/`; fix before reuse.

## 7. `motif_validation.png` / `.tsv`

**Finding.** Atlas-independent sequence test: called peaks flank genuine polyadenylation sites, but the
peak 3' end is not the cleavage site — it stops ~90–105 nt short (consistent with 10x R2 coverage
ending before the poly(A) junction). AATAAA|ATTAAA sits within +0..+100 nt downstream in 37% of sites
vs 15% in genic-shuffle nulls, with full canonical architecture; the internal-priming filter halves
A-rich downstream tracts.

**Headline numbers.**
- AATAAA|ATTAAA, prescribed −40..−5 window (assumes peak end = cleavage): real 7.7% vs null 5.4% —
  near-null, nowhere near the 70–80% literature rate for true PAS.
- Re-anchored +0..+100 nt: real 37.1% (ip 37.4%) vs null 14.7% (2.5×); any-of-12 hexamers 65.8% vs 43.7%.
- AATAAA profile: >2× null across +29..+94 nt, mode +75 nt, bimodal sub-peaks +48/+75.
- Composition: A-fraction crests at 0.426 at +98 nt then cliffs to ~0.30 background (modal cleavage
  position); T overtakes A past +99.
- Internal-priming proxy (+10..+30): real 6.3% → ip-filtered **3.0%** (628/20,605 = 3.05%) vs null 2.2%
  — 51.7% relative reduction (run-A6 6.0→2.8%, fracA70 2.4→1.2%).
- n: real 22,629, ip 20,605, null 3 seeds × ~22,620. PolyASite not used anywhere.

**Verdict: SOUND.** Strand-correct point collapse verified formally, by an independent
genomic-coordinate route reproducing every count exactly, and by manual samtools revcomp check; all
on-figure numbers match the TSV.

**Must travel with it.**
- Write **3.0% (or 3.05%), not 3.1%**, for the ip-filtered internal-priming rate.
- Soften "downstream U-rich element": T over +101..+150 is only ~1 pp above null — say "T overtakes A
  past the modal cleavage position".
- The +90–105 nt implied-cleavage offset assumes canonical 15–30 nt AATAAA-to-cleavage spacing on the
  +75 profile mode; `pasbed.bed` carries no finer cleavage estimate.
- Even re-anchored, 37% is below curated-PAS rates (70–80%); per-site offset variability and residual
  false calls are not separable here. The +48/+75 bimodality (two offset populations) is uninvestigated.
- Prescribed-window numbers (near-null) are all still in the TSV — do not quote them as evidence
  against the calls without the anchor-offset explanation.
- Methods trap documented in the script: with `bedtools getfasta -s`, asymmetric windows break
  minus-strand index mapping (phantom −525 hump); symmetric ±160 nt windows are immune.

## 8. `benchmark_headtohead.png` / `.tsv` (+ `benchmark_consolidated.tsv`, 1,066 rows)

**Finding.** Six PAS callers on two datasets describe a trade surface, not a ranking: precision,
recall, replicate reproducibility and call count trade off against each other, every rank order
changes between the two datasets, and **PeakATail ranks LAST of the five de novo tools on human PBMC**
(F1@100 0.134) while placing second on mouse testis (0.258–0.260).

**Headline numbers.**
- PBMC F1@100 (detected-gene denominator): polyApipe **0.261** > SCAPTURE 0.199 > Sierra 0.177 >
  scAPAtrap 0.150 > **PeakATail 0.134**. scUTRquant 0.290 shown but not ranked (annotation-based).
- Testis F1@100: polyApipe **0.308/0.312** > PeakATail 0.258/0.260 ≈ scAPAtrap 0.248/0.248 ≈
  SCAPTURE 0.253 (mouse1 only) > Sierra 0.177/0.195. scUTRquant 0.389, unranked.
- SCAPTURE is the precision arm among the **de novo** tools on **both** datasets — P@100 0.652
  (human, 35,759 calls) and 0.694 (mouse1, 23,645) — ahead of Sierra (0.256 human / 0.509–0.581
  mouse). The fewest-calls arm is dataset-dependent: SCAPTURE on human (35,759 vs Sierra's 106,170),
  Sierra on mouse (21.3k–22.6k vs SCAPTURE's 23,645). scAPAtrap calls the most (787,138 human) and is
  least precise (0.100).
- Call counts span 21k–787k (37-fold) across tools, so precision and recall are only interpretable
  jointly (panel b).
- Replicate reproducibility @100 bp (mouse, both directions): scAPAtrap 0.79–0.88, Sierra 0.79–0.83,
  PeakATail 0.73–0.75, polyApipe 0.49–0.50. Chance level ≤0.009 (same query vs a genic-shuffled copy
  of the other replicate).
- **Depth stratification is the mechanism**: 62% of polyApipe's mouse calls are depth-1 singletons and
  they reproduce only **0.31** of the time, against **0.81** for its depth ≥2 calls. PeakATail 2%
  singletons, Sierra 1%, scAPAtrap 0% (its own `reducePeaks(min.cells=10, min.count=10)` removes them).
- Resources, clean runs: wall time PBMC SCAPTURE 12:13 > scAPAtrap 4:04 (resumed) > PeakATail 3:27 ≈
  polyApipe 3:26 > Sierra 2:43 > scUTRquant 0:30. Peak RSS PeakATail **100.7 GiB** on the human BAM vs
  4–23 GiB for every other tool/dataset — the largest single resource gap in the benchmark.
- Every real precision clears the 3-seed genic-shuffle null (~0.022 human, ~0.014 mouse), but by very
  different margins: 4.4–36× on human (scAPAtrap 4.4×, PeakATail 5.2×, scUTRquant 36×) and 17–57× on
  mouse. Do not quote a single "×above null" figure for the benchmark.

**Verdict: VERIFIED (adversarial re-verification pass, 2026-08-19).** Nine plotted values were
recomputed from scratch — straight from the genome-filtered point BEDs with `bedtools closest`, not
from the score TSVs or the consolidated table — and all reproduce exactly: human P@100 PeakATail
0.117538 / SCAPTURE 0.651808 / Sierra 0.255854, human R@100 PeakATail 0.156189 / scAPAtrap 0.300376
(denominator 285,136), mouse P@100 SCAPTURE 0.693550 / PeakATail 0.369159 / polyApipe (m2) 0.411943 /
Sierra (m2) 0.581064, mouse R@100 scAPAtrap 0.245315 (denominator 126,686). Scoring the same mouse
point set against the *human* atlas gives 0.016 rather than 0.694, confirming the species split is
real. All four panel-(d) concordances and polyApipe's depth stratification (62.9% depth-1, 0.3079 vs
0.8093) were recomputed independently and match to 6 dp. All 1,066 consolidated rows round-trip
against the correct per-dataset score TSV (910 values + 130 nulls + 26 replicate rows, zero
mismatches). Runtimes were checked against the raw `time -v` logs, retries included. Corrected in
this pass: the SCORING footnote paired 285,136 sites with the wrong gene count (14,851 → **14,949**);
the deck called Sierra the fewest-calls arm, which is false on human; panel (e) hatched SCAPTURE's
mouse bars for an *accuracy* caveat while its own footnote promised runtime-only hatching; "35-fold"
→ 37-fold; "26–57× above null" → 4.4–57×; the panel-(d) SCAPTURE note quoted the human evaluated-peak
count inside a mouse panel. Every accuracy number is
harvested unmodified from the `score_*.tsv` files produced by the single scoring path
(`scripts/benchmark_tools/score_tool.py`, itself validated against `benchmark_curated.tsv`), and every
figure value that also appears in `manuscript/09_headtohead_results.md` or
`scripts/benchmark_tools/README_STATUS.md` matches it exactly. The reproducibility panel reuses the
same genome-filtered point sets `score_tool.py` wrote, and reproduces the committed
`replicate_concordance.tsv` values to 6 dp for PeakATail, Sierra and polyApipe.

**⚠ PENDING RE-VERIFICATION (harness bug 0g, found 2026-08-19, fix committed, figure NOT yet
re-rendered).** `benchmark_headtohead.py` built PeakATail's per-PAS depth as
`dict(zip(pasbed.bed names, annotated_matrix.mtx row sums))`. `pasbed.bed` is **coordinate**-sorted;
the matrix rows are **pas_id**-ordered (`05_annotated_matrix/<strategy>/annotated_pas_ids.tsv`). The
two are permutations of the same id set, so `assert len(names) == len(sums)` passed while
**45,921/45,921 (mouse1) and 46,672/46,672 (mouse2) PAS — 100% — were paired with another PAS's
matrix row** (the two orderings differ at every single position). The depth *value* actually changed
for 45,756/45,921 = 99.64% (mouse1), 46,475/46,672 = 99.58% (mouse2) and 276,596/277,164 = 99.80%
(PBMC); the small remainder is coincidental ties (mostly depth 0/1). Fixed in `depth_map()`, re-keyed on
`annotated_pas_ids.tsv`, with a set-equality assert replacing the length assert
(`scripts/manuscript_figures/tests/test_headtohead_depth_keying.py`). **The numbers below are the
shipped ones; they are reported here, not silently rewritten. Nothing outside PeakATail's depth is
touched — polyApipe / Sierra / scAPAtrap rows of `depth_concordance.tsv` are byte-identical, and
`conc_all` (the 0.73–0.75 quoted above and in manuscript/09) is depth-free and **unchanged**.**

| value | where | shipped | re-keyed | Δ |
|---|---|---:|---:|---:|
| PeakATail depth≤1 concordance, m1→m2 | panel (d) whisker, `depth_concordance.tsv` | 0.720856 | **0.775401** | +0.054545 |
| PeakATail depth≤1 concordance, m2→m1 | panel (d) whisker | 0.755074 | **0.701897** | −0.053177 |
| PeakATail depth≤1 **bar** (mean of the two) | panel (d) | 0.737965 | **0.738649** | +0.000684 |
| PeakATail depth≥2 concordance, m1→m2 | panel (d) whisker | 0.746420 | **0.745286** | −0.001134 |
| PeakATail depth≥2 concordance, m2→m1 | panel (d) whisker | 0.733616 | **0.734472** | +0.000856 |
| PeakATail depth≥2 **bar** (mean of the two) | panel (d) | 0.740018 | **0.739879** | −0.000139 |
| `n_depth1` / `n_depth2plus`, m2→m1 | `depth_concordance.tsv` | 739 / 45,915 | **738 / 45,916** | −1 / +1 |
| `frac_depth1`, m2→m1 | `depth_concordance.tsv` | 0.015840 | **0.015819** | −0.000021 |
| PeakATail top-22k P@100 | figure "READ WITH CARE" footnote; `manuscript/09` §"Why PeakATail underperforms" claim 1; `manuscript/10` §1 | 0.1186 (quoted 0.119) | **0.1199** (0.120) | +0.0014 (+1.15%) |
| PeakATail top-36k P@100 | same three places | 0.1198 (quoted 0.120) | **0.1181** (0.118) | −0.0017 (−1.43%) |
| PeakATail top-106k P@100 | same three places | 0.1179 (quoted 0.118) | **0.1178** (0.118) | −0.0001 (−0.11%) |

**Panel (d) bar heights barely move (≤0.0007) — but only because the two directions swapped: the
depth≤1 whisker widens from [0.721, 0.755] to [0.702, 0.775]. Do not read the small bar deltas as
"the bug did not matter"; the pairing was 100% wrong and the agreement is coincidental.**

**The top-N arms were generated by the same defect and are the more consequential hit.** They came
from an uncommitted one-off that zipped `pas.bed` (coordinate order) onto the matrix rows; that
pairing was reproduced bit-for-bit (`--emulate-bug`, differences confined to ties at the exact
boundary depth: 4/22,000, 6/36,000, 24/106,000 ids). Under correct keying the **selected PAS sets are
almost entirely different** — the shipped top-22k shares only **7.89%** of its members with the
correctly-ranked top-22k (top-36k **13.07%**, top-106k **38.22%**) — yet the precision is unchanged to
±0.002. **The manuscript's claim survives and is in fact strengthened**: PeakATail's own depth
ranking carries essentially no PAS information, which is precisely the read-level-evidence deficit
diagnosed in `manuscript/10_caller_fix_plan.md` §1. Re-scored arms (same `score_tool.py` path, same
references) are in `results/benchmark_tools/pbmc_10k_v3/peakatail/topN_rekeyed/`; generator now
committed as `scripts/manuscript_figures/peakatail_topn_rank.py`.

**Still correct, checked explicitly:** the panel-(d) footnote's zero-UMI counts (715/45,921 mouse1,
555/46,672 mouse2) are row properties of the matrix and are keying-invariant; polyApipe's 62.9%
depth-1 / 0.3079 vs 0.8093 stratification is untouched; every `conc_all` and every accuracy number in
panels (a)–(c) and (e) is depth-free and untouched.

**⚠ THE RE-KEYED COLUMN IS NOT THE FINAL ANSWER — IT IS BLOCKED ON BUG 0a.** Fix 0g makes the
harness read the matrix with the tool's *own* row index, which is the only correct thing a harness
can do; it does **not** make the counts in that matrix correct. Item **0a** of
`manuscript/10_caller_fix_plan.md` records that `annotate()` looks rows up by *position* in a
compacted `pas_ids` list against a full-height matrix, so **45,860/45,921 = 99.87% of annotated PAS
on this same testis mouse1 run carry another PAS's counts**, and that fix is still **uncommitted** —
every run produced after `04e0b3a` (2026-03-23) is affected, including the 2026-08-19 benchmark runs
behind panel (d) and the top-N arms. Cross-check on the shipped artifacts: after the 0g re-key, only
**226/45,921 = 0.49%** of mouse1 PAS have an annotated depth equal to the row `pas_id − 1` of
`filterdmatrix.mtx` (the row 0a identifies as the correct one), and `filterdmatrix.mtx` has 1,540
all-zero rows — exactly the condition that makes the positional lookup skew. **So the "re-keyed"
column above is "correctly read from a known-corrupt matrix", not "correct". It will move again when
0a lands, and the "PeakATail's own depth ranking carries no PAS information" reading of the top-N
result cannot be attributed to PeakATail's ranking until then.**

**To clear this flag:** land 0a first (otherwise the numbers below are re-verified against a matrix
that is itself mis-keyed), then regenerate with `python3 scripts/manuscript_figures/benchmark_headtohead.py`
(the `depth_concordance.tsv` cache has already been rebuilt with the fix; the pre-fix copy is kept
beside it as `depth_concordance.pre0g_buggy.tsv`), update the three top-N precisions in the figure's
hardcoded "READ WITH CARE" footnote, in `manuscript/09_headtohead_results.md` claim 1 and in
`manuscript/10_caller_fix_plan.md` §1, then re-verify.

**New in this figure (not in manuscript/09 — carry these forward):**
- scAPAtrap mouse replicate concordance was "pending" in manuscript/09; it is now computed: **0.88
  (m1→m2) / 0.79 (m2→m1)**, i.e. *above* PeakATail's 0.73–0.75 at higher recall. The sentence
  "PeakATail is best-balanced on reproducibility among the high-recall tools" must be re-scoped: it
  holds against polyApipe, not against scAPAtrap. The honest qualifier is that scAPAtrap's concordance
  is measured on a set its own min.count=10 / min.cells=10 filter has already depth-cleaned, and it
  buys that reproducibility at precision 0.24–0.25 vs PeakATail's 0.37.
- SCAPTURE is the precision leader among de novo tools on both datasets. Any framing that calls Sierra
  "the precision arm" is wrong on PBMC (Sierra 0.256 vs SCAPTURE 0.652) and second-place on testis.
- The polyApipe singleton claim is now quantified end-to-end (0.31 vs 0.81 concordance), not inferred
  from the singleton fraction alone.

**Must travel with it.**
- **Depth is not comparable across tools.** polyApipe depth = poly(A)-read `peakdepth`; PeakATail,
  Sierra and scAPAtrap depth = total UMIs from their own count matrices. Never write a cross-tool
  depth threshold; the panel compares each tool against *itself*.
- SCAPTURE's mouse arm is **mouse1 only** — the mouse2 run was truncated by disk exhaustion and
  SCAPTURE exited 0 anyway. It is hatched in panels (b) and (c) and must never be quoted as a
  two-replicate result. It has no reproducibility bar and no depth at all (its BED score column is 0
  for every evaluated peak — 162,346 human, 81,138 mouse1 — verified).
- Hatching in panel (e) marks **runtime** caveats only, never accuracy. scAPAtrap's PBMC 4:04:29 is a
  *resumed* run after an OOM kill at 8:54:21 that had already completed three stages — it understates a
  from-scratch run (12:58:50 of machine time in total) and must be quoted with that sentence attached.
  PeakATail's PBMC 3:27:06 is attempt 4 of 4; attempts 1–3 (1:33:21, 1:52:42, 4:52:58) died on three
  CellRanger-input bugs and are disclosed, not summed. scUTRquant's mouse wall time is
  BAM-conversion-dominated (64/54 min of 76/64 min) by a prep step this benchmark imposed.
- scUTRquant is annotation-based; its precision against an annotation-derived atlas and its
  cross-replicate agreement are both near-tautological. It is drawn open-marker/separated and must
  never be ranked with the de novo tools.
- SCAPTURE is annotation-**guided** (peaks called inside gene models, then filtered by the DeepPASS
  sequence classifier). Its precision lead is best read as evidence for the evidence-type argument in
  manuscript/09, not as a like-for-like coverage-caller comparison.
- Recall is always against the per-dataset detected-gene-restricted atlas (human 285,136 sites /
  14,949 PBMC-detected genes, `shared_refs_pbmc`; mouse 126,686 sites / 14,014 genes,
  `gse104556/shared_refs`). It is not comparable to any recall computed against the full atlas, which
  is also in `benchmark_consolidated.tsv` under `reference=atlas_full`.
- **One exception, disclosed in the table's `notes`:** the non-plotted `sierra_summit` diagnostic arm
  was scored against the older Laughney-derived human denominator (285,220 sites / 14,851 genes)
  rather than `shared_refs_pbmc`. Its `atlas_detected` recall/F1 are therefore not comparable to the
  other human arms (recall@100 0.064610 as scored; 0.064818 on the PBMC-native set). Every plotted
  arm is on the correct per-dataset denominator.
- Sensitivity/diagnostic arms (`polyapipe_all`, `scapture_all`, `sierra_summit`, `sierra_pooled`,
  `peakatail_top{22000,36000,106000}`) are in `benchmark_consolidated.tsv` but deliberately **not**
  plotted; do not mix them into the ranked comparison.
- Concordance is a nearest-neighbour fraction, not a symmetric statistic: the two directions differ
  only through their denominators (drawn as the vertical range on each bar).

## QUARANTINE NOTE (2026-08-20)
`spermatogenesis_control.*` and `kinnex_truth_validation.*`: adversarial verdicts PROBLEM — science
reproduces, prose overstates. Binding corrections in `11_verifier_corrections.md`; do not quote either
figure until regenerated (scheduled with the Stage 3 re-tests).

## fdr_calibration_v2 — switch-test calibration on a correctly keyed matrix (2026-08-21, verified SOUND)

Supersedes `fdr_calibration` (v1, mis-keyed matrix). 2×2: (a) null p histograms, (b) QQ, (c) q<0.05
hits per null run vs TRUE, (d) TRUE hits surviving permutation-calibrated q. Six arms: fisher reads /
fisher cells / nb_pairwise, each with default top-200 marker pre-selection (solid) and with
`--marker-top-n 0` (dashed/hollow). Headline: only fisher-cells without pre-selection gives valid
FDR control (conservative, 3.0% null p<0.05, 0/20 null runs with hits); defaults are
anti-conservative (20.3% / 13.0% / 24.7%).
**Caveats that travel with it:** marker mode also changes the Fisher gene denominator (restricted
matrix), so marker-on vs marker-off are different tests, not subsets; single mouse/tissue; TRUE hit
counts are detectability upper bounds, not precision; label-permutation null tests the global null
only; KS rejects for every arm (discrete Fisher) and is not gated. Full write-up:
[14_switch_calibration_v2.md](14_switch_calibration_v2.md).

## final_benchmark — Fig 2, final Stage-2 run vs the competitor panel (2026-08-21, numbers verified SOUND in 15)

Four panels: (a) PBMC 10k v3 and (b) testis mice — atlas-agreement precision@100 vs detected-gene recall for
every tool, PeakATail drawn as a path of operating points (shipped → both tiers → tier-1 ≥1 mol → +IP →
precision default), F1 isolines, the pre-registered gate P ≥ 0.50 and the original gate (P ≥ 0.38, F1 > 0.261),
3-seed genic-shuffle nulls; (c) call-set sizes (log); (d) PeakATail tier decomposition on the PBMC IP arm
(tier-2 / tier-1 singletons / ≥2-molecule default with each slice's P@100). Headline: precision default P@100
0.717 / 0.741 / 0.755 (gate PASS on all three), recall at or below polyApipe, F1 tied.
**Caveats that travel with it:** atlas-agreement precision is not ground-truth precision (Kinnex check in 16);
recall denominator is the detected-gene atlas (full-atlas recall 0.087 / 0.096–0.098 in caption); single PBMC
donor; the non-IP ≥2-molecule file is not the default and is not shown; the IP-filter minus-strand bug (Stage 1d)
will move the PBMC default to ~0.706 on re-run. Full caption: `figures/final_benchmark.caption.md`; write-up
[15_final_gate.md](15_final_gate.md).

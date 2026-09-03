# FIGURES_MANIFEST — PeakATail draft submission (Genome Biology, Method track)

Machine-checkable roster of every figure in the submission scheme. Generated at the
2026-09-02 rename pass; the old→new stem trace is `FIGURE_MAP.tsv` (one row per old stem).
Convention: unpadded `fig1_<slug>` / `figS1_<slug>` (supersedes 21 §4.1's zero-padding — Fig 1
shipped on the unpadded convention and D6's rename-once rule forbids repadding it).
One stem owns every artefact: `scripts/manuscript_figures/<stem>.py` →
`manuscript/figures/<stem>.{png,pdf,caption.md}` + `results/figures/manuscript/<stem>*.tsv`
(every plotted value; the caption is written by the script, never by hand). Since the
2026-09-02 submission surgery pass the sidecar is structured `## Legend` (the journal
legend — the single source of the caption, opening `Figure N | Title.`, carrying the
panel-by-panel description and every binding caveat that used to sit in the on-figure
footers) followed by `## Provenance` (sources, TSV pointers, commits); the images carry
panel letters, short titles, axis labels and data annotations only.

House rules (properties of every script, verified at this pass): `LC_ALL=C` mandatory
(tr_TR machine locale), 600 dpi PNG (submission surgery pass 2026-09-02; was 300 dpi in the
working phase) + fonttype-42 PDF with no Type 3, Okabe-Ito
colourblind-safe palette, an audit TSV row for every plotted value, the 8-pixel blank-edge
assert, `OMP_NUM_THREADS=1` for renders. Version/env-var switches keep their pre-rename
names (`FINAL_BENCHMARK_VERSION`, `TRUSTED_NOVEL_VERSION`, `LAUGHNEY_SWITCHES_VERSION`,
`FIG3_WORKDIR`; `SPERMATOGENESIS_VERSION` added 2026-09-02 with the Fig 5 v2 supersession,
default `v2`; `FIG1_STYLE` added 2026-09-02 with the Fig 1 redesign, default `simple`, with
`detailed` writing the separate `fig1_overview_detailed.*` stem): they select runs or renders,
not stems, and the reproduction commands recorded in 15/16/19/20 must keep working. Since the
2026-09-02 design pass the six mains additionally import the shared modules
`scripts/manuscript_figures/_pubstyle.py` (palette / tool identity / type scale / rcParams /
sentence case) and, for Figs 2-3, `_molsweep.py` (one molecule-support sweep implementation). Score TSVs under `results/benchmark_tools/*final*` are
read-only inputs.

Status vocabulary:
- **submission-ready** — verified on the v2 chain; rename executed; regenerate freely.
- **regenerate-on-v2** — (retired 2026-09-02: no row carries it any more; Fig 5, its last
  holder, is regenerated and verified on code `9dfdefb3`) the render was the verified v1
  record pending regeneration on the merged code.
- **to-build** — script does not exist yet; sources are verified and listed.
- **rebuild-pending** — an old script/asset exists but the figure must be rebuilt under the
  listed binding preconditions before any number from it may be quoted.
- **blocked** — a required input (a null, a gap item) does not exist yet.

## Main figures

Publication design pass, 2026-09-02/03 (`DESIGN_DIRECTIVES.md`, all six items executed and
audited): one message per panel; development history off the mains (S8/S12 carry it); exactly
one precision/recall plane across Figs 1-3 and it is Fig 2 (item 5, redundancy); sentence-case
prose labels with canonical identifiers exempt; every figure- and panel-level descriptive
subtitle in the sidecar `## Legend`, never on the image. All six mains import the single style
module `scripts/manuscript_figures/_pubstyle.py` (`PAL`, `TOOL_STYLE`, `TYPE`, `apply_rc()`,
`sentence_case()`) - no script defines a local palette or type size any more - are laid out at
the final print width 7.09 in / 180 mm (Fig 1's retained detailed render excepted at 7.35 in),
render 600 dpi PNG + fonttype-42 vector PDF, and assert a text-overlap, minimum-type-size
(6 pt) and off-canvas audit alongside the 8-pixel blank-edge check. Fig 2 and Fig 3 share the
molecule-support sweep through `scripts/manuscript_figures/_molsweep.py`. NO-LOSS RULE: every
sentence and value that left an image is in that figure's `## Legend`, every audit TSV for
retained data regenerates byte-identically (verified 2026-09-03 against the pre-redesign
scripts), a panel that was removed keeps its TSV emission under a `not plotted since
2026-09-02` header, and an EXPECT assert whose value left the image became a labelled
SIDECAR-GENERATION CHECK rather than being deleted.

| # | Title | Stem | Script | Sources of truth | Status | Regenerate |
|---|---|---|---|---|---|---|
| Fig 1 | How PeakATail finds 3' ends in single cells, and what it gives you - four plain-language SCHEMATIC panels (tail-carrying reads mark the end; ends pile into a called site; the internal-priming look-alike and the genome check; trusted sites -> per-cell counts -> calibrated comparison) with exactly two headline badges (~0.6% clip rate = 0.5730% genome-wide; the trusted set's rule). No panel plots a data series; every quantity is measured by the script, written to five audit TSVs and printed in the Legend | `fig1_overview` | `scripts/manuscript_figures/fig1_overview.py` | 19, 16, 13, 10, 23 §2; every value in the sidecar Legend + `fig1_overview*.tsv` | **submission-ready** (v2, `9dfdefb`; redesigned 2026-09-02 under directive 3, technical no-loss audit 2026-09-03). Caveat that travels: the second badge names the paper's pre-registered TRUSTED SET, not a shipped CLI default - `--polya-min-umis` is 1 at the caller and the genome check is the opt-in `--ip-filter` | `export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/fig1_overview.py` - data-rich six-panel record: prefix `FIG1_STYLE=detailed` (writes the separate `fig1_overview_detailed.*` stem; both styles write the same five TSVs byte-identically, so no restore run is needed) |
| Fig 2 | Where PeakATail sits against the field, and what it costs at a matched call budget: the paper's ONLY precision/recall plane (a PBMC, b testis) carrying the molecule-support operating curve and the two user-choosable points (>=2 precision default 0.7062/0.7450/0.7572, >=1 sensitivity arm), pre-registered gate P >= 0.50 PASS x3, plus precision at a matched call budget (c/d, 25 §2) | `fig2_accuracy` | `scripts/manuscript_figures/fig2_accuracy.py` (+ `_molsweep.py`) | 19 (FIXED) + `final_v2_verify/VERIFIED_v2.md`; 25 §2/§6/§8; v1 record 15 | **submission-ready** (v2; directive 4 executed - the "shipped (pre-fix)" point, both no-IP intermediates and the arrowed dev path are `plotted = False` rows quantified in the Legend and drawn in S12/S8; the tier-decomposition panel moved to Fig S2 and its table stays in `_tiers.tsv`; F1 isolines, the superseded 0.38 gate and every per-call-set genic-shuffle null stay as rows in `_reference_lines.tsv`) | `export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/fig2_accuracy.py` - v1 record: prefix `FINAL_BENCHMARK_VERSION=v1` |
| Fig 3 | Robustness: matching-window resolution 10-100 bp (a/b), biological-replicate and pre-registered cross-donor agreement with its call-count ceiling (c), and compute for the CURRENT caller vs the five competitors on the same box (d). The trade surface moved to Fig 2 and the v1->v2 arrow to S8/S12; the stem is retained for file identity | `fig3_tradeoff` | `scripts/manuscript_figures/fig3_tradeoff.py` (+ `_molsweep.py`) | 19 §1/§3/§4/§5; 13 §1; 12 CORRECTION; 26 §5 (FIXED) | **submission-ready** (v2; directives 4 + 5 executed - the molecule sweep is still recomputed and asserted here and written to `fig3_tradeoff.tsv` with `plotted = False`, and the v1 / uncontended / shipped-caller compute rows stay as `not_plotted_history` / `not_plotted_caveat` / `audit_only`) | `export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/fig3_tradeoff.py` (cache: `$FIG3_WORKDIR`, default `/mnt/ssd0/emaout/fig3_trade_repro_work`; delete to force recompute) |
| Fig 4 | Which `switch diff` test configurations control the false-discovery rate: of six, only Fisher/cells/no-marker-preselection does (3.0% null p<0.05, 0/20 null runs with a q<0.05 hit); every configuration keeping the default top-200 marker pre-selection is anti-conservative (20.3 / 13.0 / 24.7%), as are Fisher reads and NB pairwise without it. Three panels on one shared configuration axis: false-positive rate, false calls per null run, cost in the real run | `fig4_calibration` | `scripts/manuscript_figures/fig4_calibration.py` | 14 (verified SOUND); `results/fdr_calibration_v2/report.json` | **submission-ready** (reframed 2026-09-02 per directive 4's deliberate exception - "which configurations control FDR", not self-repair; audited 2026-09-03: the finding is NOT softened, the Legend still states the anti-conservative verdict for every default configuration verbatim. The null p-value histogram and the QQ-vs-uniform panel are retired; `_hist.tsv` is still written and the leading-bin densities, min null p, 1e-16-floor counts and KS D are asserted in the sidecar-generation check and quoted in the Legend) | `export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/fig4_calibration.py` |
| Fig 5 | Testis biology control, six one-message panels: per-gene monotone shortening 31.4% / 30.5% vs nulls 17.4% / 21.0% (z 10.5 / 6.4); composition-controlled per-cell residual falls at every step (Cliff's delta 0.55 / 0.63); cross-mouse PAS replication 6.1-11.2k per pair, 0 in all 15 null pairings, and the per-gene effect with its rho-vs-null strip; the two PI-endorsed negatives inside a dashed "not claim carriers" frame | `fig5_spermatogenesis` | `scripts/manuscript_figures/fig5_spermatogenesis.py` | 18 **v2 section** (SOUND 2026-09-02; merged code `9dfdefb3`); `results/stage3_spermatogenesis_v2/` | **submission-ready** (v2 default; the 2026-09-02 v2 re-run was adversarially verified to the digit - the LAST FREEZE BLOCKER of 21 gap 2 is cleared. Design pass removed nothing: all six panels and every plotted number stay; the title, subtitle and numbered cautions footer live in the Legend) | `export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/fig5_spermatogenesis.py` - v1 record: prefix `SPERMATOGENESIS_VERSION=v1` |
| Fig 6 | Tumour-cohort application, four panels: the Stage-3 v2 replication funnel real vs the mean of 10 patient-wise label-shuffle nulls (nothing survives under the null), the top 12 of the 47 multi-patient cell-type pairs with the K>=3 sensitivity subset, patient support, and the genomic-context honesty panel (kept on the image deliberately) | `fig6_cohort` | `scripts/manuscript_figures/fig6_cohort.py` | 20 (SOUND 2026-08-22; v2_code default) | **submission-ready** (v2_code; the effect-size histogram is retired - `fig6_cohort_effects.tsv` is still written under a `not plotted since 2026-09-02` header and its median 0.3390 / mean 0.3609 / modal bin 0.29 / 291 below floor / 730 removed by the per-patient floor are asserted in the sidecar-generation check and quoted in the Legend) | `export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/fig6_cohort.py` - v1 record: prefix `LAUGHNEY_SWITCHES_VERSION=v1_code` |

## Supplementary figures (ordered by first citation: S1–S7 with Figs 1–2 and the negative result, S8 with Fig 3, S9 with Fig 4, S10 from R6, S11–S12 the audit trail)

| # | Title | Stem | Script | Sources of truth | Status | Regenerate |
|---|---|---|---|---|---|---|
| Fig S1 | Dataset table made visual, incl. per-BAM poly(A) clip rate (corrected 0.573% genome-wide, 23 §2). Cited from Methods/R1. Replaces retired `cohort_qc` | `figS1_datasets` | `scripts/manuscript_figures/figS1_datasets.py` | v2 run_config/run_manifest JSONs; cohort manifest; 04; 10 §R2; 23 §2 | **submission-ready** (verified 2026-09-02 adversarial pass: cold re-run byte-identical; values traced to 23 §2/§3, 26, 20, run_config JSONs). Caveats: unmeasured cells drawn as em-dashes, never approximated; head-sample estimator values comparable only to each other, never to the genome-wide bar; 1.152% appears only as the SUPERSEDED reference | `export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/figS1_datasets.py` |
| Fig S2 | v2 peak-call QC: widths, per-gene PAS counts, tier composition, molecule-count distributions (arm always named, MUST-NOT-CLAIM 8), ~90–105 nt cleavage-offset census with its boxed NEGATIVE disposition (correction fails every 28 §2 criterion; nothing adopted). Cited from R1 | `figS2_peakqc` | `scripts/manuscript_figures/figS2_peakqc.py` | v2 run outputs (read-only score TSVs, `run/pasbed.bed`); 19; 12 CORRECTION; `results/paramsweep/` (arms_all.tsv, leadB/, VERDICTS.md §2–§3) | **submission-ready** (built+verified 2026-09-02: cold re-run byte-identical; tier ns/P@100 and offset values traced to 19 and VERDICTS.md). Caveat: panels e–g cite VERDICTS §2–§3, self-marked provisional pending its own verifier — stated in the caption Legend (moved off the image at the 2026-09-02 surgery pass); the disposition is negative (report inferred cleavage, do not move calls), so no manuscript claim rests on it positively | `export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/figS2_peakqc.py` |
| Fig S3 | Why the benchmark is built this way: dump vs curated regimes (spread <0.004 vs 0.137), chance = local density of the curated reference (ratio 0.97 at 100 bp), one reconciled null across every verified score TSV (defaults 33–55× their nulls). Cited with Fig 2 | `figS3_nulldesign` | `scripts/manuscript_figures/figS3_nulldesign.py` (rebuild of `null_control.py`, executed 2026-09-02) | 07 §1/§2 (SOUND); `benchmark_curated.tsv`; `fig2_accuracy.tsv` (19, FIXED); densities recomputed fresh from `polyasite2.GRCh38.96.rep_sites.bed6`; 05 R4/R5 only as the superseded record | **submission-ready** (verified 2026-09-02; both MANDATORY preconditions confirmed met: (i) the 0.753-vs-0.615 conflict is shown only as the RETIRED discrepancy over 07 §1's one reconciled null, (ii) every density is a fresh recompute on the curated point reference — no dump density carried). Caveats: two sub-panels dropped with reasons in the caption sidecar (match-class composition — only verified record is dump-derived; quantitative interval-vs-point on curated — never run); the null≈density ratio 0.97 is a consistency observation of this recompute, labelled as such; SCAPTURE mouse 2 not plotted per fig2_accuracy convention (only mouse 1 carries a verified null-bearing row) | `export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/figS3_nulldesign.py` (density cache: `results/figures/manuscript/.cache_figS3_nulldesign/`; delete to force the ~2 min recompute; `null_control.*` stays only as the superseded record) |
| Fig S4 | Second-donor validation (26): pbmc4k scored by the unchanged pre-registered default — gate PASS, P@100 0.8279 (4 libraries now) + cross-donor concordance both directions (84.2% fwd / 39.9% rev vs its 44.4% arithmetic ceiling). Cited with Fig 2 | `figS4_seconddonor` | `scripts/manuscript_figures/figS4_seconddonor.py` | 26 (verifier FIXED; same scorer/nulls); pbmc4k run score TSVs | **submission-ready** (verified 2026-09-02: cold re-run byte-identical; 0.8279/n 20,672, 84.2/81.3/39.9/36.1%, matched-N 0.886–0.907 all traced to 26). Caveats: the matched-N caveat is its own panel c and MUST travel ("precision generalises; it does not improve"); panel c's five donor-1 matched-N points exist only in 26's verifier table (no results-tree TSV) — cited constants, recorded in the TSV source column; second *library/chemistry/CellRanger*, only presumptively a second individual | `export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/figS4_seconddonor.py` |
| Fig S5 | Per-tool long-read re-ranking and the genuine / internal-priming / unsupported decomposition (incl. Sierra ~half-IP). Cited with Fig 2 | `figS5_longread_rerank` | REBUILD of `kinnex_truth_validation.py` → `figS5_longread_rerank.py` (not yet executed) | 11 (BINDING corrections); 05 QUARANTINE NOTE | **rebuild-pending (quarantined)** — rebuild under 11: never "every stringency"; support is an alignment-record count, not de-duplicated UMIs (de-duplicate or relabel the axis); numbers unquotable until rebuilt. **Disposition 2026-09-02 (verifier): KEPT pending, slot retained** — blocker confirmed real (no de-duplicated truth-set support exists; ~10.5% of records repeat a (CB,UMI) at the same terminus, 11). Unblocking run: de-duplicate the Kinnex truth-set records (or relabel the support axis), re-score `results/reliability/trusted_novel_final_v2_pbmc/` + competitor call sets at the ≥5/≥20 thresholds only, verify, then build. Content (Sierra decomposition, per-tool ρ) is not carried by Fig 2/S7, so not droppable on carries-it grounds; 21 §S5 records it strengthens R3 without blocking submission — if still unbuilt at freeze, drop it THEN and renumber S6→S5 … S12→S10 in one pass | — (rebuild first) |
| Fig S6 | Hexamer at −40..−5 of the clip-seeded cleavage point (76.7%), A-fraction profile, hexamer adds +0.046 overall and nothing among atlas-novel sites. Cited with Fig 2 | `figS6_motif` | REBUILD of `motif_validation.py` → `figS6_motif.py` (re-anchored; not yet executed) | 07; 05 §7 | **blocked** on gap 5: the shuffled-position null for the re-anchored window does not exist; without it the panel cannot ship. **Disposition 2026-09-02 (verifier): KEPT pending, slot retained** — blocker confirmed real (no re-anchored shuffled-position null artefact exists anywhere under `results/`; 21 §S6 calls it cheap to compute). Unblocking run: the gap-5 shuffled-position null for the re-anchored −40..−5 window, verified, then the re-anchored rebuild. If still unbuilt at freeze, drop and renumber with S5's pass | — (blocked) |
| Fig S7 | The pre-registered trusted-novel definition and its failure: 46,524 → 35,712 → 7,259 → 4,329; concordance 0.765 → 0.811 → 0.485 → 0.502 vs the 0.70 target; 69.4% intronic; 27.0% at Kinnex IP decoys vs 7.6%. NEGATIVE result; cited from the Kinnex corroboration paragraph AND R6; stays in the Author Summary | `figS7_novelfunnel` | `scripts/manuscript_figures/figS7_novelfunnel.py` | 16 §v2 (FIXED); 19 §4 | **submission-ready** (v2; demoted from 21's main Fig 6 by map amendment (d) — figure moved, prominence did not; all C11 caveats travel verbatim) | `export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/figS7_novelfunnel.py` — v1 record: prefix `TRUSTED_NOVEL_VERSION=v1` |
| Fig S8 | Compute detail behind Fig 3d: 293.7 GB → 12.53 GB, 3:45:53 → 34:37 (concurrency disclosed, MUST-NOT-CLAIM 5), cohort 9:06:26 → 1:06:26 at identical output (505,197 unified PAS), competitor runtimes with retry/skip caveats. Cited with Fig 3 | `figS8_compute` | `scripts/manuscript_figures/figS8_compute.py` | `fig3_tradeoff_compute.tsv` (verified); 19 §3/§5; 15 §5 | **submission-ready** (verified 2026-09-02: cold re-run byte-identical; wall/RSS values traced to `fig3_tradeoff_compute.tsv`, 19 §3/§5, 15 §5, 10 §5). Caveats: competitor retry/skip caveats carried verbatim (scAPAtrap resumed-run understates from-scratch, scUTRquant sums two attempts, SCAPTURE prebuild untimed); MUST-NOT-CLAIM 5 in force on every panel (concurrency disclosure in the caption Legend since the 2026-09-02 surgery pass; the per-bar `*` markers stay on the image); 150 GB stop-signal line drawn from 10 §5 | `export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/figS8_compute.py` |
| Fig S9 | Per-stage-pair and per-stratum null rates (2.6–4.0%), permutation-calibrated q, reads-vs-cells count mode, NB dispersion-floor diagnostic. Cited with Fig 4 | `figS9_calibration_extended` | `scripts/manuscript_figures/figS9_calibration_extended.py` | `fig4_calibration_*.tsv` (verified); 14 (SOUND); `results/fdr_calibration_v2/` (read-only; no permutation re-run) | **submission-ready** (verified 2026-09-02: cold full-recompute re-run byte-identical; B0 2.94–3.13%, 20.3/13.0/9.6/3.0%, q_perm 66,630 vs 60,332, floor 8.5%→67%, min null p 4e-182 all traced to 14 / per_pair TSV). Caveat: panel b's per-stratum bars are a deterministic recomputation from the archived B0 null p-values under a declared binning (2.5–4.0%; 14's verified range is 2.6–4.0% and its exact bin edges were not persisted) — stated in the caption Legend and the strata TSV, with the short `band: verified range 2.6-4.0%` label kept on the image (2026-09-02 surgery pass); pooled rate matches 14's 3.03% exactly | `export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/figS9_calibration_extended.py` (`FIGS9_REUSE_STRATA=1` reuses the strata TSV for layout-only re-renders) |
| Fig S10 | The dropped claim stated honestly: PAS-profile clustering recovers GEX types (PBMC AMI 0.708 / ARI 0.502; Laughney median AMI 0.662 / ARI 0.463, 17/17) AND collapsing 275,370 sites to 14,891 gene totals gives AMI 0.698 — site resolution adds nothing. Absorbs `parameter_sweep`'s Leiden finding (+112.5%). Cited from R6 | `figS10_clustering` | `scripts/manuscript_figures/figS10_clustering.py` (rebuild of `clustering_concordance.py`, executed 2026-09-02) | `pbmc_novelty.tsv` verified ablation rows (01 §S4); `B2_gex_celltyping/concordance.csv` (05 R2, FIXED; READ-ONLY); `parameter_sweep.tsv` (05 R3) | **submission-ready** (verified 2026-09-02: cold re-run byte-identical; 0.708/0.6978, Laughney medians 0.6616/0.4627 (17/17), Leiden +112.5% traced to the three named sources). The mandatory condition holds: the ablation is a co-equal panel on the SAME AMI axis as the positive result. Caveats: reference partition is marker-derived, not curated truth (in the caption Legend since the 2026-09-02 surgery pass); the isoform-only USAGE spaces are excluded as mis-keyed (01 §S4), and the ablation verdict does not depend on them; VERDICTS §1's inert-prominence label carried in the caption Legend | `export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/figS10_clustering.py` |
| Fig S11 | The pre-registration timeline as a figure: original gate, stale Stage-2 failure, disclosed post-hoc sweep, default committed `0e27b1a` 01:18:59, arms 02:59:36, v2 re-run 16:14, every gate outcome on every arm. Makes the pre-registration auditable. Cited from R2/R6 and Methods | `figS11_gatehistory` | `scripts/manuscript_figures/figS11_gatehistory.py` | 12 CORRECTION; 13; 15; 19; 24; 26; 27; table T3 | **submission-ready** (verified 2026-09-02: cold re-run byte-identical; original gate traced to 10 §5, Stage-2 FAIL values to 12, commit 01:18:59/+1:40:37 to 15 §1, v1/v2/pbmc4k gate values to 15/19/26, prime Δ 0.000000 identity to 24 Amendment 3). Caveats: docs 10/12/24/27 record dates without clock times — those markers sit at display positions flagged in the audit TSV's display_note column; Stage-2 values are cited from 12 (its call sets are superseded; no verified score TSVs survive), everything else re-read programmatically from score TSVs | `export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/figS11_gatehistory.py` |
| Fig S12 | The four-way version benchmark (shipped → v1 `4efeb125` → v2 `9dfdefb` → prime `6954082d`): prime's default byte-identical to v2's manuscript arm (pre-registered criterion FAILS at Δ 0.000000 as an identity, Amendment 3), flag-matched dominance panels, #97 compute drop. Everything prime is EXPLORATORY w.r.t. the manuscript gates (24 §3.4). Audit-trail companion to S11; cited from Methods/Discussion | `figS12_versions` | `scripts/manuscript_figures/figS12_versions.py` | 24 §3 (Amendment 3); `results/prime_bench/README.md` + `fourway.tsv`; reads `fig2_accuracy.tsv` | **submission-ready** (regenerate after any Fig 2 regeneration — it reads the Fig 2 audit TSV) | `export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/figS12_versions.py` |

## Retired outright (scripts kept for the record; do not regenerate, do not cite)

`cohort_qc`, `parameter_sweep` (carries the binding `results/paramsweep/VERDICTS.md` §1 label:
`--min-pas-prominence` PROVEN byte-identically inert on the lambda strategies), `benchmark_strategies`,
`benchmark_curated`, `benchmark_headtohead` (depth panels flagged, blocked on caller bug 0a),
`benchmark_tools_running`, `pbmc_novelty`, `fdr_calibration` (v1, mis-keyed matrix),
`spermatogenesis_control` (quarantined v1). Full reasons and replacement pointers:
`FIGURE_MAP.tsv` and the RETIRED block of `../05_figure_index.md`.

## Supplement verification pass — 2026-09-02 (adversarial)

The eight supplements built 2026-09-02 (S1, S2, S3, S4, S8, S9, S10, S11) were adversarially
verified the same day: every script cold re-run (all eight PNGs byte-identical to the builders'
outputs; figS3 additionally byte-identical from a from-scratch density recompute with its cache
deleted), all EXPECT and 8-px edge asserts passing, every PDF pdffonts-checked (TrueType
Identity-H only, zero Type 3), every caption confirmed script-written, and ≥4 plotted values per
figure traced to the named verified documents/TSVs. Targeted hunts came back clean: no
dump-derived number outside the marked SUPERSEDED/RETIRED records, no clip rate other than the
corrected 0.5730% (1.152% only as the superseded reference), no "SCAPTURE mouse-2 invalid"
phrasing, no FDR-style claim sourced to the 10-permutation shuffle nulls (figS9's calibration
claims rest on 14's 20-permutation harness, verified SOUND), and figS4's matched-N caveat is a
mandatory co-equal panel. Dispositions for the two unbuilt slots: S5 and S6 are KEPT pending
with their unblocking runs named in their rows — no slot was dropped, so the S1–S12 numbering
is unchanged and gapless. If either is still unbuilt at freeze, drop it then and renumber the
downstream stems in a single pass (gaps in supplementary numbering are a desk-reject nuisance).

## Submission surgery pass — 2026-09-02 (footer→legend, 600 dpi; adversarially verified)

All 16 built figures (Figs 1–6, S1–S4, S7–S12) were converted for submission: the on-figure
footer captions, figure-level titles/subtitles and multi-sentence in-panel descriptions moved
into each sidecar's `## Legend` (single source of the journal legend; `## Provenance` holds
sources), the canvases were tightened by the freed footer space, and every PNG is now 600 dpi
(PDFs stay vector fonttype 42). Widths: fig3 and fig4 reach the 180 mm double-column norm;
the rest stay wider because legibility outranked the width target (each script's canvas
comment records the blocking element; at 600 dpi a proportional print-scale reduction to
180 mm retains well over 300 effective dpi everywhere).

Adversarial verification of the pass (same day): every script cold re-run (exit 0, all EXPECT
asserts pass, all sixteen 8-px edge checks = 0 ink pixels); every PNG confirmed 600 dpi and
every PDF pdffonts-checked (TrueType Identity-H only, zero Type 3); every rendered PNG viewed
at print size (no orphaned footer references, no blank bands, no new overlaps). No-loss audit:
each Legend was diffed against the pre-surgery on-image strings (git HEAD scripts) and against
the binding-caveat lists of `../05_figure_index.md` — two gaps found and fixed (fig3's Legend
gained the 12-CORRECTION post-hoc-sweep disclosure clause; fig6's Legend Record paragraph
gained the v1-code numbers 14,480 / 74,954-PAS universe and the
`LAUGHNEY_SWITCHES_VERSION=v1_code` reproduction pointer — both from values already asserted
in the scripts; no number was invented or altered). Data unchanged: all 115 audit TSVs under
`results/figures/manuscript/` are byte-identical across the cold re-runs, and for eight
figures (fig1, fig2, fig3, fig5, figS3, figS4, figS8, figS12) the pre-surgery HEAD scripts
were additionally re-run and their TSVs verified byte-identical to the surgical outputs —
moving text changed no data.

## Freeze rule

The rename pass is executed (2026-09-02). Per 21 §4.4 no numbers were changed in it: every
renamed figure's audit TSVs are byte-identical to the pre-rename outputs (verified with `cmp`
at the pass). The spermatogenesis v2 re-run on `9dfdefb` (Fig 5) — the last number-changing
event before submission — landed later the same day, was adversarially verified to the digit
(18's v2 section, SOUND), and `fig5_spermatogenesis.py` now renders it by default
(`SPERMATOGENESIS_VERSION=v1` reproduces the v1 record).

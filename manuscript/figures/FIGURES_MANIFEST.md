# FIGURES_MANIFEST — PeakATail draft submission (Genome Biology, Method track)

Machine-checkable roster of every figure in the submission scheme. Generated at the
2026-09-02 rename pass; the old→new stem trace is `FIGURE_MAP.tsv` (one row per old stem).
Convention: unpadded `fig1_<slug>` / `figS1_<slug>` (supersedes 21 §4.1's zero-padding — Fig 1
shipped on the unpadded convention and D6's rename-once rule forbids repadding it).
One stem owns every artefact: `scripts/manuscript_figures/<stem>.py` →
`manuscript/figures/<stem>.{png,pdf,caption.md}` + `results/figures/manuscript/<stem>*.tsv`
(every plotted value; the caption is written by the script, never by hand).

House rules (properties of every script, verified at this pass): `LC_ALL=C` mandatory
(tr_TR machine locale), 300 dpi PNG + fonttype-42 PDF with no Type 3, Okabe-Ito
colourblind-safe palette, an audit TSV row for every plotted value, the 8-pixel blank-edge
assert, `OMP_NUM_THREADS=1` for renders. Version/env-var switches keep their pre-rename
names (`FINAL_BENCHMARK_VERSION`, `TRUSTED_NOVEL_VERSION`, `LAUGHNEY_SWITCHES_VERSION`,
`FIG3_WORKDIR`): they select runs, not stems, and the reproduction commands recorded in
15/16/19/20 must keep working. Score TSVs under `results/benchmark_tools/*final*` are
read-only inputs.

Status vocabulary:
- **submission-ready** — verified on the v2 chain; rename executed; regenerate freely.
- **regenerate-on-v2** — the current render is the verified v1 record; every number must be
  regenerated on code `9dfdefb` before submission (this is the sole remaining freeze blocker).
- **to-build** — script does not exist yet; sources are verified and listed.
- **rebuild-pending** — an old script/asset exists but the figure must be rebuilt under the
  listed binding preconditions before any number from it may be quoted.
- **blocked** — a required input (a null, a gap item) does not exist yet.

## Main figures

| # | Title | Stem | Script | Sources of truth | Status | Regenerate |
|---|---|---|---|---|---|---|
| Fig 1 | What PeakATail does, drawn on the data it was measured on — method overview: the clip-evidence channel (0.573% genome-wide clip rate, superseding 1.152%), the pre-registered default as a decision tree (167,565 → 46,524), one real locus | `fig1_overview` | `scripts/manuscript_figures/fig1_overview.py` | 19, 16, 13, 10, 23 §2; every printed number cited on-figure | **submission-ready** (v2, `9dfdefb`) | `export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/fig1_overview.py` |
| Fig 2 | Accuracy vs the field: precision default P@100 0.7062 / 0.7450 / 0.7572, pre-registered gate P ≥ 0.50 PASS ×3, matched-N framing (25 §2/§8.2), tier decomposition | `fig2_accuracy` | `scripts/manuscript_figures/fig2_accuracy.py` | 19 (FIXED) + `final_v2_verify/VERIFIED_v2.md`; 25 §2/§8; v1 record 15 | **submission-ready** (v2; stale SCAPTURE panel-b caption line fixed at rename per 21 §11) | `export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/fig2_accuracy.py` — v1 record: prefix `FINAL_BENCHMARK_VERSION=v1` |
| Fig 3 | The trade surface: molecule sweep P 0.706→0.941 vs R_det 0.175→0.083 (only ≥2 pre-registered), 10 bp window check, cross-mouse replicate agreement, compute 293.7 GB → 12.53 GB | `fig3_tradeoff` | `scripts/manuscript_figures/fig3_tradeoff.py` | 19 §1/§3/§4/§5; 13 §1; 12 CORRECTION | **submission-ready** (v2) | `export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/fig3_tradeoff.py` (cache: `$FIG3_WORKDIR`, default `/mnt/ssd0/emaout/fig3_trade_repro_work`; delete to force recompute) |
| Fig 4 | Statistical calibration: of six `switch diff` configurations only Fisher/cells/no-marker-preselection controls FDR (3.0% null p<0.05, 0/20 null runs with q<0.05 hits); shipped defaults anti-conservative (20.3 / 13.0 / 24.7%) with three named mechanisms | `fig4_calibration` | `scripts/manuscript_figures/fig4_calibration.py` | 14 (verified SOUND); `results/fdr_calibration_v2/report.json` | **submission-ready** (caption sidecar added at rename — was the one asset without one; optional sperm-v2 overlay only after gap 2) | `export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/fig4_calibration.py` |
| Fig 5 | Testis biology control: per-gene monotone shortening 31.5% / 30.2% vs nulls 17.5% / 21.0% (z 10.5 / 5.8); composition-controlled per-cell residual falls at every step (Cliff's δ 0.56 / 0.61); cross-mouse replication 6.0–11.3k PAS/pair, 0 in all 15 null pairings; negatives boxed | `fig5_spermatogenesis` | `scripts/manuscript_figures/fig5_spermatogenesis.py` | 18 (FIXED; **v1 code `4efeb125`**); `results/stage3_spermatogenesis_final/` | **regenerate-on-v2** — the `9dfdefb` re-run is the LAST FREEZE BLOCKER (21 gap 2); no spermatogenesis number in the abstract until it lands and is verified | `export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/fig5_spermatogenesis.py` (renders the v1 record until `results/stage3_spermatogenesis_final/` is regenerated on `9dfdefb`) |
| Fig 6 | Tumour-cohort application: Stage-3 v2 replication funnel vs 10 patient-wise label-shuffle nulls (0 replicated in every null), per-pair yields, effect-size floor, genomic-context honesty panel | `fig6_cohort` | `scripts/manuscript_figures/fig6_cohort.py` | 20 (SOUND 2026-08-22; v2_code default) | **submission-ready** (v2_code) | `export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/fig6_cohort.py` — v1 record: prefix `LAUGHNEY_SWITCHES_VERSION=v1_code` |

## Supplementary figures (ordered by first citation: S1–S7 with Figs 1–2 and the negative result, S8 with Fig 3, S9 with Fig 4, S10 from R6, S11–S12 the audit trail)

| # | Title | Stem | Script | Sources of truth | Status | Regenerate |
|---|---|---|---|---|---|---|
| Fig S1 | Dataset table made visual, incl. per-BAM poly(A) clip rate (corrected 0.573% genome-wide, 23 §2). Cited from Methods/R1. Replaces retired `cohort_qc` | `figS1_datasets` | NEW `figS1_datasets.py` | v2 run_config/run_manifest JSONs; cohort manifest; 04; 10 §R2; 23 §2 | **to-build** (cheap) | — (build first) |
| Fig S2 | v2 peak-call QC: widths, per-gene PAS counts, tier composition, molecule-count distributions (arm always named, MUST-NOT-CLAIM 8), ~90–105 nt cleavage-offset analysis (07 §3). Cited from R1 | `figS2_peakqc` | NEW `figS2_peakqc.py` | v2 `pas_support.tsv`; 07 §3; `scripts/paramsweep/` offset census (prime TASK E) | **to-build** | — (build first) |
| Fig S3 | Why the benchmark is built this way: point-vs-interval matching, match-class composition, and that the 18.4M-entry dump could not rank strategies (spread 0.004) where the curated reference can (0.355–0.492). Cited with Fig 2 | `figS3_nulldesign` | REBUILD of `null_control.py` → `figS3_nulldesign.py` (not yet executed) | 05 §5 (superseded record); 07; 21 §4.3-S4 | **rebuild-pending** — two MANDATORY preconditions: (i) reconcile the 0.753 vs 0.615 null conflict to `score_tool.py`'s gene-body-shuffled null (gap 7); (ii) recompute every density figure on the curated PolyASite 2.0 point reference — the 23/33/45/61/73% and "every 168 bp" numbers are dump-derived and may not be carried over | — (rebuild first; `null_control.*` stays only as the superseded record) |
| Fig S4 | Second-donor validation (26): pbmc4k scored by the unchanged pre-registered default — gate PASS, P@100 0.8279 (4 libraries now) + cross-donor concordance 84.2% within 100 bp. Cited with Fig 2 | `figS4_seconddonor` | NEW `figS4_seconddonor.py` | 26 (same scorer/nulls) | **to-build** (small; gap 8's citable asset) | — (build first) |
| Fig S5 | Per-tool long-read re-ranking and the genuine / internal-priming / unsupported decomposition (incl. Sierra ~half-IP). Cited with Fig 2 | `figS5_longread_rerank` | REBUILD of `kinnex_truth_validation.py` → `figS5_longread_rerank.py` (not yet executed) | 11 (BINDING corrections); 05 QUARANTINE NOTE | **rebuild-pending (quarantined)** — rebuild under 11: never "every stringency"; support is an alignment-record count, not de-duplicated UMIs (de-duplicate or relabel the axis); numbers unquotable until rebuilt | — (rebuild first) |
| Fig S6 | Hexamer at −40..−5 of the clip-seeded cleavage point (76.7%), A-fraction profile, hexamer adds +0.046 overall and nothing among atlas-novel sites. Cited with Fig 2 | `figS6_motif` | REBUILD of `motif_validation.py` → `figS6_motif.py` (re-anchored; not yet executed) | 07; 05 §7 | **blocked** on gap 5: the shuffled-position null for the re-anchored window does not exist; without it the panel cannot ship | — (blocked) |
| Fig S7 | The pre-registered trusted-novel definition and its failure: 46,524 → 35,712 → 7,259 → 4,329; concordance 0.765 → 0.811 → 0.485 → 0.502 vs the 0.70 target; 69.4% intronic; 27.0% at Kinnex IP decoys vs 7.6%. NEGATIVE result; cited from the Kinnex corroboration paragraph AND R6; stays in the Author Summary | `figS7_novelfunnel` | `scripts/manuscript_figures/figS7_novelfunnel.py` | 16 §v2 (FIXED); 19 §4 | **submission-ready** (v2; demoted from 21's main Fig 6 by map amendment (d) — figure moved, prominence did not; all C11 caveats travel verbatim) | `export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/figS7_novelfunnel.py` — v1 record: prefix `TRUSTED_NOVEL_VERSION=v1` |
| Fig S8 | Compute detail behind Fig 3d: 293.7 GB → 12.53 GB, 3:45:53 → 34:37 (concurrency disclosed, MUST-NOT-CLAIM 5), cohort 9:06:26 → 1:06:26 at identical output (505,197 unified PAS), competitor runtimes with retry/skip caveats. Cited with Fig 3 | `figS8_compute` | NEW `figS8_compute.py` | `fig3_tradeoff_compute.tsv` (verified); 19 §3/§5 | **to-build** (over existing verified TSVs) | — (build first) |
| Fig S9 | Per-stage-pair and per-stratum null rates (2.6–4.0%), permutation-calibrated q, reads-vs-cells count mode, NB dispersion-floor diagnostic. Cited with Fig 4 | `figS9_calibration_extended` | NEW `figS9_calibration_extended.py` | `fig4_calibration_*.tsv` (verified); 14 | **to-build** (over existing verified TSVs) | — (build first) |
| Fig S10 | The dropped claim stated honestly: PAS-profile clustering recovers GEX types (PBMC AMI 0.708 / ARI 0.502; Laughney median AMI 0.662 / ARI 0.463, 17/17) AND collapsing 275,370 sites to 14,891 gene totals gives AMI 0.698 — site resolution adds nothing. Absorbs `parameter_sweep`'s Leiden finding (+112.5%). Cited from R6 | `figS10_clustering` | REBUILD of `clustering_concordance.py` → `figS10_clustering.py` (+ ablation panel; not yet executed) | 05 §2 (superseded record); R6 ablation numbers | **rebuild-pending** — the positive result without the ablation on the same axes must never ship | — (rebuild first) |
| Fig S11 | The pre-registration timeline as a figure: original gate, stale Stage-2 failure, disclosed post-hoc sweep, default committed `0e27b1a` 01:18:59, arms 02:59:36, v2 re-run 16:14, every gate outcome on every arm. Makes the pre-registration auditable. Cited from R2/R6 and Methods | `figS11_gatehistory` | NEW `figS11_gatehistory.py` | 12 CORRECTION; 13; 15; 19; 24; 26; 27; table T3 | **to-build** | — (build first) |
| Fig S12 | The four-way version benchmark (shipped → v1 `4efeb125` → v2 `9dfdefb` → prime `6954082d`): prime's default byte-identical to v2's manuscript arm (pre-registered criterion FAILS at Δ 0.000000 as an identity, Amendment 3), flag-matched dominance panels, #97 compute drop. Everything prime is EXPLORATORY w.r.t. the manuscript gates (24 §3.4). Audit-trail companion to S11; cited from Methods/Discussion | `figS12_versions` | `scripts/manuscript_figures/figS12_versions.py` | 24 §3 (Amendment 3); `results/prime_bench/README.md` + `fourway.tsv`; reads `fig2_accuracy.tsv` | **submission-ready** (regenerate after any Fig 2 regeneration — it reads the Fig 2 audit TSV) | `export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/figS12_versions.py` |

## Retired outright (scripts kept for the record; do not regenerate, do not cite)

`cohort_qc`, `parameter_sweep` (carries the binding `results/paramsweep/VERDICTS.md` §1 label:
`--min-pas-prominence` PROVEN byte-identically inert on the lambda strategies), `benchmark_strategies`,
`benchmark_curated`, `benchmark_headtohead` (depth panels flagged, blocked on caller bug 0a),
`benchmark_tools_running`, `pbmc_novelty`, `fdr_calibration` (v1, mis-keyed matrix),
`spermatogenesis_control` (quarantined v1). Full reasons and replacement pointers:
`FIGURE_MAP.tsv` and the RETIRED block of `../05_figure_index.md`.

## Freeze rule

The rename pass is executed (2026-09-02). Per 21 §4.4 no numbers were changed in it: every
renamed figure's audit TSVs are byte-identical to the pre-rename outputs (verified with `cmp`
at the pass). The remaining number-changing event before submission is the spermatogenesis
v2 re-run on `9dfdefb` (Fig 5), after which `fig5_spermatogenesis.py` is re-run unchanged
and re-verified.

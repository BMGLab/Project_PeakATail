# Manuscript figures

Final sequential numbering for the draft submission (Genome Biology, Method track), executed
2026-09-02. Naming convention: unpadded `fig1_<slug>` / `figS1_<slug>`; one stem owns every
artefact — `scripts/manuscript_figures/<stem>.py` writes `<stem>.png` (300 dpi), `<stem>.pdf`
(vector, fonttype 42, no Type 3), `<stem>.caption.md` (the caption is script-written, never
hand-edited) here, plus `results/figures/manuscript/<stem>*.tsv` (every plotted value, for
auditing). Never edit these files by hand; edit the script and re-run it.

The machine-checkable roster (status + exact regenerate command per figure) is
[`FIGURES_MANIFEST.md`](FIGURES_MANIFEST.md); the old→new stem trace is
[`FIGURE_MAP.tsv`](FIGURE_MAP.tsv). Full findings, verification verdicts and the caveats that
must travel with each figure are in [`../05_figure_index.md`](../05_figure_index.md).

## Main figures

| # | File | Shows |
|---|---|---|
| 1 | `fig1_overview` | Method overview on real data: the clip-evidence channel (0.573% genome-wide clip rate — supersedes 1.152%), the pre-registered default as a decision tree (167,565 → 46,524), one real locus |
| 2 | `fig2_accuracy` | Final Stage-2 run vs competitors on PBMC and testis; precision default P@100 0.706/0.745/0.757 passes the pre-registered P ≥ 0.50 gate ×3; matched-N framing per 25 |
| 3 | `fig3_tradeoff` | The reliability trade surface: molecule sweep P 0.706→0.941 vs R_det 0.175→0.083 (only ≥2 pre-registered), 10 bp window check, cross-mouse replicate agreement, compute 293.7 GB → 12.53 GB |
| 4 | `fig4_calibration` | Only Fisher cells-mode without marker pre-selection controls FDR (3.0% null p<0.05, 0/20 null runs with hits); shipped defaults anti-conservative (20.3/13.0/24.7%) |
| 5 | `fig5_spermatogenesis` | Testis control: per-gene monotone shortening above shuffle nulls in both mice (31.4% / 30.5% vs 17.4% / 21.0%); per-cell residual falls at every stage; 6.1–11.2k replicated PAS/pair, 0 in all 15 null pairings; negatives boxed. **v2 record (merged code 9dfdefb3, verified SOUND 2026-09-02; `SPERMATOGENESIS_VERSION=v1` reproduces the v1 render)** |
| 6 | `fig6_cohort` | Tumour-cohort application: Stage-3 v2 replication funnel, 0 replicated in each of 10 patient-wise label-shuffle nulls, effect floor, genomic-context honesty panel |

## Supplementary figures

Built and verified (the eight missing supplements were built 2026-09-02 and adversarially
verified the same day — cold re-runs byte-identical, values traced to the verified documents,
PDFs Type-3-free, captions script-written; see FIGURES_MANIFEST.md's verification-pass note
for caveats per figure):

| # | File | Shows |
|---|---|---|
| S1 | `figS1_datasets` | Dataset table made visual + per-BAM poly(A) clip rate (corrected 0.5730% genome-wide; 1.152% only as the SUPERSEDED reference; unmeasured cells are em-dashes) |
| S2 | `figS2_peakqc` | v2 peak-call QC (tiers, widths, per-gene calls, molecule support — arm always named) + the tier-2 cleavage-offset census with its boxed NEGATIVE disposition (nothing adopted) |
| S3 | `figS3_nulldesign` | Why the benchmark is built this way: dump vs curated regimes, chance = local reference density (ratio 0.97), one reconciled null everywhere (defaults 33–55× null) |
| S4 | `figS4_seconddonor` | pbmc4k gate PASS 0.8279 (4 libraries), cross-donor concordance both directions, and the mandatory matched-N panel (precision generalises; it does not improve) |
| S7 | `figS7_novelfunnel` | The pre-registered trusted-novel NEGATIVE result, 0.485 vs the 0.70 target — demoted from a main slot, prominence kept in R6/abstract/Author Summary |
| S8 | `figS8_compute` | Compute detail behind Fig 3d: 293.7 GB → 12.53 GB, concurrency disclosed (MUST-NOT-CLAIM 5), cohort 8.2× at identical output, competitor retry/skip caveats verbatim |
| S9 | `figS9_calibration_extended` | Calibration deep-dive behind Fig 4: per-pair/per-stratum null rates, reads-vs-cells, permutation-calibrated q, NB dispersion-floor diagnostic |
| S10 | `figS10_clustering` | The dropped clustering claim stated honestly — gene-total ablation co-equal on the same AMI axis (site resolution adds nothing); Leiden +112.5% |
| S11 | `figS11_gatehistory` | The pre-registration timeline + every gate outcome on every arm (table T3 drawn), incl. both FAILs and the untested pre-registration |
| S12 | `figS12_versions` | The four-way version benchmark; everything prime is EXPLORATORY (Amendment 3: FAIL at Δ 0.000000 is an identity) |

Pending (slots KEPT at the 2026-09-02 verification pass — blockers confirmed real, unblocking
runs named in FIGURES_MANIFEST.md; numbering stays gapless because no slot was dropped):
`figS5_longread_rerank` (REBUILD of quarantined kinnex_truth_validation under 11's binding
corrections; needs the truth-set support de-duplication first) and `figS6_motif` (REBUILD,
re-anchored; BLOCKED on the gap-5 shuffled-position null). If either is still unbuilt at
freeze, drop it then and renumber the downstream stems in one pass.

## Retired (kept on disk as the superseded record; do not regenerate, do not cite)

| Old file | Why retired | Replaced by |
|---|---|---|
| `cohort_qc` | Headline is a tautology of the upstream TIER filter; the 16,500-site 17/17 intersection is used by no analysis | `figS1_datasets` |
| `parameter_sweep` | `results/paramsweep/VERDICTS.md` §1 PROVED (byte-identity, both species) that `--min-pas-prominence` is a no-op on the lambda strategies; the 13 shipped branches swept only annotation-trim/clustering knobs, and the underlying run is the pre-clip-seeded lambda_gradient cohort — the wrong caller for every current claim | Leiden finding → `figS10_clustering`; the sweep that matters is Fig 3a |
| `benchmark_strategies` | 18.4M-dump precision (0.996–1.000, spread <0.004) cannot rank strategies; its null (0.753) contradicts null_control's (0.615) | lesson → `figS3_nulldesign` |
| `null_control` (as-is) | Computed on the 18.4M-entry dump; no number from it may reach the paper | rebuilt as `figS3_nulldesign` |
| `benchmark_curated` | Round-2 intermediate; superseded by the single scoring path (`score_tool.py`) behind `fig2_accuracy` | `fig2_accuracy` |
| `benchmark_headtohead` | Superseded by `fig2_accuracy` + 25's matched-N framing; PeakATail depth panels flagged (bug 0g re-keyed, blocked on 0a) — must not be cited | `fig2_accuracy` |
| `benchmark_tools_running` | Progress/status figure; supports no claim | — |
| `pbmc_novelty` | Superseded by the pre-registered trusted-novel analysis | `figS7_novelfunnel` |
| `fdr_calibration` | v1 on the mis-keyed matrix; kept renderable for the record only | `fig4_calibration` |
| `spermatogenesis_control` | Quarantined v1 (verdict PROBLEM — prose overstates) | `fig5_spermatogenesis` |
| `kinnex_truth_validation` | QUARANTINED (11's corrections binding); numbers unquotable until rebuilt | `figS5_longread_rerank` |
| `motif_validation` | To be re-anchored on cleavage points; blocked on the re-anchored shuffled-position null (gap 5) | `figS6_motif` |
| `clustering_concordance` | Positive result must never ship without the gene-level ablation on the same axes | `figS10_clustering` |

## Regenerating

```bash
cd /mnt/ssd1/Projects/PeakATail_wd
export LC_ALL=C          # REQUIRED: machine locale is tr_TR; GNU sort otherwise misorders
                         # BED files and bedtools returns silently wrong distances
export OMP_NUM_THREADS=1 # house rule for figure renders
python3 scripts/manuscript_figures/<stem>.py
```

Version switches (kept through the rename — they select runs, not stems; v1 renders stay
reproducible): `FINAL_BENCHMARK_VERSION=v1` (fig2_accuracy), `TRUSTED_NOVEL_VERSION=v1`
(figS7_novelfunnel), `LAUGHNEY_SWITCHES_VERSION=v1_code` (fig6_cohort), `FIG3_WORKDIR`
(fig3_tradeoff cache).

PDFs are vector throughout (no rasterized panels) with subsetted TrueType fonts embedded and
no Type 3 fonts, which is what journals require. Every script asserts the 8-pixel blank-edge
check on its PNG.

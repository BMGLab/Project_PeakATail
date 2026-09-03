# Manuscript figures

Final sequential numbering for the draft submission (Genome Biology, Method track), executed
2026-09-02. Naming convention: unpadded `fig1_<slug>` / `figS1_<slug>`; one stem owns every
artefact — `scripts/manuscript_figures/<stem>.py` writes `<stem>.png` (600 dpi,
publication resolution since the 2026-09-02 submission surgery pass), `<stem>.pdf`
(vector, fonttype 42, no Type 3), `<stem>.caption.md` (the caption is script-written, never
hand-edited) here, plus `results/figures/manuscript/<stem>*.tsv` (every plotted value, for
auditing). Never edit these files by hand; edit the script and re-run it.

Legend location (submission surgery pass, 2026-09-02): the on-figure footer captions,
figure-level titles and multi-sentence in-panel notes of the working phase moved into each
sidecar's `## Legend` section, which is the single source of the journal legend (panel
letters, short titles, axis labels and data annotations stay on the images); `## Provenance`
below it holds sources and audit pointers. Every binding caveat listed in
`../05_figure_index.md` sits in the Legend section, not in Provenance.

The machine-checkable roster (status + exact regenerate command per figure) is
[`FIGURES_MANIFEST.md`](FIGURES_MANIFEST.md); the old→new stem trace is
[`FIGURE_MAP.tsv`](FIGURE_MAP.tsv). Full findings, verification verdicts and the caveats that
must travel with each figure are in [`../05_figure_index.md`](../05_figure_index.md).

## Main figures

Redesigned 2026-09-02/03 under [`DESIGN_DIRECTIVES.md`](DESIGN_DIRECTIVES.md): one message per
panel, development history off the mains, exactly one precision/recall plane across Figs 1–3
(it is Fig 2), sentence-case labels, and the shared style module
`scripts/manuscript_figures/_pubstyle.py` (palette, per-tool marker identity, type scale,
`apply_rc()`, `sentence_case()`) so all six read as one system. Every main is laid out at the
final print width — 7.09 in / 180 mm — saved at 600 dpi PNG + fonttype-42 vector PDF, and each
script asserts a text-overlap, minimum-type-size (6 pt) and off-canvas audit plus the 8-pixel
blank-edge check before it writes.

| # | File | Shows |
|---|---|---|
| 1 | `fig1_overview` | Four plain-language **schematic** panels for a reader who has never run a caller: some reads keep a piece of the poly(A) tail; tail-carrying read ends pile up at one position; the internal-priming look-alike and the genome check that removes it; trusted sites → per-cell counts → a calibrated cell-type comparison. Exactly two headline badges (~0.6% clip rate; the trusted set's rule). No panel plots a data series — every quantity they stand for is measured by the script, written to the audit TSVs and printed in the Legend. **The data-rich six-panel record stays reachable: `FIG1_STYLE=detailed` → `fig1_overview_detailed.{png,pdf,caption.md}`** (7.35 in canvas, working-phase design, deliberately not held to the design laws) |
| 2 | `fig2_accuracy` | Where PeakATail sits against the field, and what it costs at a matched call budget. Panels a/b: the paper's **only** precision/recall plane (PBMC, testis) carrying the molecule-support operating curve with the two user-choosable points (≥2 precision default, ≥1 sensitivity arm) and the pre-registered P ≥ 0.50 gate — PASS ×3 (0.706/0.745/0.757). Panels c/d: precision at a matched call budget (25 §2). **No development history** — the shipped pre-fix caller and the no-IP intermediates are `plotted = False` rows and live in S12/S8 |
| 3 | `fig3_tradeoff` | Robustness: matching-window resolution 10–100 bp, biological-replicate and pre-registered cross-donor agreement with its call-count ceiling, and compute for the **current** caller against the five competitors on the same box. The trade surface moved to Fig 2 (directive 5); the v1→v2 arrow moved to S8/S12. Stem retained for file identity |
| 4 | `fig4_calibration` | Which `switch diff` configurations control the false-discovery rate: only Fisher · cells · markers off (3.0% null p<0.05, 0/20 null runs with a q<0.05 hit); **every configuration that keeps the default top-200 marker pre-selection is anti-conservative** (20.3 / 13.0 / 24.7%), as are Fisher reads and NB pairwise without it. Three panels on one shared configuration axis: false-positive rate, false calls per null run, cost in the real run |
| 5 | `fig5_spermatogenesis` | Testis control, six one-message panels: per-gene monotone shortening above shuffle nulls in both mice (31.4% / 30.5% vs 17.4% / 21.0%); composition-controlled per-cell residual falls at every stage; cross-mouse PAS replication (6.1–11.2k per pair, 0 in all 15 null pairings) and per-gene effect; the two PI-endorsed negatives inside a dashed "not claim carriers" frame. **v2 record (merged code 9dfdefb3, verified SOUND 2026-09-02; `SPERMATOGENESIS_VERSION=v1` reproduces the v1 render)** |
| 6 | `fig6_cohort` | Tumour-cohort application, four panels: the replication funnel real vs the label-shuffle null (nothing survives under the null), the top 12 of 47 multi-patient cell-type pairs with the K≥3 subset, patient support, and the genomic-context honesty panel — kept on the image because a caveat a reader must not miss cannot live only in a caption |

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
(figS7_novelfunnel), `LAUGHNEY_SWITCHES_VERSION=v1_code` (fig6_cohort),
`SPERMATOGENESIS_VERSION=v1` (fig5_spermatogenesis), `FIG3_WORKDIR`
(fig3_tradeoff cache), and — added at the 2026-09-02 design pass — `FIG1_STYLE=detailed`
(fig1_overview), which renders the retained data-rich six-panel record to its own
`fig1_overview_detailed.*` stem. Both Fig 1 styles compute the same numbers and write the same
five audit TSVs byte-identically; the style only chooses what is drawn, so it never needs a
restore run. Note: fig2, fig5 and fig6 write both versions' outputs over the same
file stems, so after a v1 render re-run the default to restore the v2 paper files.

PDFs are vector throughout (no rasterized panels) with subsetted TrueType fonts embedded and
no Type 3 fonts, which is what journals require. Every script asserts the 8-pixel blank-edge
check on its PNG.

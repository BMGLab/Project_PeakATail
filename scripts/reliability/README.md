# scripts/reliability — Stage-3 reliability program tools

Tools that implement the operational definitions pre-registered in
`manuscript/13_reliability_positioning.md` (PI decision, 2026-08-21: reliability
first; claims are low-FP PAS detection, reliable cell-type switches, 3'UTR
length, trusted de novo PAS).

These are **tools**. The final Stage-2 caller run they are meant to consume
does not exist yet. Every number produced by running them on the July 2026
Laughney outputs or on `results/fdr_calibration/` is a **stale dry run that
proves mechanics only** and must not be reported as a result.

Environment: `tools/PeakATail/.venv/bin/python` (3.10, pandas 2.2, scipy 1.15).
Tests: `tools/PeakATail/.venv/bin/python -m pytest scripts/reliability/tests -v`.

---

## replication_filter.py — cross-sample replication filter for switches

Spec row: *"replication | none required today | a switch is reported only if
called (same direction) in >= 2 independent samples (Laughney 17 samples;
testis 2 mice) — the single strongest FP control we have"*, plus the
*effect-size floor* row (|ΔPDUI| >= 0.1 or Δproportion, reported next to q).

### Inputs (one per sample, `--sample NAME=PATH` or `--sample PATH`)

| kind | what is accepted | per-sample test statistic used |
|---|---|---|
| `ema switch diff` output | the output dir (holds `differential/`) or `differential/` itself; files `differential/<strategy>_<c1>_vs_<c2>.tsv` | `qvalue` + the signed effect column |
| `ema switch length` per-cell table | `pdui_classic.tsv` (or the dir holding it) | computed here: `dpdui = mean PDUI(c1) − mean PDUI(c2)`, two-sided Mann-Whitney U on per-cell PDUI, BH per pair (`--min-cells`, default 10) |
| pre-summarised PDUI table | TSV with `gene_id cluster1 cluster2 (dpdui|delta_pdui) qvalue` | as given |

Schema confirmed 2026-08-21 against `tools/PeakATail/docs/cli/switch-diff.md`,
`ema/switch_test/runner.py::_augment_diff_df`, `ema/switch_test/strategies/fisher.py`
and the July Laughney outputs. The fisher per-pair TSV columns are:

```
pas_id gene_id chrom start end strand cluster1 cluster2 pvalue qvalue n_cells
n_cells_cluster1 n_cells_cluster2 n_cells_expr_cluster1 n_cells_expr_cluster2
n_reads_pas_cluster1 n_reads_pas_cluster2 n_reads_gene_cluster1 n_reads_gene_cluster2
odds_ratio delta_proportion log2fc
```

(the docs page still lists the older `statistic` column; the code writes
`odds_ratio` + `delta_proportion` + `log2fc`). `nb_pairwise` writes
`pvalue qvalue log2fc dispersion n_cells test_stat`. `nb_multi` writes one
omnibus table without pairs or a signed effect and is rejected.

**Signed effect and its sign convention** (`--effect-col`, default auto =
`delta_proportion` > `dpdui` > `log2fc`):

| column | definition in code | `+` means |
|---|---|---|
| `delta_proportion` | `prop(c1) − prop(c2)` (within-gene PAS proportion, in `--count-mode` units) | higher in **cluster1** |
| `log2fc` (fisher) | `log2((prop(c2)+eps)/(prop(c1)+eps))` | higher in **cluster2** |
| `log2fc` (nb_pairwise) | cluster2-indicator coefficient / ln 2 | higher in **cluster2** |
| `dpdui` | `mean PDUI(c1) − mean PDUI(c2)` | cluster1 **longer** |

Direction consistency is evaluated on one effect column across all samples,
so the convention only matters for the human-readable `utr_direction` /
`distal_utr_direction` columns. Pair orientation is canonicalised
lexicographically (`N_vs_T`); when a sample wrote `T_vs_N` the effect sign is
flipped before comparison (`flipped` tracked internally; tested).

### Definitions

- *called* in a sample: `qvalue < --fdr` (0.05) **and** finite effect.
- `replication_count` = max(#samples called with effect > 0, #samples called
  with effect < 0) — the same-direction count.
- `direction_consistent` = at least one call and no call in the opposite
  direction.
- `passes_replication` = `replication_count >= K` (`--min-samples`, default 2);
  with `--discordant-policy exclude` (default) it also requires
  `direction_consistent`, i.e. one opposite-direction call in any sample
  vetoes the switch. `allow` relaxes that.
- **Unit of replication** (`--sample-group SAMPLE=GROUP`, repeatable, and/or
  `--sample-group-tsv FILE` with header columns `sample`, `group`; added
  2026-08-21 for manuscript/13 addendum item 2, *"unit of replication =
  patient, not GSM ... for a patient with two GSMs, either GSM counts once"*):
  when a map is given, `n_pos` / `n_neg` (and the floor counts) are the numbers
  of **distinct groups** with >= 1 same-direction call, so
  `replication_count` counts patients. A group whose own samples disagree in
  sign counts on both sides and is therefore vetoed under `exclude`
  (conservative). Unlisted samples are singleton groups; unknown samples in
  the TSV are ignored, in an explicit spec they are an error. New columns:
  `n_groups_tested`, `n_groups_called`; `summary.json["params"]` records
  `unit_of_replication`, `sample_groups`, `n_groups`; per-pair
  `n_samples_tested` / `n_groups_tested` in `real.per_pair`. Without a map the
  table is identical to the per-sample one (tested).
- Effect-size floor (`--effect-floor`, default 0.1): the same counts are
  recomputed over calls with `|effect| >= floor` →
  `replication_count_floor`, `direction_consistent_floor`,
  `effect_floor_pass` (the flag), `passes_replication_floor` (= both).
  The floor is a flag and a second column set, never a hidden filter, so
  both "q-only" and "q + floor" numbers are always reported.
- `--level gene` (switch-diff input): a gene replicates iff one of its PAS
  replicates (same PAS, same direction, >= K samples); the row carries the
  best PAS (`representative_pas_id`), `n_pas_tested`, `n_pas_replicated`,
  and `distal_pas_replicated` / `distal_utr_direction` (`c1_longer` /
  `c2_longer`) from the strand-aware 3'-most PAS of the gene (needs
  coordinates; `+`: largest `end`, `−`: smallest `start`). PDUI input is
  already gene/transcript-level (`feature_id = gene_id[:transcript_id]`).

### Permutation-style control (`--null-dirs`, `--null-perms`)

Expected replication under independent per-sample label shuffles:

- `--null-dirs [NAME=]GLOB` (repeatable): per-sample null switch outputs,
  each a switch-diff output dir such as the existing
  `results/fdr_calibration/null/perm_XX/` layout (the null driver permutes
  `obs[label]` across cells and re-runs the FULL `ema switch diff`, markers
  included). Forms: `NAME=GLOB`, a GLOB containing `{sample}`, or one bare
  GLOB per sample in `--sample` order.
- `--null-perms N` (per-cell PDUI input only): the tool permutes the
  cell -> `cluster` map inside each sample itself (one label per cell,
  broadcast to all of that cell's genes, label multiset preserved — the same
  operation the switch-diff null driver applies to `obs[label]`). Per-cell
  PDUI does not depend on the labels, so this is an exact label-shuffle null
  for the Mann-Whitney/dPDUI summary. The `n_features` of a null combination
  is close to, not identical with, the real one because the `--min-cells`
  filter is re-applied to the shuffled split.
- `--null-universe real` (default): every null combination is scored only on
  the (pair, feature) set that was tested in the real run, so real and null
  replication counts refer to the same universe. A label shuffle changes
  which features pass the per-sample filters (markers / `--min-cells`); on
  the stale July PDUI tables, where a feature's finite-PDUI cells are ~91%
  concentrated in one stage, the unrestricted null universe was 3x the real
  one. `all` keeps every null feature (reported as `n_features` per combo
  either way, plus `n_features_real` / `n_features_null_mean`).
- Combinations: by default perm *i* of every sample is paired (index-aligned,
  `min(n_perms)` combos); `--null-combos N --seed S` draws N random
  independent indices per sample instead. Each combination runs the
  identical filter (K, fdr, floor, pairs, level).
- `null_control.tsv` also carries per-pair columns `rep__<c1>_vs_<c2>` /
  `repf__<c1>_vs_<c2>` (0 when a pair is absent from a combination), and
  `summary.json["null"]["per_pair"]` holds the same mean/sd/empirical-p block
  per pair.
- Reported in `summary.json["null"]`: per-combo counts, mean/sd/median/max,
  `empirical_p_ge_observed = (1 + #{null >= observed}) / (1 + n_combos)` and
  `expected_false_replicated_fraction = mean(null replicated) / observed`
  (an empirical FDR estimate for the replicated set, valid to the extent the
  per-sample null mimics the real pipeline). Samples without null outputs
  are dropped from the null with a warning (and listed in the JSON).

### Outputs (`--out DIR`)

| file | content |
|---|---|
| `all_features.tsv` | one row per (canonical pair, feature) seen in >= 1 sample; counts/flags above; `q__<sample>`, `effect__<sample>`, `called__<sample>`; `min_q_called`, `mean_effect_consensus`, coordinates |
| `replicated.tsv` | rows with `passes_replication == True` |
| `null_control.tsv` | one row per null combination (`perm__<sample>` indices, counts) |
| `summary.json` | params, sample paths + mtimes, `pairs_seen`, `per_sample_called`, `real` counts (total and per pair), `null` block, `note` |

### Example (STALE dry run, mechanics only)

```bash
export LC_ALL=C
PY=/mnt/ssd1/Projects/PeakATail_wd/tools/PeakATail/.venv/bin/python
R=/mnt/ssd2/Laugney_Aligned/peakatail_experiments/runs/A2_base_FIXED_reclustered/B3_stage_switch/diff
$PY scripts/reliability/replication_filter.py \
  --sample MAST=$R/CELL_TYPES_WANSLEEBEN_HOGAN_2013_MAST/fisher \
  --sample DC=$R/CELL_TYPES_WANSLEEBEN_HOGAN_2013_DC/fisher \
  --sample TREG=$R/CELL_TYPES_WANSLEEBEN_HOGAN_2013_TREG/fisher \
  --min-samples 2 --fdr 0.05 --effect-floor 0.1 \
  --out <scratch>/replication_dryrun \
  --note "STALE July-2026 Stage-2 per-celltype outputs; mechanics only"
```

Those July dirs are per-**cell-type** stage contrasts (Normal / StageIA /
StageIVprimary) from the pre-fix caller, not per-patient outputs, and their
Fisher p-values are from the uncalibrated `--count-mode reads` path
(manuscript/13: 38.6% null p<0.05). Using three cell types as "samples" only
exercises parsing, orientation canonicalisation and the flags. The real
Stage-3 input is one `ema switch diff` run per Laughney patient (17) / per
testis mouse (2) on the final caller's matrices, plus their per-sample null
perms.

### Limitations / decisions to make before the real run

- BH q-values are taken as written by `ema switch diff` (per pair). If the
  Stage-3 calibration ships permutation-calibrated q-values, point
  `--effect-col`/fdr at those tables — the filter is agnostic to how q was made.
- Gene-level direction from `min q` across PAS of a gene is a screen; use
  `distal_utr_direction` (or PDUI input) for lengthening/shortening claims.
- Per-isoform PDUI tables (`--isoform-agg per_isoform`) yield one feature per
  transcript (`gene:transcript`); transcripts sharing the same proximal/distal
  pair are counted separately. Collapse them upstream (`per_gene`) if one row
  per gene is wanted.
- The per-cell PDUI summary is a Mann-Whitney on cells (pseudo-replication
  across cells within a sample is exactly what the cross-sample replication
  requirement compensates for); it is not a substitute for the Stage-3
  calibrated test.

---

## trusted_novel_pas.py — pre-registered "trusted de novo PAS" + Kinnex validation

Spec: `manuscript/13_reliability_positioning.md` §3. Two sub-commands, stdlib + `bedtools`
(2.30 on PATH) only; every sort is ASCII/`LC_ALL=C`.

### Definition implemented by `call` (all five must hold)

| # | criterion | how it is evaluated | CLI |
|---|---|---|---|
| 1 | clip-supported | site name (col 4) present in the clip-supported tier BED | `--clip-bed` (omit only if `--pas-bed` is already clip-only; warns) |
| 2 | ≥ 2 distinct molecules | BED col 5 ≥ `--min-support` | `--support-unit molecules\|reads` is **required**; `reads` (pre-Stage-1c outputs) is accepted with a loud warning and `stale_input_dry_run: true` in the manifest |
| 3 | not internal-priming | `--ip-from-annot ANNOT_BED:COL` (tool's internal-priming flag = col 10 of the 11-column per-dataset `run/03_gtf_annotation/<dataset>/annotatedpas.bed`; NOTE the top-level `run/annotatedpas.bed` of the 4efeb125 final run has only 9 columns and NO IP column — the tool raises on a too-short line rather than guessing) or recomputed: ≥ 6 consecutive A **or** ≥ 70 % A in the transcript-strand window **+10..+30** downstream of the cleavage site (same window as `motif_validation.py` panel c; the in-tool filter uses −10..+30, reproducible with `--ip-lo -10 --ip-hi 29`) | `--ip-lo/--ip-hi/--ip-a-stretch/--ip-a-frac` |
| 4 | canonical hexamer | AATAAA/ATTAAA = **strong** tier; any of the 12 canonical variants = **weak** tier (superset); hexamer must lie fully inside **−40..−5** of the cleavage site on the transcript strand | `--hex-lo/--hex-hi` |
| 5 | atlas-novel | `bedtools closest -s -d` distance to the union of `--atlas` files ≥ 100 bp (strand-matched; `--atlas-ignore-strand` for strand-agnostic). Any atlas path containing `hg19` is refused (PolyA_DB 3.2 is hg19-only; see `data/references/atlases/README.md`) | `--atlas` (repeatable), `--atlas-min-dist` |

Cleavage base c (0-based) = `end-1` on `+`, `start` on `-`. Sequence comes from ONE
`bedtools getfasta -s -bedOut` call on a **symmetric** window `[c-W, c+W+1)`, W = 40, so
index i ↔ transcript-relative r = i − W on both strands — the asymmetric-window minus-strand
trap documented in `scripts/manuscript_figures/motif_validation.py` cannot occur. Sites whose
window leaves the contig, or whose contig is absent from the FASTA, are `sequence_testable = 0`
and fail (conservative). Funnel order is the pre-registered order: contigs (optional) →
clip → support → sequence testable → not IP → hexamer (any of 12) → atlas-novel =
**trusted_novel**; → AATAAA/ATTAAA = **trusted_novel_strong**.

Outputs (`--outdir`):

| file | content |
|---|---|
| `trusted_novel.bed` | BED6 cleavage points (name, support, strand) + `orig_start orig_end hex_tier hex_hits ip_flag atlas_dist_bp`; header line starts with `#` |
| `trusted_novel.strong.bed` | the AATAAA/ATTAAA subset, same layout |
| `funnel.tsv` | sequential survival per stage (`n_pass`, `n_fail_at_stage`, fractions) |
| `funnel_marginal.tsv` | each criterion alone over the whole input (which filter bites hardest) |
| `site_annotations.tsv` | one row per input site with every criterion, `ip_source` (annot/recomputed/untestable), `fail_stage` |
| `stages/NN_<stage>.bed` | cleavage-point BED6 of the survivors of every stage — the inputs for `validate` |
| `run_manifest.json` | argv, params, bedtools version, input md5s, funnel, warnings, `stale_input_dry_run` |

### Validation (`validate`) — atlas-independent, vs shuffled null

For each `--query LABEL=BED` × `--truth LABEL=BED` (Kinnex truth point BEDs, e.g.
`results/benchmark_tools/kinnex_truth/x3p/x3p_truth_t{5,20,100,500}.point.bed`; the
`≥ 5/20/100 UMI` sets are the ones `kinnex_truth_plan.md` documents): fraction of query
sites with a same-strand truth 3' end within `--window` (25) bp, Wilson 95 % CI, and a
null from `--seeds` (10) replicates of `bedtools shuffle -chrom -noOverlapping` (width,
chromosome and strand of every record preserved; `--incl` optionally restricts the landing
space, e.g. `results/benchmark_tools/shared_refs/genebodies.merged.bed`). Reported:
`null_mean/sd/min/max`, `enrichment`, `z_vs_null`, `empirical_p = (1+#null≥obs)/(1+seeds)`,
`meets_target_0.70` (13 §3 target). Both query and truth are restricted to the contigs in
`--chrom-sizes` so denominators match. Files: `validation.tsv`, `validation_null_seeds.tsv`,
`run_manifest.json`. Feeding every `stages/*.bed` as queries shows how each filter moves
the long-read concordance.

### Stale dry-run (mechanics only — numbers are NOT results)

`dryrun_stale_pbmc.sh [outdir]` runs both steps on the STALE Stage-2 PBMC output
`results/benchmark_tools/pbmc_10k_v3/peakatail_clipseeded/` (`pas.bed` + `pas_tier1.bed`,
col 5 = reads → `--support-unit reads`; that run had no IP annotation, so IP is recomputed)
against PolyASite 2.0 GRCh38 and the x3p Kinnex truth at 5/20/100/500 UMI. Default outdir
`results/reliability/_dryrun_stale_stage2_pbmc/`, which carries a
`DRYRUN_STALE_INPUTS_NOT_RESULTS.txt` marker. For the real Stage-3 run: point `--pas-bed`
at the final caller's clip-supported default output, `--support-unit molecules`, and
`--ip-from-annot <run>/03_gtf_annotation/<dataset>/annotatedpas.bed:10` (NOT the 9-column top-level `<run>/annotatedpas.bed`; or let it recompute). Executed 2026-08-21 on the final IP-filter arm: `results/reliability/trusted_novel_final_pbmc/` (verifier-checked; see its `REPORT_PROVISIONAL.md` §9).

### Tests

`tests/test_trusted_novel_pas.py` — synthetic C/G-only genome (no hexamer or A-run can
occur by accident on either strand) with implanted signals: hexamer window boundaries
(−40 and −5 inclusive, −41 and −4 excluded), IP run/fraction rule and window edges,
minus-strand sequence mapping through `getfasta -s` (plus-strand site at the same base
sees `TTTATT`), edge/missing-contig handling, `closest` point-distance semantics, full
`call` funnel counts and membership, `--atlas-ignore-strand`, `--restrict-chroms`,
`--ip-from-annot` override + fallback, `--support-unit reads` warning, hg19 refusal,
`validate` observed fraction / CI / per-seed null rows / strand handling.

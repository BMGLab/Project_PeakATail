# Stage-3 Laughney cohort run: tested-PAS universe (addendum item 3), fixed BEFORE any switch statistic is viewed

**Status:** policy document, written 2026-08-21 (cohort run PID 1994653 still peak-calling dataset 2/17; no
`unified/` files, no clusters.h5ad and no switch test exist yet). Implements
`manuscript/13_reliability_positioning.md` Addendum 2026-08-21 05:30, item 3:
*"Tested PAS universe = the pre-registered precision default (tier-1 AND IP-pass AND >= 2 molecules;
TIER_FILTER=tier1_ge2). Coverage-only tier-2 PAS are never switch-tested."*
The pre-registration fixes WHAT the universe is; this document fixes HOW it is computed in cohort mode
(design A, one `ema run -c` over 17 datasets), where the tool's own files do not carry a cohort support value.
Every number produced downstream is PROVISIONAL until its verifier passes.

## 0. What the frozen snapshot (tools/pa-polya-run-4efeb125) actually writes in multi-dataset mode (verified in code + on the chr21 smoke run)

| fact | where verified |
|---|---|
| Per-dataset caller BED col 5 = `clip_umis` = distinct (CB, UMI) poly(A)-clip molecules of that caller PAS; sidecar `peakcalling/<ds>_0.{pos,neg}.support.tsv` carries `pas_id clip_reads clip_umis clip_reads_f3844 clip_umis_f3844 window_reads tier` keyed by the OLD per-dataset pasnumber | `ema/countmatrix/paswrite.py` SUPPORT_COLUMNS |
| Caller `tier`: 1 = clip-seeded cluster, 2 = coverage-only candidate. Tier-2 rows have `clip_umis == 0` by construction (smoke: all 1,955 tier-2 rows; cohort dataset 1: checked again by the builder, asserted) | `support.tsv` of smoke + GSM3516662 |
| `unified/multi_sample_merged.bed` col 5 = score of the FIRST member of the merged interval only (`bedtools merge -s -d 100 -c 4,5,6 -o collapse,first,first`) -> it is ONE member's clip_umis, NOT a cohort support value | `ema/datasets/pas_merge.py` step 3 |
| `unified/multi_sample_pas_mapping.tsv` (`dataset_id strand old_pasnumber new_pas_id`) gives every member of every unified PAS | same |
| `--ip-filter --ip-filter-mode filter` runs on the UNIFIED set; survivors = `run/posbed.bed` + `run/negbed.bed`; per-dataset matrices are annotated from that post-filter set, so an IP-flagged unified PAS never enters any `07_clustering/<ds>/clusters.h5ad` (smoke: 483/483 and 66/66 var pas_id in posbed+negbed; re-asserted per GSM at injection) | `COHORT_LAYOUT_NOTES.txt`, `04_pas_gene_assignment/peak_filters_stats.json` |
| `run/annotatedpas.bed` col 9 `TIER_1/TIER_2` is the UTR-DISTANCE tier of the gene annotation (`ema/annotate/find_close.py`), unrelated to the caller's clip tier. It is NEVER used for the universe | `find_close.py` lines 6-31 |
| `clusters.h5ad` var_names == var['pas_id'] == new unified id; `layers['counts']` integral | smoke verify |

## 1. Universe rule (cohort-level; ONE universe shared by all 17 GSMs)

For every unified PAS `u` (row of `unified/multi_sample_merged.bed`), with members M(u) = all
(dataset, strand, old_pasnumber) rows of the mapping that point to `u`:

* `clip_umis_sum(u)` = sum over M(u) of the member's sidecar `clip_umis`. Members are disjoint caller
  clusters (within a dataset) or different datasets (different cells), so the sum is a count of DISTINCT
  molecules supporting `u` in this run.
* `any_tier1(u)` = at least one member has caller tier 1 (a merged interval that contains a clip-seeded
  cluster is clip-supported). Because tier-2 members carry 0 clip molecules, `clip_umis_sum >= 2`
  implies `any_tier1`; both are still computed and recorded.
* `ip_pass(u)` = `u` present in `run/posbed.bed` or `run/negbed.bed` (post-filter unified set).

**`in_universe(u) = ip_pass(u) AND any_tier1(u) AND clip_umis_sum(u) >= 2`.**

The per-GSM tested set is `in_universe` intersected with the var of that GSM's `clusters.h5ad` (the tool's
own per-dataset matrix filters -- min_cells 3, min PAS per cell 50, gene-annotated -- still apply; the
universe never ADDS a PAS to a GSM). `n_pas_after_tier_filter` in each GSM's `label_report.json` is the
reported "universe size per GSM".

Rationale (product reading of item 3 + s1): the precision-first default is defined on the OUTPUT OF A RUN
("clip-supported tier only, IP filter on, >= 2 distinct molecules of clip support"). In design A the run
is the cohort run, its PAS are the unified PAS, and the run's clip evidence for a unified PAS is the pooled
evidence of its members. Using one cohort-level universe also keeps the tested family identical across
GSMs, which is what makes the cross-patient replication count (item 2) a comparison of like with like.

## 2. Recorded diagnostics that do NOT enter the rule (for the verifier; any switch to them is a post-hoc change that would require re-running BH)

* `umis__<dataset_id>` per unified PAS (that GSM's own clip molecules) and, per GSM at injection,
  `n_pas_after_tier_filter_with_own_umis_ge2` (the stricter "supported with >= 2 molecules in THIS GSM"
  variant) and `n_pas_after_tier_filter_with_own_umis_0` (PAS in the GSM's matrix purely through the
  unified counting, with no clip call of its own in that GSM).
* `all_tier1(u)` (every member clip-seeded), `clip_umis_max(u)`, `n_datasets(u)`, `n_members(u)`,
  `clip_umis_f3844_sum(u)`.
* `in_annotated(u)` = present in `run/annotatedpas.bed` (gene-assigned post-IP set).
* The fraction of multi-member unified PAS whose bed col 5 differs from `clip_umis_sum` (documents why
  col 5 cannot be used directly).

## 3. Everything else in the switch stage is as pre-registered (addendum item 4, manuscript/14)

* `ema switch diff --strategy fisher --count-mode cells --marker-top-n 0 --counts-layer counts
  --cluster-key celltype --fdr 0.05 --no-plots --gtf GRCh38.99 --pasbed <cohort_run>/run/unified/multi_sample_merged.bed`
  (BED6; switch diff reads only chrom/start/end/name/strand from it; it must cover 100 % of the tested
  pas_ids -- asserted at injection, abort otherwise). `--min-cells-per-group` = tool default 10, as in
  manuscript/14 and the design-B driver. `--isoform-agg` = tool default `per_gene`.
* Labels: `final_label` from `labels/confirmed_labels.tsv` (LABEL_POLICY.md), full-id join, unconfirmed
  dropped, types with < 20 confirmed cells in the GSM dropped (`stage3_label_inject.py --min-cells 20`).
* Nulls: N_PERM = 10 label shuffles per GSM, `numpy.random.default_rng(seed=k)` k = 1..10, permutation of
  the cell order applied to the TRUE `celltype` vector, identical `switch diff` call into `<id>/null/perm_XX/`.
* Calibration check per GSM (task rule, anti-conservative direction only, as manuscript/14 arm B0):
  FLAG if pooled null p<0.05 rate > 7 % OR if > 25 % of perms carry any q<0.05 hit (>= 3 of 10). Null
  q<0.05 rate and hits per perm are reported alongside.
* A GSM with < 2 retained types (GSM3516671-MetBrain per `per_gsm_retained.tsv`) has no pair and is
  recorded as `NOPAIRS`, not as a failure.

## 4. Files

`scripts/stage3/stage3_build_universe.py` -> `<cohort_run>/universe/{universe.tsv, universe_ids.txt,
universe.bed, universe_summary.json}` (input md5s, this document's md5, all counts). The injection step
(`stage3_label_inject.py --universe-tsv`) refuses to run without it in cohort mode.

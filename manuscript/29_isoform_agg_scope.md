# Aggregation scope of the switch test: `per_gene` vs `within_utr` vs `between_utr`

Record for the addition to MAIN.md's Methods and the switch-test caveats. Run 2026-09-06.

## What the paper used

Every switch result in the manuscript ran at `--isoform-agg per_gene`, the shipped default.
It was never passed explicitly; the run logs of `stage3_spermatogenesis_v2/*/diff/true`
carry no `--isoform-agg` flag, and `RunConfig.isoform_agg` (config_schema.py:916) defaults
to `"per_gene"`.

Under `per_gene` each PAS is tested against the remainder of its **gene**. A significant
result therefore means "this site's share of its gene changed", which is not the same
statement as "this gene switched between two ends of one 3'UTR".

## The other scopes

- `within_utr` — each PAS against the other PAS sharing its 3'UTR isoform. This is the
  scope that corresponds to tandem 3'UTR APA.
- `between_utr` — PAS collapsed to 3'UTR-level counts, testing which 3'UTR isoform a group
  prefers. Genes with >= 2 annotated UTRs only. A different unit and a different question.

## `within_utr` was unusable in the released tool

Two defects, stacked, on a path with **no test coverage** (`grep -rl within_utr tests/`
returns nothing):

1. `NameError` at `runner.py:526`. `_build_diff_isoform_groups` reads `gene_fallback_map`,
   a name bound only in `run_length`. Every invocation died here. The comment 20 lines
   above documents the deliberate switch to the rank-free `_build_gene_id_map`; the stale
   reference was left behind by that change.
2. `TypeError` in `strategies/fisher.py:173`, reached once (1) is fixed. A spliced 3'UTR
   contributes several `three_prime_utr` records per transcript, so the same PAS was
   appended repeatedly to one UTR group; duplicate columns make `int(agg1[p])` receive a
   Series.

Both fixed in `tools/PeakATail` branch `fix/within-utr-nameerror` (commit `9ceafc6`).
No released PeakATail version can run `within_utr`.

## Agreement, after the repair

Same cells, same labels, same test, changing only the flag. Restricted to the PAS x
stage-pair tested by both scopes, q < 0.05.

| | mouse 1 | mouse 2 |
|---|---|---|
| tested by both | 161,874 | 165,964 |
| `per_gene` significant | 47,246 | 48,945 |
| `within_utr` significant | 44,219 | 46,319 |
| both | 42,953 | 44,855 |
| **% of `per_gene` reproduced** | **90.9** | **91.6** |
| **% of `within_utr` reproduced** | **97.1** | **96.8** |
| `per_gene` only | 4,293 | 4,090 |
| `within_utr` only | 1,266 | 1,464 |

`within_utr` emits more rows overall (202,895 / 209,428 against 175,947 / 179,573) because
a PAS lying in several transcripts' UTRs gets one row per isoform.

Reading: the scope changes **under a tenth** of the calls, but it changes what the call
*means*. Both numbers belong in the paper.

## A comparison that was wrong, and why

An earlier pass compared `per_gene` against `between_utr` at gene level and read 68%
agreement. That number is not a like-for-like comparison and must not be quoted as one:
`between_utr` tests UTR-isoform preference on genes with >= 2 UTRs, a different unit and
question. Worse, an intermediate version read **0%** agreement, because `between_utr`
output leaves `gene_id`, `chrom`, `start`, `end` and `strand` empty and carries gene
identity only in `diff_group_id` — the join was on a blank column. That schema
inconsistency is a third issue to report.

## Deliberately not done

The pre-registered analysis was **not** re-run at `within_utr`. The 20 label-permutation
nulls per mouse that establish the switch test's calibration were all computed at
`per_gene`; changing scope invalidates that evidence and would require the full null suite
re-run plus a pre-registration amendment. The tandem-3'UTR-scoped analysis belongs in the
follow-up project (`/home/biolab/PaperCodes/APA_Reanalysis_Followup`), where the 130-gene
tandem set already lives.

## Artefacts

- `results/figures/manuscript/utr_mode_pas_agreement.tsv` — the table above
- `results/figures/manuscript/utr_mode_comparison.tsv` — the `between_utr` gene-level counts
- runs under the session scratchpad `utrmode/{mouse1,mouse2}/{within_utr_FIX2,between_utr}_*`

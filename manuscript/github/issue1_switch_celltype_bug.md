All 24 `SWITCH_CELLTYPE` tasks in the `RERUN_2026-08_fixed` sweep failed (exit 2, retried once each), so **no per-celltype `switch diff/length/trend` output exists** — the cohort's core APA result is missing. Everything upstream succeeded and is cached.

## Mechanism

`pipeline/main.nf` line 215:

```groovy
def pb  = "${h5.parent}/pasbed.bed"
```

`h5` is a Nextflow `path` input, staged into the task dir under its **bare basename**. A single-segment relative `Path` has no parent, so `h5.parent` is `null` and the command becomes:

```
ema switch diff -i CELL_TYPES_WANSLEEBEN_HOGAN_2013_TREG.h5ad --pasbed null/pasbed.bed ...
→ File 'null/pasbed.bed' does not exist
```

## The fix — and the trap in it

```diff
-    def pb  = "${h5.parent}/pasbed.bed"
+    def pb  = "${cohort}/unified/multi_sample_merged.bed"
```

**Do not "fix" this to `${cohort}/pasbed.bed`** — that fails *silently* instead of crashing. The per-celltype h5ads from `SWITCH_COMBINE` live in the **unified** multi-sample PAS space:

| BED | entries | covers the combined h5ad's 115,450 PAS |
|---|---|---|
| `unified/multi_sample_merged.bed` | 301,860 | **115,450 / 115,450 (100%)** |
| `pasbed.bed` (cohort root) | 16,500 | 16,405 / 115,450 (**14%**) |

The unified ID space itself is sound: two samples' `clusters.h5ad` share 8,851 PAS IDs with **100%** `gene_id` agreement.

## Rerun

One-line edit, then `nextflow run pipeline/main.nf -resume` (same launch flags as `logs/launch.out`). Only the ~24 short `SWITCH_CELLTYPE` tasks re-execute.

## Same-pass secondary fix

The three `ema switch length` calls (lines 226–230) pass no `--gtf`, so per-isoform aggregation silently degrades to `per_gene`. If per-isoform PDUI is wanted for the manuscript, add `--gtf ${params.gtf} --isoform-agg per_isoform`.

Full diagnosis with verification commands: `PeakATail_wd/scripts/pipeline/rerun_switch_celltype_fix.md` on the BioLab box.

@TRextabat — the sweep dir is yours, so the edit + `-resume` is your call; everything above was verified read-only from our side.

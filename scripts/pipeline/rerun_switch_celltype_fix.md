# Fix for the failed SWITCH_CELLTYPE tasks (RERUN_2026-08_fixed)

**Status:** all 24 `SWITCH_CELLTYPE` tasks failed (exit 2, retried once each) in the
2026-08-11/12 sweep. Everything upstream succeeded, so only these ~24 short tasks need to
rerun. Diagnosed 2026-08-12 from `/mnt/ssd2/.../RERUN_2026-08_fixed/` (read-only inspection).

## The bug

`pipeline/main.nf` line 215:

```groovy
def pb  = "${h5.parent}/pasbed.bed"
```

`h5` is a Nextflow `path` input, so it is staged into the task work directory under its
**bare basename** (`CELL_TYPES_WANSLEEBEN_HOGAN_2013_TREG.h5ad`). A single-segment relative
`Path` has **no parent**, so `h5.parent` is `null` and the string interpolates to the literal
`null/pasbed.bed`. Confirmed in the failed task's `.command.sh`:

```
ema switch diff -i CELL_TYPES_WANSLEEBEN_HOGAN_2013_TREG.h5ad --pasbed null/pasbed.bed ...
→ File 'null/pasbed.bed' does not exist
```

## The fix — and the trap in it

```diff
-    def pb  = "${h5.parent}/pasbed.bed"
+    def pb  = "${cohort}/unified/multi_sample_merged.bed"
```

`cohort` is already an absolute path (the sibling `def out = "${cohort}/B3_switch"` resolved
correctly and SWITCH_MATCH succeeded using it).

**Do NOT "fix" this to `${cohort}/pasbed.bed`.** That is the obvious guess and it is wrong —
it fails silently instead of crashing, which is worse. The per-celltype h5ads produced by
`SWITCH_COMBINE` live in the **unified** multi-sample PAS space, verified by ID coverage:

| BED file | entries | covers h5ad's 115,450 PAS |
|---|---|---|
| `unified/multi_sample_merged.bed` | 301,860 | **115,450 / 115,450 (100%)** |
| `pasbed.bed` (cohort root) | 16,500 | 16,405 / 115,450 (14%) |

Using the cohort `pasbed.bed` would silently drop ~86% of PAS and produce a plausible-looking
but badly truncated switch result.

The unified space itself is sound — checked `07_clustering/<GSM>/clusters.h5ad` for two
samples: 8,851 shared PAS IDs, **100%** agree on `gene_id`, so `new_pas_id` is genuinely
consistent across datasets and `SWITCH_COMBINE` merged correctly.

## Rerun

```bash
cd /mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed
# edit pipeline/main.nf line 215, then:
nextflow run pipeline/main.nf -resume -c pipeline/nextflow.config   # match the original launch flags in logs/launch.out
```

Everything upstream is cached; only the 24 `SWITCH_CELLTYPE` tasks should re-execute.

## Secondary issue worth fixing in the same pass

Lines 226–230: the three `ema switch length` calls do **not** pass `--gtf`, so per-isoform
aggregation is unavailable and the run silently uses `per_gene`. July's
`check_pas_in_utr.py` already flagged this. If per-isoform PDUI is wanted for the manuscript,
add `--gtf ${params.gtf}` (and `--isoform-agg per_isoform`) to those calls.

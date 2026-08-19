## What this fixes

Count-matrix rows are keyed by **list position** while `pas_ids` is **compacted** (empty rows removed), so the two indexings drift by the number of empty rows. Measured on testis mouse1 (98,574×9,564, 21.7M nnz): **99.87% of annotated PAS carry another PAS's counts**; Spearman(own depth, shipped counts) = **0.004**. Regression bisected to `04e0b3a` (2026-03-23) — every count-derived output since then is affected (annotated_matrix.mtx, preprocessed/clusters.h5ad, differential tests, PDUI). Coordinate-only outputs (pasbed.bed, F1 accuracy scoring) are NOT affected.

## The fix

`make_dataframe()` subsets the matrix to `pas_ids` so row *i* **is** `pas_ids[i]` (verified on the real 274 MB matrix: 97,034/97,034 rows now carry their own counts, +3.2 s). Every boundary hardened to make the bug class impossible: `annotate()`, `preprocessing()`, `write_annotated_matrix()`, `write_pas_gene_artifacts()` all raise on shape/length mismatch. Note: shape checks alone could NOT catch this — the shipped matrix had the right row *count* with the wrong row *content*.

## Tests

`tests/test_matrix_pas_id_row_alignment.py`: 8 tests, **6 fail on develop / 8 pass here**; suite otherwise identical to develop (7 pre-existing failures, zero new). Adversarially re-verified by an independent pass (SOUND).

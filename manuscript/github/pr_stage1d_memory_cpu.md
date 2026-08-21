# perf: 294 GB → 12 GB and 3h46 → 28 min on the PBMC 10k run, byte-identical output

**Base: `feat/polya-evidence` (merge after it).** Branch: `perf/clip-memory`, 8 commits on top of `4efeb12` (5 changes + tests/docs + one verifier repair).

This PR changes **no output**. Every file the pipeline writes is byte-identical
at the default settings, on three real datasets, verified by `cmp` against the
runs this PR is supposed to reproduce. What changes is that the final PBMC 10k
benchmark run now needs **12.4 GB instead of 293.7 GB** and **27 min instead of
3 h 46 min**, and `--threads` finally does something for peak calling.

---

## 1. The problem

The final PBMC 10k v3 run (CellRanger BAM, `--threads 16`, `clip_seeded`,
issue #10 §2 / `15_final_gate.md` §5) peaked at **293.7 GB RSS** with a wall of
**3 h 46 min** at ~143 % CPU — i.e. effectively single-threaded despite
`--threads 16`. The shipped caller used 105.6 GB / 3 h 27 min on the same BAM,
so Stage 1b/1c looked like it had doubled memory. Mouse (STARsolo, 10k cells)
ran at 23.1 GB / 1 h 01 min.

Two things had to be true at the end: the same numbers out, and a footprint
back at or below the pre-1b 140 GB.

## 2. What the profile actually found

**The 293.7 GB was never in the caller.** It was one function in the clustering
stage, `ema/clustering/strategies/leiden_tfidf.py::_tfidf_signac_method1`, which
did `X.toarray()`, `.astype(float64)`, `tf = X / cell_totals`,
`np.log1p(tf * idf * scale)` on the cells × PAS matrix: four concurrent dense
float64 copies. The peak RSS of **every** recorded run equals
`4 × n_cells × n_PAS × 8 B` plus a small baseline:

| run | matrix | predicted | measured |
|---|---|---|---|
| shipped caller | 11,836 × 275,370 | 104.3 GB | 105.6 GB |
| pre-1b clip-seeded | 20,737 × 212,082 | 140.7 GB | 140.2 GB |
| final (this PR's target) | 23,303 × 390,493 | 291.2 GB | 293.7 GB |
| final, IP arm | 23,302 × 317,574 | 236.8 GB | 239.5 GB |
| mouse1 | 10,339 × 65,075 | 21.5 GB | 23.1 GB |

Confirmed by replaying `clustering()` alone on the final run's
`preprocessed.h5ad`: 287.8 GB in 903 s.

So the "1b/1c doubled memory" was **indirect**: the tier-1 counting fix puts all
read ends (not the ~1 % clip reads) into the matrix, so 23,328 instead of 20,737
cells pass `min_read=1500` and 390 k instead of 212 k PAS survive `min_cells=3`
→ matrix area ×2.07 → RSS ×2.1. Inside the caller itself 1b/1c cost ~0.4 GB on
a gene-dense chromosome.

**The wall time was the two full BAM scans.** 3 h 45 min 53 s decomposed as
`+` pass 5626 s + `-` pass 5221 s (80 %; each iterates all 734 M records at
7–8 µs/record) + CB filter 738 s + annotation 73 s + clustering 1875 s. Peak
calling is one pure-Python main thread; `--threads` was never consulted on the
default path (it only fed `ResourceManager.get_n_jobs` for `--tiles` and the
downstream pool), and the 143 % CPU was that thread plus htslib BGZF threads.

## 3. Changes

| # | Commit | What |
|---|---|---|
| 1 | `a95ea97` | **Sparse TF-IDF.** Same scalar sequence (`/ cell_total`, `* idf`, `* scale`, `log1p`) applied to the stored entries; `~3 × nnz × 8 B` instead of `4 × cells × PAS × 8 B`. |
| 2 | `a6261b5` | **Per-(contig, strand) parallel peak calling** with a deterministic merge (`ema/countmatrix/chrom_parallel.py`, `--peak-workers`). |
| 3 | `4c86442` | **`read_check` reorder** (strand test before the CB tag lookup and the CIGAR walk) + interned `"<RG>_<CB>"` composite. |
| 4 | `5420398` | **CB filter streams integers** (12 B/non-zero) instead of a pandas frame with a per-row `cb_str` object column (~95 B/non-zero). |
| 5 | `5c0e891`, `0a3bed0`, `bf9310d` | Docs (performance section, `--threads` / `--peak-workers` semantics), CHANGELOG, strategy-kwargs test. |
| 6 | `4dff7c3` | **Verifier repair:** the parallel merge now seeds `pas_id` from `Peak.pasnumber` and writes it back, so the second BAM of a multi-dataset run continues the first one's numbering exactly as the legacy loop does (it restarted at 1; harmless downstream, which keys on `(dataset_id, strand, pasnumber)`, but not byte-identical). Regression test runs a second "dataset" both ways. Single-BAM output unchanged (re-verified on the slice and mouse1 after the change). |

### Why each is bit-exact, not approximately equal

**1. Sparse TF-IDF.** Every structural zero maps to `log1p(0 × idf × scale) == 0`
and was dropped by the old trailing `csr_matrix(tfidf)` anyway, so the sparsity
pattern is identical after `eliminate_zeros()`. The per-entry arithmetic is the
same operations in the same order. The row sums and per-PAS cell counts are sums
of integer counts, exact in float64 in any summation order. Checked against the
old implementation on the real `preprocessed.h5ad` of two runs — CSR structure,
`data`, `X_lsi` and Leiden labels all identical (0.14 s / 0.56 GB vs 3.1 s /
4.5 GB on the slice; 0.76 s / 1.55 GB vs 17.4 s / 22.2 GB on mouse1).

**2. Parallel peak calling.** The streaming state machine is reset at every
chromosome change and the clip-seeded emitter flushes per chromosome, so one
job = `peak_calling(direction, region=(contig, 0, length))` sees exactly the
read sequence the sequential pass saw for that contig. The merge then reproduces
the sequential emission order exactly:

* `+` strand first, contigs in BAM header order, then `-` strand;
  `pas_id` renumbered 1..N across both strands;
* the shared barcode index is rebuilt by walking each job's local `cb.tsv` in
  local first-write order and appending unseen barcodes — the sequential
  `BarcodeIndex` singleton also assigns ids at write time, so the column
  assignment is the same;
* matrices and support sidecars are re-keyed through those two maps.

Unlike `--tiles`, every config value the caller reads (`barcode_tag`, `cb_len`,
`seqlen`, `ignore_chro`, thresholds, strategy kwargs) travels in the job spec
rather than being expected in the spawned child's module state.

**3. `read_check`.** Every rejection returns the same `(0,0,0,0,0)` sentinel, so
the checks may run in any order; only *which* check fires first changes. The
interned composite is the same string value.

**4. CB filter.** Row semantics (blank/comment skipping, <3-token rows dropped,
per-file dimension-header detection, stable sort, contiguous re-index, output
format) are unchanged. The Pass-1 grouping moves from the raw token to its
integer — the same partition while every barcode token is canonical — so a
non-canonical token (`"007"`, `"+7"`, `"7.0"`), a >3-token row, a non-integer
token or an absurd column index defers the whole call to `_filter_cb_legacy`,
the byte-for-byte copy of the original algorithm. On such pathological input the
result is now the reference's, which the previous fast path only approximated
(it truncated float-valued tokens instead of skipping the row).

## 4. Before / after

`/usr/bin/time -v`, same machine, same inputs, same flags as the runs being
reproduced; peak RSS is the largest single process.

| Run | Wall before | Wall after | CPU before | CPU after | Peak RSS before | Peak RSS after |
|---|---|---|---|---|---|---|
| PBMC chr19+21 slice, `--threads 16` | 14 min 41 s | **6 min 11 s** | 148 % | 270 % | 5.56 GB | **1.16 GB** |
| GSE104556 mouse1 testis, `--threads 12 --ip-filter` | 1 h 01 min 09 s | **9 min 03 s** | 140 % | 594 % | 23.12 GB | **3.68 GB** |
| PBMC 10k v3 full BAM, `--threads 16` | 3 h 45 min 53 s | **27 min 43 s** | 143 % | 689 % | 293.74 GB | **12.45 GB** |

The gains are not only parallelism. Re-running the slice with
`--peak-workers 1` — the legacy single-process caller, same code — gives
11 min 27 s / 1.71 GB against 14 min 41 s / 5.56 GB: peak calling 733 s → 613 s
(`+` 439 → 391 s, `-` 294 → 222 s, from the `read_check` reorder) and the CB
filter 64 s → 25 s.

Stage breakdown of the full PBMC run, before → after:

| stage | before | after |
|---|---|---|
| peak calling (both strands) | 10,847 s, 1 thread | 674 s in 252 jobs on 16 workers + 174 s merge |
| CB filter | 738 s | 490 s |
| annotation | 73 s | 54 s |
| clustering | 1,875 s | 249 s |
| **peak RSS / where** | 293.7 GB, dense TF-IDF | 12.4 GB, clustering (no caller worker above 3.0 GB) |

The 120 GB acceptance target is met by a factor of ~10; the run is also below
the 105.6 GB of the *shipped* caller and the 140 GB of the pre-1b clip-seeded
one, on a matrix twice their area.

## 5. Identical-output proof

Byte comparison (`cmp`) of the new run against the run each reproduces:

* **Full PBMC 10k v3** vs `results/benchmark_tools/pbmc_10k_v3/peakatail_clipseeded_final/run`
  — **18/18 identical**: `peakcalling/default_0.{pos,neg}.{bed,mtx}`,
  `default_0.{pos,neg}.support.tsv`, `default_0.cb.tsv`, `posbed.bed`,
  `negbed.bed`, `pasbed.bed`, `pas_support.tsv`, `posmatrix.mtx`,
  `negmatrix.mtx`, `filterdmatrix.mtx`, `filtered_cb.tsv`, `annotatedpas.bed`,
  `annotated_matrix.mtx`, `pas_gene.tsv`. Plus `preprocessed.h5ad` and
  `clusters.h5ad`: `X` (CSR `data`/`indices`/`indptr`), `X_lsi` and the Leiden
  labels all equal.
* **Mouse1 testis (`--ip-filter`)** vs
  `results/benchmark_tools/gse104556/peakatail_clipseeded_final/mouse1/run`
  — **18/18 identical**, same list, plus both h5ads (27 clusters, identical
  labels). This arm also exercises the internal-priming filter
  (`n_ip_flagged=37118`, same as the reference).
* **PBMC chr19+21 slice** vs the pre-change run of the same code
  — **21/21 identical** (the 18 above plus `01_peak_calling/default/raw/pas.bed`,
  `annotated_cells.tsv`, `annotated_pas_ids.tsv`), plus `clusters.h5ad`
  including the `counts` layer.

The same slice run with `--peak-workers 1` (legacy caller through `main.py`)
also reproduces the reference bytes on the 12 files checked — the fallback
path is not a second implementation of the outputs.

Unit tests (`pytest tests`: **1257 passed**, 2 pre-existing environment failures
in `tests/test_pyproject_install.py`, unchanged from the baseline's 1215):

* `tests/test_chrom_parallel_identity.py` — a synthetic 4-contig BAM (both
  strands, clipped and bare loci, an ignored contig, reads without RG,
  per-contig barcode orders) called sequentially and in parallel with 2–3
  workers, for `lambda_gradient` and `clip_seeded`: the five pipeline files and
  both sidecars compared as bytes; the merge helpers pinned on hand-made parts;
  `get_mapping()` order after the merge; the no-poly(A) case.
* `tests/test_tfidf_sparse_identity.py` — sparse vs the old dense formula on
  random count matrices (explicit zeros, duplicate/unsorted indices, empty rows
  and columns, dense and non-CSR input, three scale factors), and a guard that
  the input is never densified.
* `tests/test_read_check_fastpath.py` — the full grid of read shapes against a
  verbatim copy of the previous implementation, both directions.
* `tests/test_matrixfilter_vectorized.py` — randomised multi-file matrices ×
  `min_read ∈ {0, 1, 40}` against the reference implementation (bytes of both
  outputs and `filtered_cb_list`), a 7-row chunk size, and each deferral case.

## 6. Limits and what is deliberately not in this PR

* **The merge is single-threaded** (174 s on PBMC: it rewrites 2.8 GB of
  MatrixMarket text). It is the obvious next target if peak calling needs to go
  below 10 min.
* **Everything after peak calling is still single-process.** `make_dataframe`'s
  `mmread` → CSC → CSR → CSC chain and `preprocessing`'s transpose are the
  remaining multi-GB copies (they are inside the 12.4 GB peak, together with
  clustering).
* **One pass for both strands** (two state machines per contig job) would halve
  the decode + `read_check` work again. It is a larger refactor of
  `peak_calling` and is left for a follow-up.
* **`--tiles` is untouched and still unusable for a real run**: its workers
  crash under spawn (config not propagated — also on `develop`), it builds jobs
  for all 194 header contigs, and `merge_tiles` sorts contigs as strings and
  unions barcodes per strand, leaving `neg.mtx` columns mismatched against the
  combined `cb.tsv`. `--pipeline` remains unreachable from the CLI. Both deserve
  their own issue: fix or remove.
* **`--peak-workers` needs a BAM index.** Without a `.bai` (or with
  `--peak-workers 1`) the legacy single-process caller runs, unchanged, and logs
  why.
* Per-worker memory scales with the deepest contig, not with `--threads`:
  budget ~2.5 GB per worker (measured max 3.0 GB on PBMC chr11 (+)).

Refs #95 (§2 memory / CPU). Base `feat/polya-evidence` landed in `develop` as #93; this branch applies cleanly on `develop` and is independent of #96 (IP-filter strand fix).

_Whole-process-tree memory (verifier, PSS/RSS sampler): PBMC 19.3 GB RSS-sum / 18.1 GB PSS-sum with the parent at 12.2 GB; mouse1 7.3 GB; slice 3.1 GB — the per-process `ru_maxrss` numbers above are the parent only._

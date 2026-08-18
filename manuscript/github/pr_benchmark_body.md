## Motivation

The current benchmark (`compute_metrics`) scores the **whole predicted interval** (mean peak width ~821 bp) against the reference with **no strand requirement**, via `predicted.window(reference)`. A null-control analysis on the `RERUN_2026-08_fixed` grid showed the consequence: against the 18.4M-entry PolyASite reference (which tiles 33% of the genome within 100 bp), **shuffled peaks score 0.62–0.75 precision** under the identical procedure — the reported 0.9986 sits on that floor. Scored at the peak's 3′ base with strand matching, real PAS hold **0.994** while the null drops to **0.35–0.46**, so the discriminating benchmark is point-mode + strand-matched, and it must ship with its null. Closes #__BENCH_ISSUE__.

## What changed

- `compute_metrics(...)` gains opt-in kwargs, **all defaulting to legacy behavior**: `point_strand=`, `restricted_reference_bed=`, `null_shuffle=`, `null_n_seeds=`, `null_base_seed=`, `null_gene_body_bed=`.
- New `"point_strand"` JSON block: strand-aware 3′-most-base scoring (BED end exclusive: `+` → end−1, `−` → start), boundary-inclusive distance ≤ cutoff, same-chrom + same-strand matching; the recall-like quantity is named `reference_coverage`, with optional `recall_restricted` against a caller-supplied restricted reference.
- New `"null"` block: N-seed width/chromosome/strand-preserving shuffled baseline (numpy RNG, deterministic per seed), optionally constrained to gene bodies, scored identically; per-seed + mean.
- Standalone exports `compute_point_strand_metrics(...)`, `compute_null_baseline(...)`; `run_benchmark(...)` passes the kwargs through (its `summary.csv` unchanged unless opted in).
- Pure Python/numpy for the new path — no new dependency; the legacy pybedtools path is untouched.

## Back-compat

Legacy top-level keys (`n_predicted`, `n_reference`, `cutoffs`) are **byte-identical** with all new options enabled (asserted by test via sorted-JSON equality). `sweep_analysis.py` readers are unaffected.

## Tests

`tests/test_benchmark_point_strand.py`: 17 new tests — ± strand point reduction, the end-exclusive off-by-one, strand-mismatch rejection at distance 0, exact cutoff boundary, null determinism/width preservation, JSON serializability, legacy byte-compat (bedtools integration).

```
42 passed, 0 failed  (17 new + 25 existing benchmark tests; module path asserted to resolve to this branch before every run)
```

An independent adversarial review pass (separate agent, own fixtures) confirmed the strand/off-by-one behavior at the 258/259 bp boundary and null placement never escaping gene bodies.

## How the manuscript uses it

The sweep's benchmark step passes `point_strand=True, null_shuffle=True, null_gene_body_bed=<gene_end.bed>` so every published precision number carries its null alongside — reviewers get the honest comparison by construction.

---
*Opened by Claude (AI assistant) for the manuscript effort; the gh session on this machine authenticates as @yasinkaymaz.*

🤖 Generated with [Claude Code](https://claude.com/claude-code)

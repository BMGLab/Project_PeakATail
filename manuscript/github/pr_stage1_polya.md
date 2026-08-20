## What this adds

Poly(A) read evidence in the caller — the fix the benchmark demanded (PeakATail was last among de novo tools on PBMC; the winners all use orthogonal evidence).

**Phase 1**: `polya.py::clip_site()` — strand-aware terminal soft-clip A/T detector (92× wrong-end specificity, reproduced on 365k real reads with 0 mismatches vs an independent re-implementation). Per-PAS clip-read counts land in the previously-unused BED score column; `--polya-*` flags mirror the `ip_*` block; wired identically across the monolithic, pipeline, and tile paths (equality-tested).

**Phase 2**: new registered strategy `clip_seeded` — clip clusters (25 bp single-linkage, read-weighted mode) are primary candidates; non-overlapping coverage peaks from the delegated default strategy are an explicit **second tier** (score==0 is an exact tier tag; `--polya-mode filter` drops exactly that tier).

## The measured result (chr19+21, identical scoring path, control reproduces the published numbers exactly)

| arm | n | P@100 | R@100 | F1@100 |
|---|---:|---:|---:|---:|
| shipped (lambda_gradient) | 12,928 | 0.161 | 0.196 | 0.177 |
| **clip_seeded, both tiers** | 19,354 | 0.237 | 0.378 | **0.291** |
| **tier 1 (clip-supported)** | 11,318 | 0.359 | 0.338 | **0.348** |
| tier 2 (coverage-only) | 8,036 | 0.065 | 0.040 | 0.050 |

Tier-1 recall sits at the measured evidence ceiling (~34% of atlas sites carry clip reads on this slice) — the caller extracts essentially all the recall this evidence type provides. Chemistry gate passed on STARsolo/mouse too (1.65% clip rate, 540× specificity), so this is not CellRanger-specific.

## Honest limits (carried into any manuscript text)

Internal-priming filter not yet applied on human (tier-1 precision is not final); tier-2 F1 is poor and the tiers must be reported separately downstream; UMI dedup of clip stacks is partial. Full-genome verdict comes from the Stage 2 re-run (pre-registered gate: F1 > 0.261 on full PBMC).

## Rebased onto current develop (post #81–#91)

Rebased cleanly onto `a3ddf4e` (only CHANGELOG needed a manual merge). One deliberate interaction with #90's `--auto-cleavage-offset`: that correction shifts every PAS 3′ end downstream by the estimated R2-read-length offset, but **clip-seeded tier-1 PAS are already placed at the observed cleavage site** — shifting them would push them past it. `rewrite_bed_3prime_offset` therefore gains `skip_supported=`, and `main.py` sets it whenever the strategy is `clip_seeded`, so only the coverage-only tier (score 0) is shifted. Covered by `tests/test_cleavage_offset_skips_clip_supported.py` (both tiers, both strands, zero-offset no-op). Base: develop. +74 tests.

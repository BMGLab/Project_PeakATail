# Issue: Stage 1d — IP-filter minus-strand window, memory footprint (280 GiB on PBMC), and three output-bookkeeping defects

**Labels:** bug, performance, Stage-1d. Found by the independent verifier of the final benchmark run (code `4efeb125`, 2026-08-21).

## 1. Internal-priming filter tests the wrong side on the minus strand
`ema/experimental/internal_priming.py` uses one forward-coordinate window `[pos-10, pos+30)` for both
strands. On `+` that is 10 bp upstream / 30 bp downstream of the cleavage site (correct: the
genomic A-stretch that primes oligo-dT sits downstream). On `-` the same genomic window is 30 bp
*upstream* / 10 bp downstream in transcript orientation, so the filter mostly tests upstream T-runs.
Re-applying a strand-symmetric rule to the survivors of the final run flags 3.5% (PBMC), 4.2% / 3.6%
(mice) of minus-strand sites that passed, and the current rule removed ~2,500 PBMC sites the
symmetric rule would keep. Benchmark effect is negligible (P@100 0.7167 → 0.7182) but the rule is
wrong. Fix: reflect the window on `-` (`[pos-30, pos+10)` genomic) and scan the reverse complement;
add a minus-strand unit test with an implanted downstream (transcript-wise) A-run.

## 2. Memory / CPU
Final PBMC 10k run: peak RSS **293.7 GB (280.1 GiB)** without IP filter, 239.5 GB with, wall 3:45 h,
~1.4 cores used despite `--threads 16`. That is 2.1× the pre-1b run (140 GB) and 2.8× the shipped
caller (105.6 GB); the tier-1 quantification changes (1b/1c) roughly doubled memory. Mice: 23 GB.
Please profile `ema/countmatrix/polya.py` ClipStream/ClipAccumulator (+ `{end1: {cb: n}}` maps) and
the candidate-partition path; likely wins: stream per chromosome, drop per-end1 maps after flush,
use array-backed counters. Also: peak calling is effectively serial — `--threads` only helps
downstream.

## 3. Output bookkeeping
- In `--ip-filter-mode filter`, `annotatedpas.bed` has no per-site IP column (9 columns); only
  `04_pas_gene_assignment/peak_filters_stats.json` counts. Downstream "trusted novel PAS" needs a
  per-site flag: write it in both modes (or always emit `ip_flags.tsv`).
- `pasbed.bed` silently drops gene-assigned PAS with zero counts after the CB filter
  (PBMC 1,814 = 0.45%; mouse1 2,809 = 3.7% of `annotatedpas.bed`). Document, or record them with
  count 0.
- `ema/cli/run.py` logs "--seq-len / yaml seqlen not set; defaulting to 150" before the Click value
  is merged, although the resolved config has 91/98. Misleading.
- In multi-dataset runs `unified/multi_sample_merged.bed` col 5 is the *first* merged member's
  clip-molecule count (`bedtools merge -o first`), not cohort support; a unified PAS can merge a
  tier-1 member with a tier-2 member. Please emit cohort-level `clip_umis_sum` / `any_tier1` per
  unified PAS (we currently recompute from the per-dataset `*.support.tsv` via
  `multi_sample_pas_mapping.tsv`).

## Evidence
`manuscript/15_final_gate.md` §5–6; verifier outputs `results/benchmark_tools/final_verify/`.

## 4. Addendum (2026-08-21, Kinnex validation) — the IP rule is too loose
Applying the pre-registered "trusted novel PAS" funnel to the final PBMC default output and validating
against Kinnex long reads: 24.5% of atlas-novel, hexamer-bearing, IP-filter-passing sites lie within
25 bp of a Kinnex *internal-priming decoy* terminus (6.6% for atlas-known sites). The current rule
(window −10..+30, ≥6 consecutive A or ≥70% A) misses a large internal-priming class that the Kinnex
rule (+1..+18 downstream, ≥12 of 18 A) catches. Proposal: add `--ip-rule {legacy,kinnex}` (default
unchanged for now; to be pre-registered and re-benchmarked before becoming default), implemented on
the strand-corrected window from §1. Details: `manuscript/16_trusted_novel_kinnex.md`.

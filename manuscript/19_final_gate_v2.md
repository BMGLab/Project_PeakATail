# v2 benchmark — verified gate results on the merged fixes (code `9dfdefb`, run 2026-08-21 16:14–16:50)

**Verifier verdict: FIXED** (independent recomputation of every number; zero unexplained movement; both
Stage-1d PRs did exactly what they promised). Full verified table:
`results/benchmark_tools/final_v2_verify/VERIFIED_v2.md`. Code = `develop` after #93 (Stage 1b+1c),
#96 (IP-filter strand fix), #97 (performance); frozen worktree `tools/pa-polya-run-9dfdefb3`; suite
1,277 passed. Gates unchanged from their pre-registrations ([13](13_reliability_positioning.md) §1,
[10](10_caller_fix_plan.md) §5).

## 1. Pre-registered precision default — gate P@100 ≥ 0.50: **PASS on all three**

| dataset | n scored | P@100 | R_det@100 | F1_det | v1 → v2 |
|---|---:|---:|---:|---:|---|
| PBMC 10k v3 (IP arm) | 46,524 | **0.7062** | 0.1754 | 0.2811 | 0.7167 → (−1.1 pp P, +0.5 pp R, +2,130 sites) |
| testis mouse 1 | 26,255 | **0.7450** | 0.2048 | 0.3213 | 0.7414 → (+0.4 pp) |
| testis mouse 2 | 26,526 | **0.7572** | 0.2080 | 0.3264 | 0.7554 → (+0.2 pp) |

Verified to 6 decimals by an independent bedtools pipeline; nulls 33–55× below; every movement equals the
measured #96 impact (residual 0: '+' strand 0 changes; '−' strand 5,599 newly removed / 11,217 newly kept
on PBMC — the verifier's own re-implementation reconstructs the v2 call set key-for-key). The IP filter on
mouse 1 buys +5.3 pp precision for −1.2 pp recall (v1 wording "+4.9/−1.4" superseded).

## 2. Original two-sided gate (PBMC ≥1-molecule outputs) — unchanged, FAIL on P

Tier-1 ≥1 mol: F1_det 0.3001 pass / P 0.3032 fail (no-IP; byte-identical to v1 per #97's identity promise —
19/19 files cmp-equal). IP arm tier-1: P 0.3520 / R_det 0.2685 / F1 0.3046. The ≥1-molecule outputs remain
the sensitivity arm.

## 3. Compute (production confirmation of #97)

All four arms **concurrently**: PBMC 34:37 / 12.53 GB and 32:59 / 11.12 GB; mice 14:10 / 3.07 GB and
16:00 / 3.65 GB (v1: 3:45:53 / 293.7 GB; 1:01 / 23 GB). Uncontended single-run figures from the PR:
27m43s PBMC. Quote peak RSS as the production number and disclose the concurrency when quoting wall time.

## 4. Verifier notes carried forward

- Each arm emits two identical ≥2-molecule files (POSTHOC + PRESPEC names) — harmless, tidy later.
- The cosmetic "--seq-len defaulting to 150" log line persists (issue 10 §3 item; resolved config correct).
- §4-of-15's "72.5% single-molecule" is 72.1% by the verifier's recount (definition nit).
- Kinnex v2 (builder-crosschecked, adversarial check running): trusted-novel negative result stands and
  slightly strengthens — concordance 0.485 (v1 0.521), decoy fraction 27.0% (24.5%); the corrected IP
  filter flags fewer raw peaks, confirming the leniency of the *rule* (not the strand bug) as the cause.
  `--ip-rule kinnex` (issue #95 addendum) remains the path to any future trusted-novel v2 definition.

**Paper numbers now come from this run.** 15 remains the v1 record; every figure/text number moves to v2.

## 5. Second production confirmation of #97 — the 17-sample Laughney cohort

The Stage-3 cohort run (17 datasets, ~224 GB of BAM, identical YAML and flags, one unified PAS space)
was executed on both codes on the same box:

| | v1 (`4efeb125`) | v2 (`9dfdefb`, `--peak-workers 16`) |
|---|---|---|
| wall | 9 h 06 m 26 s | **1 h 06 m 26 s** (8.2×) |
| peak RSS | 13.63 GB | 13.53 GB |
| unified PAS | 505,197 | 505,197 |

Same output scale at 1/8 the wall time. The smoke test isolates the other fix cleanly: on a chr21
2-dataset slice the plus strand is byte-identical to v1 while minus-strand internal-priming removal
falls 21.9% → 18.4% — i.e. exactly the #96 strand correction and nothing else.

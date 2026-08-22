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
- Kinnex v2 (builder-crosschecked, adversarially verified): trusted-novel negative result stands and
  slightly strengthens — concordance 0.485 (v1 0.521), decoy fraction 27.0% (24.5%); the corrected IP
  filter flags fewer raw peaks than the buggy one did. **CORRECTED 2026-08-21:** the final clause of this
  bullet previously read *"`--ip-rule kinnex` (issue #95 addendum) remains the path to any future
  trusted-novel v2 definition"*. **That is withdrawn — it is measured dead.** On the no-IP arm at
  ≥2 molecules the *shipped* rule removes 15,586 sites and lifts P@100 0.5907 → **0.7062**, while the
  Kinnex rule removes 900 and lifts it only to **0.5944**; genome-wide the shipped rule flags **68,855 of
  402,860 peaks (17.09 %)** against the Kinnex rule's **12,186 (3.02 %)** — the shipped rule is
  **5.65× stricter**, so the Kinnex rule is a no-op on top of it (31 of 46,524 sites) and a regression as
  a replacement. The residual decoy proximity is therefore **not** attributable to rule leniency in the
  sense that swapping rules would fix. `results/perf_gap/D2_rule_sweep/` §2(5) and
  `results/perf_gap/VERIFY/README.md` correction 1 (the ratio itself: [22](22_performance_roadmap.md)
  §6 item 5 — do **not** quote the "~17×" that appears in the verifier's own prose);
  [16](16_trusted_novel_kinnex.md) §"IP-rule note", [25](25_competitive_position.md) §10.

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

---

## 6. Competitive framing — ADDED 2026-08-21 (verified; see [25](25_competitive_position.md))

The gate numbers above are unchanged. What changes is how they may be *compared*.

The standing summary line built on this table — *"the precision default is the most atlas-concordant de
novo call set in the panel, at polyApipe-or-lower recall (−12 % PBMC, −18 % / −17 % mouse), F1 essentially
tied"* — is **misleading as a statement about the tool**: it sets **46,524** of our calls against
**120,916** of polyApipe's, and both precision and recall move with call budget.

Verified replacement (`results/perf_gap/D3_head_to_head/`, `D2_rule_sweep/`, `D1_recall_decomposition/`;
verifier verdict **FIXED**, `results/perf_gap/VERIFY/README.md`):

* **At matched call count PeakATail leads every de novo competitor on recall at every N tested, and on
  precision at every N from 20,320 upward, on both datasets.** At polyApipe's own N = 120,916:
  **0.4036 / 0.2336** vs **0.3800 / 0.1988**. At our own N = 46,524: **0.7062 / 0.1754** vs
  **0.6049 / 0.1375**. Mouse 1 at N = 71,983: **0.4488 / 0.3034** vs **0.4353 / 0.2231**. Below
  N ≈ 15,000 polyApipe's precision is the higher of the two (0.9671 vs 0.9527 at N = 9,456), so the
  precision half of this claim carries a lower bound on N and the recall half does not
  ([25](25_competitive_position.md) §2.3.3).
* **On both mice the ≥1-molecule arm in §2 above already beats polyApipe on precision, recall and F1
  simultaneously** — 0.5686 / 0.2806 / 0.3758 (mouse 1) and 0.5880 / 0.2826 / 0.3818 (mouse 2) against
  0.4005 / 0.2499 / 0.3078 and 0.4119 / 0.2508 / 0.3118.
* **The ≥2-molecule default is not the F1 optimum** on any of the three datasets (the ≥1-molecule arm
  reaches 0.3046 / 0.3758 / 0.3818 against the default's 0.2811 / 0.3213 / 0.3264). It is a
  reliability-first choice and the paper must say so rather than let a reader infer 0.2811 is the tool's
  best.
* **The recall gap is mostly unreachable:** 79.96 % of the 235,111 missed PBMC atlas sites have no peak of
  any kind within 100 bp (mouse 84.53 %), so only 20.04 % / 15.47 % is reachable by any threshold, filter
  or tier policy, and the hard ceiling downstream of peak calling is R_det **0.3406** / **0.3278**.
* **The denominator is inflated:** only 44.6 % of the detected-gene atlas is corroborated by PBMC Kinnex
  long reads at ≥5 UMI — a **donor-mismatched, site-level** panel ([16](16_trusted_novel_kinnex.md)), so
  "not corroborated" is not "not real"; against the x3p-corroborated subset the same default recalls
  **0.4372**, not 0.1754. The 44.6 % is the **union** count (127,188 / 285,136) and 0.4372 is against the
  **x3p** subset (107,260): two separate facts, never divided across each other
  ([25](25_competitive_position.md) §5).

**Honesty rule added with this section: never quote a single-point head-to-head without the matched-N
context** — see [25](25_competitive_position.md) §9. Qualifiers: the claim is among **de novo** tools
(catalog-based scUTRquant is above PeakATail on both axes at several matched N); the precision half is
not established for 20,320 ≤ N ≤ 35,759 because of competitor-side ties, and is **false** below
N ≈ 15,000 ([25](25_competitive_position.md) §2.3.1, §2.3.3).

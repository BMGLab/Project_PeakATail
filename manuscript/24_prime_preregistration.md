# 24 — `peakAtail-prime`: pre-registration of the algorithmic-improvement comparison

**Written 2026-08-21, BEFORE any accuracy number from the `peakAtail-prime` branch exists.**
At the time of writing the branch is `9dfdefb` with zero behavioural commits on top; the only
artefacts that exist are the dev harness (`scripts/prime/`), the baseline test-suite counts and a
byte-identity reference run of the *unmodified* branch. Nothing in this file was informed by a prime
result, because there is none.

This document exists because the branch is allowed to change the caller's defaults. Once defaults
move, "did it get better?" stops being answerable after the fact — so the answer is fixed here,
first.

---

## 1. The four code versions being compared

| name | commit | what it is | where its outputs live |
|---|---|---|---|
| **shipped** | pre-Stage-1 `develop` | the caller as distributed before the Stage-1 poly(A)-evidence work: default coverage strategy (`lambda_gradient`), and the first `clip_seeded` arm | `results/benchmark_tools/pbmc_10k_v3/peakatail/` (lambda_gradient, n=277,111) and `.../peakatail_clipseeded/` (clip_seeded, n=377,236); `results/benchmark_tools/gse104556/peakatail{,_clipseeded}/` |
| **v1** | `4efeb125` | Stage-1a/1b caller; the numbers in manuscript [15](15_final_gate.md) | `results/benchmark_tools/pbmc_10k_v3/peakatail_clipseeded_final{,_ipfilt}/`; `results/benchmark_tools/gse104556/peakatail_clipseeded_final/mouse{1,2}/`; frozen tree `tools/pa-polya-run-4efeb125` |
| **v2** | `9dfdefb` | Stage 1b+1c + IP-strand fix (#96) + performance (#97). **The current manuscript numbers** ([19](19_final_gate_v2.md)) | `results/benchmark_tools/pbmc_10k_v3/peakatail_clipseeded_final_v2{,_ipfilt}/`; `results/benchmark_tools/gse104556/peakatail_clipseeded_final_v2/mouse{1,2}/`; frozen tree `tools/pa-polya-run-9dfdefb3` |
| **prime** | branch `peakAtail-prime`, cut from `9dfdefb` | this work | code: `tools/pa-prime`; runs: `/mnt/ssd0/emaout/peakatail_benchmark/prime/<label>/`; small artefacts: `results/prime/` |

The PI wants a **three-way** comparison (shipped / v2 / prime), with v1 available as the intermediate
step. That comparison is only meaningful if the baselines can be regenerated, so:

> **Reproducibility rule.** The `peakAtail-prime` branch must retain a compatibility mode that
> reproduces v2 output **byte-for-byte** on the same inputs, for the lifetime of the branch. This is
> enforced by a regression test on a fixture and re-validated on the real PBMC chr19+21 slice before
> any prime accuracy number is quoted. If byte-identity is ever lost, the three-way comparison is
> void until it is restored.

## 2. Datasets, scoring protocol and denominators — all unchanged from v2

**Datasets (the same three arms the v2 gate used).**

| dataset | BAM | `--seq-len` |
|---|---|---|
| PBMC 10k v3 (human, CellRanger, CB/UB) | `data/benchmark/pbmc_10k_v3/pbmc_10k_v3_possorted_genome_bam.bam` | 91 |
| GSE104556 testis mouse 1 (STARsolo) | `data/benchmark/gse104556/starsolo/Mouse1_scRNAseq/Aligned.sortedByCoord.out.bam` | 98 |
| GSE104556 testis mouse 2 (STARsolo) | `data/benchmark/gse104556/starsolo/Mouse2_scRNAseq/Aligned.sortedByCoord.out.bam` | 98 |

A PBMC **chr19+chr21 slice** (`/mnt/ssd0/emaout/peakatail_benchmark/perf/slice/pbmc_chr19_21.bam`,
41,608,052 + 9,290,404 records, md5 `7b94a6684e32449bde928b59a4e0bea5`) is the development loop only.
**No headline claim may rest on the slice**; every number that enters the manuscript comes from a
full-BAM run on the three datasets above.

**Scorer.** `scripts/benchmark_tools/score_tool.py` with the reference arguments of
`scripts/benchmark_tools/stage2_final_launch.sh`, invoked through `scripts/prime/score_pas.sh`. Not
one line of the metric is re-implemented: the wrapper reads the values out of the scorer's own TSV.
Matching is strand-matched `bedtools closest -s -d -t first` on 1-bp points, `LC_ALL=C` throughout.

**Metrics and their denominators (frozen):**

| symbol | scorer row | denominator |
|---|---|---|
| P@10 / P@25 / P@100 | `precision / real / atlas_full / {10,25,100}` | number of scored calls `n` |
| R_det@100 | `recall / real / atlas_detected / 100` | human **285,136** (`results/benchmark_tools/shared_refs_pbmc/pas2.in_detected_genes.bed`); mouse **126,686** (`results/benchmark_tools/gse104556/shared_refs/pas2.in_detected_genes.bed`) |
| F1_det@100 | `f1 / real / atlas_detected / 100` | harmonic mean of the two above |
| null P@100 | `precision / null_genic / atlas_full / 100`, 3 shuffle seeds | same `n`, positions reshuffled inside merged detected gene bodies |
| Kinnex P@25 | fraction of calls within 25 bp, same strand, of an x3p long-read 3' end | `results/benchmark_tools/kinnex_truth/x3p/x3p_truth_t{5,20}.point.bed` |
| Kinnex decoy@25 | same against long-read internal-priming termini | `.../x3p_decoy.point.bed` |

"Precision" in this document always means **agreement with the curated atlas (PolyASite 2.0 rep sites
∪ protein-coding TES) at 100 bp**, which is what the manuscript calls precision. It is not a claim
about truth; the atlas-independent read-out is the Kinnex column, and manuscript
[16](16_trusted_novel_kinnex.md) states its limits (donor-mismatched, site-level).

**The v2 default operating point that prime must beat** (from `19_final_gate_v2.md`, arm
`pas_PRESPEC_precision_default.bed` = tier-1 ∩ IP-pass ∩ ≥2 molecules):

| dataset | n | P@10 | P@25 | P@100 | R_det@100 | F1_det@100 |
|---|---:|---:|---:|---:|---:|---:|
| PBMC 10k v3 | 46,524 | 0.520871 | 0.638875 | **0.706195** | **0.175443** | 0.281060 |
| testis mouse 1 | 26,255 | 0.579661 | 0.682042 | **0.745001** | **0.204790** | 0.321268 |
| testis mouse 2 | 26,526 | 0.594473 | 0.692980 | **0.757219** | **0.208002** | 0.326357 |

PBMC Kinnex read-out for the same arm: t5 P@25 0.764659, t20 P@25 0.558378, decoy@25 0.129997
(`results/algo_headroom/A3_scoring_model/tables/headline_arms.tsv`, which reproduces
`results/reliability/trusted_novel_final_v2_pbmc/validate_incl_genebodies/`).

## 3. What counts as success for prime — decided now

Prime's purpose is to **lift** the precision/recall curve, not slide along it. Filtering already
slides along it, so a change that trades precision for recall (or the reverse) is not a success no
matter how good the trade looks.

**Comparison discipline.** The headline comparison is **default vs default**: v2's default output as
defined in [13](13_reliability_positioning.md) §1 against whatever prime's default emits, each run
with its own defaults on the same BAM with the same `--seq-len`, `--threads` free to differ.

### 3.1 PRIMARY criterion (the one that decides adoption)

Let `ΔP = P@100(prime default) − P@100(v2 default)` and
`ΔR = R_det@100(prime default) − R_det@100(v2 default)`, per dataset.

> **Prime's default PASSES iff, on all three datasets (PBMC, mouse 1, mouse 2) simultaneously:**
> **(i) `ΔP ≥ −0.005`** (absolute; no material loss of atlas-agreement precision), **and**
> **(ii) `ΔR ≥ +0.010`** (absolute; a real gain in detected-gene recall — on PBMC that is
> 0.175443 → ≥ 0.185443, i.e. ≥ +5.7 % relative), **and**
> **(iii) `ΔF1_det@100 > 0`.**

The caller is deterministic — the same code on the same BAM produces byte-identical output — so
these are not noisy estimates and no tolerance for run-to-run variation is needed. The `−0.005` in
(i) is a *materiality* allowance, fixed now at roughly one third of the smallest precision movement
the project has treated as meaningful (the 1.1 pp v1→v2 PBMC precision change,
[19](19_final_gate_v2.md) §1), not a noise allowance.

Two strictness labels are reported alongside, so the PI can see which was achieved without the
criterion moving:
* **STRICT up-and-right** — `ΔP ≥ 0` and `ΔR ≥ +0.010` on all three.
* **RIGHT-AT-EQUAL-PRECISION** — `|ΔP| ≤ 0.005` and `ΔR ≥ +0.010` on all three.

**FAIL ⇒ the change does not become the default.** If the criterion is missed on any of the three
datasets, the new behaviour stays behind a non-default flag and the branch's default remains v2
behaviour. A miss is reported as a negative result in the manuscript record, in the same way
[16](16_trusted_novel_kinnex.md)'s trusted-novel negative result was.

**Supporting (non-deciding) views, always reported with the primary:**
* **matched precision** — prime's score threshold moved so that `|P@100(prime) − P@100(v2)| ≤ 0.002`;
  report `R_det@100`. This is the cleanest statement of "the curve moved", and is how
  `results/algo_headroom/` measured headroom.
* **matched call count** — prime's threshold moved so that `n(prime) = n(v2) ± 1 %`; report P and R.
* the full trade curve (n, P@100, R_det@100, F1) over the prime score, so a reader can see the whole
  curve rather than one point on it.

### 3.2 SECONDARY criteria (reported for every prime arm; none of them can rescue a primary FAIL)

1. **Resolution.** P@10 and P@25, and the retention ratios P@10/P@100 and P@25/P@100. Target:
   no loss vs v2 (PBMC v2 retention P@10/P@100 = 0.7376, P@25/P@100 = 0.9047).
2. **Atlas-independent concordance.** Kinnex x3p P@25 against t5 and t20, and the internal-priming
   **decoy@25** rate, PBMC only. Target: t20 P@25 not below v2's 0.558378 **and** decoy@25 not above
   v2's 0.129997. A recall gain bought in long-read internal-priming decoys is not a recall gain —
   `results/algo_headroom/A1_end_pileup/` measured exactly that failure mode for the end-position
   channel.
3. **Cross-mouse reproducibility.** The two testis mice are biological replicates on the same
   protocol. Report per-mouse P/R and the mouse1↔mouse2 call-set concordance (fraction of mouse-1
   calls with a mouse-2 call within 25 bp, same strand, and the reverse). Target: prime's
   cross-mouse concordance ≥ v2's, and `|ΔP|` between mice no larger than v2's 0.0122.
4. **Compute.** Wall time and peak RSS from `/usr/bin/time -v` on the full PBMC BAM and mouse 1,
   uncontended, compared with v2 (PBMC 27 m 43 s / 12.45 GB single-run; mouse 1 9 m 03 s / 3.68 GB —
   [19](19_final_gate_v2.md) §3 and the #97 PR). Guard rail: prime must not exceed **2× v2 wall time
   or 1.5× v2 peak RSS**. Exceeding the guard rail does not by itself fail the primary criterion but
   must be declared prominently and needs the PI's explicit acceptance.
5. **Test suite.** `1,277 passed, 9 skipped, 4 xfailed, 1 xpassed` with exactly the 2 known
   environment failures in `tests/test_pyproject_install.py`. Any new failure is a stop signal.
6. **Byte-identity of the compatibility mode**, re-validated on the slice at the end of the branch.

### 3.3 Threshold-selection rule — no threshold may be tuned on the metric it is reported against

Prime is expected to introduce at least one per-site score with a threshold that sets the default
operating point. That threshold is a free parameter and is therefore the single easiest way to
manufacture a result. Fixed now:

> **Rule T.** Any per-site score threshold shipped as a **default** is chosen on data **disjoint**
> from the data it is evaluated on, and never by reading the metric being reported.
>
> * **Primary protocol — cross-dataset.** The default threshold is fixed on **GSE104556 testis
>   mouse 1** and then evaluated, unchanged, on **PBMC 10k v3** and **testis mouse 2**. Mouse 1's own
>   post-hoc numbers are reported and labelled `THRESHOLD-FITTING SET — not evidence`.
> * **Secondary protocol — chromosome-held-out**, for anything that must be fixed within one dataset:
>   fit on odd-numbered autosomes, evaluate on even + X + Y (the fold definition already used by
>   `results/algo_headroom/A3_scoring_model/`, whose chromosome-grouped 2-fold and 5-fold AUCs agree
>   to 0.001). Both folds reported.
> * The threshold-selection objective is declared **before** it is run and is stated in the commit
>   that introduces the threshold.
> * **Kinnex long-read truth is never used to select anything** — not a threshold, not a feature, not
>   a model. It is a read-out only. The same holds for the PBMC atlas numbers in §2.
> * Any model coefficients are trained **offline** and shipped as constants evaluated with numpy;
>   scikit-learn does not enter the tool's runtime import path. The training script, its inputs and
>   its fold definition are committed with the coefficients.
> * A threshold chosen any other way ships **flag-off** and is labelled exploratory, permanently.

### 3.4 Status of every prime number with respect to the manuscript's gates

> **Every number produced on the `peakAtail-prime` branch is EXPLORATORY** with respect to the
> pre-registered gates in [13](13_reliability_positioning.md) (and [10](10_caller_fix_plan.md) §5),
> which continue to describe the **v2** default and are unchanged by this work.
>
> Prime cannot retroactively become "the pre-registered default". If prime passes §3.1 and the PI
> decides to adopt it, adoption requires (a) its **own** pre-registration document, written before
> the adoption run, restating the default output, its gates and its thresholds; (b) a fresh
> verification pass by the independent verifier on the adoption run; and (c) an explicit statement in
> the manuscript of what changed and which figures moved. Until all three exist, the manuscript's
> reported default is v2 and prime appears — if at all — as a clearly labelled exploratory arm.

The pre-registered gates in 13 §1 (P@100 ≥ 0.50 on PBMC and both mice) are *not* prime's success
criterion; §3.1 is strictly harder than them, and passing 13 §1 while failing §3.1 is a FAIL.

## 4. Flag policy

1. **Every behavioural change is behind a CLI flag / config option.** No behavioural change is
   unconditional.
2. **New behaviour may be the DEFAULT on this branch** — that is the point of "prime" — but the v2
   behaviour must remain reachable through explicit flags.
3. **A single compatibility switch restores v2 behaviour in one step**, so a reviewer never has to
   reconstruct a flag list. Reaching v2 output must not require knowing which changes exist.
4. **Byte-identity is a test, not a claim.** A regression test asserts it on a committed fixture and
   runs in the normal suite; `scripts/prime/identity_check.py` re-validates it on the real slice.
5. Each flag's default, and the evidence for that default (which measurement in
   `results/algo_headroom/` or which measurement made on this branch), is recorded in
   `tools/pa-prime/PRIME_PLAN.md` and in the commit that introduces it. Where a measurement says an
   idea does not pay, the flag ships **defaulted OFF** and the plan says so.

## 5. Harness (fixed before any result)

| script | what it does |
|---|---|
| `scripts/prime/common.sh` | paths, the reference-argument sets copied from the launcher, `assert_code` (refuses to run unless `ema` resolves inside the intended worktree — five worktrees of this package exist), `guard_load` |
| `scripts/prime/make_slice.sh` | verifies/creates the PBMC chr19+21 slice BAM |
| `scripts/prime/run_slice.sh` | runs the caller on the slice with arbitrary extra flags into `/mnt/ssd0/emaout/peakatail_benchmark/prime/<label>/` |
| `scripts/prime/score_pas.sh` | scores any pas.bed with `score_tool.py` + the launcher's reference arguments; prints n / P@100 / P@25 / P@10 / R_det@100 / F1 and the Kinnex 25 bp concordance |
| `scripts/prime/identity_check.py` | byte-identity of two run trees; the compatibility-mode check |

Operational constraints in force: ≤ 8 threads for this work; the Stage-3 replication chain under
`/mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v3/` is never touched; nothing is written to
`/mnt/ssd2`, to another worktree under `tools/`, or to `results/benchmark_tools/*final*`.

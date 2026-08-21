# Issue: cohort mode — promote 1-molecule sites that a second sample independently calls (PROVISIONAL: verifier first)

**Labels:** enhancement, accuracy, cohort, research, **provisional**. Measured 2026-08-21 in
`results/algo_headroom/A5_cross_sample/` on the frozen v2 outputs of GSE104556 mouse 1 and mouse 2
(code `9dfdefb`). **This measurement has NOT been through the independent verifier.** Everything
below is provisional and must not enter the manuscript before it has.

Related: **#94**; `manuscript/20_stage3_replication.md` (the cohort PAS space and patient-level
replication filter this would live next to); `manuscript/23_algorithm_roadmap.md` §2 row G and §3 Step 4.

## 1. Symptom

The caller's largest single discard is the 1-molecule tier-1 site: PBMC 72.1 % of tier-1 sites, mouse
1 26,544 of 52,807 (50.3 %). In a **single** library there is nothing to distinguish a real site with
one molecule from a spurious one. In a **cohort** — which is how the tool is actually used for
switch analysis — there is: another sample may have called the same position with full support. The
caller does not use that, even in cohort mode.

## 2. Measured diagnosis (PROVISIONAL)

Design (no atlas used to select anything): pool = mouse-1 tier-1 sites with exactly 1 clip molecule
(IP-pass by construction, both mice ran `--ip-filter`); signal = a mouse-2 **default** call within
W bp on the same strand; evaluation = the scorer's own conventions against PolyASite 2.0 GRCm38 and
the 126,686-site detected-gene denominator. The pipeline reproduces the published mouse arms exactly
(mouse 1 P@100 0.7449 / R_det 0.2048 vs `19_final_gate_v2.md`'s 0.745001 / 0.204790).

**The added calls are 1.6× better than an expression-matched control:**

| arm | n | P@10 | P@100 | R_det@100 |
|---|---|---|---|---|
| whole 1-molecule pool | 26,544 | 0.2430 | 0.3941 | 0.1013 |
| **corroborated** (m2 default >= 2 mol within 10 bp) | 2,816 | **0.5163** | **0.6683** | 0.0206 |
| corroborated, stricter (m2 >= 3 mol, 10 bp) | 1,365 | 0.5678 | **0.7143** | 0.0108 |
| corroborated, strictest (m2 >= 5 mol, 5 bp) | 275 | 0.6073 | **0.7491** | 0.0023 |
| CONTROL random singletons | 2,816 | 0.2447 | 0.3967 | 0.0114 |
| CONTROL top-2,816 by `window_reads` | 2,816 | 0.1491 | 0.4798 | 0.0163 |
| **CONTROL expression-matched** (`window_reads` decile) | 2,816 | 0.2383 | **0.4144** | 0.0128 |

**At matched call count the corroborated union beats every control on both axes:**

| arm | n | P@100 | R_det@100 |
|---|---|---|---|
| shipped default (mouse 1) | 26,263 | 0.7449 | 0.2048 |
| **default + corroborated (m2 >= 2 mol, 10 bp)** | 29,079 | **0.7374** | **0.2207** |
| default + expression-matched singletons | 29,079 | 0.7129 | 0.2143 |
| default + random singletons | 29,079 | 0.7111 | 0.2134 |
| default + top-`window_reads` singletons | 29,079 | 0.7192 | 0.2114 |

It **replicates in the reverse direction** (mouse-2 singletons corroborated by mouse-1: added-call
P@100 0.6310 at m1 >= 3 mol / 10 bp against a 0.4106 pool; union n 27,747, P@100 0.7516,
R_det 0.2151 against the m2 default's 0.7572 / 0.2080).

**Three honest limits, stated up front.**
1. Against the shipped default this is **+7.8 % relative recall for −0.0075 precision** — it
   **fails** the pre-registered materiality allowance of `24_prime_preregistration.md` §3.1(i)
   (ΔP >= −0.005). At the precision-safe setting (m2 >= 3 mol) ΔR is only +0.0085, below §3.1(ii)'s
   +0.010. It sits exactly on the boundary and does not clear it on one mouse.
2. **Internal priming replicates across samples too** — the genome is the same in both animals — so
   corroboration cannot suppress the residual-IP class the way sequence information or an independent
   truth can. This must be checked against Kinnex decoys before anyone believes the precision figure.
   (Mouse has no Kinnex truth, so this check needs a multi-library human dataset.)
3. It **cannot help a single-library dataset**, i.e. it cannot move the PBMC headline at all.

## 3. Proposed change

- [ ] **Send `results/algo_headroom/A5_cross_sample/` to the independent verifier first.** Nothing
      below starts until it comes back.
- [ ] Ask the verifier specifically for: a decoy-class check on a multi-library human cohort; whether
      corroboration survives after the calibrated score of **issue14** is applied (they may be
      measuring the same thing — the score also ranks singletons, and the two recovered-site sets have
      not been intersected).
- [ ] If it survives: implement in **cohort mode only**, as a promotion rule on the unified PAS space
      that `20_stage3_replication.md` already builds — a 1-molecule site is promoted when >= N other
      samples carry a default-quality call within W bp on the same strand. Defaults to propose:
      W = 10 bp, corroborator >= 3 molecules, N = 1 (the precision-safe setting).
- [ ] Never in single-sample mode; never as a way to reach a precision target in a benchmark table.
- [ ] Report the promoted set as a separate column/flag, not merged silently into the default arm.

## 4. How to reproduce

    export LC_ALL=C
    cd /mnt/ssd1/Projects/PeakATail_wd/results/algo_headroom/A5_cross_sample
    cat README.md                      # design, controls, fidelity check
    column -t -s$'\t' A5_cross_sample_summary.tsv
    ./score.sh union_k2_w10.bed UNION  # any arm can be re-scored in ~10 s

## 5. What it is expected to buy

On a two-sample mouse cohort, **+7.8 % relative recall for −0.75 precision points**, or +4.2 % recall
at −0.16 precision points at the precision-safe setting — measured, controlled for expression, and
replicated in both directions, but on the boundary of the pre-registered criterion and on **one
species with two samples**. It buys nothing on PBMC. Its real interest is that it is an **atlas-free
selection signal on a dataset that is not PBMC**, which is exactly the axis on which every other
positive result in `23_algorithm_roadmap.md` is untested.

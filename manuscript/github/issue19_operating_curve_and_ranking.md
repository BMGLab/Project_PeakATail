# Issue: PeakATail has an operating *curve* and we ship one point — document it, name the ranking column, and stop letting `>=2 molecules` be read as the tool's best

**Labels:** documentation, enhancement, benchmark. Measured by the competitive-gap study
(`results/perf_gap/D2_rule_sweep/`, `D3_head_to_head/`) and confirmed by its independent verifier
(`results/perf_gap/VERIFY/`, verdict FIXED, `out/VERIFICATION_RESULTS.tsv`), 2026-08-21, against the
frozen v2 tree `9dfdefb`.

Related: **#94** (Stage 1b/1c poly(A) evidence — this is the arm structure that work produced);
`manuscript/22_performance_roadmap.md` §4.1; sibling drafts **issue13** (`--ip-filter` should be
default-on) and **issue14** (a calibrated per-site score, which would replace the integer knob
described here). **Zero caller behaviour changes in this issue.**

## 1. Symptom

The caller emits `pas.bed` with `col5 = distinct clip molecules`, and the shipped default keeps
`col5 >= 2`. That threshold is an *operating point on a curve*, but nothing in the tool, the README
or the manuscript says so, and one consequence is a claim about our own tool that is measurably
false:

> "the precision default is the most atlas-concordant de novo call set in the panel on both datasets,
> **at polyApipe-or-lower recall** (PBMC −12 %, mouse −18 % / −17 %)"
> — `manuscript/01_outline_and_journals.md`

That sentence compares **46,524 of our calls with 120,916 of polyApipe's**. It is true of the
`>=2`-molecule default and false of the tool.

Second symptom: `col5` is a perfectly good ranking column and we never say so, so no user and no
reviewer can construct the curve themselves.

## 2. Measured diagnosis

**(a) On both mice, an arm we already publish wins outright — no ranking, no new code.**

| dataset | arm | n | P@100 | R_det@100 | F1_det |
|---|---|---:|---:|---:|---:|
| mouse 1 | PeakATail tier-1, IP-pass, **>=1 molecule** | 52,792 | **0.5686** | **0.2806** | **0.3758** |
| mouse 1 | polyApipe (full set) | 89,527 | 0.4005 | 0.2499 | 0.3078 |
| mouse 2 | PeakATail tier-1, IP-pass, **>=1 molecule** | 51,842 | **0.5880** | **0.2826** | **0.3818** |
| mouse 2 | polyApipe (full set) | 87,983 | 0.4119 | 0.2508 | 0.3118 |

Higher on precision, recall and F1 simultaneously, on both animals. (Verifier: reproduced to six
decimals from the competitors' own frozen call sets.)

**(b) On PBMC, ranking by `col5` and cutting at polyApipe's own call count wins both axes.**

| N | PeakATail P@100 / R_det | polyApipe P@100 / R_det |
|---:|---|---|
| 46,524 (ours) | **0.7062 / 0.1754** | 0.6049 / 0.1375 |
| 106,170 (Sierra's) | **0.4269 / 0.2204** | Sierra 0.2559 / 0.1350 |
| **120,916 (polyApipe's own N)** | **0.4036 / 0.2336** | 0.3800 / 0.1988 |
| 333,920 (our max) | 0.2044 / **0.3043** | scAPAtrap at its full 787,138: 0.0999 / 0.3004 |

**(c) The measured operating curve** (default + a random x % of single-molecule tier-1 calls — the
*no-information* control, so this is a lower bound on what a ranked cut achieves):

| x | n | P@100 | R_det | F1_det | Kinnex x3p <=25 bp |
|---:|---:|---:|---:|---:|---:|
| 0 (default) | 46,524 | 0.7062 | 0.1754 | 0.2811 | 0.7647 |
| 30 % | 82,834 | 0.4908 | 0.2044 | 0.2886 | 0.5592 |
| 60 % | 119,144 | 0.4080 | 0.2332 | 0.2968 | 0.4785 |
| 80 % | 143,356 | 0.3752 | 0.2510 | 0.3008 | 0.4483 |
| 100 % (>=1 mol) | 167,565 | 0.3520 | 0.2685 | 0.3046 | 0.4256 |

x = 30 is the first point that beats polyApipe on both axes; at x = 60 we are at polyApipe's own call
count and dominate it on precision, recall, F1 **and** long-read concordance (0.4785 against
polyApipe's **0.3895** — recomputed independently while writing this draft, see
`results/perf_gap/CHECK22/README.md` §2, which also reproduces the default's 0.7647).

**(d) The support sweep, for the documentation table** (PBMC, n / P@25 / P@100 / R_det / F1 / Kinnex):
`>=1` 167,565 / 0.2880 / 0.3520 / 0.2685 / 0.3046 / 0.4256 — `>=2` 46,524 / 0.6389 / 0.7062 / 0.1754 /
0.2811 / 0.7647 — `>=3` 29,814 / 0.7738 / 0.8358 / 0.1419 / 0.2427 / 0.8822 — `>=5` 20,320 / 0.8489 /
0.9074 / 0.1122 / 0.1997 / 0.9432 — `>=10` 13,560 / 0.8882 / 0.9412 / 0.0828 / 0.1522 / 0.9673.

**(e) The default is not the F1 optimum on either dataset.** PBMC `>=1` F1 0.3046 vs the default's
0.2811; mouse 1 0.3758 vs 0.3213. The default is defensible on reliability-first grounds and on
nothing else.

**Two honesty constraints that must travel with these numbers.**
1. The **precision** half of the matched-N claim is **not robust to the competitor's tie-break** in
   the band 20,320 <= N <= 35,759: polyApipe's `peakdepth` is more heavily tied than our molecule
   counts (18,857 sites at depth 2, 6,161 at depth 3), and under an oracle tie-break polyApipe reaches
   0.9193 / 0.8611 / 0.7620 at N = 20,320 / 29,814 / 35,759 against our 0.9074 / 0.8358 / 0.7579 (our
   own optimistic bracket at 35,759 is 0.7917 — the intervals overlap). It **is** unambiguous at
   N = 46,524 (polyApipe best case 0.6636 vs 0.7062) and at N = 120,916, where polyApipe is its own
   full set with no ties at all. **The recall half is robust at every N.**
2. Ranking by **window read depth** instead of clip molecules gives P 0.3900 at N = 46,524 against
   0.7062 — depth is a much worse ranking and is mildly *anti*-correlated with correctness at fixed
   molecule support. `col5` is the right column; say which one it is.

## 3. Proposed change

**Tool / docs (no behaviour change):**
- [ ] Document `pas.bed` **col5 as the ranking column** ("distinct clip molecules supporting the
      cleavage site; monotone in reliability") in the README, the CLI help and the output schema doc.
- [ ] Ship the support sweep of §2(d) as a table in the README, with the explicit sentence that
      **`>=2` is a reliability-first default, not the F1 optimum**, and that `>=1` is the sensitivity
      arm.
- [ ] Add a convenience so a user gets the curve from one run instead of five: either
      `--emit-support-arms 1,2,3,5,10` writing one BED per cut, or document that
      `awk '$5>=k' pas.bed` is the supported way to move the operating point. **Do not** change the
      default.
- [ ] Cross-check with **issue14**: if the calibrated per-site score ships, the curve becomes a score
      threshold rather than an integer knob and this documentation must be rewritten, not duplicated.

**Manuscript (owned by `21_paper_architecture.md` / `05_figure_index.md`, listed here so the
programme is visible in one place):**
- [ ] Replace the single-row head-to-head table with a **matched-N figure**: P@100 and R_det versus N,
      one curve per tool, our curve produced by ranking on col5.
- [ ] Restate the claim as **"best precision and best recall among de novo tools at polyApipe's own
      call budget, on both datasets"**, carrying the tie-break caveat of §2 constraint 1.
- [ ] Correct the "recall at or below polyApipe" line in `01_outline_and_journals.md`.
- [ ] Say in C6 (`21_paper_architecture.md`) that the default is **not** the F1 optimum and why we
      keep it anyway.
- [ ] Methods disclosure required by the matched-N analysis: **SCAPTURE emits no non-circular
      per-peak ranking column** (its only numeric per-peak field is atlas-derived and its own
      quantification table is empty), so no matched-N statement below its own N is possible for it;
      **Sierra's BED score column is all zeros**, so its ranking was reconstructed from its count
      matrix.
- [ ] Also correct the standing summary line "both tiers 333,920 / 0.1932 / 0.2979 / 0.2415": the
      measured v2 values are **IP arm both tiers 333,920 / 0.2044 / 0.3043 / 0.2445** and **no-IP arm
      both tiers 402,765 / 0.1932 / 0.3363 / 0.2455**; 0.2979 / 0.2415 are the **pre-v2** row of
      `manuscript/15_final_gate.md` line 36 (n 328,302).

## 4. How to reproduce

    export LC_ALL=C
    cd /mnt/ssd1/Projects/PeakATail_wd
    column -t -s$'\t' results/perf_gap/D3_head_to_head/tsv/01d_matched_n_headline_pbmc.tsv
    column -t -s$'\t' results/perf_gap/D3_head_to_head/tsv/01e_matched_n_headline_mouse1.tsv
    column -t -s$'\t' results/perf_gap/D3_head_to_head/tsv/01b_tiebreak_sensitivity.tsv
    awk -F'\t' 'NR==1 || $1 ~ /support sweep|operating curve/' \
        results/perf_gap/D2_rule_sweep/D2_arm_sweep.tsv | column -t -s$'\t'
    sed -n '1,60p' results/perf_gap/VERIFY/README.md          # correction 3: the tie-break qualification

Any arm can be re-scored in seconds with `results/perf_gap/D2_rule_sweep/03_score.sh`, which passes
the launcher's reference arguments verbatim (verifier-confirmed).

## 5. What it is expected to buy

**Nothing in the tool and everything in the paper.** No published number moves: the pre-registered
default stays 0.7062 / 0.1754 (PBMC) and 0.7450 / 0.2048 / 0.7572 / 0.2080 (mice). What changes is
that the head-to-head stops being a single row that reads as "high precision, weak recall" and becomes
a curve that shows PeakATail above polyApipe over a wide range on PBMC and outright dominant on both
mice. This is the cheapest item on the performance roadmap and the only one with **zero**
pre-registration exposure — the `>=1`-molecule arm is already the pre-registered sensitivity arm
(`13_reliability_positioning.md` §1), and the ranked-prefix curve is a descriptive analysis that must
simply be labelled post-hoc.

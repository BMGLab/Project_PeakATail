# Issue: `--ip-filter` is the biggest measured accuracy lift in the caller and is OFF by default; and the poly(A) clip-rate constant we ship is 2× wrong

**Labels:** enhancement, accuracy, docs. Measured by the algorithmic-headroom study
(`results/algo_headroom/`) and its independent verification (`results/algo_headroom/VERIFY/`,
verdict FIXED), 2026-08-21, against the frozen v2 tree `9dfdefb`.

Related: **#95** (IP rule — strand window, `--ip-rule kinnex` addendum), **#99 §2** (clip-rate warning
samples only the head of chr1). This issue supplies the *constant* that #99 §2 is missing, and the
*evidence* that #95's rule is worth keeping at all.

## 1. Symptom

`--ip-filter` is an opt-in flag (`ema/cli/config_schema.py:529`, `is_flag=True`, default off). Every
benchmark arm we report runs with it (`scripts/benchmark_tools/stage2_final_launch.sh:72-77`), but a
user who runs `ema` out of the box does not get it, and no document states what it is worth.

Separately, `ema/countmatrix/polya.py` (docstring) and `manuscript/10_caller_fix_plan.md` state that
**1.15 %** of cell-barcoded reads carry a poly(A) soft clip. That number is not reproducible.

## 2. Measured diagnosis

**(a) The IP filter is the largest measured lift in the whole caller, and it is a *lift*, not a
trade.** Sweeping `tier1 & >=k clip molecules` with and without the veto and interpolating recall at
**matched precision** on each truth (`VERIFY/tables/v8_ipveto_value.tsv`, PBMC 10k v3):

| matched precision | R_det@100 with veto | without | relative gain |
|---|---|---|---|
| atlas P@100 = 0.60 | 0.2033 | 0.1887 | **+7.7 %** |
| atlas P@100 = 0.7062 (the shipped default) | 0.1754 | 0.1612 | **+8.8 %** |
| atlas P@100 = 0.80 | 0.1512 | 0.1366 | **+10.7 %** |
| atlas P@100 = 0.8358 | 0.1419 | 0.1258 | **+12.8 %** |
| long-read Kinnex P@25(t5) = 0.7647 | 0.2722 | 0.2313 | **+17.7 %** |
| long-read Kinnex P@25(t5) = 0.85 | 0.2422 | 0.1971 | **+22.9 %** |

It lifts the curve rather than sliding along it because it brings in information the molecule
thresholds do not have — genomic sequence. For comparison, the entire set of detector changes we
tested (relaxing `--polya-min-clip`, purity, boundary run, seeded extension, an end-position/pileup
channel) is worth **−10 % to +0.4 %** on the same axis.

**(b) The continuous A-content covariate the filter already computes is a stronger discriminator
than the binary flag, and stronger than a 49-feature model.** Separating candidates near a genuine
long-read terminus from candidates near a long-read *internal-priming* terminus
(`VERIFY/tables/v3_truth_vs_decoy_auc.tsv`): `ip_tool_afrac` inverted **AUC 0.7790**; a full
gradient-boosted model 0.7655; raw clip molecule count 0.6785; the **binary** IP flag alone 0.534.
The float is already computed and thrown away.

**(c) The clip-rate constant.** Genome-wide, with the tool's own `read_check` + `clip_site`
criteria: **0.5730 %** (3,195,067 / 557,564,408 accepted CB reads); 0.5364 % on
`check_clip_rate`'s own denominator; 0.5030 % of all valid-CB mapped reads. Per chromosome
0.3669 % (chr21) to 0.8179 % (chr19). The verifier reimplemented the criteria from scratch in gawk
over `samtools view` and got **7,567,263 accepted CB reads / 27,768 qualifying clips / 0.003669 on
chr21 — identical to six decimals** from a separate implementation. `check_clip_rate()`'s
head-of-BAM sample of 200,000 CB reads returns **2.2565 %**, i.e. the QC estimator is ~4× high on a
coordinate-sorted BAM (this is #99 §2's mechanism, with the number attached).

## 3. Proposed change

- [ ] Make `--ip-filter` **default-on** whenever `--genome-fasta` is supplied; add `--no-ip-filter`.
      Warn (do not fail) when no FASTA is given, naming the measured cost.
- [ ] Emit the **continuous** downstream A-fraction per site (`ip_afrac`) into `pas_support.tsv` and
      `annotatedpas.bed` in both `--ip-filter-mode` settings — it is already computed in
      `ema/experimental/internal_priming.py`, and #95's bullet 3 asks for the per-site flag anyway.
      Do this on the strand-corrected window from **#95 §1**.
- [ ] Correct the clip-rate figure in `ema/countmatrix/polya.py`'s docstring and in
      `manuscript/10_caller_fix_plan.md`: **~0.57 % genome-wide** (0.37–0.82 % per chromosome), not
      1.15 %.
- [ ] Fix `check_clip_rate()` to sample across the BAM (every Nth record, or per-chromosome strata),
      or compute the rate during the main pass and warn at the end — **this is #99 §2**; close them
      together.
- [ ] Regression test: the QC estimator on a fixture whose head is unrepresentative must not differ
      from the full-file rate by more than a stated factor.

## 4. How to reproduce

    export LC_ALL=C
    cd /mnt/ssd1/Projects/PeakATail_wd
    cat results/algo_headroom/VERIFY/tables/v8_ipveto_value.tsv        # the IP sweep
    cat results/algo_headroom/VERIFY/tables/v3_truth_vs_decoy_auc.tsv  # afrac vs the model
    cat results/algo_headroom/A2_clip_sensitivity/tables/T0*  T3d*     # the clip-rate census
    # verifier's from-scratch chr21 check: VERIFY/VERIFICATION_REPORT.md §5

## 5. What it is expected to buy

Nothing on our published numbers — the pre-registered arm already runs with `--ip-filter`, so
P@100 0.7062 / 0.7450 / 0.7572 and R_det 0.1754 / 0.2048 / 0.2080 (`manuscript/19_final_gate_v2.md`)
are unchanged. It buys **users** the 8–23 % relative recall at matched precision that we currently
keep behind a flag, it removes a stated fact that is 2× wrong, and it retires a QC estimator that is
4× off — which matters because that estimator could equally *mask* a genuinely destroyed evidence
channel (it already caused a wrong pre-registered data-quality exclusion, #99 §2).

**Caveat to state in the changelog:** the like-for-like *no*-IP mouse arm at v2 code does not exist
(only `results/benchmark_tools/gse104556/peakatail_clipseeded_v3/mouse1`, Stage-1c code), so the
sweep above is PBMC-only. A mouse re-run without `--ip-filter` costs ~9 min / 3.7 GB and should be
done before the number goes into the paper.

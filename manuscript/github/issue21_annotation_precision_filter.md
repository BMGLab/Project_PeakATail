# Issue: add a `--strict` post-filter that drops non-3'UTR calls with no canonical hexamer — and do NOT ship the PBMC-tuned intronic variant, which does not transfer to mouse

**Labels:** enhancement, accuracy, needs-pre-registration. Measured by
`results/perf_gap/D3_head_to_head/` (17-candidate filter search, PBMC and mouse 1) and confirmed by
the independent verifier (`results/perf_gap/VERIFY/`, verdict FIXED; the whole filter table plus a
200-random-site spot check of the underlying per-call annotations, **200/200 sites, 0 mismatches on
every one of 15 columns**), 2026-08-21, frozen tree `9dfdefb`.

Related: **#95** (the IP veto — this issue shows internal priming is *not* the residual error mode of
the default, which is the argument for keeping #95's rule as-is), **#99** (same GTF machinery).
Sibling drafts: **issue20** (the mirror-image promotion rule on the same annotation axis — read both
before building either), **issue19** (the operating curve this sits on).
Roadmap context: `manuscript/22_performance_roadmap.md` §4.3.

## 1. Symptom

The precision default's residual error has never been characterised, so we cannot tell a user what
kind of mistake the tool makes, and we ship no way to trade a little recall for more precision other
than raising the molecule threshold (which is coarse: `>=2` -> `>=3` costs 36 % of the call set).

## 2. Measured diagnosis

**(a) The residual error is thin and nameable.** Using "hard FP" = no strand-matched atlas site within
100 bp **and** no Kinnex x3p >=5 UMI 3' end within 25 bp, the PBMC default has **7,027 / 46,524 =
15.10 %** hard FPs, concentrated as follows:

| stratum | FP rate within stratum | share of all FP mass |
|---|---:|---:|
| intronic | **0.3489** | 77.64 % |
| 3'UTR | 0.0179 | — |
| support = 2 molecules | **0.2980** | 70.86 % |
| support >= 10 molecules | 0.0165 | — |
| no canonical hexamer | — | 46.83 % |

**Internal priming is not the residual error mode and should stop being discussed as if it were.**
Two stricter re-tests flag **36** and **31** of the 46,524 default calls (the +10..+30 rule and the
Kinnex +1..+18 rule respectively); the tool's own rule flags 0 by construction. The residual is
low-support intronic peaks with no poly(A) signal.

**(b) polyApipe's error is bulk, ours is thin — and that is the reliability-first claim, measured.**
polyApipe: **59,856 / 120,916 = 49.50 %** hard FPs (8.5x our count, 3.3x our rate); 64.5 % of its
whole output is single-read peaks with a 62.6 % FP rate there; 8,083 of its calls (6.7 %) are flagged
by *our* IP rule and 70.7 % of those are false. Mouse 1 (atlas-only FP definition, **not** comparable
to the human number): PeakATail 6,695 / 26,255 = 25.5 %; polyApipe 53,670 / 89,527 = 59.95 %.

**(c) The transferable filter.** `drop calls that are not in a 3'UTR and carry no canonical hexamer`:

| dataset | n kept | P@100 | R_det@100 | F1_det | hard FPs removed | evidenced calls removed |
|---|---:|---:|---:|---:|---:|---:|
| PBMC default (baseline) | 46,524 | 0.7062 | 0.1754 | 0.2811 | — | — |
| PBMC, filter applied | **38,824** | **0.7993** | 0.1683 | 0.2781 | **44.3 %** | 11.6 % |
| mouse 1 default (baseline) | 26,255 | 0.7450 | 0.2048 | 0.3213 | — | — |
| mouse 1, filter applied | **22,623** | **0.8033** | 0.1942 | 0.3127 | **33.5 %** | 7.1 % |

**Be honest about the shape of this: F1 falls on both datasets.** It is a precision-for-recall trade,
which is exactly what a `--strict` mode is for, and exactly why it must not become the default.

**(d) The PBMC-tuned variant must not ship.** `drop intronic AND no hexamer` is the best single rule
on PBMC — it removes 2,495 hard FPs (35.5 %) for 2,988 evidenced calls (7.6 %), P 0.7062 -> **0.7801**,
and it is one of only **two of 17** candidates that raise precision *and* F1 at once (F1 0.2811 ->
0.2825). It removes **8.45 %** of mouse FPs. The mouse-testis FP mass sits in **"other exon"**
(FP rate 0.5395), not introns. A single tuned filter would not survive replication; the wider
non-3'UTR form would.

**(e) For contrast, polyApipe's analogous fix is strictly dominated by our default.** `peakdepth >= 3`
removes 96.3 % of its FPs and takes P 0.380 -> 0.828, but leaves 24,043 calls at R_det 0.1077 —
against our default's 0.7062 / 0.1754.

## 3. Proposed change

- [ ] Pre-register the rule and the adoption criterion under `24_prime_preregistration.md` §3 **before**
      re-running, and fix the variant in writing now: the proposal is the **non-3'UTR AND no-hexamer**
      form. The intronic variant is explicitly excluded, on the mouse evidence above.
- [ ] Implement as `--strict` (or `--post-filter annotation`), **default OFF**. With the flag off,
      output must be **byte-identical** to v2 (`24` §1 reproducibility rule).
- [ ] Reuse the GTF already loaded for gene assignment (**#99**'s path) and the canonical 12-hexamer
      search in `scripts/reliability/trusted_novel_pas.py` (-40..-5, correct minus-strand handling).
      **Do not** add a second annotation source or a second sequence scan — **issue20** wants the
      identical two predicates and they must share one implementation.
- [ ] Emit the two predicates as **columns** (`feature_class`, `hexamer`) in `pas_support.tsv`
      regardless of whether the filter is on, so a user can build their own rule and so downstream
      analyses can stratify. This is the cheap half of the issue and can land first.
- [ ] Document the measured trade in the CLI help: "raises P@100 to ~0.80 on both benchmark datasets
      at a cost of ~4-7 % of recall and ~0.003-0.009 F1".
- [ ] **Measure the interaction with issue20 before either ships.** They act on the same annotation
      axis in opposite directions (promote 3'UTR singletons; drop non-3'UTR hexamer-fail doubletons)
      and have **never been scored together**. Their combined effect is unmeasured and **must not be
      estimated by addition**.
- [ ] Regression tests: an intronic call with no hexamer is dropped under `--strict` and kept without
      it; a 3'UTR call with no hexamer is kept in both; the hexamer search is strand-correct.

## 4. How to reproduce

    export LC_ALL=C
    cd /mnt/ssd1/Projects/PeakATail_wd
    column -t -s$'\t' results/perf_gap/D3_head_to_head/tsv/07_filter_candidates_pbmc.tsv
    column -t -s$'\t' results/perf_gap/D3_head_to_head/tsv/07b_filter_candidates_mouse1.tsv
    column -t -s$'\t' results/perf_gap/D3_head_to_head/tsv/05_fp_forensics_pa_default_pbmc.tsv
    column -t -s$'\t' results/perf_gap/D3_head_to_head/tsv/06_fp_forensics_polyapipe_pbmc.tsv
    head -1 results/perf_gap/D3_head_to_head/tsv/03_calls_pa_default_pbmc.tsv   # per-call master table

Filter evaluator: `results/perf_gap/D3_head_to_head/d3_filter_eval.py`; annotation:
`d3_feature_class{,_mouse}.sh`, `d3_annotate_sites.py`.

## 5. What it is expected to buy

**For users:** a documented, cross-species-validated way to reach P@100 ~0.80 without the blunt
instrument of `>=3` molecules, which costs 36 % of the call set for a similar precision gain
(0.8358 at R_det 0.1419, against 0.7993 at R_det 0.1683 here — the filter is the better trade at that
precision).

**For the paper, and this is the bigger prize:** a *qualitative* error-mode claim that reliability-first
positioning needs and that we have never been able to make. PeakATail makes a thin, characterisable
error — 35.5 % of its FPs live in a single cell of the cross-tab (intronic AND no hexamer) — while
polyApipe makes a bulk error (64.5 % of its output is single-read peaks, 62.6 % of them false).
Both statements come from identical per-call machinery applied to both tools, and both survived the
verifier's 200-site column-by-column spot check.

**No published number changes** while the flag is off. **Tissue caveat that must go into the paper:**
the dominant FP class is tissue-dependent (intronic in PBMC, other-exon in mouse testis), so any
recommended post-filter must be the wider non-3'UTR form, not the PBMC-tuned intronic one.

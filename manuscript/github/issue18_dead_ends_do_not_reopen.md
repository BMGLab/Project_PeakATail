# Issue (documentation): four measured algorithmic dead ends — do not re-propose without new evidence

**Labels:** documentation, research, wontfix. Record of measurements made 2026-08-21
(`results/algo_headroom/A1_end_pileup`, `A2_clip_sensitivity`, `A3_scoring_model`, `A4_resolution`)
and independently verified (`results/algo_headroom/VERIFY/VERIFICATION_REPORT.md`, verdict FIXED).

This issue exists so that the next person who has one of these ideas — and they are all reasonable
ideas — finds the measurement instead of spending a week rediscovering it. Related: **#94**, **#95**,
**#98**, **#99**, and `manuscript/23_algorithm_roadmap.md` §2. **None of these is a bug. Each is a
plausible idea that was measured and does not pay.**

## 1. An end-position / read-3'-end pileup caller — DEAD

*Idea:* call sites from where reads terminate rather than from coverage shape.

*Measured:* at the shipped default's precision on chr19+21, the best pileup-only rule reaches
P@100 0.793 at **R_det 0.0638** against the default's 0.7385 / 0.2085 — a **2.6–3.3× loss of recall**,
held out per chromosome (chr19 3.1×, chr21 4.3×) and reproduced under the long-read truth (Kinnex
R@100 0.494 vs 0.173 on chr19). As an *added* channel on top of the default, the extra calls cap at
**36.1 % precision** at every bar tried against the default's 73.9 %, so union F1 stays flat while
precision falls.

*And the verifier made it worse:* A1's null was drawn uniformly inside gene bodies (median local
depth **2** molecule-ends per ±100 bp against **8** at missed atlas sites and **243.5** at recovered
ones). Standardised to the missed-atlas depth distribution, every pileup lift collapses:

| bar | lift over RAW null | lift over **depth-matched** null |
|---|---|---|
| mol5 >= 2 / 5 / 10 / 20 | 1.99× / 2.38× / 3.02× / 3.82× | **0.99× / 0.95× / 1.00× / 1.01×** |
| sharp + multi-cell + edge | 5.39× / 7.86× | **1.41× / 1.28×** |
| **clip >= 1 / >= 2** | 25.8× / 15.4× | **11.36× / 4.69×** |

**A molecule-end pileup at a true missed PAS is exactly as likely as at a random position of the same
local read depth.** The mechanism is in the code: 96.1–96.3 % of accepted reads carry no 3'-side soft
clip, so their 3' end is `start + 91` — a shifted copy of coverage (see **issue15**).

*Reproduce:* `VERIFY/tables/v4_a1_key_number_depthmatched.tsv`, `VERIFY/code/v4_a1_null.py`.

## 2. Relaxing `--polya-min-clip` (6 → 3/2/1) — DEAD

*Idea:* a more sensitive clip detector finds more sites.

*Measured:* k = 6 → 3 buys **+28.9 % more clip reads** and costs **−0.7 % to −10.0 % recall at
matched atlas precision** (−3.3 % to +2.6 % on Kinnex). k = 2 and k = 1 are −10 % to −45 %. Purity
and boundary-run relaxations move nothing. Site-level: clusters 16,338 (k6) → 22,542 (k3) → 66,436
(k1) while P@100 falls 0.5526 → 0.4706 → 0.2590 in lockstep.

*Why:* the added evidence is mostly internal priming. Increment truth:decoy against the long-read
truth is **0.83–0.89 : 1** for every relaxation, against **3.67 : 1** for the evidence we already
have — i.e. more than half of what relaxation adds is IP. The verifier notes the truth class is 4.2×
deeper than the decoy class in that comparison, so depth-matching pushes the increment ratio
**further below parity**: the conclusion is conservative.

*The exception is the span rule*, whose increment is 3.0 : 1 — same quality as existing evidence.
That is **issue15**, and it is the only relaxation on the recommended list.

*Also dead in the same family:* a "corroborated read-end" gate. 75.2 % of read-ends at clip sites are
uncorroborated, but ranking the shipped tier-1 PAS by read-end mass instead of clip molecules
collapses precision at every operating point (top-2,000: P@100 **0.502 vs 0.808**).

*Reproduce:* `A2_clip_sensitivity/tables/T3d*, T4*, T7b*, T10*`; `VERIFY` §2.2.

## 3. Deconvolving merged calls into multiple PAS — DEAD (recovers exactly zero sites)

*Idea:* one call covering two nearby atlas sites is costing us recall.

*Measured:* **0** of the 11,606 chr19+21 atlas sites the default misses are recovered by perfect
deconvolution. Realistic gain **R_det@25 +0.0010**, R_det@10 +0.0005.

*Why (verified two ways):* `score_tool.py` matches with `bedtools closest -s -d -t first` in both
directions, which is **many-to-one**. Empirically, 3,058 slice atlas sites are recalled at 100 bp
using only **2,050 distinct calls** — 33 % of current recall is already many-to-one — and replacing
every call with three copies at (mode−1, mode, mode+1) leaves **P@100 exactly unchanged at 0.7385**.
The scorer already credits both members of a doublet, and it does not penalise duplication. A4's own
permuted-offset and ±1 bp-duplication controls reproduce the **entire** apparent precision gain
(0.7617 real vs 0.7629 permuted vs 0.7626 duplicated).

*If it is ever re-proposed,* the honest framing is: it does not increase the matched count under the
real scorer, and any gain it appears to show is the scorer's permissiveness.

*Reproduce:* `A4_resolution/tables/T0*, T2c*, T4*, T4b*, T4d*`; `VERIFY` §4.

## 4. The `edge100` 3'-boundary "resolution refinement" — WITHDRAWN

*Idea:* gate calls on a sharp 3' coverage boundary to improve base-pair resolution.

*Measured:* gating the 2,895 slice default calls at `edge100 >= 0.85` (keeping 1,736) raises Kinnex
P@10/P@100 retention 0.8636 → 0.9271. But at the **same call count**, simply keeping the 1,736 calls
with the most clip molecules — a number the caller already writes to `pas_support.tsv` — gives
retention 0.9233 with **better** independent-truth precision: Kinnex P@10 **0.7765 vs 0.6959**,
P@100 **0.8410 vs 0.7506**. `edge100` is confounded with clip support (median clip molecules 4 → 9
across the gate) and scores identically at internal-priming decoys (0.839) and true termini (0.836),
so it cannot be used as a precision filter either. **It is a more expensive way to raise the molecule
threshold.**

*Reproduce:* `VERIFY/tables/v6_edge_gate_control.tsv`.

## 5. A second BAM pass for cross-cell statistics — DEAD

Adding distinct clip cells, top-cell share, 3'-end counts at ±5/25/100, distinct cells among those
ends, pileup sharpness ratios and end-position entropy changes held-out AUC by **0.000 to 0.002** on
both truths (chr19+21, 19,242 candidates, 50.9 M BAM records). Not worth a second pass.

*Reproduce:* `A3_scoring_model/tables/` (BAM-channel ablation).

## 6. Two claims that must not be quoted from the A1–A4 reports

* **"+39 % / +48 % recall from a scoring model" (atlas-trained).** Circular. On every empirical
  independent axis the atlas-trained model is **worse** than the current rule: −8 % at matched call
  count, −8 % at matched long-read precision, −11 % at matched high-confidence long-read precision.
  The honest figure is **+6.6 % at matched call count and matched expression** (with +7 precision
  points) from a **Kinnex-trained** score — see **issue14**.
* **"+16 % on the atlas-independent curve."** A matched-odds artefact: at odds 5.88 the model sits at
  n ≈ 84,700 with long-read precision 0.587 against the rule's 0.765; it matches the odds only by
  halving the decoy denominator.

## 7. The ceiling that bounds all of the above

Of the **11,606** chr19+21 detected-gene atlas sites the shipped default does not recall at 100 bp,
only **2,478 (21.4 %)** have *any* candidate of any tier within 100 bp and only **1,505 (13.0 %)**
have a tier-1 & IP-pass candidate. The union of the **most aggressive** variant of all four ideas is
**2,360 sites (20.3 %)** — against a naive sum of 3,584, i.e. **34 % double-counted**, and already
95 % of the ceiling. **Do not add these recall estimates together.** Four in five sites we miss have
nothing to promote, re-score, split or relax into existence.

*Reproduce:* `VERIFY/tables/v5_recoverable_sets.tsv`, `v5_union_matrix.tsv`.

## What would reopen any of these

New **evidence**, not new geometry: a channel the library does not currently expose (long-read-guided
or hexamer-anchored candidate generation, paired bulk 3'-seq), or a different candidate generator.
Re-proposing the same statistic on the same reads should come with a depth-matched null and a
matched-precision comparison on an atlas-independent truth — the three controls that killed these
four ideas.

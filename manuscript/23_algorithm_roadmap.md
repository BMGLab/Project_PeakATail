# 23 — Algorithm roadmap: what to change in the *caller*, ranked by measured headroom

**Status: PROVISIONAL — written 2026-08-21 for the PI, from the four algorithmic-headroom
measurements (`results/algo_headroom/A1_end_pileup`, `A2_clip_sensitivity`, `A3_scoring_model`,
`A4_resolution`) and the adversarial verification of all four
(`results/algo_headroom/VERIFY/VERIFICATION_REPORT.md`, verdict **FIXED**).**
No code was changed, no pipeline was re-run, nothing was committed and nothing was filed on GitHub.

## Number policy for this file

Every number below is one of:

* **[V]** — confirmed by the independent verifier (`VERIFY/VERIFICATION_REPORT.md`,
  `VERIFY/tables/*`). Where the verifier *corrected* an A1–A4 number, the corrected value is used and
  the correction is stated in the same sentence.
* **[19]/[24]** — quoted from `19_final_gate_v2.md` (verifier: FIXED) or `24_prime_preregistration.md`.
* **[A*]** — measured by A1–A4 and **not** independently re-verified; always labelled as such.
* **[A5, new here]** — measured while writing this document, in
  `results/algo_headroom/A5_cross_sample/` (README + `score.sh` + every scored BED). It is
  **PROVISIONAL and has not been seen by the verifier.** It exists because no A1–A4 measurement
  covered the cross-sample row of the table and I would otherwise have had to speculate.

Three A1–A4 headline claims are **withdrawn** by the verifier and are not used anywhere in this
document: A3's "+16 %/+39 % on the atlas-independent curve" (a matched-odds artefact), A3's
"truth-vs-decoy AUC is the non-circular core" (beaten by a covariate the caller already computes),
and A1's `edge100` boundary-resolution refinement (beaten by ranking on clip molecules).

## Relationship to the other planning documents

| document | owns |
|---|---|
| `13_reliability_positioning.md`, `19_final_gate_v2.md` | the **shipped v2 default** and its pre-registered gates. Unchanged by this file. |
| `22_performance_roadmap.md` | **does not exist yet.** If it lands, it owns the *filtering* axis — thresholds, tiers, IP-rule variants, post-calling cell/PAS filters. This file deliberately does **not** cover those; where a finding is a filtering finding I say so and hand it over. |
| **this file (23)** | the **detection algorithm**: what evidence is extracted from the BAM, how candidates are generated, positioned and ranked. |
| `24_prime_preregistration.md` | how any change from this file is *allowed to be measured* on the `peakAtail-prime` branch, and what would have to happen before it could become a default. §3.1 (the primary criterion), §3.3 (Rule T, threshold selection) and §3.4 (everything is exploratory) govern every recommendation below. |

---

## 0. The decisions, on one screen

| # | Decision | Why |
|---|---|---|
| **A1** | **One change dominates: a calibrated per-site score, trained on an *independent* truth, used as a re-ranker inside the existing tier-1 and IP gates, replacing only the `>=2 clip molecules` threshold.** Everything else can wait. | It is the only measured change that moves precision **and** recall at the same call count on a truth the model never saw: **P@100 0.7766 vs 0.7062, R_det@100 0.1870 vs 0.1754, decoy@25 0.0954 vs 0.1300 at n = 46,524 and matched expression** [V]. |
| **A2** | **Its first deliverable is not the score. It is the cross-dataset transfer test.** Fit on GSE104556 mouse 1, evaluate unchanged on PBMC and mouse 2, per [24] §3.3. | Everything measured so far is **one PBMC library** with chromosome-disjoint folds. Chromosome-disjoint CV tests positional generalisation, **not** dataset or species transfer [V]. That is the first question a reviewer asks. |
| **A3** | **The single biggest measured lift is already in the tool and is not on by default: `--ip-filter`.** Default it and *say so in the paper*. | Swept head-to-head at matched precision it buys **+7.7 % / +8.8 % / +12.8 % relative recall** at atlas P@100 = 0.60 / 0.7062 / 0.8358 and **+13.4 % to +22.9 %** at matched long-read precision [V, `v8_ipveto_value.tsv`]. This is a *writing* and *default* change, not an engineering one. |
| **A4** | **Four ideas are dead. Do not spend a week on any of them:** an end-position / read-3'-end pileup caller; relaxing `--polya-min-clip`; deconvolving merged calls; the `edge100` boundary refinement. | All four verified negative, two of them *more* negative after verification (§2 rows A, B, D and §2.1). |
| **A5** | **One cheap, zero-risk change can land before submission and changes no headline number: emit an `inferred_cleavage` column at a calibrated −2 bp, and stop recommending `cleavage_offset = 95`.** | The −2 bp constant is agreed in sign by **both** truths at base-pair matching (atlas P@1 0.2119 → 0.2960, Kinnex P@1 0.1782 → 0.2642) and moves **nothing** at W >= 25 [V, `VERDICT_FOR_PI.tsv`]. The shipped docstring's "try 95" is destructive for the clip-seeded caller (P@10 0.5209 → 0.0551) [A4]. |
| **A6** | **The honest headline of this whole exercise is a ceiling, not a gain.** Of the 11,606 chr19+21 detected-gene atlas sites the shipped default misses at 100 bp, only **2,478 (21.4 %)** have *any* candidate of any tier within 100 bp, and only **1,505 (13.0 %)** have a tier-1 & IP-pass candidate [V]. | Four in five sites we miss have nothing to promote, re-score, split or relax into existence. That sentence belongs in the paper. |

---

## 1. The argument in a paragraph

Everything the team has fixed so far — internal-priming windows, the `>=2` clip-molecule threshold,
hexamer gating, tier splits — is a *filter*, and a filter can only move the tool **along** its
precision/recall trade curve. The evidence for that is arithmetic and it is now measured: on PBMC
the same candidate pool gives (P@100 0.3520, R_det 0.2685) at `>=1` molecule and (0.7062, 0.1754) at
`>=2` [19], and sweeping the threshold traces one curve whose endpoints are fixed by the candidate
set the caller generates. **Lifting** that curve means one of exactly three things: (i) extracting
*more evidence per read* from the same BAM, (ii) *ranking the same candidates better*, so that at
any call count the calls kept are more often real, or (iii) *positioning the same call more
accurately*, which buys resolution rather than recall. The four measurements say, without much
ambiguity, that (i) is nearly exhausted and (ii) is where the remaining headroom lives — and that
both are bounded by a hard ceiling in the library itself. On (i): 96.1–96.3 % of accepted reads carry
no 3'-side terminal soft clip at all, so their 3' end is `start + 91`, a shifted copy of coverage
[A1]; the genome-wide qualifying poly(A) clip rate is **0.5730 %**, not the 1.152 % the code and
`10_caller_fix_plan.md` claim (independently reproduced from scratch on chr21 at 0.003669 by the
verifier) [V]; and at *verified* PAS the detector already captures ~73 % of the tail-bearing reads,
with the remaining short-A class ~9× less specific — relaxing `--polya-min-clip` from 6 to 3 buys
+28.9 % more clip reads and **loses** 0.7–10.0 % recall at matched precision, because more than half
of what it adds is internal priming (increment truth:decoy 0.83 : 1 against 3.67 : 1 for the evidence
we already have) [A2, and the verifier notes depth-matching makes this *worse*, not better]. On the
pileup idea specifically, the verifier killed it outright: A1's null was drawn uniformly inside gene
bodies (median local depth 2 molecule-ends per ±100 bp against 8 at missed atlas sites and 243.5 at
recovered ones), and after direct standardisation to the missed-atlas depth distribution **a
molecule-end pileup at a true missed PAS is exactly as likely as at a random position of the same
local depth — 0.95×, 0.95×, 1.00×, 1.01× at the mol5 >= 2/5/10/20 bars — against 11.36× for the clip
channel** [V]. On (ii), the one unambiguous positive in the whole exercise: a gradient-boosted score
trained **only on the Kinnex long-read truth**, with all downstream-A features dropped, on
chromosome-disjoint folds, and then evaluated on the PolyASite atlas **it has never seen**, beats the
`>=2` rule at identical call count *and* identical window-read decile on every axis at once —
P@10 0.5567 vs 0.5209, P@25 0.6835 vs 0.6389, P@100 0.7766 vs 0.7062, R_det@100 0.1870 vs 0.1754,
internal-priming decoy rate 0.0954 vs 0.1300 — or +14.3 % relative recall at matched atlas precision
[V]. That is a lift, it is small (+7 % recall with +7 precision points, **not** the +39 %/+48 % the
atlas-trained model advertised), and it is bounded: the candidate pool itself caps R_det@100 at
0.269 for tier-1 & IP-pass, the aggressive union of *all four* ideas reaches only 20.3 % of the
missed atlas sites against a naive sum of 30.9 % (34 % double-counted), and 78.6 % of the sites we
miss have no candidate of any tier within 100 bp at all [V]. **The next real gain after the score is
not in the detection geometry. It is in the candidate generator or in the chemistry.**

---

## 2. The headroom table

One row per algorithmic idea. "Upside" is always **recall at matched precision** or **resolution**,
never recall alone — a recall number at a lower precision is a slide along the curve, not a lift.
"Pre-reg?" means: would adopting this change the pre-registered default of `13_reliability_positioning.md`
§1, and therefore require its own pre-registration under [24] §3.4 before any manuscript number moves?

| # | Idea | Measured upside | Measurement | Risk it is illusory | Difficulty | Pre-reg? |
|---|---|---|---|---|---|---|
| **A** | **End-position / read-3'-end pileup modelling** (call from where reads end, not from coverage) | **None, and negative.** −2.6× to −4.3× recall at matched precision on both chromosomes and both truths; as an added channel the extra calls cap at 36.1 % precision against the default's 73.9 % | A1 stages 0–7; verifier `v4_a1_null.py` | **Already realised.** The published 2.0–7.9× lifts were a *local-read-depth artefact*; depth-standardised they are **0.95–1.01×** [V] | — (do not build) | n/a |
| **B** | **More sensitive clip detection** (`--polya-min-clip` 6 → 3/2/1, purity, boundary run) | **None.** −0.7 % to −10.0 % atlas recall at matched precision; increment truth:decoy **0.83 : 1** vs 3.67 : 1 for existing evidence. *Sub-item:* dropping the `read.py` span rule is **+1.76 % clip reads genome-wide, projected ~+0.4 % recall** — the only variant positive on both truths at every operating point | A2 T3d/T4/T7/T7b; verifier §2.2, §5 | Low risk of being *illusory*; the verifier notes the truth class is 4.2× deeper than the decoy class, so depth-matching pushes the increment ratio **further below parity**. The span-rule sub-item is real but A2 already discounted its own slice estimate 2.6× | Trivial (a constant) / Small (one predicate) | Yes for the span rule (it changes the candidate set **and** the count matrix) |
| **C** | **Calibrated per-site score replacing the hard `>=2` gate** (tier-1 and IP veto retained as gates) | **+6.6 % relative recall AND +7.0 precision points at matched call count and matched expression** (R_det 0.1870 vs 0.1754, P@100 0.7766 vs 0.7062, decoy 0.0954 vs 0.1300); **+14.3 % relative recall at matched atlas precision** (0.2005 vs 0.1754) | Verifier §1.5 / §7.2, `v9_expr_matched_kinnex_trained.tsv`, `v3_decontamination.tsv` — Kinnex-trained, downstream-A features dropped, chromosome-disjoint folds, evaluated on the atlas it never saw | **The real risk is generalisation, not circularity.** Circularity was tested four ways and the transfer result survives; but the model was fitted and tested on **one PBMC library**. Chromosome-disjoint CV is not dataset or species transfer [V]. Also: **P@10 does not improve at matched precision** (0.4934 vs 0.5209) — a re-ranker cannot move a call | Medium (offline training, numpy-only inference, ~50 features all already computable at call time) | **Yes — this replaces the pre-registered default.** [24] §3.4 applies in full |
| **D** | **Multi-PAS deconvolution** (split merged calls into their component sites) | **Zero.** 0 of the 11,606 missed atlas sites recovered; realistic gain R_det@25 **+0.0010**, R_det@10 +0.0005 | A4 T0/T2c/T4/T4b/T4d; verifier §4 | **Already realised.** The scorer's matching is many-to-one: 3,058 slice atlas sites are recalled at 100 bp using only **2,050 distinct calls** (33 % many-to-one), and tripling the call set at (mode−1, mode, mode+1) leaves P@100 **exactly** unchanged at 0.7385 [V]. A4's own permuted-offset and ±1 bp-duplication controls reproduce the entire apparent precision gain | — (do not build) | n/a |
| **E** | **Cleavage-offset calibration** (constant −2 bp, or a support-conditional offset) | **Resolution only, at base-pair matching.** Atlas P@1 0.2119 → **0.2960**, Kinnex P@1 0.1782 → **0.2642**; both truths agree in sign. **Nothing moves at W >= 25** (P@100 0.7064 → 0.7064) | A4 T5/T5b/T5e, endorsed in `VERDICT_FOR_PI.tsv` [V] | Low, but **not universal**: on mouse the P@1 optimum is −1 bp, not −2 (0.3682 vs 0.3166 unshifted) [A4], so the constant is dataset/chemistry-dependent. It also costs a little at intermediate windows on the atlas (P@10 0.5209 → 0.5093) while gaining on Kinnex (0.4979 → 0.5108) [A4] | Trivial as a **column**; a shift of the reported coordinate is not trivial | **No, if shipped as an extra column.** Yes, if it moves the reported coordinate |
| **F** | **Internal priming as a covariate rather than a veto** | **As a gate: the largest measured lift in the whole exercise — +7.7 % / +8.8 % / +10.7 % / +12.8 % relative recall at atlas P@100 = 0.60 / 0.7062 / 0.80 / 0.8358, and +13.4 % / +17.7 % / +22.9 % at matched long-read precision.** As a *feature*: the binary flag alone is worth ~nothing (AUC 0.534), but the tool's own **continuous** `ip_tool_afrac`, inverted, scores AUC **0.7790** at separating genuine long-read termini from IP decoys — beating the entire A3 model's 0.7655 | Verifier §7.1 (`v8_ipveto_value.tsv`) and §6(c) (`v3_truth_vs_decoy_auc.tsv`) | Low. The obvious contamination hypothesis (the Kinnex decoy class *is* a downstream-A rule, and the features contain its ingredients) was **tested and rejected**: refitting with all 12 downstream-A features removed changes atlas R_det 0.2445 → 0.2432 and decoy 0.0593 → 0.0603 [V]. Caveat both truths carry an IP sequence rule, so neither is neutral on this question | Trivial (flip a default) + Small (expose one float) | **Yes** — `--ip-filter` becoming default changes the shipped default, though *not* the pre-registered arm, which already has it on |
| **G** | **Cross-cell and cross-sample evidence** | **Cross-cell: none.** A second BAM pass for distinct clip cells, top-cell share, end counts at ±5/25/100, pileup sharpness and end-position entropy moves held-out AUC by **0.000–0.002** [V]. **Cross-sample: real but modest, and only in multi-sample designs.** Mouse-1 1-molecule sites corroborated by a mouse-2 default call within 10 bp score **P@100 0.6683** against **0.4144** for expression-matched singletons and 0.3941 for the whole singleton pool; at matched call count (n = 29,079) the corroborated union beats both controls on **both** axes (P@100 0.7374 / R_det 0.2207 vs 0.7129 / 0.2143 expression-matched) — but against the shipped default it is +7.8 % relative recall for **−0.0075 precision**, which **fails** [24] §3.1(i) | [A5, new here] `results/algo_headroom/A5_cross_sample/` | **Untested by the verifier**, and there is a mechanism that must be checked before anyone believes it: **internal priming replicates across samples too** — the genome is the same in both mice — so corroboration cannot suppress the residual-IP class the way an independent truth can. It also **cannot help PBMC**, which is a single library and the paper's headline dataset | Medium (cohort-mode plumbing already exists; single-sample mode gains nothing) | **Yes**, and only for cohort mode |

### 2.1 The two positive claims that were withdrawn, so nobody re-proposes them

* **A1's `edge100` 3'-boundary "resolution refinement".** Applied as written (gate the 2,895 slice
  default calls at `edge100 >= 0.85`, keeping 1,736) it gives Kinnex retention 0.9271 — but at the
  same call count, **simply keeping the 1,736 calls with the most clip molecules**, a number the
  caller already writes to `pas_support.tsv`, gives retention 0.9233 with *better* independent-truth
  precision (Kinnex P@10 **0.7765 vs 0.6959**, P@100 **0.8410 vs 0.7506**) [V, `v6_edge_gate_control.tsv`].
  `edge100` is also confounded with clip support (median clip molecules 4 → 9 across the gate). It is
  a more expensive way to raise the molecule threshold.
* **A3's atlas-trained headline.** "+39 % atlas / +16 % on the atlas-independent curve" is not
  available. The atlas figure is circular; the "independent" figure is a **matched-odds artefact** —
  at odds 5.88 the model sits at n ≈ 84,700 with **KinP_t5 0.587 against the rule's 0.765**, i.e. it
  matches the odds only by halving the decoy denominator while being 18 precision points worse on
  real long-read termini. On matched call count it is −8 %, on matched KinP_t5 −8 %, on matched
  KinP_t20 −11 % [V]. **Row C uses the Kinnex-trained transfer result instead**, which is one third
  the size and survives every axis.

### 2.2 Do not add these numbers up

On one common reference set — the 11,606 chr19+21 detected-gene atlas sites the shipped default misses
at 100 bp — the four ideas in their headline configurations recover 162 (A1), 608 (A2 k>=3), 474 (A3
Kinnex-trained at matched precision) and 0 (A4). **Naive sum 1,244; true union 1,040 (16 % double-
counted).** In their most aggressive configurations: naive sum 3,584, **true union 2,360 = 20.3 %**
(34 % double-counted) [V, `v5_union_matrix.tsv`]. Pairwise: A1∩A2 = 88, A2∩A3 = 109, A1∩A3 = 11,
anything∩A4 = 0. The aggressive union is already **95 % of the ceiling** (2,478 sites with any
candidate at all) — the ideas are competing for the same small pool, not adding to each other.

### 2.3 Slice representativeness (applies to every chr19+21 number above)

chr19+21 is **27 % clip-richer** (0.7254 % vs 0.5730 % genome-wide), **18 % tier-1-richer**
(tier1:tier2 1.19 vs 1.01) and **19 % easier on recall** (R_det@100 0.2085 vs 0.1754) than the genome
[V, `v7_representativeness.tsv`]. Every bias runs the same way, so **slice negatives are conservative
and any slice-only positive is optimistic.** Rows A, B and D are negatives measured on the slice and
are therefore safe. Row C's headline is genome-wide (all 402,765 PBMC candidates). Row G is
genome-wide on mouse.

---

## 3. The recommended programme

**One change dominates. Say it plainly: build the score, and do not spend engineering time on
anything else in this document until the score has survived a second dataset.**

### Step 0 — free, this week, no code: stop under-selling the filter we already ship

`--ip-filter` is an opt-in flag (`ema/cli/config_schema.py:529`, `is_flag=True`), and no document in
the manuscript set quantifies what it is worth on a like-for-like curve. It is worth **more than every
detector change tested, combined**: +7.7 % to +12.8 % relative recall at matched atlas precision and
+13.4 % to +22.9 % at matched long-read precision [V]. It lifts rather than slides because it brings
in information the thresholds do not have — genomic sequence. Three actions:

1. Make it **default-on** when `--genome-fasta` is available (Issue 13). This does not change the
   pre-registered arm, which already runs with it.
2. Add the measured curve to the paper (§5). Ideally replicate the sweep on mouse first — the
   like-for-like no-IP v2 mouse arm **does not exist** (only `peakatail_clipseeded_v3/mouse1`, which
   is Stage-1c code), and a mouse re-run without `--ip-filter` costs ~9 min and 3.7 GB [19] §3.
3. Correct the clip-rate claim everywhere it appears (`ema/countmatrix/polya.py` docstring,
   `10_caller_fix_plan.md`): **0.5730 % genome-wide, not 1.152 %**, and the head-sample QC estimator
   returns 2.2565 % where the truth is 0.5364 % [V, and the chr21 leg independently reproduced].
   The sampling half of this is already open as **#99 §2**; add the corrected constant to it.

### Step 1 — the score, but the transfer test first (Issue 14)

The deliverable order is deliberately inverted, because the result that decides everything is not the
model, it is whether the model transfers.

1. **Transfer test before any implementation.** Refit the verifier's exact recipe
   (`VERIFY/code/v3_a3_decontaminate.py`: `HistGradientBoostingClassifier(max_iter=400,
   learning_rate=0.06, max_leaf_nodes=31, min_samples_leaf=100, l2_regularization=1.0)`, Kinnex-t5
   label, all 12 downstream-A features dropped) on **GSE104556 mouse 1** and evaluate unchanged on
   **PBMC** and **mouse 2** — the primary protocol of [24] §3.3. Mouse has no Kinnex truth, so the
   mouse-side label has to be the mouse atlas, which makes the mouse→PBMC direction *atlas-trained*
   and therefore only a generalisation test, not a headline. **If it does not transfer, stop here and
   report that as the result.**
2. Only if it transfers: implement the score as **numpy constants evaluated inline** — scikit-learn
   must not enter the tool's runtime import path [24] §4/§3.3 — behind `--score-model`, defaulted OFF
   on the prime branch until the §3.1 criterion is met on all three datasets.
3. **Keep tier-1 and the IP veto as hard gates.** Every configuration in which the score was allowed
   to override the IP veto looked spectacular on the atlas and no better than the rule on long reads.
   Applying the veto at *selection* time on top of the model still cuts the decoy rate from 0.088 to
   0.060 [A3, consistent with V].
4. **Ship the probability as a column** regardless of what the default threshold ends up being. A
   calibrated per-site confidence is a better paper claim than a better set, and it turns `>=2` from a
   commitment into a default. Calibration against the training label is essentially perfect
   (ECE 0.0038); against the *other* truth the same probability is optimistic in the mid-range
   (predicted 0.6 → observed 0.48) [A3] — which is itself a publishable demonstration of what
   curated-atlas benchmarking rewards.
5. **Do not build**: a second BAM pass for clip-cell counts or end-pileup statistics (0.000–0.002 AUC
   [V]); a model that merely re-weights the existing clip counts (+1.9 % to +5.8 %, inside the noise of
   a threshold change [A3]); any claim of improved base-pair resolution from the score (P@10 does not
   improve at matched precision).

### Step 2 — the `read.py` span rule (Issue 15)

`ema/countmatrix/read.py:106` — `if span > seq_len: return 0,0,0,0,0` — discards **13.74 % of all
valid-CB reads genome-wide, 96 % of them spliced**, before the clip detector and before the count
matrix [A2]. For *detection* it is small (+1.76 % clip reads, projected ~+0.4 % recall) but it is the
only relaxation positive on **both** truths at every operating point, and its increment has the same
truth:decoy quality as the evidence we already have (3.0 : 1 vs 3.67 : 1) [A2]. **The larger prize is
quantification: 88.8 M reads never enter the matrix** because a spliced alignment's *reference* span
exceeds `--seq-len`. This should be a span-vs-**query**-length test. It is Step 2 rather than Step 1
only because it changes both the candidate set and the matrix, so it needs the full three-dataset
re-run and a re-verification, and its detection payoff is ~2 % of the score's.

### Step 3 — the cleavage-offset column, and killing a bad recommendation (Issue 16)

Cheapest item in the document and the only one that can land before submission without touching a
headline number: emit `inferred_cleavage = cleavage − 2 bp` (transcript orientation) as a **column**,
never as a coordinate move. Both truths agree in sign at base-pair matching; nothing at W >= 25 moves.
At the same time, **retract the "try 95" recommendation** in `ema/cli/config_schema.py:503` and in
the `issue6_cleavage_offset.md` draft: the +95 bp offset was estimated for the *coverage* caller and
is destructive for the clip-seeded one (P@10 0.5209 → 0.0551, P@100 0.7062 → 0.6655) [A4].

### Step 4 — cross-sample corroboration, **cohort mode only**, and only after verification (Issue 17)

[A5, new here] shows the signal is real and is not expression: corroborated singletons score P@100
0.6683 against 0.4144 for an expression-matched control, replicated in the reverse direction. But it
fails the pre-registered materiality allowance, it cannot help a single-library dataset (i.e. it
cannot move the PBMC headline at all), and residual internal priming replicates across samples by
construction. It belongs in the cohort/Stage-3 path where a cohort PAS space already exists, not in
the caller's single-sample default. **Send [A5] to the verifier before anyone builds on it.**

### Never (closed by measurement — record them so they are not re-proposed)

An end-position/pileup caller (row A); relaxing `--polya-min-clip` (row B); deconvolving merged calls
(row D); the `edge100` boundary gate (§2.1); a second BAM pass for cross-cell statistics (row G);
selling the atlas-trained model's +39 %/+48 % (§2.1).

---

## 4. What is actually a ceiling — chemistry, not code

These are the parts of the gap that no algorithm reaches, with the measured proportion.

| the gap | measured size | source |
|---|---|---|
| Accepted reads with **no 3'-side terminal soft clip at all** — their 3' end is `start + 91`, a shifted copy of coverage | **96.14 % (chr21) / 96.27 % (chr19)** | [A1] stage 0 |
| Genome-wide **qualifying poly(A) clip rate** (corrects the published 1.152 %) | **0.5730 %** (3,195,067 / 557,564,408 accepted CB reads); chr21 leg independently reproduced at 0.003669 from a from-scratch reimplementation | [V] §5 |
| Reads terminating at a **verified** PAS that carry **no tail in the read at all** — the fragment ends at or before the poly(A) junction | **40.46 %** | [A2] T12 |
| …of the same reads, the class the detector **already captures** | **38.44 %** (~73 % of the tail-bearing reads) | [A2] T12 |
| …the short (1–5 nt) A tail the criteria reject — recoverable in principle, **~9× less specific** | **13.93 %** | [A2] T12 |
| Detected-gene atlas sites with **no poly(A) clip read at all** within 25 bp | **76 %** | [A4] T6 |
| Ceiling of an **infinitely sensitive** clip detector (`min_clip = 1`): fraction of detected-atlas sites with >= 1 / >= 2 clip molecules within ±25 bp | **34.9 % / 22.7 %** (current criteria 24.0 % / 15.3 %) | [A2] |
| Missed atlas sites with **any candidate of any tier** within 100 bp — the bound on everything in §2 | **21.4 %** (2,478 / 11,606); tier-1 & IP-pass **13.0 %** | [V] §3 |
| Hard clips hiding evidence | **0** in 663,383,342 records | [A2], [V] §5 |

**The sentence the paper should use** (Discussion / Limitations, and it should appear almost verbatim):

> The dominant limit on de novo 3'-end calling from 10x Genomics 3' data is the library, not the
> algorithm: only 0.57 % of cell-barcoded reads carry a qualifying non-templated poly(A) soft clip,
> 96 % of accepted reads carry no 3'-side clip at all — so their inferred 3' end is a fixed offset
> from their alignment start and therefore a shifted copy of coverage — and at cleavage sites we
> independently verify with long reads, 40 % of the reads that terminate there contain no tail in the
> read whatsoever. We measured the consequence directly: of the curated atlas sites in detected genes
> that our precision-first default does not recall, only 21 % have any candidate of any kind within
> 100 bp, and a molecule-end pileup at a true missed site is no more likely than at a random position
> of the same local read depth (0.95–1.01×), while a poly(A) clip is 11× more likely. Better peak
> geometry cannot recover what the chemistry did not record; the remaining headroom is in ranking the
> candidates that do exist, and in orthogonal evidence.

---

## 5. Impact on the manuscript

**This work belongs in THIS paper, as one Results paragraph, one Discussion paragraph and one
supplementary panel — not as a follow-up paper.** Reasons, in order of weight:

1. It is the **negative result the paper's own thesis needs.** [21] D1 fixes the paper as a
   *calibrated-reliability* methods paper whose fifth honest result is that "the discriminator is
   evidence type, not peak-calling sophistication". §1 and §4 above are the measurement of exactly
   that claim, and without them it is an assertion.
2. The **`--ip-filter` curve (§3 Step 0) is a Results-grade addition that needs no code change** and
   directly strengthens L1: it shows the precision-first default's advantage comes from *sequence
   information*, not from a threshold. This is the one item I would actively push into the paper.
3. The **depth-matched null** (§1) is a methodological point worth making once, in Methods or a
   supplementary note: "lift over a genic-shuffle null" is not a safe statistic for anything
   3'-biased, because a naive null sits at median local depth 2 against 8 at missed atlas sites and
   243.5 at recovered ones, and that alone manufactures 2–8× lifts [V].
4. The **atlas-vs-long-read decomposition** (an atlas-trained score wins 16 precision points on the
   atlas while losing 8–11 % on every long-read axis, and the sites it swaps in are ~4× better by the
   atlas and ~2× worse by long reads) is a publishable methods result about *benchmarking*, and it
   fits the paper's Fig 6 negative-result slot without needing a new figure.

**What is cheap enough to land before submission, and what it does to the headline numbers:**

| change | code cost | effect on the headline numbers |
|---|---|---|
| `--ip-filter` documented as the largest measured lift; default-on | one line + docs | **None.** The pre-registered arm already has it on; P@100 0.7062 / 0.7450 / 0.7572 and R_det 0.1754 / 0.2048 / 0.2080 [19] are unchanged |
| clip-rate correction 1.152 % → 0.5730 %, + fix the head-sample estimator | docstring + `10_caller_fix_plan.md` + #99 §2 | **None** on P/R. It corrects a *stated fact* that is currently 2× wrong, and removes a QC estimator that is 4× off and could equally mask a genuinely destroyed evidence channel |
| `inferred_cleavage` column at −2 bp | small, additive | **None at W >= 25** (P@100 0.7064 → 0.7064). Enables a bp-exact resolution sentence: atlas P@1 0.2119 → 0.2960, Kinnex P@1 0.1782 → 0.2642 |
| retract `cleavage_offset = 95` from docs / `issue6` | docs only | **None** — but it prevents a reader reproducing our run with a flag that costs 46 points of P@10 |
| the score (row C) | medium | **Would move every headline number** ⇒ needs its own pre-registration, a fresh three-dataset run and a verifier pass [24] §3.4. **Do not attempt before submission.** |
| the span-rule fix (row B) | small code, large blast radius | Changes candidates **and** the matrix ⇒ same as above. **Follow-up.** |

**Follow-up paper material** (if the score transfers): "a calibrated per-site confidence for
single-cell 3'-end calls, trained on long reads and validated on a curated atlas", with the
transfer test as its central experiment and the atlas-vs-long-read decomposition as its methods
contribution. That is a real paper; it is not this one.

---

## 6. Issue drafts

Written to `manuscript/github/`, continuing the local numbering after `issue12_*` (GitHub #99).
**Nothing filed, nothing committed.** Each cross-references the open issues so Amir sees one
programme rather than five requests.

| draft | title | maps to | cross-refs |
|---|---|---|---|
| `issue13_ip_filter_default_and_cliprate.md` | Make `--ip-filter` default-on; correct the poly(A) clip-rate constant (1.15 % → 0.57 %) | §3 Step 0 | #95 (IP rule), **#99 §2** (clip-rate sampling — this is its missing constant) |
| `issue14_calibrated_site_score.md` | Replace the `>=2 clip molecules` gate with a calibrated per-site score — **transfer test first** | §3 Step 1 | #94, #95, [24] §3.1/§3.3/§3.4 |
| `issue15_readpy_span_rule.md` | `read.py` span rule discards 13.7 % of valid-CB reads (96 % spliced) before both the clip detector and the count matrix | §3 Step 2 | #98 (per-isoform), #99 (gene assignment) |
| `issue16_cleavage_offset_column.md` | Emit an `inferred_cleavage` column at −2 bp; retract the `cleavage_offset = 95` recommendation | §3 Step 3 | **issue6 draft (#72 lineage)** — this partially retracts it |
| `issue17_cross_sample_corroboration.md` | Cohort mode: corroborate 1-molecule sites across samples (PROVISIONAL, verifier first) | §3 Step 4 | #94, `20_stage3_replication.md` |
| `issue18_dead_ends_do_not_reopen.md` | Record of four measured dead ends, so they are not re-proposed | §2 rows A/B/D, §2.1 | all of the above |

---

## 7. Provenance and how to reproduce

| claim block | where it lives |
|---|---|
| A1–A4 primary measurements | `results/algo_headroom/{A1_end_pileup,A2_clip_sensitivity,A3_scoring_model,A4_resolution}/` (each has `tables/`, its scripts, and a `VERIFIER_ADDENDUM.md`) |
| every corrected number in this file | `results/algo_headroom/VERIFY/VERIFICATION_REPORT.md`, `VERIFY/tables/v2…v9`, `VERIFY/code/v1…v9` |
| the one-line verdict per idea | `results/algo_headroom/VERIFY/tables/VERDICT_FOR_PI.tsv` |
| row G / [A5] | `results/algo_headroom/A5_cross_sample/` — `README.md` (design + controls), `score.sh`, `A5_cross_sample_summary.tsv`, every scored BED. Reproduce with `./score.sh union_k2_w10.bed UNION` |
| shipped baselines | `19_final_gate_v2.md`; scorer `scripts/benchmark_tools/score_tool.py` with the reference arguments of `scripts/benchmark_tools/stage2_final_launch.sh` |

[A5] fidelity check: its 12-line scorer reproduces the published mouse arms exactly — mouse 1
n 26,263 / P@10 0.5796 / P@25 0.6819 / P@100 0.7449 / R_det 0.2048 against [19]'s 0.579661 /
0.682042 / 0.745001 / 0.204790, and mouse 2 0.5944 / 0.6929 / 0.7572 / 0.2080 against 0.594473 /
0.692980 / 0.757219 / 0.208002.

## 8. What this document does **not** answer

1. **Does the score transfer across datasets or species?** Untested and untestable from the existing
   artefacts. It is Step 1's first deliverable and the reviewer's first question.
2. **Is there a better candidate generator?** Every idea in §2 re-ranks, re-positions or relaxes the
   *existing* candidate set. The ceiling (§4) says that set is the binding constraint: 78.6 % of the
   atlas sites we miss have no candidate of any tier within 100 bp, and tier-2 coverage peaks are not
   a reservoir (P@100 0.0568 [19]; no end statistic rescues them above ~0.09 [A1]). **Nobody has
   measured what a different candidate generator would find** — that is the next measurement, and it
   is a bigger question than anything in this file.
3. **Orthogonal evidence.** Nanopore/Kinnex-guided calling, a hexamer-anchored generator, or paired
   bulk 3'-seq would each add a channel rather than re-rank one. Out of scope here; the tool is
   de novo by design ([21] D1) and the Kinnex truth is a read-out only ([24] §3.3).
4. **The filtering axis.** Post-calling cell/PAS filters drop 31.2 % of tier-1 clip clusters on
   chr19+21 (5,095 of 16,338), including 217 with >10 clip molecules; the dropped set is genuinely
   worse (P@100 0.192 vs 0.360) so the filter is doing sensible work [A2], and un-dropping them
   recovers 317 missed atlas sites (2.73 %) [V]. That seam belongs to the filtering roadmap
   (`22_performance_roadmap.md` when it exists), not here.

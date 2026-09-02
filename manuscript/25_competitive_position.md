# 25 — The competitive position (verified; the quotable source for every other document)

**Status: VERIFIED-SOURCED, 2026-08-21.** Every number in this document is either (i) recomputed in
`results/perf_gap/` or `results/algo_headroom/` and confirmed by the adversarial verifier
(`results/perf_gap/VERIFY/README.md`, verdict **FIXED**, 69 checked quantities in
`out/VERIFICATION_RESULTS.tsv`; `results/algo_headroom/VERIFY/VERIFICATION_REPORT.md`), or (ii)
quoted from an already-verified manuscript document ([19](19_final_gate_v2.md) §1,
[15](15_final_gate.md) §3, [16](16_trusted_novel_kinnex.md)). **No number here is new, modelled or
estimated.** Where a claim is post-hoc, it is labelled post-hoc and may not be used as evidence.

This file exists because the standing one-liner — *"the precision default is the most atlas-concordant
de novo call set in the panel, at polyApipe-or-lower recall (−12 % PBMC, −18 % / −17 % mouse), F1
essentially tied"* — is **misleading as written**. It is a true statement about one operating point and
a false statement about the tool. Every document that compares PeakATail with a competitor should quote
**this** file. Roadmap context, ranked fixes and issue drafts live in
[22](22_performance_roadmap.md) and [23](23_algorithm_roadmap.md); this file is the *claim* layer, they
are the *action* layer.

---

## 0. Definitions, denominators, and what makes the comparison legitimate

All numbers below use **one scorer** and **identical denominators** for every tool
([19](19_final_gate_v2.md); `scripts/benchmark_tools/score_tool.py`):

| quantity | definition |
|---|---|
| `P@100` | fraction of called 1-bp points with a **strand-matched** PolyASite 2.0 representative site within 100 bp (`bedtools closest -s -d -t first`); atlases 569,005 (GRCh38) / 301,006 (GRCm38) |
| `R_det@100` | fraction of the **detected-gene-restricted** atlas recovered within 100 bp; denominators **285,136** (PBMC 10k v3) and **126,686** (GSE104556 mouse) |
| `F1_det` | harmonic mean of the two |
| **hard FP** | a call with **no** strand-matched atlas site within 100 bp **and** no Kinnex x3p ≥5-UMI long-read 3′ end within 25 bp (human only; the mouse variant is atlas-only and **not** numerically comparable) |

Precision here is **atlas-agreement** precision, not ground-truth precision: an atlas-novel true site
counts as a false positive. That caveat is unchanged and still travels with every number.

**The machinery was validated three ways before any comparison was made**
(`results/perf_gap/D3_head_to_head/README.md` §0, all confirmed by the verifier): the independent
scorer reproduces the published `score_*.tsv` values to six decimals (PeakATail default
0.706195 / 0.175443; polyApipe 0.379966 / 0.198796); the hexamer annotation matches the shipped
`site_annotations.tsv` on 46,524/46,524 sites; and the re-applied internal-priming rule agrees with the
caller's own per-peak flag on 651,196/651,196 raw peaks.

---

## 1. The correction, in one paragraph

The single-operating-point head-to-head compares **46,524 PeakATail calls with 120,916 polyApipe
calls**. Precision and recall both move with call budget, so a comparison at unequal budgets measures
the budget, not the tool. When the budgets are equalised — each tool truncated to the same N under
**its own shipped ranking** — PeakATail leads polyApipe on **recall at every N tested**, on both
datasets, and on **precision at every N from 20,320 upward** (with the tie-break qualification of §2.3
inside 20,320 ≤ N ≤ 35,759; below N ≈ 15,000 polyApipe's precision is the higher of the two — §2.3.3);
and on **both** mice an arm that is *already published*
(tier-1 ∩ IP-pass ∩ ≥1 molecule) beats polyApipe on precision, recall **and** F1 simultaneously. The
pre-registered ≥2-molecule default is a deliberately conservative point on a curve that dominates
polyApipe's single point; it is **not** the tool's frontier and it is **not** its F1 optimum.

---

## 2. Matched call count — the like-for-like comparison

Source: `results/perf_gap/D3_head_to_head/tsv/01d_matched_n_headline_pbmc.tsv`,
`01e_matched_n_headline_mouse1.tsv`, full grids `01_matched_n_pbmc.tsv` / `01_matched_n_mouse1.tsv`.
Verifier: *"D3's whole matched-N grid … reproduced exactly"*, with one qualification (§2.3 below).

### 2.1 PBMC 10k v3 — P@100 / R_det@100

| N | PeakATail | polyApipe | scAPAtrap | Sierra | SCAPTURE\* | scUTRquant† |
|---:|---|---|---|---|---|---|
| 9,456 | 0.9527 / **0.0609** | **0.9671** / 0.0587 | 0.4236 / 0.0224 | 0.6412 / 0.0354 | — | 0.8754 / 0.0529 |
| 13,560 | 0.9412 / **0.0828** | **0.9483** / 0.0777 | 0.4324 / 0.0327 | 0.5788 / 0.0456 | — | 0.8774 / 0.0756 |
| 14,939 | 0.9358 / **0.0896** | **0.9395** / 0.0832 | 0.4247 / 0.0353 | 0.5613 / 0.0486 | — | 0.8773 / 0.0830 |
| 20,320 | **0.9074 / 0.1122** | 0.8708 / 0.0988 | 0.4157 / 0.0472 | 0.5069 / 0.0590 | — | 0.8698 / 0.1095 |
| 24,043 | **0.8797 / 0.1248** | 0.8277 / 0.1077 | 0.3908 / 0.0523 | 0.4743 / 0.0648 | — | 0.8631 / 0.1264 |
| 29,814 | **0.8358 / 0.1419** | 0.7438 / 0.1160 | 0.3629 / 0.0597 | 0.4361 / 0.0729 | — | 0.8446 / 0.1493 |
| 35,759 (SCAPTURE) | **0.7579 / 0.1507** | 0.6861 / 0.1245 | 0.3383 / 0.0660 | 0.4066 / 0.0804 | 0.6518 / 0.1176 | 0.8165 / 0.1676 |
| 40,519 (scUTRquant) | **0.7263 / 0.1613** | 0.6496 / 0.1310 | 0.3233 / 0.0708 | 0.3855 / 0.0854 | — | 0.7926 / 0.1777 |
| **46,524 (our default)** | **0.7062 / 0.1754** | 0.6049 / 0.1375 | 0.3066 / 0.0763 | 0.3636 / 0.0915 | — | set too small |
| 106,170 (Sierra) | **0.4269 / 0.2204** | 0.4002 / 0.1871 | 0.2216 / 0.1159 | 0.2559 / 0.1350 | — | — |
| **120,916 (polyApipe)** | **0.4036 / 0.2336** | 0.3800 / 0.1988 | 0.2101 / 0.1233 | — | — | — |
| 167,565 (our ≥1-mol arm) | 0.3520 / 0.2685 | — | 0.1858 / 0.1448 | — | — | — |
| 333,920 (our max) | **0.2044 / 0.3043** | — | 0.1448 / 0.2052 | — | — | — |
| 787,138 (scAPAtrap) | not reachable | — | 0.0999 / 0.3004 | — | — | — |

\* SCAPTURE has **no usable ranking column** (BED score all zeros; its only numeric per-peak field is
atlas-derived and therefore circular), so it cannot be truncated: its row is a random subset, not a
top-N, and no matched-N statement about SCAPTURE below its own N is possible.
† scUTRquant is **catalog-based**, not de novo — see §2.4.

### 2.2 GSE104556 testis mouse 1 — P@100 / R_det@100

| N | PeakATail | polyApipe | scAPAtrap | Sierra | SCAPTURE\* | scUTRquant† |
|---:|---|---|---|---|---|---|
| 22,550 (Sierra) | **0.7735 / 0.1865** | 0.7388 / 0.1339 | 0.4851 / 0.0965 | 0.5087 / 0.1071 | 0.6970 / 0.1496 | 0.8352 / 0.2043 |
| 23,645 (SCAPTURE) | **0.7677 / 0.1935** | 0.7260 / 0.1374 | 0.4793 / 0.0993 | — | 0.6935 / 0.1551 | 0.8321 / 0.2109 |
| **26,255 (our default)** | **0.7450 / 0.2048** | 0.6950 / 0.1443 | 0.4687 / 0.1060 | — | — | 0.8228 / 0.2245 |
| 35,732 (scUTRquant) | **0.6728 / 0.2368** | 0.6132 / 0.1677 | 0.4305 / 0.1270 | — | — | 0.7825 / 0.2586 |
| 52,792 (our ≥1-mol arm) | **0.5686 / 0.2806** | 0.4964 / 0.1925 | 0.3802 / 0.1575 | — | — | — |
| **71,983 (our max)** | **0.4488 / 0.3034** | 0.4353 / 0.2231 | 0.3372 / 0.1843 | — | — | — |
| 89,527 (polyApipe) | not reachable | 0.4005 / 0.2499 | 0.3069 / 0.2044 | — | — | — |
| 135,538 (scAPAtrap) | not reachable | — | 0.2502 / 0.2453 | — | — | — |

### 2.3 The three things that must be said with these tables

1. **The tie-break caveat (verifier qualification, `VERIFY/README.md` correction 3).** Both rankings
   are heavily tied — polyApipe's `peakdepth` more so than our clip-molecule counts (64.5 % of its
   calls at depth 1, 18,857 at depth 2, 6,161 at depth 3; ours: 121,041 singletons, 16,710 at exactly
   2). Where a matched N cuts *inside* a tie group the result depends on the secondary key. Under a
   favourable **oracle** tie-break polyApipe reaches P 0.9193 / 0.8611 / 0.7620 at
   N = 20,320 / 29,814 / 35,759 against our 0.9074 / 0.8358 / 0.7579 (our own optimistic bracket at
   35,759 is 0.7917 — the intervals overlap). **So: the recall half of the claim is robust to
   tie-breaking on both sides at every N; the precision half is not established in the band
   20,320 ≤ N ≤ 35,759.** It *is* unambiguous at **N = 46,524** (polyApipe's best case 0.6636 vs our
   0.7062) and at **N = 120,916**, where polyApipe is its own full set with no ties at all
   (0.3800 vs 0.4036). polyApipe cannot realise an oracle tie-break, so the as-shipped comparison
   stands — but quote the matched-N claim at 46,524 and 120,916, not in the ambiguous band.
   (`tsv/01b_tiebreak_sensitivity.tsv`.)
2. **scAPAtrap's three smallest-N rows are not a real ranking** (19,410 of its sites sit at its
   `min(1000, …)` cap, so any N ≤ 19,410 cuts arbitrarily inside the capped group), and
   *"PeakATail at scAPAtrap's N"* is **not reachable** — 787,138 exceeds our largest possible output
   (333,920). Say "not reachable", never extrapolate. (Mouse 1: scAPAtrap's own N = 135,538 likewise
   exceeds our maximum 71,983; `01e` row 135,538, scAPAtrap 0.2502 / 0.2453.)
3. **The precision lead reverses below N ≈ 15,000, and the tables above show it.** At the three
   smallest matched N the as-shipped ranking puts **polyApipe's precision above ours**: 0.9671 vs
   0.9527 at N = 9,456, 0.9483 vs 0.9412 at 13,560, 0.9395 vs 0.9358 at 14,939 (`01d`, column
   `denovo_winner_P` = polyApipe on all three rows). Our recall is higher at each of them, so this is
   a genuine precision/recall exchange at very small budgets, not a tie artefact. **Therefore: the
   precision half of the matched-N claim is stated as "at every N from 20,320 upward", never as "at
   every matched N".** The recall half carries no such restriction.

### 2.4 What the matched-N result does and does not license

* **Licensed:** *"Among **de novo** tools, at matched call count PeakATail leads on recall at every N
  tested and on precision at every N from 20,320 upward, on both datasets."* At polyApipe's own N:
  **0.4036 / 0.2336 vs 0.3800 / 0.1988** (+6.2 %
  relative precision, **+17.5 % relative recall**). At our own N: **0.7062 / 0.1754 vs
  0.6049 / 0.1375** (+16.7 % / +27.6 % relative). Mouse 1 at our max N = 71,983:
  **0.4488 / 0.3034 vs 0.4353 / 0.2231** (+3.1 % / +36.0 % relative).
* **Licensed:** *"scAPAtrap's recall advantage is a call-count artefact."* At 333,920 calls PeakATail
  reaches R_det **0.3043**, above scAPAtrap's full-set **0.3004** at 787,138 calls, with **2.05×** the
  precision (0.2044 vs 0.0999) and 2.4× fewer calls.
* **Licensed, with the §2.3.1 tie qualifier:** *"polyApipe's own errors are removable by a trivial
  filter."* `peakdepth ≥ 3` removes 96.3 % of its hard FPs (57,623 of 59,856) and lifts its P@100
  0.380 → **0.828**, at 24,043 calls and R_det **0.1077** — where our own matched-N point is
  **0.8797 / 0.1248**. The **recall** comparison there is robust; the **precision** comparison sits
  inside the tie-ambiguous band 20,320 ≤ N ≤ 35,759, so do **not** write "and it is then dominated"
  without that qualifier.
* **NOT licensed:** any claim over **scUTRquant**. It is catalog-based and is shown unranked; at
  matched N it is above PeakATail on both axes at several points (PBMC N = 35,759: 0.8165 / 0.1676 vs
  0.7579 / 0.1507; every mouse row). The de novo restriction is not a rhetorical hedge — it is the
  boundary of the claim, and dropping it makes the claim false.
* **NOT licensed:** *"PeakATail leads polyApipe on precision at **every** matched N."* **False** — at
  N = 9,456 / 13,560 / 14,939 polyApipe's precision is the higher of the two (§2.3.3). The recall lead
  holds at every N; the precision lead is bounded below at N ≈ 20,320.
* **NOT licensed:** any matched-N statement about **SCAPTURE** below its own N (§2.1 footnote), or any
  claim in the tie-ambiguous precision band (§2.3.1).

---

## 3. On both mice the competitive target is already met by a published arm

`results/perf_gap/D2_rule_sweep/D2_arm_sweep.tsv` rows `M1_sup_ge1` / `M2_sup_ge1`, competitor rows
recomputed there from the competitors' own frozen score files; the same values are published in
[19](19_final_gate_v2.md) §2 and [01](01_outline_and_journals.md) §A.

| arm (tier-1 ∩ IP-pass ∩ ≥1 molecule) | n | P@100 | R_det@100 | F1_det |
|---|---:|---:|---:|---:|
| **mouse 1 PeakATail** | 52,792 | **0.5686** | **0.2806** | **0.3758** |
| mouse 1 polyApipe | 89,527 | 0.4005 | 0.2499 | 0.3078 |
| **mouse 2 PeakATail** | 51,842 | **0.5880** | **0.2826** | **0.3818** |
| mouse 2 polyApipe | — | 0.4119 | 0.2508 | 0.3118 |

**Higher precision, higher recall and higher F1 simultaneously, on both replicate mice, with no new
code and no matched-N argument at all.** This is the cleanest single refutation of "recall at or below
polyApipe": the arm is in the frozen outputs and is already in the paper as the sensitivity arm.

On PBMC the same thing is true of the *curve* rather than of a shipped arm: the measured operating
curve (default + a random x % of singletons, `D2_arm_sweep.tsv` group 8) crosses polyApipe's point
cleanly — `PI_curve_rand60pct`, n 119,144 (essentially polyApipe's call-set size), scores
**0.4080 / 0.2332 / 0.2968** with Kinnex concordance **0.4785** against polyApipe's
**0.3800 / 0.1988 / 0.2610** and **0.3895** — dominant on all four, and the random top-up is the
*conservative* control, not a tuned rule.

---

## 4. The recall budget — most of the gap is unreachable, and that is measured

Source: `results/perf_gap/D1_recall_decomposition/tables/decomposition_*.tsv` and
`recovery_ceiling_*.tsv`; the class partition is exact (classes sum to 651,196 calls and 285,136 sites)
and was re-derived independently by the verifier. Summary table also in [22](22_performance_roadmap.md) §2.1.

Of the **235,111** detected-gene atlas sites the PBMC default misses:

| what could reach it | share of the misses |
|---|---:|
| (a) a tier-1 IP-pass call exists but has only **1 clip molecule** | 11.28 % |
| (b) a tier-1 call exists but was **removed by the IP filter** | 3.45 % |
| (c) only a **tier-2** (coverage-only) call exists | 4.78 % |
| (d) only a **raw pre-call peak** exists | 0.52 % |
| **(a)+(b)+(c)+(d) = everything any threshold, filter or tier policy can ever reach** | **20.04 %** |
| | *(the four component percentages are rounded and sum to 20.03; 20.04 % is the value computed from the counts)* |
| **(e) no peak of any kind within 100 bp** | **79.96 %** (188,006 sites) |

Mouse 1 (100,742 misses): reachable **15.47 %**, class (e) **84.53 %** (85,153 sites).

**Hard ceiling.** Keeping *every one* of the peaks the caller ever produced reaches **R_det 0.3406**
(PBMC: 651,196 peaks, 14× the default's output, at P@100 0.1493) and **0.3278** (mouse 1: 195,498
peaks, 7.4×, at P@100 0.2799). An oracle that
recovered all of it with no new false positives reaches the same 0.3406 at P 0.854. **No scoring rule,
threshold, IP rule or tier policy downstream of peak calling can exceed those numbers**, because there
is nothing there to promote. (Consequence for [07](07_curated_benchmark_report.md)'s "Roadmap
publication bars": the R ≥ 0.60 bar is unreachable without changing peak calling itself — see
[22](22_performance_roadmap.md) §5.)

**Class (e) is not an assay-reach failure.** Only **1,579 sites (0.67 % of the misses)** have no
same-strand cell-barcoded read at all (mouse: 6,628 = 6.58 %); the median class-(e) site carries 71
reads in ±100 bp, and **35,224 sites (12.35 % of the whole denominator)** are class (e) *and*
Kinnex-confirmed at ≥5 UMI *and* carry ≥100 reads. It is an evidence-*generation* failure at the
seeding step — which is why it is an algorithm question ([23](23_algorithm_roadmap.md)), not a
threshold question.

---

## 5. The denominator is inflated — state it precisely or not at all

`results/perf_gap/D1_recall_decomposition/tables/denominator_audit_pbmc_10k_v3.tsv`.

| denominator | sites | R_det (default) | ceiling (any raw peak) |
|---|---:|---:|---:|
| detected-gene atlas, as scored today | 285,136 | 0.1754 | 0.3406 |
| ∩ Kinnex **x3p** ≥5 UMI within 100 bp | 107,260 | **0.4372** | 0.6736 |
| ∩ Kinnex **GEM-X** ≥5 UMI | 114,576 | 0.4134 | 0.6461 |
| ∩ Kinnex **x3p ∪ GEM-X** ≥5 UMI | 127,188 | 0.3796 | 0.6188 |
| ∪-corroborated **and** in a gene with ≥10 k raw reads | 94,023 | 0.3695 | 0.6135 |

**Only 44.6 % of the detected-gene atlas is corroborated by PBMC long reads at ≥5 UMI**
(127,188 / 285,136, the **union** count). Against the corroborated subset the shipped default recalls
**0.4372** — this figure is recall against the **x3p** subset (107,260 sites, 37.6 % of the
denominator), not against the union. **Rows must not be divided across each other** (verifier
correction 4): quote *"44.6 % corroborated (union)"* and *"0.4372 against the x3p-corroborated subset"*
as two separate facts, or quote the union pair (127,188 / 0.3796).

Caveats that travel with this: the Kinnex panel is **donor-mismatched** and site-level
([16](16_trusted_novel_kinnex.md)), and ≥5 UMI has its own sensitivity limit — *"not confirmed"* is not
*"not real"*. The atlas is a union over hundreds of tissues; one PBMC library uses a small subset of it
(the atlas places 20.0 sites in an average detected human gene — 9.7 in mouse — while the default emits
4.0 per gene; `results/perf_gap/D1_recall_decomposition/README.md` §4.1).

---

## 6. The ≥2-molecule default is not the F1 optimum, and that is on purpose

Measured, both datasets (`D1 recovery_ceiling_*.tsv`; `D2_arm_sweep.tsv` group 1):

| dataset | default (≥2 mol) P / R / **F1** | ≥1-molecule arm P / R / **F1** |
|---|---|---|
| PBMC 10k v3 | 0.7062 / 0.1754 / **0.2811** | 0.3520 / 0.2685 / **0.3046** |
| testis mouse 1 | 0.7450 / 0.2048 / **0.3213** | 0.5686 / 0.2806 / **0.3758** |
| testis mouse 2 | 0.7572 / 0.2080 / **0.3264** | 0.5880 / 0.2826 / **0.3818** |

The best PBMC F1 measured anywhere in the 69-arm sweep is **0.3253**
(`PI_union_sup2_or_1hex12`), also above the default. **The default is therefore chosen on
reliability-first grounds — a call set a biologist can act on without re-validating each site — and not
on F1.** The paper must say this in so many words; otherwise a reader reasonably infers that 0.2811 is
the best the tool can do, which is false. (Corollary: the standing "F1 essentially tied with polyApipe"
line understates the tool twice over — it compares at unequal N *and* at a point deliberately off the
F1 optimum.)

---

## 7. False-positive forensics — the two tools fail in different ways

`results/perf_gap/D3_head_to_head/tsv/05_fp_forensics_pa_default_pbmc.tsv` and
`06_fp_forensics_polyapipe_pbmc.tsv` (PBMC, hard-FP definition of §0);
`05b` / `06b` for mouse 1 (atlas-only definition, **not** comparable with the human numbers).

| | PeakATail default | polyApipe |
|---|---:|---:|
| calls | 46,524 | 120,916 |
| **hard FPs** | **7,027 = 15.10 %** | **59,856 = 49.50 %** |
| FPs that are minimum-support calls | 70.9 % (support = 2 molecules) | **81.6 % (peakdepth = 1)** |
| FP rate at minimum support | 29.8 % | **62.6 %** |
| FP rate at ≥10 support | 1.7 % | 0.3 % |
| calls in 3′UTR | **51.8 %** | 18.8 % |
| FPs that are intronic | 77.6 % | 69.3 % |
| FPs that are intergenic | 8.7 % | **23.2 %** |
| FPs with no canonical hexamer | 46.8 % | 67.8 % |
| calls flagged internal-priming by the caller's own rule | **0 (by construction)** | **8,083 (6.7 %); 70.7 % of those are FPs** |
| Kinnex x3p ≥5 UMI concordance ≤25 bp | **0.7647** | 0.3895 |

**Composition, in one sentence each:**

* **PeakATail default — a thin error.** Everything internal-priming-like is already gone (re-testing the
  survivors with two stricter rules flags 36 and 31 calls of 46,524); what remains is **low-support
  intronic peaks without a poly(A) signal** — the single cell `intronic ∧ no hexamer` holds 35.5 % of
  all its FPs at a 45.5 % FP rate within that cell, while 3′UTR calls (51.8 % of output) carry a 1.8 %
  FP rate.
* **polyApipe — a bulk error.** It emits **78,016 single-read peaks (64.5 % of its output)** of which
  62.6 % are hard FPs; it applies no effective internal-priming control (6.7 % of its calls are
  IP-flagged by our rule despite its own `misprime` exclusion, and 70.7 % of those are FPs); nearly a
  quarter of its FPs are intergenic; and 82.2 % of the atlas sites it recovers and we do not are backed
  by a **single poly(A) read**, with 44.4 % of them intronic and 9.1 % IP-flagged by our rule.

Mouse 1, atlas-only definition and therefore a *weaker, more pessimistic* measure on both sides:
PeakATail 6,695 / 26,255 = **25.50 %**, polyApipe 53,670 / 89,527 = **59.95 %**. The mouse FP mass sits
in a different place (`other exon` 53.95 % FP rate; introns only 13.8 % of calls), which is why the
PBMC-optimal post-filter does not transfer — see [22](22_performance_roadmap.md) §4.3.

**Set algebra behind it** (`tsv/02_set_algebra_summary.txt`): of the 17,332 PBMC atlas sites polyApipe
recovers and the default does not, **81.5 % already have a PeakATail tier-1 call with exactly one clip
molecule** and only **5.4 % have no PeakATail evidence of any kind** (mouse 1: 67.6 % / 14.1 %).
polyApipe's exclusive recall is a **scoring** difference, not a detection difference.

---

## 8. How to state this in the paper — the exact sentences

Use these verbatim, or an edit that preserves every qualifier.

**8.1 The replacement one-liner (Results, and the standing summary in [01](01_outline_and_journals.md)):**

> At its pre-registered operating point PeakATail is the most atlas-concordant de novo call set in the
> panel on both datasets (P@100 0.7062 PBMC; 0.7450 / 0.7572 mice). Its recall at that point
> (R_det 0.1754 / 0.2048 / 0.2080) is below polyApipe's, but that comparison sets 46,524 of our calls
> against 120,916 of polyApipe's: at matched call count PeakATail leads polyApipe on recall at every call
> budget tested and on precision at every budget above 20,000 calls, on both datasets — at polyApipe's
> own N of 120,916, 0.4036 / 0.2336 against 0.3800 / 0.1988 —
> and on both mice the ≥1-molecule arm exceeds polyApipe on precision, recall and F1 simultaneously
> (0.5686 / 0.2806 / 0.3758 and 0.5880 / 0.2826 / 0.3818 against 0.4005 / 0.2499 / 0.3078 and
> 0.4119 / 0.2508 / 0.3118). The default is a deliberately conservative point on a curve that dominates
> polyApipe's operating point, not the tool's frontier.

**8.2 Whenever a single-point head-to-head is shown (figure captions, tables):**

> Call-set sizes in this panel span 35,759 to 787,138 (22×) and both precision and recall move with call
> budget, so a single-point comparison measures the operating point, not the tool. The matched-call-count
> comparison is in [25](25_competitive_position.md) §2.

**8.3 On the default's F1 (Results or Discussion):**

> The ≥2-molecule default is not the F1 optimum on either dataset — the ≥1-molecule arm reaches
> F1_det 0.3046 (PBMC) and 0.3758 / 0.3818 (mice) against the default's 0.2811 and 0.3213 / 0.3264. The
> default is chosen for reliability, not for F1: it is the arm whose calls a biologist can act on
> without re-validating each site. We report both arms as labelled operating points.

**8.4 On recall, in the Discussion:**

> 79.96 % of the detected-gene atlas sites the PBMC default misses have no peak of any kind within
> 100 bp (84.53 % in mouse), so at most 20.04 % / 15.47 % of the gap is reachable by any threshold,
> filter or tier policy, and the hard ceiling for anything downstream of peak calling is R_det 0.3406 /
> 0.3278. The residual is an evidence-generation limit at the seeding step, not a thresholding choice.

**8.5 On the denominator, in the Discussion:**

> The recall denominator is itself an upper bound on what any tool can find in one library: only 44.6 %
> of the detected-gene atlas is corroborated by matched-tissue Kinnex long reads at ≥5 UMI — from a
> **different donor** than the short-read library, so "not corroborated" is not "not real". Against the
> long-read-corroborated subset the same default recalls 0.4372 (x3p subset), not 0.1754. Absolute
> recall against a pan-tissue atlas should not be read as sensitivity.

**8.6 On error composition (Results):**

> 15.10 % of the default's PBMC calls are false by both truths (no atlas site within 100 bp and no
> long-read 3′ end within 25 bp), against 49.50 % for polyApipe. The residual error is not internal
> priming — it is low-support intronic peaks without a poly(A) signal, 35.5 % of all false positives
> falling in that single class.

**8.7 Language to avoid** (each of these is now known to be wrong or unqualified):

| do not write | why |
|---|---|
| "recall at or below polyApipe (−12 % / −18 %)" | property of one operating point, at 2.6× fewer calls |
| "the F1 lead is not a meaningful margin" *(alone)* | true at that point, but the point is off the F1 optimum by choice; say why |
| "PeakATail wins at matched N" *(unqualified)* | true **among de novo tools**; false against catalog-based scUTRquant |
| "PeakATail wins on precision at every matched N" | **false** below N ≈ 15,000, where polyApipe is higher (§2.3.3); and not established for 20,320 ≤ N ≤ 35,759 (§2.3.1). Write "at every N from 20,320 upward" |
| "PeakATail's recall lead only holds above 20,000 calls" | over-correction: the **recall** lead holds at every matched N on both datasets; only the precision lead is bounded below |
| "PeakATail at scAPAtrap's N" | not reachable (787,138 > our maximum 333,920) |
| any matched-N claim about SCAPTURE | it has no non-circular ranking column |

---

## 9. Honesty rules — the standing ones that still hold, plus the new one

Unchanged and still binding: precision is **atlas-agreement** precision, never ground-truth precision;
no single cross-dataset accuracy number; F1 never quoted without naming the truth set; the non-IP
≥2-molecule file is not the default; no "trusted novel PAS" claim ([16](16_trusted_novel_kinnex.md));
reproducibility leadership is claimed against polyApipe only.

**NEW RULE (2026-08-21): never quote a single-point head-to-head without the matched-N context.**
Any sentence, table row, caption or slide that compares a PeakATail number with a competitor number
must either (i) be at matched call count, or (ii) carry the call counts of both tools and a pointer to
§2 of this document. A head-to-head table row with two Ns that differ by more than ~10 % and no
matched-N companion is a defect, not a result.

**Two corollary rules.** (a) State whether a comparison is restricted to **de novo** tools every time
the restriction is load-bearing. (b) State the tie-break status of any matched-N claim inside
20,320 ≤ N ≤ 35,759 (§2.3.1), and never state the precision half without its lower bound N ≈ 20,320
(§2.3.3).

**Documents corrected on 2026-08-21:** [15](15_final_gate.md) §3, [19](19_final_gate_v2.md) §4 and new
§6, [01](01_outline_and_journals.md) (§A one-liner, R3 key claim, two honesty bullets, one new honesty
bullet), [05](05_figure_index.md) Fig 2 headline, [16](16_trusted_novel_kinnex.md) (IP-rule note),
`github/issue10_stage1d_ipfilter_memory.md` §4, `figures/final_benchmark.caption.md` (now `fig2_accuracy.caption.md`) and its generator.

**Still carrying the retired framing — OUTSTANDING, not fixed here:**
[21](21_paper_architecture.md) in several places — lines 42, 45, 184 and §8 "Gap 0" (line 488) — states
the recall deficit without the matched-N context (its lines 569 and 661 already say the gap "is the price
of the pre-registered ≥2-molecule threshold, a chosen operating point", which is the right framing and
should be propagated). 21 is owned by another work stream and the correction is already on the checklists
in [22](22_performance_roadmap.md) §6 item 2,
`github/issue19_operating_curve_and_ranking.md` and `github/issue22_perfgap_negative_results.md`. Until
21 is corrected, the manuscript set is internally inconsistent on this point; quote **this** file.
[09](09_headtohead_results.md) line 82 ("whose F1 lead rides on depth-1 singletons") is about the
*shipped pre-fix* caller, a different and still-true claim, and is deliberately left alone.

---

## 10. Where every number comes from

| block | primary artefact | verification |
|---|---|---|
| §2 matched-N grids, tie-break brackets | `results/perf_gap/D3_head_to_head/tsv/01*.tsv` | `results/perf_gap/VERIFY/README.md` (grid reproduced exactly; correction 3 adds the competitor-side tie bracket) |
| §3 mouse ≥1-mol arms, PBMC operating curve | `results/perf_gap/D2_rule_sweep/D2_arm_sweep.tsv` (`M1_sup_ge1`, `M2_sup_ge1`, group 8 = the measured operating curve) | verifier re-derived 11 D2 arms from the genome; "D2's mouse *already wins* headline on both mice" reproduced |
| §4 recall budget, ceiling, class (e) | `results/perf_gap/D1_recall_decomposition/tables/decomposition_*.tsv`, `recovery_ceiling_*.tsv`, `read_evidence_by_class_*.tsv` | whole PBMC decomposition and ceiling reproduced; `results/perf_gap/CHECK22/README.md` for the 35,224 hard core |
| §5 denominator audit | `.../tables/denominator_audit_pbmc_10k_v3.tsv` | reproduced; verifier correction 4 fixes how the two columns may be quoted |
| §6 F1 arms | `recovery_ceiling_*.tsv`, `D2_arm_sweep.tsv` group 1 | reproduced |
| §7 FP forensics, set algebra | `.../D3_head_to_head/tsv/02,05,05b,06,06b,10_*.tsv` | reproduced, incl. a 200-site / 15-column spot check of the per-call master table |
| §7 Kinnex concordance 0.7647 / 0.3895 | `results/perf_gap/D2_rule_sweep/kinnex/kinnex_competitors.tsv` | independently recomputed in `results/perf_gap/CHECK22/README.md` §2 |
| gate numbers (0.7062 / 0.7450 / 0.7572) | [19](19_final_gate_v2.md) §1 | verifier verdict FIXED; `results/benchmark_tools/final_v2_verify/VERIFIED_v2.md` |

**Two arithmetic traps recorded so they are not repeated** ([22](22_performance_roadmap.md) §6):
the standing "both tiers 333,920 / 0.1932 / 0.2979 / 0.2415" line is a mash-up of v2 and pre-v2 rows —
the measured v2 values are **IP arm 333,920 / 0.2044 / 0.3043 / 0.2445** and **no-IP arm
402,765 / 0.1932 / 0.3363 / 0.2455**; and the "~17×" in `results/perf_gap/VERIFY/README.md` is a slip
in that document's own prose — 68,855 / 12,186 = **5.65×** (§ [16](16_trusted_novel_kinnex.md) note).

---

**Known stale artefact.** `manuscript/figures/final_benchmark.png` / `.pdf` were rendered before this
correction and still carry the retired in-figure clause. The generator
(`scripts/manuscript_figures/final_benchmark.py`) and the sidecar caption have been corrected; the
figure must be re-rendered before any version of it is circulated. (Now resolved: the figure was
re-rendered from the corrected generator — commit `606968b`'s re-render — and renamed to
`fig2_accuracy` at the 2026-09-02 rename pass.)

*Everything in this file is PROVISIONAL as a manuscript claim until a verifier passes on this document
itself; every number in it is already verified at source.*

# 22 — Performance roadmap: what to fix next, ranked by measured gain

**Status: every number below is either (i) confirmed by the independent adversarial verifier of the
perf-gap study (`results/perf_gap/VERIFY/`, verdict **FIXED**, 69 checked quantities in
`VERIFY/out/VERIFICATION_RESULTS.tsv`), (ii) confirmed by the verifier of the parallel
algorithmic-headroom study (`results/algo_headroom/VERIFY/`, verdict FIXED) and labelled as such, or
(iii) recomputed in this document from the primary per-site table, with the command given.
Where the verifier corrected a diagnostic, the corrected value is used and the correction is named.**

Sources: `results/perf_gap/D1_recall_decomposition/`, `D2_rule_sweep/`, `D3_head_to_head/`, all three
adversarially verified on 2026-08-21 against the frozen v2 tree `9dfdefb`. Scorer, references and
denominators are the launcher's throughout (`scripts/benchmark_tools/stage2_final_launch.sh`);
recall denominators 285,136 (human) / 126,686 (mouse); "P@100" and "R_det" as defined in
[19](19_final_gate_v2.md).

This document is about **engineering priorities**. The complementary document
[23](23_algorithm_roadmap.md) covers the *algorithmic* headroom study (A1–A5) and owns GitHub issue
drafts 13–18; this one covers the *competitive-gap* study (D1–D3) and owns drafts 19–22. They share
two conclusions and contradict each other nowhere; the overlaps are cross-referenced explicitly in §7.

---

## 1. The competitive position, in three sentences

1. At its pre-registered operating point PeakATail is the most atlas-concordant de novo call set in
   the panel on both datasets — P@100 **0.7062** (PBMC 10k v3), **0.7450** / **0.7572** (mice) — but
   it recalls only R_det **0.1754** / **0.2048** / **0.2080**, which is **below polyApipe's**
   0.1988 / 0.2499 / 0.2508 ([19](19_final_gate_v2.md) §1; [15](15_final_gate.md) §3).
2. That recall deficit is, however, an artefact of comparing **46,524 calls with 120,916**: at
   polyApipe's own call count — where polyApipe is its own full set, so no tie-break bracket applies —
   PeakATail scores **0.4036 / 0.2336** against polyApipe's **0.3800 / 0.1988**, winning both axes,
   and on **both** mice the already-published ≥1-molecule arm beats polyApipe on precision, recall and
   F1 *simultaneously* (**0.5686 / 0.2806 / 0.3758** and **0.5880 / 0.2826 / 0.3818** against
   **0.4005 / 0.2499 / 0.3078** and **0.4119 / 0.2508 / 0.3118**).
3. What is **not** an artefact is the ceiling: **79.96 %** of the 235,111 detected-gene atlas sites
   the PBMC default misses have **no peak of any kind** within 100 bp, so no threshold, filter or
   re-scoring change can take R_det past **0.3406** (PBMC) / **0.3278** (mouse) — the roadmap's
   published bar of R ≥ 0.60 ([07](07_curated_benchmark_report.md), "Roadmap publication bars") is unreachable without
   changing peak calling itself.

> **The one-liner that must be retired.** "Recall at or below polyApipe (−12 % PBMC, −18 % mouse)"
> ([01](01_outline_and_journals.md) §A/§Key claim) is true of the ≥2-molecule *default* and false of
> the *tool*. The honest replacement is in §4.1.

---

## 2. The recall budget: how much is recoverable, and from where

### 2.1 The decomposition (PBMC 10k v3, 235,111 missed of 285,136 = 82.46 %)

Priority-ordered, mutually exclusive, exhaustive; **the partition is exact** — the verifier
re-derived every class from the raw 652,665-peak universe with its own code and reproduced all eight
counts, with the class call counts summing to 651,196 and the site counts to 285,136.

| class | what it is | missed sites | % of missed | calls in class | class P@100 | Kinnex-union confirmed |
|---|---|---:|---:|---:|---:|---:|
| — | **recovered by the default** | 50,025 | — | 46,524 | **0.7062** | 0.9651 |
| (a) | tier-1 IP-pass call exists, **1 clip molecule only** | 26,526 | 11.28 % | 121,041 | 0.2158 | 0.7359 |
| (b1) | tier-1 ≥2 molecules, **removed by the IP filter** | 3,814 | 1.62 % | 15,586 | 0.2461 | 0.6841 |
| (b2) | tier-1 1 molecule, removed by the IP filter | 4,310 | 1.83 % | 39,804 | 0.1205 | 0.4749 |
| (c1) | only a **tier-2** (coverage-only) IP-pass call | 10,198 | 4.34 % | 166,355 | 0.0558 | 0.4873 |
| (c2) | only a tier-2 call, IP-flagged | 1,029 | 0.44 % | 13,455 | 0.0697 | 0.3537 |
| (d) | only a **raw pre-call peak** | 1,228 | 0.52 % | 248,431 | 0.0782 | 0.7394 |
| (e) | **no peak of any kind within 100 bp** | 188,006 | **79.96 %** | — | — | 0.2579 |

Mouse 1 (100,742 missed of 126,686 = 79.52 %) has the same shape: (a) 9,609 (9.54 %), (b1) 1,269,
(b2) 1,067, (c1) 2,877, (c2) 161, (d) 606, **(e) 85,153 (84.53 %)**.

**Read it as a budget.** Everything a scoring, thresholding or filtering change can ever reach is
classes (a)+(b)+(c)+(d) = **20.04 %** of the PBMC misses and **15.47 %** of the mouse misses.
The remaining four fifths is not reachable by any post-peak-calling rule.

### 2.2 The measured ceiling (not a model — these are scored arms)

| arm | n calls | P@100 | R_det@100 | F1_det |
|---|---:|---:|---:|---:|
| **A0** precision default (tier-1 ∧ IP-pass ∧ ≥2 mol) | 46,524 | 0.7062 | 0.1754 | 0.2811 |
| **A1** + class (a) = drop the ≥2-molecule rule | 167,565 | 0.3520 | 0.2685 | **0.3046** |
| **A2** + class (b) = drop the IP filter too | 222,955 | 0.3032 | 0.2970 | 0.3001 |
| **A3** + class (c) = trust tier-2 = every called PAS | 402,765 | 0.1932 | 0.3363 | 0.2455 |
| **A4** + class (d) = every raw peak the caller ever made | 651,196 | 0.1493 | **0.3406** | 0.2077 |
| **ORACLE** perfect recovery of (a)–(d), zero new FPs | 93,629 | 0.8540 | **0.3406** | 0.4855 |

Mouse: 0.7450/0.2048 → 0.5686/0.2806 → 0.4992/0.2991 → 0.4093/0.3231 → 0.2799/**0.3278**;
oracle 41,844 / 0.8400 / 0.3278 / 0.4716.

**The hard ceiling is R_det 0.3406 (PBMC) / 0.3278 (mouse).** Keeping 14× the default's output at
precision 0.149 still leaves 65.94 % of the detected-gene atlas untouched, and even a perfect oracle
recovery reaches only 0.3406 at F1 0.4855.

### 2.3 The marginal cost of recall (PBMC)

| block | calls added | block P@100 | atlas sites gained | sites per call | calls per +0.01 R_det |
|---|---:|---:|---:|---:|---:|
| the default itself | 46,524 | 0.7062 | 50,025 | **1.075** | 2,652 |
| (a) singletons | 121,041 | 0.2158 | 26,526 | 0.219 | 13,011 |
| (b) IP-removed | 55,390 | 0.1559 | 8,124 | 0.147 | 19,441 |
| (c) tier-2 | 179,810 | 0.0568 | 11,227 | 0.062 | 45,667 |
| (d) raw uncalled | 248,431 | 0.0782 | 1,228 | 0.005 | 576,845 |

The default is already at its own call-budget ceiling; each subsequent block is 5–200× less
efficient. **Any engineering effort aimed at the scoring/filtering layer is bidding for at most
+0.165 R_det at a precision cost from 0.71 to 0.15.**

### 2.4 Class (e) is a *seeding* failure, not a reach failure

Measured from the BAMs (`-F 2820 -d CB`, same-strand; the verifier re-measured chromosome 20
independently and matched D1's per-site column on **6,310/6,310 sites, zero mismatches**):

| | reads = 0 | ≥1 | ≥5 | ≥20 | ≥100 | median |
|---|---:|---:|---:|---:|---:|---:|
| class (e), PBMC (n 188,006) | 1,579 | 186,427 | 173,862 | 144,002 | 77,321 | **71** |
| recovered by the default (n 50,025) | 0 | 50,025 | 50,023 | 49,945 | 46,280 | 638 |

Only **0.67 %** of the misses are "nothing was sequenced there". And R_det by gene expression
saturates far below what long reads confirm is real:

| gene raw reads | atlas sites | R_det | Kinnex-confirmed fraction |
|---|---:|---:|---:|
| <100 | 16,179 | 0.0068 | 0.0355 |
| 100–1k | 34,280 | 0.0418 | 0.1145 |
| 1k–10k | 75,097 | 0.1675 | 0.3817 |
| 10k–100k | 138,228 | 0.2232 | 0.5785 |
| ≥100k | 21,352 | **0.2368** | **0.6584** |

**The sharpest single statement of the deficit:** 35,224 sites — **12.35 % of the whole denominator**
— are simultaneously long-read-confirmed at ≥5 UMI, covered by ≥100 same-strand cell-barcoded reads
in ±100 bp, and produce **no peak of either tier**.
*(Recomputed here, not taken from D1: `zcat results/perf_gap/D1_recall_decomposition/tables/per_site_classification_pbmc_10k_v3.tsv.gz | awk -F'\t' 'NR>1 && $2=="e_nopeak" && $5!="NA" && $5+0<=100 && $10+0>=100' | wc -l` → 35,224. The verifier separately confirmed the class-(e) Kinnex-confirmed total of 48,489 post-contig-filter.)*

### 2.5 The denominator is half the story, and it is free to fix

Only **127,188 of 285,136** detected-gene atlas sites (44.6 %) are corroborated by PBMC long reads at
≥5 UMI (union of x3p and GEM-X; x3p alone 107,260, GEM-X alone 114,576). Against the corroborated
subsets the shipped default's recall is:

| denominator | n sites | recovered | R | ceiling (any raw peak) |
|---|---:|---:|---:|---:|
| detected-gene atlas (as published) | 285,136 | 50,025 | **0.1754** | 0.3406 |
| ∧ Kinnex x3p ≥5 UMI | 107,260 | 46,892 | **0.4372** | 0.6736 |
| ∧ Kinnex GEM-X ≥5 UMI | 114,576 | 47,369 | 0.4134 | 0.6461 |
| ∧ Kinnex either | 127,188 | 48,281 | 0.3796 | 0.6188 |

PolyASite 2.0 is a union over hundreds of tissues and conditions and puts ~20 sites in an average
human detected gene; one PBMC library uses a small subset of that. Reporting a long-read-restricted
recall *alongside* the atlas recall costs a second denominator and no code.
**Verifier note (do not misread this table):** "44.6 % corroborated" is the *union* count (127,188)
while 0.4372 is the recall against the *x3p* count (107,260). Each row is internally consistent; the
numbers must not be divided across rows.

---

## 3. Ranked candidate fixes

Every ΔP/ΔR below is **measured** on a scored call set unless the row says EXTRAPOLATED or UNMEASURED.
"Pre-reg?" means: would this have to be pre-registered under
[24](24_prime_preregistration.md) §3.1/§3.4 before it may change the paper's default?

| # | fix | what it changes | measured effect (PBMC → / mouse 1 →) | difficulty | risk to existing results | pre-reg? |
|---|---|---|---|---|---|---|
| **1** | **Report an operating curve, not one point** | nothing in the caller; the head-to-head figure and one paragraph of Methods | **no change to any published number.** Adds: PBMC at N = 120,916 **0.4036 / 0.2336** vs polyApipe 0.3800 / 0.1988; mouse ≥1-mol **0.5686 / 0.2806 / 0.3758** and **0.5880 / 0.2826 / 0.3818** vs 0.4005 / 0.2499 / 0.3078 and 0.4119 / 0.2508 / 0.3118 | **zero code** | **none** to numbers; presentational only (the ranked-prefix curve is a new descriptive analysis and must be labelled post-hoc) | **No** — the ≥1-molecule arm is already pre-registered as the sensitivity arm ([13](13_reliability_positioning.md) §1) |
| **2** | **3′UTR promotion of single-molecule tier-1 calls** — `≥2 mol OR (1 mol ∧ annotated 3′UTR ∧ hexamer)` | promotes a named subset of class (a) | conservative (∧ any-12 hexamer): 46,524 → **54,220** calls, P 0.7062 → **0.7194**, R_det 0.1754 → **0.2015**, F1 0.2811 → **0.3148**, hard-FP rate 0.1510 → **0.1462**, Kinnex ≤25 bp 0.7647 → 0.7361. Mouse 1: n **30,974**, **0.7549 / 0.2402 / 0.3644**. Aggressive (3′UTR only): PBMC 61,232 / **0.7034 / 0.2150 / 0.3294**, mouse 34,865 / **0.7359 / 0.2545 / 0.3782** | **trivial** — the GTF is already loaded | **medium**: it makes the *default* annotation-dependent, which weakens the "de novo" claim; and it is post-hoc on gate data | **YES** |
| **3** | **Precision post-filter: drop non-3′UTR calls with no canonical hexamer** | removes a named FP class from the default | PBMC 46,524 → **38,824**, P 0.7062 → **0.7993**, R_det 0.1754 → **0.1683**, F1 0.2811 → 0.2781, hard FPs −44.3 % for −11.6 % evidenced calls. Mouse 1 → **22,623**, P **0.8033**, R_det 0.1942, F1 0.3127, FPs −33.5 % for −7.1 % | **trivial** | **medium**: same annotation dependence as #2, and **F1 falls on both datasets** — this is a reliability trade, not a free win | **YES** |
| **4** | **Cohort / cross-replicate borrowing** (`--min-samples`-style promotion) | promotes class (a) sites a second sample independently calls | mouse 1: replicated singletons P **0.5848** vs size-matched random **0.3882**; union arm 33,621 / **0.7099 / 0.2394 / 0.3580**. Mouse 2 reproduces: 33,722 / **0.7144 / 0.2404 / 0.3597**. **PBMC value UNMEASURED — the flagship human dataset is one library** | moderate (needs cohort mode; a Stage-3 driver exists) | low on published numbers (cohort-only); **unmeasured** effect on the Stage-3 unified PAS space | **YES** |
| **5** | **Hexamer union as a labelled sensitivity arm** (`≥2 OR 1-mol strong hexamer`) | a presentational arm only | PBMC 71,074 / **0.6136 / 0.2167 / 0.3202**; mouse 1 34,139 / 0.7034 / 0.2410 / 0.3590. **But**: against a size-matched *random* singleton top-up it gains **+0.078 atlas precision and only +0.007 long-read concordance** (0.6090 vs 0.6021) | trivial | **high — the gain is largely circular** (see §5.3) | **YES**, and it should not become a default |
| **6** | `--ip-rule kinnex` (issue **#95** addendum) | replaces or supplements the IP window rule | **DO NOT BUILD.** On top of the shipped filter it removes **31 of 46,524** calls (P 0.7062 → 0.7065). As a *replacement* on the no-IP arm the shipped rule lifts P 0.5907 → 0.7062 while the Kinnex rule reaches only **0.5944**; genome-wide the shipped rule flags **68,855 of 402,860 peaks (17.09 %)** against the Kinnex rule's **12,186 (3.02 %)** | — | — | closes #95's addendum with a **negative** |
| **7** | Cleavage-offset **coordinate correction** | moves the reported cleavage point | **DO NOT MOVE THE COORDINATE.** −1/−2/−3 bp shifts move PBMC P@100 0.7062 → 0.7062 / 0.7062 / 0.7063 and P@25 0.6389 → 0.6386 / 0.6378 / 0.6363, while Kinnex internal-priming decoy proximity rises 0.1300 → 0.1360 / 0.1419 / 0.1477 | — | — | see §7 — an *additional column* is a separate, accepted proposal (issue **16** draft) |
| **8** | Tier-2 rescue at any depth | promotes class (c) | **DO NOT BUILD.** Class (c1) P@100 **0.0558**, (c2) 0.0697; best gated variant costs −0.163 P@100 and −0.178 Kinnex concordance for +0.005 R_det | — | — | — |
| **9** | `-F 3844 ≥ 2` strict-filter support | changes which clip molecules count | **DO NOT ADOPT.** Dominated on PBMC (ΔP −0.035 against the plain support knob at equal recall); and PBMC singletons whose only molecule **fails** `-F 3844` have **higher** Kinnex concordance (0.3863) than those that pass (0.2858) | — | — | — |

**Unmeasured combinations, flagged so nobody assumes them.** #2 and #3 act on the same annotation
axis in opposite directions (promote 3′UTR singletons; drop non-3′UTR hexamer-fail doubletons) and
have **never been scored together**. Their combined effect is **UNMEASURED and must not be estimated
by addition** — the algorithmic-headroom study's own union analysis found 34 % double-counting when
four recall ideas were summed ([23](23_algorithm_roadmap.md); issue18 draft §7).

---

## 4. What to do before submission, and what to defer

### 4.1 DO — #1, publish an operating curve (zero code, highest value)

This is the single highest-value action in the document and it costs nothing but a figure and a
paragraph. Three things follow from it:

* **A matched-N figure** (P and R_det versus N, one curve per tool) replaces the single-row
  head-to-head table. The single-row table is what created the false impression in the first place.
* **The claim upgrades** from "best precision among de novo tools, at polyApipe-or-lower recall" to
  **"best precision *and* best recall among de novo tools at polyApipe's own call budget, on both
  datasets"** — verified at N = 120,916 (0.4036 / 0.2336 vs 0.3800 / 0.1988) and, without any
  ranking at all, on both mice from the already-published ≥1-molecule arm.
* **Two honesty obligations come with it.** (i) The precision half of the matched-N claim is **not**
  robust to the competitor's tie-break in the band 20,320 ≤ N ≤ 35,759 — under an oracle tie-break
  polyApipe reaches 0.9193 / 0.8611 / 0.7620 at N = 20,320 / 29,814 / 35,759 against our 0.9074 /
  0.8358 / 0.7579 (our own optimistic bracket at 35,759 is 0.7917; the intervals overlap). It **is**
  unambiguous at N = 46,524 (polyApipe best case 0.6636 vs 0.7062) and at N = 120,916, where
  polyApipe is its own full set with no ties at all. The recall half is robust everywhere.
  (ii) The shipped ≥2-molecule default is **not the F1 optimum** on either dataset (A1 gives 0.3046
  vs 0.2811 on PBMC and 0.3758 vs 0.3213 on mouse). Say explicitly that the default is chosen on
  reliability-first grounds, not on F1, rather than leave a reader to infer 0.2811 is the tool's best.

Draft: **`manuscript/github/issue19_operating_curve_and_ranking.md`**.

### 4.2 DO — #2, 3′UTR promotion, but pre-register it first and ship it as a flag

This is the best *measured* buy in the study, and it is the **only** singleton discriminator that
moves the atlas-independent truth by more than ×1.1:

| singleton stratum | n | atlas P@100 | Kinnex x3p ≤25 bp |
|---|---:|---:|---:|
| all singletons | 121,041 | 0.2158 | 0.2952 |
| **in an annotated 3′UTR** | 14,708 | **0.6947** (×3.22) | **0.5143** (×1.74) |
| 3′UTR ∧ any-12 hexamer | 7,696 | 0.7994 | 0.5635 |
| 3′UTR ∧ strong hexamer | 4,164 | 0.8489 | 0.5865 |
| strong hexamer, any location | 24,550 | 0.4380 (×2.03) | 0.3139 (**×1.06**) |
| not in a 3′UTR | 106,333 | 0.1496 | 0.2650 |

This table is the verifier's own — it was measured with independently parsed Ensembl 99/102
references, and it **corrected** D2's claim that "the singleton class is irreducible with the
available features" (true for sequence and depth features; false for genomic annotation, which D2's
69-arm sweep never tested).

**Do the conservative variant** (`3′UTR ∧ any-12 hexamer`): PBMC P **up** 1.3 pp *and* R_det up
2.6 pp *and* hard-FP rate **down** (0.1510 → 0.1462) for only −2.9 pp of long-read concordance;
mouse P up 1.0 pp and R_det up 3.5 pp. It is the only candidate in the study that improves precision
and recall at once on both species.

**Three conditions, non-negotiable.** (a) It is **post-hoc on gate data** and must be pre-registered
and re-run on all three datasets before it can enter the paper as anything but a labelled
exploratory arm ([24](24_prime_preregistration.md) §3.4). (b) It must ship **behind a flag with the
default unchanged** until that is done. (c) The **de-novo claim** must be protected: a promotion rule
that reads the GTF makes the *default* annotation-dependent. The residual caveat the verifier itself
flagged — annotated 3′UTRs derive from transcript ends, so a weak prior is shared with long-read
3′ ends — is far weaker than the motif/PolyASite one, but it is not zero.

Draft: **`manuscript/github/issue20_threeutr_singleton_promotion.md`**.

### 4.3 CONSIDER — #3, the transferable precision post-filter, as a `--strict` mode only

`drop non-3′UTR ∧ no canonical hexamer` removes 44.3 % of PBMC hard FPs for 11.6 % of evidenced
calls (P → 0.7993) and 33.5 % of mouse FPs for 7.1 % (P → 0.8033). **It lowers F1 on both datasets**
(0.2811 → 0.2781; 0.3213 → 0.3127), so it is a reliability trade, not an improvement, and it should
be offered rather than imposed.

The PBMC-tuned variant (`intronic ∧ no hexamer`, P → 0.7801, F1 0.2825 — one of only two of 17
candidates that raise precision and F1 at once) **must not ship**: it removes only 8.45 % of mouse
FPs, because the mouse-testis FP mass sits in "other exon" (FP rate 0.5395), not introns. This is a
concrete, measured tissue-dependence caveat the paper should carry.

Draft: **`manuscript/github/issue21_annotation_precision_filter.md`**.

### 4.4 DEFER — #4, cohort borrowing

The mouse evidence is good and replicates in both directions (singleton P 0.3942 → **0.5848**
against a size-matched random draw's 0.3882; union arms 0.7099 / 0.2394 / 0.3580 and 0.7144 /
0.2404 / 0.3597), but **its human value is entirely unmeasured** because the flagship human dataset is
a single library, and there is **no long-read truth for mouse** to check whether cross-sample
corroboration is picking up replicated internal priming — the genome is the same in both animals.

**Do not file a new ticket for this.** A draft already exists from the algorithmic-headroom
workstream: `manuscript/github/issue17_cross_sample_corroboration.md` (measured on
`results/algo_headroom/A5_cross_sample/`, **not yet verifier-passed**). The perf-gap sweep's mouse
cross-replicate arms are an **independent, verifier-confirmed corroboration** of that draft's
direction, on a different rule (25 bp window, tier-1 corroborator) and a different pipeline:
append them to issue17 as an addendum rather than opening a second ticket.

### 4.5 DEFER — #5, and never as a default

The hexamer union is trivial and beats polyApipe on every metric, but §5.3 explains why it must be
budgeted as a presentational arm and not sold as a quality fix.

### 4.6 RECORD — the negatives, so nobody rebuilds them

`--ip-rule kinnex`, tier-2 rescue, cleavage-offset coordinate correction and `-F 3844 ≥ 2` are all
measured negatives. They should be written down once, with numbers, so they do not come back.

Draft: **`manuscript/github/issue22_perfgap_negative_results.md`**.

---

## 5. What engineering cannot fix, and how the paper should say it

### 5.1 The chemistry budget

The binding constraint is upstream of every rule in §3: a 3′ 10x library exposes a poly(A) soft clip
on a small minority of reads, and clip-seeded calling can only seed where one exists.

* **The clip rate is ~0.57 %, not 1.15 %.** Genome-wide, with the tool's own `read_check` +
  `clip_site` criteria: **0.5730 %** of accepted CB reads (3,195,067 / 557,564,408), range 0.3669 %
  (chr21) to 0.8179 % (chr19). The algorithmic-headroom verifier re-implemented the criteria from
  scratch in gawk over `samtools view` and reproduced the chr21 figure **to six decimals** from a
  separate implementation (`results/algo_headroom/VERIFY/VERIFICATION_REPORT.md` §5).
  **This halves the number currently stated in [10](10_caller_fix_plan.md) and in
  `ema/countmatrix/polya.py`'s docstring — the sensitivity budget is tighter than the manuscript
  says, not looser.** (Correcting it is issue13's checklist; see §7.)
* **96.1–96.3 % of accepted reads carry no 3′-side soft clip at all**, so their reported 3′ end is a
  fixed offset from the alignment start — which is why an end-position/pileup channel carries no
  information beyond local depth after depth matching (issue18 draft §1; verified).

### 5.2 The consequence, measured

* **65.94 %** of the PBMC detected-gene atlas (**67.22 %** mouse) produces no peak of any kind, at
  loci carrying a **median of 71** same-strand cell-barcoded reads.
* **35,224 sites (12.35 % of the denominator)** are long-read-confirmed, ≥100-read-covered, and
  peakless.
* The ceiling for everything downstream of peak calling is **R_det 0.3406 / 0.3278**; the oracle
  reaches F1 **0.4855 / 0.4716**.

### 5.3 The methodological trap that outlives this study

**The precision reference is motif-selected.** Hexamer-pass rates of the *truth sets themselves*:
PolyASite 2.0 representative sites **0.6511** (any-12) / 0.4312 (strong); Kinnex ≥5 UMI 3′ ends
**0.2869** / 0.1273; random genic null 0.1894 / 0.0527. So any rule involving the poly(A) signal will
look better against the atlas than it is. Measured on the singleton class: a strong-hexamer gate moves
atlas precision **×2.03** (0.2158 → 0.4380) and long-read concordance **×1.06** (0.2952 → 0.3139).
Against a size-matched *random* singleton top-up the strong-hexamer union arm gains **+0.078 atlas
precision and +0.007 long-read concordance** — an 11:1 ratio (the any-12 union: +0.061 vs +0.013,
4.6:1). **Any recall-recovery proposal must be scored against
the Kinnex ≤25 bp metric before it is scored against PolyASite**, or the paper will ship a gain that
a reviewer with long reads can dismantle. The 3′UTR gate is the counter-example that passes this test
(×3.22 atlas *and* ×1.74 long-read), which is exactly why it is ranked #2 and the hexamer union is #5.

### 5.4 Proposed wording for the Discussion (adopt or edit, but say something equivalent)

> PeakATail's recall against PolyASite 2.0 is bounded by the evidence a 3′ 10x library exposes:
> ~0.57 % of cell-barcoded reads carry a poly(A) soft clip in this dataset, and 65.9 % of
> detected-gene atlas sites produce no candidate peak of any tier — at loci carrying a median of 71
> same-strand cell-barcoded reads, so this is a limit of clip-seeded *seeding*, not of sequencing
> depth or of the assay's reach. Perfect recovery of every site for which the caller does produce a
> candidate would reach R_det 0.341 (human) / 0.328 (mouse); no threshold, filter or re-scoring
> change can exceed that. We also note that PolyASite 2.0 is a union over hundreds of tissues and
> conditions: only 44.6 % of the sites in our detected-gene denominator are corroborated by
> long reads from the same tissue, and against that corroborated subset the same call set recalls
> 0.44 rather than 0.18. We report both denominators rather than choosing the flattering one.

**And retire the R ≥ 0.60 bar.** [07](07_curated_benchmark_report.md) (“Roadmap publication bars”, and its §3 "bars need a decision") already flags it as
arithmetically unreachable; §2.2 now makes that a measured statement rather than a suspicion. The bar
should be restated against a defensible denominator or dropped.

---

## 6. Corrections this study forces on standing text

1. **The standing summary line "both tiers 333,920 / 0.1932 / 0.2979 / 0.2415" is a mash-up**, and
   the verifier traced it: 0.2979 / 0.2415 come from [15](15_final_gate.md) line 36, the **pre-v2**
   row at n 328,302. Measured v2 values: **IP arm, both tiers 333,920 / 0.2044 / 0.3043 / 0.2445**;
   **no-IP arm, both tiers 402,765 / 0.1932 / 0.3363 / 0.2455**.
2. **"Recall at or below polyApipe (−12 % / −18 %)"** ([01](01_outline_and_journals.md)) is true of
   the default and false of the tool — see §1 and §4.1.
3. **The default is not the F1 optimum** on either dataset — [21](21_paper_architecture.md) C6 should
   say so explicitly.
4. **[16](16_trusted_novel_kinnex.md)'s premise is inverted.** It states the caller's IP rule is
   "looser than Kinnex's". At the level of what each rule *removes from a call set*, the shipped rule
   flags **68,855 of 402,860 peaks (17.09 %)** and the Kinnex rule **12,186 (3.02 %)** — the shipped
   rule is **5.7× stricter**. 16's own finding (27.0 % of trusted-novel sites lie near a decoy
   terminus) stands; that residual is simply not fixable by swapping in the Kinnex genomic rule.
5. **Do not quote "~17×" from `results/perf_gap/VERIFY/README.md`.** That line is an arithmetic slip
   in the verifier's own prose: 68,855 / 12,186 = **5.65**, and the verifier's own bottom-line
   summary says 5.7×. The two counts are both confirmed; only the ratio in that one sentence is wrong.
6. **The (b1)/(b2) split is priority-order-dependent** — reversing the within-group order moves 203
   sites (3,814 / 4,310 → 3,611 / 4,513). The class (b) total 8,124 and every headline number are
   invariant. Quote the total, or state the order.
7. **The clip rate is ~0.57 %, not 1.15 %** — §5.1; the fix is on issue13's checklist.

---

## 7. Relationship to the algorithmic-headroom roadmap ([23](23_algorithm_roadmap.md), issues 13–18)

The two studies were run independently, on the same frozen tree, and verified separately. Where they
touch, they agree:

| topic | 23 / issues 13–18 | this document | reconciliation |
|---|---|---|---|
| the ceiling | "of 11,606 chr19+21 missed sites only 21.4 % have any candidate within 100 bp" (issue18 §7) | 20.04 % of PBMC misses genome-wide are in classes (a)–(d) | **same finding, two scales.** Genome-wide is the number to quote |
| cohort borrowing | issue17 draft, mouse A5, **not yet verified** | verifier-confirmed mouse cross-replicate arms, a different rule and pipeline | **independent corroboration.** Append to issue17; do not file twice (§4.4) |
| cleavage offset | issue16: emit an `inferred_cleavage` **column** at −2 bp; do not move the coordinate | moving the coordinate does not help at any window ≥ 25 bp and worsens decoy proximity ×1.7 at W = 10 | **compatible.** Both say: report the shift, do not apply it |
| hexamer / motif features | issue14: the atlas-trained score is circular; the Kinnex-trained one is honest | the hexamer gate is circular (atlas ×2.03, long-read ×1.06) | **same mechanism, two applications** |
| singleton discrimination | issue14: a calibrated score ranks singletons; issue18: "no post-hoc feature separates them" | **genomic annotation does** (3′UTR: atlas ×3.22, long-read ×1.74) | issue14's Phase 0 transfer test and this document's #2 may be **measuring the same thing** — intersect the recovered-site sets before building both |
| clip rate | issue13: 0.5730 %, verified | adopted here as the chemistry budget (§5.1) | **one number, one fix, issue13's checklist** |
| `--ip-filter` | issue13: default it on; it is the largest measured lift in the caller | the IP filter removes 5.7× more than the Kinnex rule and lifts P 0.5907 → 0.7062 | **mutually reinforcing** — #6 here is the negative that keeps #95's addendum from displacing it |

**Issue-number note.** This task specified drafts numbered "13 onward"; 13–18 were claimed by the
algorithmic-headroom workstream between the task being set and this document being written, so these
drafts start at **19**. Nothing is renumbered.

---

## 8. Issue drafts produced (drafts only — nothing filed, nothing committed)

| file | fix | rank | needs pre-registration |
|---|---|---|---|
| `manuscript/github/issue19_operating_curve_and_ranking.md` | report an operating curve; expose the ranking column | #1 | no |
| `manuscript/github/issue20_threeutr_singleton_promotion.md` | 3′UTR promotion of single-molecule tier-1 calls | #2 | **yes** |
| `manuscript/github/issue21_annotation_precision_filter.md` | `--strict` non-3′UTR ∧ no-hexamer post-filter | #3 | **yes** |
| `manuscript/github/issue22_perfgap_negative_results.md` | four measured negatives + four standing-text corrections | #6–#9 | n/a (closes #95's addendum with a negative) |

All four cross-reference the open issues **#94** (Stage 1b/1c poly(A) evidence), **#95** (IP rule),
**#98** (per-isoform degenerate pairs) and **#99** (PAS→gene assignment, clip-rate sampling), and the
sibling drafts 13–18, so the programme reads as one sequence rather than a pile of tickets.

---

## 9. Primary artefacts

* `results/perf_gap/D1_recall_decomposition/` — decomposition, ceiling, marginal cost, read evidence,
  denominator audit, per-site classification of all 285,136 / 126,686 sites, and
  `work/verify/reverify.sh` (one-command re-derivation of every class count from the shipped BEDs).
* `results/perf_gap/D2_rule_sweep/` — 69 scored arms (`D2_arm_sweep.tsv`; pre-verification values
  preserved in `D2_arm_sweep.tsv.preverify`), the Kinnex concordance table, the size-matched random
  controls, and per-site master joins (`annot/master_*.tsv`) from which any further arm is one awk line.
* `results/perf_gap/D3_head_to_head/` — matched-N grids in both directions, set algebra, per-call
  master tables for PeakATail and polyApipe, FP forensics, the 17-candidate filter table, and the
  promotion arms.
* `results/perf_gap/VERIFY/` — the adversarial pass: `out/VERIFICATION_RESULTS.tsv` (69 quantities),
  `out/vscores.tsv`, the independently reconstructed class BEDs, and the four corrections.

# Issue (documentation): four measured negatives from the competitive-gap study, and seven standing statements they correct — including the close-out of #95's `--ip-rule kinnex` addendum

**Labels:** documentation, research, wontfix. Record of measurements made 2026-08-21
(`results/perf_gap/D1_recall_decomposition`, `D2_rule_sweep`, `D3_head_to_head`) and independently
verified (`results/perf_gap/VERIFY/`, verdict **FIXED**, 69 re-derived quantities in
`out/VERIFICATION_RESULTS.tsv`), against the frozen v2 tree `9dfdefb`.

This issue exists so nobody rebuilds a thing we already measured, and so four numbers we currently
have in writing get corrected before they reach a reviewer. Related: **#94**, **#95** (this closes its
addendum), **#98**, **#99**; sibling draft **issue18** (four *different* dead ends, from the
algorithmic-headroom study — the two lists do not overlap), **issue16** (cleavage offset: compatible,
see §3), **issue13** (`--ip-filter`, which §1 here reinforces).
Roadmap context: `manuscript/22_performance_roadmap.md` §4.6 and §6.

**None of these is a bug. Each is a reasonable idea that was measured and does not pay.**

---

## 1. `--ip-rule kinnex` as specified (#95 addendum) — DO NOT BUILD

*Idea:* `manuscript/16_trusted_novel_kinnex.md` found that 27.0 % of trusted-novel sites lie within
25 bp of a Kinnex internal-priming decoy terminus, concluded the caller's IP rule is "looser than
Kinnex's", and proposed `--ip-rule {legacy,kinnex}` (the Kinnex genomic rule being +1..+18 downstream,
>=12 of 18 A).

*Measured — as an addition on top of the shipped filter:* it removes **31 of 46,524** default calls.
P@100 0.7062 -> **0.7065**, Kinnex concordance 0.7647 -> 0.7649, R_det unchanged at 0.1754. A no-op.
(Independently rebuilt from the genome by the verifier: 46,493 / 0.706536.)

*Measured — as a replacement, applied to the no-IP arm (PBMC, >=2 molecules):*

| arm | n | P@100 | R_det | Kinnex t5 <=25 | Kinnex **decoy** <=25 |
|---|---:|---:|---:|---:|---:|
| no IP rule at all | 62,110 | 0.5907 | 0.1911 | 0.6081 | 0.3122 |
| **shipped rule** (-9..+30, >=6 consecutive A or >=70 % A) | **46,524** | **0.7062** | 0.1754 | **0.7647** | **0.1300** |
| Kinnex rule (+1..+18, >=12 A) | 61,210 | 0.5944 | 0.1899 | 0.6120 | 0.3040 |

The shipped rule removes **15,586** peaks from that arm and lifts precision by 11.6 points; the Kinnex
rule removes **900** and lifts it by 0.4. Genome-wide the shipped rule flags **68,855 of 402,860
peaks (17.09 %)** against the Kinnex rule's **12,186 (3.02 %)** — **5.7x more, not less**.

> **`manuscript/16`'s premise is inverted and must be corrected.** At the level of what each rule
> removes from a call set, the shipped rule is 5.7x *stricter* than the Kinnex rule, not looser.
> 16's own finding (27.0 % of trusted-novel sites near a decoy terminus) **stands** — that residual is
> simply not fixable by swapping in the Kinnex genomic rule.

*Verifier note:* D2's first pass computed the caller's IP window as transcript-relative -10..+30
(41 nt). `ema/experimental/internal_priming.ip_window` anchors on the BED `end` for `+` and the BED
`start` for `-`, so the true window is **-9..+30 (40 nt)**. On the correct window the recomputed rule
reproduces the caller's own per-peak flag on **402,860 / 402,860 called peaks, both directions, zero
discrepancies**. All affected rows are corrected in `D2_arm_sweep.tsv` (pre-correction values kept in
`D2_arm_sweep.tsv.preverify`). **No conclusion changed.** *Also: do not quote "~17x" from
`results/perf_gap/VERIFY/README.md` correction 1 — that is an arithmetic slip in the verifier's own
prose; 68,855 / 12,186 = 5.65, and the verifier's own bottom-line summary says 5.7x.*

- [ ] Close #95's addendum with this negative result; keep #95 §1 (the strand fix, already shipped as
      #96) and #95 §3 (per-site IP output) open on their own merits.
- [ ] Add the correction to `manuscript/16_trusted_novel_kinnex.md`.
- [ ] Note in #95 that the shipped rule's *value* is measured elsewhere and is large: it is the single
      biggest accuracy lift in the caller (**issue13** §2a), which is the argument for defaulting it
      on rather than replacing it.

*Reproduce:* `awk -F'\t' 'NR==1||$1~/^4 Kinnex/' results/perf_gap/D2_rule_sweep/D2_arm_sweep.tsv | column -t -s$'\t'`;
verifier rebuild in `results/perf_gap/VERIFY/myarms/`.

---

## 2. Tier-2 (coverage-only) rescue at any depth — DEAD

*Idea:* tier-2 calls are 166,355 sites we throw away; gate them on a hexamer, an IP check and enough
local coverage and some must be real.

*Measured:* class (c1) tier-2 IP-pass has **P@100 0.0558** and class (c2) 0.0697 (both verifier-
confirmed), against the default's 0.7062. Gating hard does not save it:

| gated tier-2 arm | n | P@100 | R_det | Kinnex t5 <=25 |
|---|---:|---:|---:|---:|
| tier-2 AND hexamer AND Kinnex-IP-pass AND window reads >=50 | 24,388 | 0.0742 | 0.0071 | 0.0691 |
| ... >=100 | 16,292 | 0.0773 | 0.0049 | 0.0789 |
| ... >=200 | 9,128 | 0.0793 | 0.0028 | 0.0864 |
| ... >=100, **strong** hexamer | 4,127 | 0.1037 | 0.0016 | 0.0880 |
| **default + the >=100 arm** | 62,816 | **0.5431** | 0.1803 | **0.5868** |

Adding the best variant to the default costs **-0.163 P@100 and -0.178 long-read concordance for
+0.005 R_det**. The whole tier-2 block buys 0.062 atlas sites per call added, against the default's
1.075 (`22_performance_roadmap.md` §2.3).

*Reproduce:* `awk -F'\t' 'NR==1||$1~/^5 tier-2/' results/perf_gap/D2_rule_sweep/D2_arm_sweep.tsv | column -t -s$'\t'`
(the two class precisions are verifier-confirmed; the four gated arms are D2-measured and were not
among the verifier's 69 re-derived quantities).

---

## 3. Moving the reported cleavage coordinate — DO NOT MOVE IT

*Idea:* the clip-seeded caller's cleavage point sits a small constant distance from the true 3' end;
shift it and every metric improves.

*Measured (PBMC, default call set, transcript-orientation shifts):*

| shift | atlas P@100 | atlas P@25 | Kinnex t5 <=25 | Kinnex **decoy** <=25 |
|---|---:|---:|---:|---:|
| 0 (shipped) | 0.7062 | **0.6389** | 0.7647 | **0.1300** |
| -1 bp | 0.7062 | 0.6386 | 0.7666 | 0.1360 |
| -2 bp | 0.7062 | 0.6378 | 0.7682 | 0.1419 |
| -3 bp | 0.7063 | 0.6363 | 0.7693 | 0.1477 |

Nothing at a 100 bp window responds; P@25 gets slightly *worse*; long-read concordance rises by
0.3-0.5 pp; and internal-priming decoy proximity rises by up to **x1.14** at 25 bp (and x1.7 at
W = 10, D2 §8). The two truths do not agree on what to optimise, which is itself the argument against
tuning on either.

**This does not contradict `issue16`, it agrees with it.** Issue16 measures the same bias at
**base-pair** resolution (atlas P@1 0.2119 -> 0.2960 at -2 bp) and its proposal is to emit
`inferred_cleavage` as an **additional column** while leaving the primary coordinate alone. That is
the right shape. The negative recorded here is specifically against **moving the reported
coordinate**, and it adds one fact issue16 does not have: at the windows the benchmark actually
reports, a shift buys nothing and costs decoy proximity.

- [ ] No action beyond linking this to `issue16` and to `issue6_cleavage_offset.md` so the "+95 bp"
      recommendation and a "-2 bp coordinate move" are both closed out in one place.

---

## 4. `-F 3844 >= 2` (strict-alignment clip support) — DO NOT ADOPT

*Idea:* count only clip molecules whose reads pass `-F 3844` (drop secondary, supplementary, QC-fail,
duplicate, unmapped); stricter evidence, better calls.

*Measured (PBMC):* `f3844 >= 2` gives n 42,015 / P@100 **0.7188** / R_det 0.1632 / Kinnex 0.7768 —
but the plain support knob reaches ~0.753 precision at that recall, so the strict filter is
**dominated** (dP **-0.0346** against the support knob at matched recall). On mouse 1 it is marginal
in the other direction (n 24,940 / 0.7737 / 0.2030, dP +0.026), i.e. it does not replicate.

*And a counter-intuitive result worth someone's attention:* PBMC singletons whose **only** clip
molecule **fails** `-F 3844` have **higher** long-read concordance (**0.3863**) than those
that pass (**0.2858**) — verifier-confirmed on the contig-filtered singleton set. The strict filter appears to be discarding
informative evidence. That is a separate investigation, not a reason to adopt the filter.

*Reproduce:* `awk -F'\t' 'NR==1||$1~/^3 -F/' results/perf_gap/D2_rule_sweep/D2_arm_sweep.tsv | column -t -s$'\t'`.

---

## 5. The methodological rule that outlives all four

**The precision reference is motif-selected, so every poly(A)-signal rule scores better against it
than it deserves.** Hexamer-pass rates of the truth sets themselves: PolyASite 2.0 representative
sites **0.6511** (any-12) / 0.4312 (strong); Kinnex >=5 UMI 3' ends **0.2869** / 0.1273; a random
genic null 0.1894 / 0.0527.

Measured on the singleton class: a strong-hexamer gate moves atlas precision **x2.03**
(0.2158 -> 0.4380) and long-read concordance **x1.06** (0.2952 -> 0.3139). Against a size-matched
**random** singleton top-up, the strong-hexamer union (`>=2 OR 1-mol strong hexamer`) gains
**+0.078 atlas precision and +0.007 long-read concordance** — an 11:1 ratio; the any-12 union gains
+0.061 atlas and +0.013 long-read, 4.6:1. The any-12 union is also *worse* on internal-priming decoy proximity
than its random control (0.1505 vs 0.1274).

- [ ] Adopt as a standing rule in `24_prime_preregistration.md`: **any recall-recovery proposal is
      scored against the Kinnex <=25 bp metric before it is scored against PolyASite**, and both
      numbers are reported.
- [ ] The hexamer union (`>=2 OR 1-mol strong hexamer`, PBMC 71,074 / 0.6136 / 0.2167 / 0.3202) may
      still be shown as a labelled sensitivity arm — it does beat polyApipe on every metric — but it
      must never be sold as a quality fix and must never become a default.
- [ ] The counter-example that **passes** this test is the 3'UTR gate (atlas x3.22, long-read x1.74):
      see **issue20**.

---

## 6. Standing statements that must be corrected

- [ ] **The "both tiers" summary line is a mash-up.** "333,920 / 0.1932 / 0.2979 / 0.2415" combines
      the v2 IP arm's n, the v2 **no-IP** arm's precision and a **pre-v2** recall/F1 from
      `manuscript/15_final_gate.md` line 36 (n 328,302). Measured v2 values: **IP arm, both tiers
      333,920 / 0.2044 / 0.3043 / 0.2445**; **no-IP arm, both tiers 402,765 / 0.1932 / 0.3363 /
      0.2455**. (Verifier-confirmed, both.)
- [ ] **"Recall at or below polyApipe (-12 % PBMC, -18 % mouse)"** (`01_outline_and_journals.md`) is
      true of the >=2-molecule default and false of the tool — see **issue19**.
- [ ] **The shipped default is not the F1 optimum** on either dataset (>=1 molecule gives 0.3046 vs
      0.2811 on PBMC and 0.3758 vs 0.3213 on mouse). Say so, and defend the default on
      reliability-first grounds instead: `21_paper_architecture.md` C6.
- [ ] **`manuscript/16`'s "looser than Kinnex" premise is inverted** — §1 above.
- [ ] **Do not quote "~17x"** from `results/perf_gap/VERIFY/README.md` correction 1; the ratio is
      **5.7x** (68,855 / 12,186). Both counts are correct; only that one sentence's arithmetic is not.
- [ ] **The (b1)/(b2) split in the recall decomposition is priority-order-dependent** — reversing the
      within-group order moves 203 sites (3,814 / 4,310 -> 3,611 / 4,513). The class (b) total
      **8,124** and every headline number are invariant. Quote the total, or state the order.
- [ ] **Do not divide across rows of the denominator audit.** "44.6 % corroborated" is the **union**
      count (127,188) while recall 0.4372 is against the **x3p** count (107,260). Each row is
      internally consistent; the cross-row division is not.
- [ ] **The clip rate is ~0.57 %, not 1.15 %** — measured and independently re-implemented by the
      algorithmic-headroom verifier; the fix is already on **issue13**'s checklist. Listed here only
      so that all the numeric corrections are enumerated in one place.

---

## 7. What would reopen any of these

New **evidence**, not new geometry. Specifically: a seeding channel the library does not currently
expose (hexamer-anchored or long-read-guided candidate generation, paired bulk 3'-seq), or a
multi-library human cohort that would let the internal-priming question in **issue17** be settled
against long-read decoys. Re-proposing any rule above on the same reads should come with (i) a
size-matched random control, (ii) a matched-precision rather than matched-count comparison, and
(iii) the atlas-independent Kinnex arm — the three controls that decided all four negatives.

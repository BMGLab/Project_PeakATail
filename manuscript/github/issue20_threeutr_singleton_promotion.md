# Issue: promote single-clip-molecule tier-1 calls that lie in an annotated 3'UTR — the only singleton discriminator that moves the long-read truth, but pre-register it before it touches the default

**Labels:** enhancement, accuracy, research, needs-pre-registration. Measured by
`results/perf_gap/D3_head_to_head/` and **added to the recommendation list by the independent
verifier** (`results/perf_gap/VERIFY/README.md`, correction 2), 2026-08-21, on the frozen v2 tree
`9dfdefb`. Verdict FIXED.

Related: **#94** (the tier-1 / clip-molecule evidence this promotes), **#95** (the IP veto, which
stays a hard gate here), **#99** (PAS -> gene assignment: the same GTF machinery). Sibling drafts:
**issue14** (a calibrated per-site score, which may recover the *same* sites — see §3 Phase 0),
**issue17** (cohort corroboration, the other singleton-promotion route), **issue19** (the operating
curve this sits on). Governed by `manuscript/24_prime_preregistration.md` §3.1 (adoption criterion),
§3.3 (Rule T) and §3.4 (post-hoc results are exploratory until re-pre-registered).
Roadmap context: `manuscript/22_performance_roadmap.md` §4.2.

## 1. Symptom

**72.2 % of PBMC tier-1 sites carry exactly one clip molecule** (121,041 of 167,629 after the contig
filter) and the default discards all of them. That single class is the caller's largest discard and
the largest single component of the recall gap that is *reachable at all*: it accounts for
**26,526 of the 235,111 missed PBMC atlas sites (11.28 %)** and **81.5 % of the sites polyApipe
recovers and we do not** (14,132 of 17,332).

The caller currently has no way to tell a real single-molecule site from a spurious one. Two
candidate discriminators exist in data the caller already has — the poly(A) hexamer and the GTF.
Only one of them is real.

## 2. Measured diagnosis

**(a) The 3'UTR gate moves both truths; the hexamer gate moves only the atlas.** Measured by the
verifier on the same 121,041 contig-filtered singletons, with its own Ensembl 99 GTF-derived
references (it did not reuse ours):

| singleton stratum | n | atlas P@100 | Kinnex x3p <=25 bp |
|---|---:|---:|---:|
| all singletons | 121,041 | 0.2158 | 0.2952 |
| **in an annotated 3'UTR** | **14,708** | **0.6947** (x3.22) | **0.5143** (x1.74) |
| 3'UTR AND any-12 hexamer | 7,696 | 0.7994 | 0.5635 |
| 3'UTR AND strong hexamer (AATAAA/ATTAAA) | 4,164 | 0.8489 | 0.5865 |
| strong hexamer, any location | 24,550 | 0.4380 (**x2.03**) | 0.3139 (**x1.06**) |
| not in a 3'UTR | 106,333 | 0.1496 | 0.2650 |
| 3'UTR promoted subset, mouse 1 | 8,610 | 0.7080 | (no mouse long-read truth) |

The bottom-but-one row is the trap: a strong-hexamer gate **doubles** atlas precision and moves
long-read concordance by **6 %**. That is a motif–atlas correlation — the atlas's own hexamer rate is
0.6511 (any-12) against 0.2869 in the long-read truth and 0.1894 in a random genic null. The 3'UTR
gate is the only candidate in a 69-arm sweep plus a 17-filter search that moves the atlas-independent
truth by more than x1.1.

**(b) The combined arms. The conservative variant improves precision AND recall AND the hard-FP rate,
on both species.**

| arm | n | P@100 | R_det@100 | F1_det | Kinnex x3p <=25 | hard-FP rate |
|---|---:|---:|---:|---:|---:|---:|
| PBMC default (shipped) | 46,524 | 0.7062 | 0.1754 | 0.2811 | 0.7647 | 0.1510 |
| PBMC + 1-mol AND 3'UTR AND any-12 hexamer | **54,220** | **0.7194** | **0.2015** | **0.3148** | 0.7361 | **0.1462** |
| PBMC + 1-mol AND 3'UTR (no hexamer requirement) | 61,232 | 0.7034 | **0.2150** | **0.3294** | 0.7045 | 0.1592 |
| PBMC + all 1-mol (the >=1 arm, for scale) | 167,565 | 0.3520 | 0.2685 | 0.3046 | 0.4256 | 0.4686 |
| mouse 1 default (shipped) | 26,255 | 0.7450 | 0.2048 | 0.3213 | — | — |
| mouse 1 + 1-mol AND 3'UTR AND any-12 hexamer | 30,974 | **0.7549** | **0.2402** | **0.3644** | — | — |
| mouse 1 + 1-mol AND 3'UTR | 34,865 | 0.7359 | **0.2545** | **0.3782** | — | — |

("hard FP" = no strand-matched atlas site within 100 bp **and** no Kinnex x3p >=5 UMI 3' end within
25 bp.) Both variants take recall past **polyApipe's full-set recall** (0.1988 PBMC, 0.2499 mouse) —
the conservative one only just (0.2015), the aggressive one comfortably (0.2150).

**(c) Why this is not simply "more calls".** The promoted subsets score 0.6947 / 0.7994 (PBMC) and
0.7080 / 0.8099 (mouse) *on their own* — at or above the shipped default's own precision in the
hexamer-gated case. For contrast, a size-matched **random** singleton top-up scores 0.5358 atlas /
0.6021 Kinnex where the hexamer union scores 0.6136 / 0.6090 — i.e. the hexamer union's whole
long-read gain over random is **+0.007**. There is no equivalent random-control failure for the 3'UTR
arm: its long-read concordance is 0.5143 against the singleton pool's 0.2952.

**(d) What it does not do.** It does not move the ceiling. Classes (a)–(d) together are only 20.0 % of
the PBMC misses; 79.96 % of missed sites have no peak of any kind
(`manuscript/22_performance_roadmap.md` §2). This buys at most +0.040 R_det on PBMC.

## 3. Proposed change

**Phase 0 — the pre-registration and the transfer test, before any caller code.**
- [ ] Write the pre-registration under `24_prime_preregistration.md` §3: state the rule
      (`tier-1 AND IP-pass AND (clip_molecules >= 2 OR (clip_molecules == 1 AND 3'UTR AND hexamer))`),
      the reference GTFs, the adoption criterion and the metrics **before** re-running.
- [ ] Fix the variant **now**, in writing: the **conservative** (hexamer-gated) one is the proposal;
      the aggressive one is a labelled sensitivity arm. Choosing between them after seeing the new
      run's numbers would be exactly the tuning Rule T forbids.
- [ ] Run it on **all three datasets** at v2/prime code and evaluate against **both** truths. The
      numbers above are post-hoc on the gate data and cannot be the adoption evidence for themselves.
- [ ] **Intersect the recovered-site set with issue14's Kinnex-trained score** before building either.
      Both promote singletons; if they recover the same sites, only one should ship, and the score is
      the more general mechanism. This intersection has **not** been done.
- [ ] Decide and record the **de-novo question** (see §5): does the promotion run by default, or only
      under `--use-annotation`?

**Phase 1 — implementation, only if Phase 0 passes.**
- [ ] Implement as a promotion in the PAS-selection step, reusing the GTF already loaded for gene
      assignment (**#99**'s code path) — 3'UTR overlap on the **same strand**, exon-level, from the
      user's own `--gtf`; no new reference file, no new download.
- [ ] Keep **tier-1** and the **IP veto** as hard gates on the promoted calls, exactly as for the
      default. (`--ip-filter` behaviour per **issue13**/**#95**.)
- [ ] Hexamer = the canonical 12-hexamer list in `scripts/reliability/trusted_novel_pas.py`, searched
      in **-40..-5** with correct minus-strand handling; reuse that code, do not re-derive it.
- [ ] Emit `promoted_by` per site (`support` / `3utr_hexamer`) in `pas_support.tsv` and as a BED
      column, so the promoted class is separable downstream and never silently merged into the
      default arm.
- [ ] Ship behind `--promote-3utr-singletons`, **default OFF**, until §3.1's adoption criterion is met
      on all three datasets.
- [ ] Regression tests: a 1-molecule call inside a same-strand 3'UTR with an AATAAA at -20 is
      promoted; the same call on the opposite strand, or with the hexamer at +20, is not; and with the
      flag off the output is **byte-identical** to v2 (the §1 reproducibility rule of `24`).

## 4. How to reproduce

    export LC_ALL=C
    cd /mnt/ssd1/Projects/PeakATail_wd
    column -t -s$'\t' results/perf_gap/D3_head_to_head/tsv/08_promote_singlemolecule_pbmc.tsv
    column -t -s$'\t' results/perf_gap/D3_head_to_head/tsv/08b_promote_singlemolecule_mouse1.tsv
    column -t -s$'\t' results/perf_gap/D3_head_to_head/tsv/11_promotion_atlas_independent_check_pbmc.tsv
    sed -n '/^2\. \*\*D2/,/^3\./p' results/perf_gap/VERIFY/README.md   # the verifier's own 3'UTR table
    # per-call master table (46,524 rows: feature class, hexamer, three IP rules, all truth distances)
    head -1 results/perf_gap/D3_head_to_head/tsv/03_calls_pa_default_pbmc.tsv

Rebuild scripts: `results/perf_gap/D3_head_to_head/03_promote_singlemol.sh` (PBMC),
`03b_promote_singlemol_mouse1.sh`, `d3_feature_class{,_mouse}.sh`.

## 5. What it is expected to buy, and the one thing it costs

**Buys (measured, post-hoc):** PBMC 0.7062 / 0.1754 / 0.2811 -> **0.7194 / 0.2015 / 0.3148** and
mouse 1 0.7450 / 0.2048 / 0.3213 -> **0.7549 / 0.2402 / 0.3644** — precision up, recall up, F1 up,
hard-FP rate down (0.1510 -> 0.1462), on both species, at a cost of 2.9 pp of long-read concordance.
It is the only fix in the perf-gap study that improves precision and recall simultaneously on both
datasets, and it is trivial: the GTF is already loaded.

**Costs — state this in the changelog and in the paper, do not bury it.** A promotion rule that reads
the GTF makes the *default* annotation-dependent, and "de novo" is a load-bearing word in our
head-to-head (it is the whole reason scUTRquant is shown but unranked). Two acceptable resolutions:
ship it **off by default** so the de-novo default is preserved and the promotion is an opt-in
"annotation-assisted" mode; or ship it on and change the paper's positioning to
"de novo detection with optional annotation-assisted recovery, reported separately". **Do not** ship
it on by default and keep calling the default de novo.

Residual caveat the verifier itself flagged: annotated 3'UTRs derive from transcript ends, so a weak
prior is shared with long-read 3' ends. It is far weaker than the motif/PolyASite one (x1.74 vs
x1.06 on the independent truth), but it is not zero, and the pre-registered run should report the
promoted subset's long-read concordance as a primary number, not a footnote.

# Issue: emit an `inferred_cleavage` column at −2 bp (base-pair resolution), and retract the `cleavage_offset = 95` recommendation for the clip-seeded caller

**Labels:** enhancement, accuracy, docs. Measured by `results/algo_headroom/A4_resolution/`
(tables T5, T5b, T5d, T5e, T7) on the frozen v2 tree `9dfdefb`; the −2 bp row is endorsed in the
independent verifier's consolidated verdict (`VERIFY/tables/VERDICT_FOR_PI.tsv`).

Related: this **partially retracts** the recommendation in `manuscript/github/issue6_cleavage_offset.md`
(the issue-#72 lineage). That analysis was done on the **coverage** caller (`lg_annotate`); it does
not transfer to the clip-seeded caller, which anchors on the soft-clip boundary rather than on where
R2 coverage runs out.

## 1. Symptom

Two separate things:

**(a)** `ema/config.py:353` ships `cleavage_offset = 0` but `ema/cli/config_schema.py:503` tells the
user *"A sane data-driven constant is ~90-100 (try 95)"*, and `issue6` recommends estimating and
applying an offset of that size. For the clip-seeded caller that advice is **destructive**.

**(b)** The clip-seeded caller's reported cleavage point sits a small, constant, **measurable**
distance from the true 3' end, and we report no corrected coordinate — so every bp-exact comparison
(ours and anyone else's) is penalised for a bias we can name.

## 2. Measured diagnosis

**(a) The +95 bp offset costs 46 points of P@10 on the clip-seeded default** (PBMC, atlas):
P@10 0.5209 → **0.0551**, P@100 0.7062 → 0.6655 (`A4 tables/T7`). The coverage caller's peak 3' end
stops ~90–105 nt short of cleavage because R2 coverage runs out; the clip-seeded caller does not have
that problem, because it anchors on the soft-clip boundary.

**(b) The real bias is −2 bp, and both truths agree in sign.** Signed distance from call to nearest
long-read 3' end, 40,052 default calls with a Kinnex t5 site within 200 bp: **median −2 bp on BOTH
strands** (+: −2, −: −3 at t5; −2/−2 at t20), modal offset −1 bp, per-chromosome median −2 on every
chromosome checked. Applying a constant −2 bp shift (transcript orientation):

| metric | unshifted | −2 bp |
|---|---|---|
| atlas P@1 | 0.2119 | **0.2960** |
| Kinnex t20 P@1 | 0.1782 | **0.2642** |
| atlas P@10 | 0.5209 | 0.5093 |
| Kinnex t20 P@10 | 0.4979 | 0.5108 |
| atlas P@100 | 0.7064 | 0.7064 |

So it is a **base-pair-resolution** result and nothing else: **no window >= 25 bp responds to any
shift in −6..+5**, and at W = 10 the two truths disagree in sign (atlas prefers +3, Kinnex prefers
−5), which is itself a warning against tuning on either.

**The bias is support-dependent**: median offset −5 bp at 2 molecules, −3 at 3–4, **−2 at >= 5**, and
the fraction of calls within ±2 bp rises 0.157 → 0.630 across the same range (`T5b`). A
support-conditional shift was tested and is **not** better than the constant (`T5e`).

**It is not universal across datasets.** On mouse testis the P@1 optimum is **−1 bp**, not −2
(P@1 0.3682 at −1 vs 0.3166 unshifted, mouse 1; `T5d`). Any shipped constant must therefore be
per-chemistry, or — better — the shift must be reported rather than applied.

## 3. Proposed change

- [ ] Emit `inferred_cleavage` as an **additional column** in `pas_support.tsv` / `annotatedpas.bed`
      = called cleavage − 2 bp in transcript orientation. **Do not move the reported coordinate**;
      a coordinate move changes the pre-registered output and would need re-pre-registration
      (24 §3.4), while a column does not.
- [ ] Also emit the per-site clip-molecule count next to it (already available) so a user can apply
      the support-dependent correction themselves if they want; ship the constant, document the
      gradient (−5 / −3 / −2 bp at 2 / 3–4 / >= 5 molecules).
- [ ] Make the constant a config value (`cleavage_offset_bp`, default −2 for `clip_seeded`) and
      document that it is chemistry-dependent — mouse testis measures −1.
- [ ] **Change the `--cleavage-offset` help text**: remove *"A sane data-driven constant is ~90-100
      (try 95)"* or scope it explicitly to the coverage strategies, and add the measured cost for
      `clip_seeded` (P@10 0.5209 → 0.0551).
- [ ] Add a note to `issue6` (or its GitHub descendant) recording the retraction, so the +95
      recommendation does not resurface.
- [ ] Regression test: on the committed fixture, `inferred_cleavage` = cleavage − 2 on `+` and
      cleavage + 2 on `−` (transcript orientation), and the primary coordinate is byte-identical to
      v2.

## 4. How to reproduce

    export LC_ALL=C
    cd /mnt/ssd1/Projects/PeakATail_wd/results/algo_headroom/A4_resolution
    cat tables/T5_offset_kinnex.txt      # signed offsets, both strands, three UMI tiers
    cat tables/T5b_offset_by_support.txt # the molecule-count gradient + per-chromosome constancy
    cat tables/T5e_support_conditional_shift.txt
    cat tables/T5d_offset_mouse.txt      # mouse optimum is -1, not -2
    grep -n "auto_cleavage\|try 95" ../../../tools/pa-polya-run-9dfdefb3/ema/cli/config_schema.py

## 5. What it is expected to buy

**No change to any published number** — P@100 0.7064 → 0.7064, and every window >= 25 bp is flat. It
buys (i) a defensible base-pair-resolution statement for the paper (atlas P@1 0.2119 → 0.2960,
Kinnex P@1 0.1782 → 0.2642, i.e. ~40–48 % more exactly-correct calls), (ii) removal of a live
recommendation in our own CLI help that would cost a user 46 points of P@10, and (iii) an honest
calibration constant to report rather than a bias to be penalised for. It is the cheapest item on the
algorithm roadmap (`manuscript/23_algorithm_roadmap.md` §3 Step 3) and the only one that can land
before manuscript submission.

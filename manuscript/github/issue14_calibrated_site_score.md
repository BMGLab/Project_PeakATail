# Issue: replace the `>=2 clip molecules` gate with a calibrated per-site score — but run the cross-dataset transfer test FIRST

**Labels:** enhancement, accuracy, research, needs-pre-registration. Measured by
`results/algo_headroom/A3_scoring_model/` and **corrected** by the independent verifier
(`results/algo_headroom/VERIFY/VERIFICATION_REPORT.md` §1.4–1.5, §7.2), 2026-08-21, on the frozen v2
tree `9dfdefb`.

Related: **#94**, **#95**. Governed by `manuscript/24_prime_preregistration.md` §3.1 (adoption
criterion), §3.3 (Rule T — no threshold tuned on the metric it is reported against) and §3.4 (every
prime number is exploratory until re-pre-registered).

## 1. Symptom

The precision-first default is a hard conjunction: `tier-1 AND IP-pass AND >= 2 distinct clip
molecules`. The molecule count is an integer with a very coarse grain — on PBMC the threshold jumps
the call set from 167,565 (k>=1) to 46,524 (k>=2) — and **72.1 % of tier-1 sites carry exactly one
clip molecule**, so the single most consequential decision in the caller is a binary applied to the
largest class it has. The other things the caller already knows about a site at call time (hexamer
presence and offset, downstream A-content as a float, 3'UTR/exon overlap, distance to the nearest
annotated 3' end, local candidate context, reads-per-cell concentration) are used for gating or not
used at all, never for ranking.

## 2. Measured diagnosis

**The rule is close to optimal *given the evidence it uses*.** Soft-combining only what it already
uses (`clip_reads`, `clip_umis`, the two `-F 3844` variants, tier, the IP flag) is worth +1.9 % to
+5.8 % — inside the noise of a threshold change.

**A score that adds the cheap covariates does lift the curve, and the size is +7 %, not +39 %.**
The honest configuration — gradient-boosted model trained **only on the PacBio Kinnex long-read
truth**, with all 12 downstream-A features dropped, on chromosome-disjoint folds (fold 0 = odd
autosomes, fold 1 = even + X + Y; verified: 0 contigs split, 0 duplicate sites, no feature intersects
any truth column), then evaluated on the **PolyASite atlas it has never seen**:

| arm | n | P@10 | P@25 | P@100 | R_det@100 | Kinnex decoy@25 |
|---|---|---|---|---|---|---|
| RULE `tier1 & IP-pass & >=2 mol` (shipped default) | 46,524 | 0.5209 | 0.6389 | 0.7062 | 0.1754 | 0.1300 |
| **score, matched call count AND matched window-read decile** | 46,524 | **0.5567** | **0.6835** | **0.7766** | **0.1870** | **0.0954** |
| score, at matched atlas precision | 56,286 | 0.4934 | — | 0.7059 | **0.2005** | 0.1048 |

That is **+6.6 % relative recall together with +7.0 precision points at identical call count and
identical expression**, or **+14.3 % relative recall at matched precision**, with the
internal-priming decoy rate down by a quarter — on a truth the model never saw
(`VERIFY/tables/v9_expr_matched_kinnex_trained.tsv`, `v3_decontamination.tsv`).

**Three things that are NOT true and must not be quoted:**
1. The **atlas-trained** model's "+39 % / +48 %" is circular. On every empirical independent axis it
   is *worse* than the rule: −8 % at matched call count, −8 % at matched long-read precision, −11 %
   at matched high-confidence long-read precision.
2. A3's "+16 % on the atlas-independent curve" is a **matched-odds artefact**: at odds 5.88 the model
   sits at n ≈ 84,700 with long-read precision 0.587 against the rule's 0.765 — it matches the odds
   only by halving the decoy denominator.
3. The score **cannot improve base-pair resolution**. At matched precision P@10 *falls* (0.4934 vs
   0.5209). A re-ranker re-orders candidates; it cannot move a call.

**Feature evidence** (held-out AUC / average precision, atlas label): all groups 0.9486 / 0.8048;
minus sequence 0.9188 / 0.7391 (largest single loss); minus GTF annotation 0.9331 / 0.7692;
**annotation-free performs the same on the independent curve**, so the score does not need the GTF.
Calibration against the training label is essentially perfect (ECE 0.0038); against the *other* truth
the same probability is optimistic in the mid-range (predicted 0.6 → observed 0.48).

**Do not build:** a second BAM pass for clip-cell counts, top-cell share, end-position entropy or
pileup sharpness — measured on chr19+21, worth **0.000–0.002 AUC**.

## 3. Proposed change

**Phase 0 — the experiment that decides everything (do this before writing any caller code).**
- [ ] Refit the verifier's exact recipe (`VERIFY/code/v3_a3_decontaminate.py`:
      `HistGradientBoostingClassifier(max_iter=400, learning_rate=0.06, max_leaf_nodes=31,
      min_samples_leaf=100, l2_regularization=1.0, early_stopping=True)`, downstream-A features
      dropped) on **GSE104556 mouse 1** and evaluate it **unchanged** on **PBMC 10k v3** and **mouse
      2** — the primary protocol of 24 §3.3.
- [ ] Report both folds and label mouse 1's own numbers `THRESHOLD-FITTING SET — not evidence`.
- [ ] **If it does not transfer, close this issue as a negative result.** Everything measured so far
      is one PBMC library; chromosome-disjoint CV tests positional generalisation, not dataset or
      species transfer.

**Phase 1 — implementation, only if Phase 0 passes.**
- [ ] Score computed at call time from features already available: clip evidence, coverage/width/cells,
      hexamer presence + offset in −40..−5, **continuous** downstream A-content (see #94/**issue13**),
      local candidate context. Ship the **annotation-free** variant so the caller stays de novo.
- [ ] Coefficients trained **offline** and shipped as numpy constants — scikit-learn must not enter
      the tool's runtime import path (24 §4). Commit the training script, its inputs and the fold
      definition alongside the constants.
- [ ] **Keep `tier-1` and the IP veto as hard gates.** Every configuration in which the score was
      allowed to override the IP veto looked spectacular on the atlas and no better than the rule on
      long reads. Applying the veto at *selection* time on top of the score still cuts the decoy rate
      from 0.088 to 0.060.
- [ ] Emit `score` (calibrated probability) as a column in `pas_support.tsv` and the PAS BEDs
      **regardless of the default threshold**, so users can pick their own operating point.
- [ ] Ships behind `--score-model`, **defaulted OFF**, until 24 §3.1 is met on all three datasets.

## 4. How to reproduce

    export LC_ALL=C
    cd /mnt/ssd1/Projects/PeakATail_wd/results/algo_headroom
    cat VERIFY/tables/v9_expr_matched_kinnex_trained.tsv   # the headline, matched count + expression
    cat VERIFY/tables/v3_decontamination.tsv               # atlas/independent curves, 4 model arms
    sed -n '1,140p' VERIFY/VERIFICATION_REPORT.md          # §1.1-1.5: leakage audit and the correction
    python3 VERIFY/code/v3_a3_decontaminate.py             # refits from A3_scoring_model/data/*.parquet

Feature table + schema: `A3_scoring_model/data/candidates.parquet`, `FEATURE_TABLE_SCHEMA.tsv`.

## 5. What it is expected to buy

At the shipped operating point: **+6.6 % relative recall and +7 precision points at the same number
of calls**, with ~26 % fewer internal-priming decoy hits — or ~+14 % recall if the operating point is
moved to matched precision. On PBMC that is R_det@100 0.1754 → 0.1870 at n = 46,524, or 0.1754 →
0.2005 at P@100 ≈ 0.706. **It changes the pre-registered default, so nothing from it may enter the
manuscript until it has its own pre-registration, a fresh three-dataset run and a verifier pass
(24 §3.4).** The secondary prize is arguably bigger than the primary one: a calibrated per-site
confidence makes `>=2` a default rather than a commitment, and "the tool emits a confidence, the
competitors emit a set" is a cleaner claim than a two-point precision win.

**Ceiling, so nobody over-scopes this.** The tier-1 & IP-pass candidate pool caps R_det@100 at
**0.269** (`19_final_gate_v2.md`, the `>=1` molecule arm), and the Kinnex-trained score at matched
precision already reaches 0.2005 = **75 % of it** (the withdrawn atlas-trained model reached 91 %,
but only against the atlas it was trained on). Of the 11,606
chr19+21 atlas sites the default misses, only 21.4 % have *any* candidate within 100 bp. After this
issue, the next gain has to come from the candidate generator, not from scoring.

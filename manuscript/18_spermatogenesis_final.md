# Spermatogenesis 3′-UTR control on the final caller — verified (2026-08-21)

**Verdict: FIXED** (all headlines reproduced by an independent verifier from the raw run outputs;
7 small report defects repaired, none load-bearing). Code `4efeb125` (final caller); both testis mice;
GEX-derived stage labels (exact-join, multiset-preserving permutations verified); PROVISIONAL only in
the sense that the post-#96/#97 v2 re-run will redo it. Full detail:
`results/stage3_spermatogenesis_final/REPORT_PROVISIONAL.md` (+ `verify/`).

> **2026-09-02: the promised v2 re-run is done and verified — the [v2 section](#v2-section-merged-caller-record-verified-sound-2026-09-02)
> below is now the record.** Everything from here to that section is the v1-code (`4efeb125`) result,
> retained as the labelled comparison exactly as [20](20_stage3_replication.md) retains its v1 row.

## The claim, in the exact form that survives (v1 code — superseded, see the v2 section)

> Per-gene 3′UTR usage shifts progressively proximal along spermatocyte → round → elongating spermatid
> in both mice: 31.5% and 30.2% of depth-guarded genes shorten monotonically across the three stages
> against label-shuffle nulls of 17.5% and 21.0% (z = 10.5 and 5.8, 20 shuffles), exceeding monotone
> lengthening (489 vs 325 and 361 vs 296 genes), and the composition-controlled per-cell distal-usage
> residual falls at every step (Cliff's δ SPC vs ES = 0.56 and 0.61, outside the entire label-shuffle null range). The
> gradient is a per-gene, equal-weight statement: the per-gene medians are not monotone (RS marginally
> above SPC) and the across-gene UMI-weighted distal index reverses at RS→ES, where protamine
> transcripts alone carry ~26% of the UMIs at ceiling PDUI.

Cautions bound to the claim: monotone lengthening is also above its null (ordered structure), so the
direction-specific evidence is the shortening-vs-lengthening excess — strong in mouse 1 (binomial
p 9.9e-9), modest in mouse 2 (p 0.0125; 1 of 40 null draws reached the same excess in absolute value (0/40 one-sided)). Wilcoxon p-values are quoted only
next to their shuffle-null column. SPG→SPC *lengthens* (kept out of the claim per [11](11_verifier_corrections.md)).

## Switch test + cross-mouse replication (the reliability payload)

- Arm B0 (fisher, cells mode, no pre-selection) **re-validates on the final caller**: null p<0.05
  3.169% / 3.166%, 0 q<0.05 hits in all 30 null BH families, 0/10 null runs with a hit. TRUE hits
  50,415 / 52,043 over the three stage pairs.
- Cross-mouse replication (gene+strand, ≤100 bp matching): **6,049 / 9,469 / 11,263** same-direction
  replicated PAS per pair; sign agreement 99.7–99.8%; enrichment 3.2–5.1× over independence;
  **0 replicated in all 15 null pairings**; per-gene effects ρ = 0.655 (n = 917; null 0.001 ± 0.032).
  This is the two-replicate version of the manuscript's replication filter and it behaves exactly as
  pre-registered.

## Negative findings (reported, not hidden)

1. **The across-gene UMI-weighted distal index reverses at RS→ES** (Cliff −0.95 / −0.86): driven by
   protamine ceiling-PDUI transcripts dominating ES UMIs. Different quantity from manuscript/11's
   `wdi` (name collision documented). Fig 5 must not use across-gene UMI weighting.
2. **The 16-gene literature panel does not reproduce**: of 16 measurable, 4 shorten in both mice
   (Prm3, Ybx2, Ppp1cc, Nsun7), 5 lengthen in both (Prm1, Tnp2, Odf1, Spata19, Smcp), 5 strictly
   discordant, 2 uninformative. The quarantined "6/6" came from a different metric and gene set. The
   literature-panel leg of Fig 5 is dropped or shown with these numbers.
3. **Tool defect** (issue draft [github/issue11_per_isoform_degenerate.md](github/issue11_per_isoform_degenerate.md)):
   `switch length --isoform-agg per_isoform --utr-unmatched drop` fails on both mice — it emits
   degenerate proximal==distal pairs (14 (gene,transcript) units × all cells; every read-bearing row
   PDUI ≡ 0.5; verifier-reproduced), and the guard at `ema/switch_test/runner.py:397` classifies the
   equality as a strand inversion (`~(prox < dist)`), blaming the wrong cause. The per_gene path passes
   the same guard cleanly on both mice — the strand-correctness re-validation of 13 §2 **passes**.

Figure: `figures/spermatogenesis_final.{png,pdf}` (+ `.caption.md`; script `scripts/manuscript_figures/spermatogenesis_final.py`; every plotted value in `results/figures/manuscript/spermatogenesis_final*.tsv`) (now `fig5_spermatogenesis` throughout: `figures/fig5_spermatogenesis.{png,pdf,caption.md}`, `scripts/manuscript_figures/fig5_spermatogenesis.py`, `results/figures/manuscript/fig5_spermatogenesis*.tsv`).

## v2 section: merged caller record, verified SOUND (2026-09-02)

**Verifier verdict: SOUND. This section is the result of record; the v1 numbers above are the
labelled comparison** (the structure [20](20_stage3_replication.md) uses). Code = the frozen MERGED
tree `9dfdefb3` (commit `9dfdefb3eb353b0817ef79c4eb9ace6d6c8aab53`; PR #96 IP-filter/strand +
PR #97 clip-memory), snapshot `WD/tools/pa-polya-run-9dfdefb3` — all 16 `DONE.ok` files carry that
sha, exit 0, and the queue logs assert `ema.__file__` inside that tree. Inputs = the pre-existing v2
benchmark arms (`results/benchmark_tools/gse104556/peakatail_clipseeded_final_v2/mouse{1,2}/run`,
h5ad sha256 verified against each run manifest; no caller re-run). Labels are code-independent and
were reused unchanged: the GEX label TSVs (md5 `37f9bbafb82d917d5a0a374e13dd2fd0` /
`c518a4c7a6c0e16fb41a3a7aad3a7bc4`, mtime 2026-08-19) are the same files the v1 record consumed, and
the verifier confirmed the cell→stage maps of the v1 and v2 input trees are identical
(order-independent digest match, both mice; stage counts 370/504/356 and 268/757/271 on SPC/RS/ES,
SPG 64/68 excluded from the claim as before). Same test configuration (arm B0:
`fisher --count-mode cells --marker-top-n 0 --counts-layer counts`), same null designs (20 label
shuffles for the PDUI arm, 5/mouse through switch diff, 15 null pairings), same summarisers. Full
detail: `results/stage3_spermatogenesis_v2/REPORT_PROVISIONAL.md` (§8 has every headline side by
side with v1; zero movements beyond its "small" tier).

**Adversarial verification (2026-09-02).** An independent verifier recounted, from the raw v2 run
outputs along its own code path (own pandas/awk, not the builder's summarisers, same null seeds):
the guarded-gene counts and monotone fractions with their 20-shuffle nulls (both mice, fixed and
free-guard variants), the Wilcoxon and residual-Cliff statistics, the switch TRUE and null hit
counts and pooled null rates (both mice, all 30 BH families), the cross-mouse replication for **all
three** stage pairs end-to-end (own greedy matcher), the 15 null-pairing zeros, the per-gene effect
correlation, the direction-excess binomials, and the pasbed churn by strand. **Every number matched
the builder's to the digit.** Label provenance, stage restriction and the frozen-tree assertion were
verified from the logs and label files as above.

### The claim, in the exact form that survives on v2 (the record)

> Per-gene 3′UTR usage shifts progressively proximal along spermatocyte → round → elongating spermatid
> in both mice: 31.4% and 30.5% of depth-guarded genes shorten monotonically across the three stages
> against label-shuffle nulls of 17.4% and 21.0% (z = 10.5 and 6.4, 20 shuffles), exceeding monotone
> lengthening (486 vs 318 and 366 vs 291 genes), and the composition-controlled per-cell distal-usage
> residual falls at every step (Cliff's δ SPC vs ES = 0.55 and 0.63, outside the entire label-shuffle null range). The
> gradient is a per-gene, equal-weight statement: the per-gene medians are not monotone (RS marginally
> above SPC) and the across-gene UMI-weighted distal index reverses at RS→ES, where protamine
> transcripts alone carry ~25% of the UMIs at ceiling PDUI.

Cautions bound to the claim, as in v1: monotone lengthening is also above its null (ordered
structure), so the direction-specific evidence is the shortening-vs-lengthening excess — strong in
mouse 1 (binomial
p 3.4e-9), and on v2 clearer in mouse 2 than it was on v1 (p 0.0039; 0 of 40 null draws reached the
excess in absolute value, where v1 had 1 of 40). Wilcoxon p-values are quoted only next to their
shuffle-null column. SPG→SPC still *lengthens* and stays out of the claim per
[11](11_verifier_corrections.md).

### v1 → v2, every claim-bearing number (verified to the digit on both sides)

| quantity | v1 code `4efeb125` | v2 code `9dfdefb3` |
|---|---:|---:|
| guarded genes (≥50 UMI/stage) | 1,553 / 1,194 | 1,548 / 1,200 |
| monotone shortening | 489 (31.5%) / 361 (30.2%) | 486 (31.4%) / 366 (30.5%) |
| shuffle null (fixed set) | 17.5% / 21.0% | 17.4% / 21.0% |
| z vs null | 10.49 / 5.83 | 10.53 / 6.42 |
| monotone lengthening | 325 / 296 | 318 / 291 |
| excess binomial p | 9.9e-9 / 0.0125 | 3.4e-9 / 0.0039 |
| null draws reaching the m2 excess | 1/40 | 0/40 |
| per-cell residual Cliff δ SPC vs ES | 0.562 / 0.608 | 0.554 / 0.632 |
| arm B0 TRUE q<0.05 (3 pairs) | 50,415 / 52,043 | 50,354 / 52,111 |
| pooled null p<0.05 | 3.169% / 3.166% | 3.169% / 3.156% |
| null q<0.05 hits (30 BH families) | 0; 0/10 runs | 0; 0/10 runs |
| replicated same-direction PAS (ES_RS / ES_SPC / RS_SPC) | 6,049 / 9,469 / 11,263 | 6,070 / 9,480 / 11,219 |
| replicated in the 15 null pairings | 0 | 0 |
| sign agreement among both-hit | 99.7–99.8% | 99.7–99.8% |
| enrichment over independence | 5.06× / 3.73× / 3.23× | 5.05× / 3.75× / 3.25× |
| per-gene effect ρ (n) | 0.655 (917) | 0.641 (923) |
| both-shorten vs both-lengthen genes | 338 vs 278 (p 0.0174) | 339 vs 273 (p 0.0085) |

**Why the numbers moved (verified cause).** The v1→v2 pasbed churn is **100% minus-strand** —
mouse1 −1,505/+1,248 PAS, mouse2 −1,302/+1,279, **zero** plus-strand changes in either mouse
(recounted by the verifier from the pasbeds themselves) — exactly the footprint of the #96
IP-filter/strand fix. Labels, stage counts, cell sets and test configuration are identical to v1, so
every movement above is bounded by that ~2% minus-strand PAS-universe change; the largest single
movement is mouse 2's z (5.83 → 6.42, its null sd shrank 0.0158 → 0.0148).

### Negative findings — all three re-checked on v2, all three still hold

1. **The across-gene UMI-weighted distal index still reverses at RS→ES** (Cliff −0.93 / −0.87).
   Protamine domination re-measured on v2: Prm1+Prm2 = 25.2% of mouse-1 guarded ES UMI (14.5% in
   mouse 2), both at PDUI ~1.00. Fig 5 must not use across-gene UMI weighting.
2. **The 16-gene literature panel still does not reproduce**, with the identical classification:
   4 shorten in both (Prm3, Ybx2, Ppp1cc, Nsun7), 5 lengthen in both (Prm1, Tnp2, Odf1, Spata19,
   Smcp), 5 discordant, 2 uninformative.
3. **Issue #98 reproduces bit-for-bit** (`per_isoform` degenerate pairs; both mice rc=1; same
   offending-row counts as v1 — m1 18,116 rows / 7 genes, m2 15,004 / 8; guard-intercept probe:
   every flagged row degenerate, 0 strand inversions). Neither #96 nor #97 touched `per_isoform`,
   so the 3′UTR-scoped "by UMI signal" arm of manuscript/11 remains unproducible, as in v1.

### Disclosures (immaterial, but they travel with the record)

- Two v1 prose-level signs flip at trivial magnitude, both outside the claim: the gene-median arm's
  SPC→RS rank-biserial (−0.034 → +0.024; the arm both reports declare **unusable at this depth** —
  guard leaves 215/194 genes, every median delta 0.0 to 5 dp) and mouse 1's supplementary SPG→SPC
  Wilcoxon crossing 0.05 (0.0354 → 0.0509; SPG is excluded from the claim).
- The v2 report's gate table quotes mouse 2's null p<0.05 as 3.15% (mean of per-run fractions); the
  pooled-by-rows convention v1's record used gives 3.156% (29,017/919,545). Both conventions clear
  the ≤7% gate; the pooled figure is quoted in the table above for like-for-like comparison.
- The mouse-1 `switch length` output was produced hours earlier in the same day by another session
  with the byte-identical command on the same frozen commit (CMD line and `DONE.ok` sha verified);
  every downstream consistency assert (cell set, cluster≡stage per row, strand convention per row,
  v2 PAS universe) passed on it, and the verifier's recount ran on it directly.

Figure: `fig5_spermatogenesis` regenerated from this v2 tree (`SPERMATOGENESIS_VERSION=v2`, now the
script default; `v1` still reproduces the v1-code render). Sources:
`results/stage3_spermatogenesis_v2/{mouse{1,2}/summary,summary}/`.

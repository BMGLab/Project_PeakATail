# Spermatogenesis 3′-UTR control on the final caller — verified (2026-08-21)

**Verdict: FIXED** (all headlines reproduced by an independent verifier from the raw run outputs;
7 small report defects repaired, none load-bearing). Code `4efeb125` (final caller); both testis mice;
GEX-derived stage labels (exact-join, multiset-preserving permutations verified); PROVISIONAL only in
the sense that the post-#96/#97 v2 re-run will redo it. Full detail:
`results/stage3_spermatogenesis_final/REPORT_PROVISIONAL.md` (+ `verify/`).

## The claim, in the exact form that survives

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

Figure: `figures/spermatogenesis_final.{png,pdf}` (+ `.caption.md`; script `scripts/manuscript_figures/spermatogenesis_final.py`; every plotted value in `results/figures/manuscript/spermatogenesis_final*.tsv`).

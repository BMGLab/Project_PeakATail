Atlas-independent motif analysis of the `lg_annotate` calls (22,629 PAS, Laughney cohort) shows the called **peak 3′ end systematically stops ~90–105 nt short of the actual cleavage site**:

- AATAAA positional density peaks at **+75 nt downstream** of the reported peak end (canonical AATAAA→cleavage spacing is 15–30 nt ⇒ implied cleavage ≈ +90..+105).
- Genomic A-fraction crests at 43% at **+98 nt** then cliffs to background — the signature of the poly(A) junction — exactly where 10x R2 coverage runs out.
- Windows that assume peak end = cleavage site contain the AATAAA signal at near-null rates (7.7% vs 5.4% shuffled); re-anchored +0..+100 the rate is 37% vs 15% null.

This has direct consequences:

1. **Benchmarks at tight cutoffs punish the offset, not the calls.** Under strand-matched point-mode scoring vs curated PolyASite 2.0, precision@100bp is 0.36–0.49 (≈20× null) but would rise substantially if the inferred cleavage position were corrected by the observed offset — @200bp already reaches 0.56.
2. **`switch geneview` coordinates and any PAS-to-atlas annotation** inherit the same shift.

## Proposal

Estimate the per-run (or per-peak) cleavage offset from the data itself — e.g., the mode of the AATAAA positional profile + canonical spacing, or the A-fraction cliff — and report an `inferred_cleavage_site` column alongside the raw peak end in `pasbed.bed`/`annotatedpas.bed`. Downstream matching/benchmarking then uses the corrected coordinate. The offset is chemistry-dependent (R2 length), so a data-driven estimate beats a constant.

Full analysis + verified numbers: `PeakATail_wd/manuscript/07_curated_benchmark_report.md` and `manuscript/figures/motif_validation.*` (BioLab box).

@TRextabat — this one is genuinely interesting rather than a bug report: it likely recovers a large chunk of benchmark precision for free, for every strategy at once.

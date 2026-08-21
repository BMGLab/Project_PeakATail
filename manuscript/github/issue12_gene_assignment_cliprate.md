# Issue: PAS→gene assignment picks the readthrough/neighbour model in overlapping loci (CD68 gets 0 PAS cohort-wide); clip-rate warning samples only the head of chr1

**Labels:** bug, annotation. Found by the Stage-3 replication verifier on the 17-sample Laughney cohort run (code `4efeb12`).

## 1. Gene assignment in overlapping loci
In loci where a readthrough/overlapping model spans a neighbour's 3′ UTR, `annotatedpas.bed` assigns
every PAS to the spanning model and the neighbour gets **zero** PAS cohort-wide, even though
`gene_end.bed` contains it. Verified examples (GRCh38.99): all PAS in chr17:7,579,636–7,582,386 are
labelled SENP3-EIF4A1 while they lie in **CD68**'s 3′ UTR (CD68: 0 assigned PAS in the whole cohort);
PTPRCAP under CORO1B; FKBP11 under AC073610.2; also GPX1←RHOA, MDK←DGKZ, STARD10←ARAP1, ARPC1A←AC004922.1,
FCER1G←NDUFS2, KRTCAP3←NRBP1. Impact: 33% of the top-30 replicated "gene switches" in the Stage-3
analysis carry the wrong gene name (21% of top-100, 12% of top-500, 7% overall) — the strongest hits are
enriched because marker genes (CD68, PTPRCAP, FKBP11) live in exactly such loci. Proposal: prefer the
gene whose annotated 3′ UTR (or last exon) contains the PAS over a model that merely spans it; report
ties; add a regression test on the CD68/SENP3-EIF4A1 locus.

## 2. Poly(A) clip-rate warning is computed on an unrepresentative sample
The caller's low-clip-rate warning samples "the first 200,000 CB reads" of the BAM. On a
coordinate-sorted BAM that is the head of chr1, not the library: GSM3516664 (MetBone) was flagged at
0.015% while its genome-wide rate is mid-cohort (306,202 clip molecules; 31.1% of tested PAS with ≥2
own molecules). This misled a pre-registered data-quality exclusion. Proposal: sample uniformly across
the BAM (e.g. every Nth read or per-chromosome strata) or compute the rate during the main pass and
warn at the end.

Evidence: `manuscript/20_stage3_replication.md` §disclosures; verifier outputs
`/mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v2/replication/verify/` on biolab.

# Figure — `fig_switch_tracks` caption

## Legend

**Replicated poly(A)-site switches shown as read coverage at the locus.** Four genes
selected from the 1,077 genes whose spermatogenesis switch replicated in both mice and
that carry exactly two replicated significant sites, so that proximal and distal are
unambiguous; genes with three or more significant sites were excluded from this figure
because no two of them can be labelled proximal and distal without ambiguity. Three
shortening examples and one lengthening example are drawn, against the 65%
shortening / 35% lengthening split of all 1,077 such genes, so the panel is
representative of the direction balance rather than of shortening alone. Each column is one
mouse, scored independently end to end. Within a column, the three tracks are the three
spermatogenic stages (SPC, spermatocyte; RS, round spermatid; ES, elongating spermatid) and
share one vertical scale, printed beneath them, because the switch is a change in the shape
of the profile and a per-track autoscale would conceal it. Coverage is the mean number of
reads per base per cell, so tracks are comparable between stages containing different
numbers of cells (mouse 1 370/504/356 SPC/RS/ES, mouse 2 268/757/271 SPC/RS/ES). Read acceptance matches the caller: SAM flag filter 3844, the
tool's own CLIP_EXCLUDE_FLAGS, which drops unmapped, secondary, QC-failed, duplicate and
supplementary records, and only reads on the gene's own strand are drawn, which is how the
caller counts; 12.1-99.8% of reads over these four loci are on that
strand in any case. Dashed lines mark the two called sites, green proximal and orange
distal; the arrow gives the direction of transcription and the grey bar the interval
between the two sites. The percentages on each track are the quantity the switch test
actually used, the proximal and distal share of that gene's reads in that stage
(n_reads_pas / n_reads_gene from the differential tables), and are not measured off the
coverage drawn here.

Per gene: *Rragc* (chr4, + strand, shortening; proximal 123,935,899, distal 123,936,992; max q = 3 x 10^-56 across both sites and both mice): mouse 1 SPC 4/96, RS 33/64, ES 92/2; mouse 2 SPC 3/94, RS 33/64, ES 94/2. *Pttg1ip* (chr10, + strand, shortening; proximal 77,597,270, distal 77,598,731; max q = 3 x 10^-61 across both sites and both mice): mouse 1 SPC 12/87, RS 87/11, ES 92/2; mouse 2 SPC 11/87, RS 88/10, ES 97/1. *Pom121* (chr5, - strand, shortening; proximal 135,377,832, distal 135,376,139; max q = 3 x 10^-20 across both sites and both mice): mouse 1 SPC 2/55, RS 16/35, ES 57/12; mouse 2 SPC 1/87, RS 12/52, ES 61/4. *Map3k11* (chr19, + strand, lengthening; proximal 5,702,861, distal 5,707,100; max q = 5 x 10^-33 across both sites and both mice): mouse 1 SPC 92/8, RS 65/32, ES 0/95; mouse 2 SPC 93/7, RS 40/57, ES 2/98.

**Caveats that travel with this figure.** These are four loci chosen to be legible, not a
random sample, and they are drawn from the replicated set, so they show what a replicated
switch looks like rather than establishing how often one occurs; the frequency claims are
in Fig 5 and manuscript/20. Coverage height reflects both site usage and gene expression,
which is why the tested usage fractions are printed rather than inferred from peak height.
Stage labels are expression-derived cluster assignments, not sorted populations. The two
mice are one study and one chemistry.

## Provenance

Coverage `results/figures/manuscript/switch_tracks_coverage.tsv.gz`, track metadata
`switch_tracks_meta.tsv`, tested usage `switch_tracks_usage.tsv`. Built by
`scripts/manuscript_figures/_switch_track_extract.py` (cached BAM pass) and
`scripts/manuscript_figures/fig_switch_tracks.py`. Switch calls from
`results/stage3_spermatogenesis_v2` (code 9dfdefb); BAMs
`data/benchmark/gse104556/starsolo/Mouse{1,2}_scRNAseq`. Canvas 7.09 x 7.5 in;
PNG 600 dpi; PDF vector, fonttype 42. The build asserts pairwise text-overlap freedom, a
6 pt minimum type size and a blank 8 px margin.

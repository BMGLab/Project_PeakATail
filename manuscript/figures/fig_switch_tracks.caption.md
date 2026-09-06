# Figure — `fig_switch_tracks` caption

## Legend

**Replicated poly(A)-site switches shown as read coverage at the locus.** Four genes
selected under three rules. First, the gene's replicated significant sites must number exactly
two across all three stage comparisons, so that proximal and distal are unambiguous:
722 of the genes whose switch replicated in both mice qualify, and a gene with three
or more sites is excluded because no two of them can be labelled proximal and distal without
ambiguity. Second, both sites must lie inside the assigned gene with no other annotated gene
overlapping the drawn window, which leaves 259. Third, both sites must fall inside the
SAME 3'UTR of the SAME transcript, which is what makes an event tandem 3'UTR polyadenylation
rather than a choice between two last exons: 130 genes qualify. That third rule removed
Prdm15, drawn in an earlier version of this figure as the lengthening example, whose two sites
lie in no annotated 3'UTR at all. That second rule is not cosmetic:
poly(A)-site-to-gene assignment in overlapping loci is a known open defect (tool issue #99),
and it removed two genes that had passed every statistical filter, Map3k11, whose apparent
distal site is really the 3' end of the neighbouring Kcnk7, and Pom121, which overlaps Nsun5.
Either would have been drawn as a 3'UTR switch that is actually a between-gene artefact.
Applying the rule also moves the direction balance from 65% shortening to 98%,
so overlapping loci were inflating the apparent lengthening class. Three
shortening examples and one lengthening example are drawn, against the 98% shortening
/ 2% lengthening split of those 130 genes (127 against 3). Restricting to one
shared 3'UTR is what produces that asymmetry: before it the clean set ran 73% / 27%, so most
apparent lengthening was not tandem 3'UTR polyadenylation. Each column is one
mouse, scored independently end to end. Within a column, the three tracks are the three
spermatogenic stages (SPC, spermatocyte; RS, round spermatid; ES, elongating spermatid) and
share one vertical scale, printed beneath them, because the switch is a change in the shape
of the profile and a per-track autoscale would conceal it. Coverage is the mean number of
reads per base per cell, so tracks are comparable between stages containing different
numbers of cells (mouse 1 370/504/356 SPC/RS/ES, mouse 2 268/757/271 SPC/RS/ES). Read acceptance matches the caller: SAM flag filter 3844, the
tool's own CLIP_EXCLUDE_FLAGS, which drops unmapped, secondary, QC-failed, duplicate and
supplementary records, and only reads on the gene's own strand are drawn, which is how the
caller counts; 99.3-99.9% of reads over these four loci lie on that
strand in any case. Dashed lines mark the two called sites, green proximal and orange
distal. Beneath each column is the gene model: tall blocks are coding exons, low blocks
untranslated, and the line with chevrons is intronic and points in the direction of
transcription. The transcript carrying both sites is drawn dark with its shared 3'UTR in
green; other isoforms of the same gene are grey, so a reader can see directly that the two
sites sit within one 3'UTR rather than in the UTRs of competing isoforms. The percentages on each track are the quantity the switch test
actually used, the proximal and distal share of that gene's reads in that stage
(n_reads_pas / n_reads_gene from the differential tables), and are not measured off the
coverage drawn here.

Per gene: *Rragc* (chr4, + strand, shortening; proximal 123,935,899, distal 123,936,992; max q = 4 x 10^-05 across both sites and both mice): mouse 1 SPC 4/96, RS 33/64, ES 92/2; mouse 2 SPC 3/94, RS 33/64, ES 94/2. *Bpgm* (chr6, + strand, shortening; proximal 34,504,465, distal 34,505,126; max q = 4 x 10^-25 across both sites and both mice): mouse 1 SPC 9/60, RS 40/26, ES 56/6; mouse 2 SPC 13/82, RS 43/34, ES 65/5. *Atg12* (chr18, - strand, shortening; proximal 46,734,255, distal 46,732,424; max q = 1 x 10^-11 across both sites and both mice): mouse 1 SPC 7/70, RS 73/16, ES 61/28; mouse 2 SPC 9/46, RS 13/1, ES 13/2. *D1Ertd622e* (chr1, - strand, lengthening; proximal 97,643,900, distal 97,636,109; max q = 5 x 10^-04 across both sites and both mice): mouse 1 SPC 70/0, RS 56/3, ES 4/8; mouse 2 SPC 61/0, RS 74/2, ES 42/33.

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
`data/benchmark/gse104556/starsolo/Mouse{1,2}_scRNAseq`. Canvas 7.09 x 8.92 in;
PNG 600 dpi; PDF vector, fonttype 42. The build asserts pairwise text-overlap freedom, a
6 pt minimum type size and a blank 8 px margin.

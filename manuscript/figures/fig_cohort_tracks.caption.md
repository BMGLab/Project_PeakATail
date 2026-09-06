# Figure — `fig_cohort_tracks` caption

## Legend

**Cell-type poly(A)-site switches in a lung-adenocarcinoma cohort.** Four genes carrying a
switch between tumour epithelium and T cells that replicated across patients, selected by
the same two rules as the spermatogenesis companion figure: the gene's replicated sites
must number exactly two, so proximal and distal are unambiguous, and both must lie inside
the assigned gene with no other annotated gene in the drawn window. Of 248 genes with two
replicated sites moving in opposite directions in this cell-type pair, 52 pass the
annotation rule. **Left**, mean read depth per cell on the gene's strand, pooled over the
patients in which the gene was tested, one track per cell type; the three stacked elements
are the two tracks and a strip marking the two called sites, the direction of transcription
and the interval between the sites. Read acceptance matches the caller, SAM flag filter
3844, and only gene-strand reads are drawn. **Right**, the same switch resolved by patient:
each row is one patient, the two markers are that patient's proximal-site usage in tumour
epithelium and in T cells, and the connecting line is the difference the test scored.
Filled markers are patients that called the switch at q < 0.05; open markers are patients
that tested the site without calling it. Under the replication rule a tested-but-not-called
patient is not a disagreement, which is why the two are drawn apart rather than pooled, and
why the counts read 4 / 3 / 2 / 2 of 9, 9, 9 and 6 patients.

Per gene: *ZBTB38* (chr3, + strand, distal-up in tumour; proximal 141,442,485, distal 141,447,567): called by 4 of 9 patients tested, and all 4 of the called patients with both cell types measured agree in direction. *YWHAB* (chr20, + strand, distal-up in tumour; proximal 44,906,531, distal 44,908,386): called by 3 of 9 patients tested, and all 3 of the called patients with both cell types measured agree in direction. *CGGBP1* (chr3, - strand, distal-up in tumour; proximal 88,057,326, distal 88,051,949): called by 2 of 9 patients tested, and all 2 of the called patients with both cell types measured agree in direction. *SPAG1* (chr8, + strand, proximal-up in tumour; proximal 100,241,436, distal 100,252,087): called by 2 of 6 patients tested, and all 2 of the called patients with both cell types measured agree in direction.

**Caveats that travel with this figure.** These are four loci chosen to be legible, not a
random sample; the frequency claims are in Fig 6 and manuscript/20. The direction is
strongly one-sided: 49 of the 52 annotation-clean genes have the distal site
higher in tumour and only 3 the reverse, so *SPAG1* is drawn as the rarity it is and
not as half of a balance. Coverage height reflects both site usage and gene expression,
which is why the per-patient tested usage is drawn rather than inferred from peak height.
Pooled tracks weight patients by their read depth; the per-patient panel is the check on
that. Cell-type labels are expression-derived and curated, not sorted populations, and
patient counts differ per gene because a gene must be expressed in both cell types in a
patient to be tested there.

## Provenance

Coverage `results/figures/manuscript/cohort_tracks_coverage.tsv.gz`, metadata
`cohort_tracks_meta.tsv`, per-patient usage `cohort_tracks_usage.tsv`. Built by
`scripts/manuscript_figures/_cohort_track_extract.py` and
`scripts/manuscript_figures/fig_cohort_tracks.py`. Switch calls from
`stage3_laughney_v3` (code 9dfdefb, the v2-code record of manuscript/20); BAMs
`/mnt/ssd2/Laugney_Aligned`. Canvas 7.09 x 7.15 in; PNG 600 dpi; PDF vector.
The build asserts pairwise text-overlap freedom, a 6 pt minimum
type size and a blank 8 px margin.

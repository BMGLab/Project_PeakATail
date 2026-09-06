# Genome-browser plots of poly(A)-site switches — how to build them

How the spermatogenesis and lung-cohort browser figures were built, and how to point the
same pipeline at a different switch set. **The selection rules matter more than the
plotting code.**

All paths below are absolute. The reference project lives on the same machine, so a
companion project can read these scripts, the reference GTFs and the example inputs
directly — nothing needs copying.

- Reference project root: `/mnt/ssd1/Projects/PeakATail_wd`
- Scripts: `/mnt/ssd1/Projects/PeakATail_wd/scripts/manuscript_figures/`
- This guide: `/mnt/ssd1/Projects/PeakATail_wd/scripts/manuscript_figures/SWITCH_TRACK_PLOTS.md`

---

## 1. Shape: two scripts, one cached table between them

Reading a BAM is slow; laying out a figure takes twenty iterations. Split at the point
where cost changes.

```
# slow, run once per gene set
_switch_track_extract.py  ->  switch_tracks_coverage.tsv.gz   per-base depth, long form
                              switch_tracks_meta.tsv          one row per track: window, cells
                              switch_tracks_usage.tsv         the tested usage fractions

# fast, run every time you change the layout
fig_switch_tracks.py      ->  fig_switch_tracks.{png,pdf,caption.md}
```

Concretely, on this machine:

```bash
export LC_ALL=C OMP_NUM_THREADS=1
python3 /mnt/ssd1/Projects/PeakATail_wd/scripts/manuscript_figures/_switch_track_extract.py
python3 /mnt/ssd1/Projects/PeakATail_wd/scripts/manuscript_figures/fig_switch_tracks.py
```

Extractors write TSVs to `/mnt/ssd1/Projects/PeakATail_wd/results/figures/manuscript/`.
Plotters write PNG, PDF and a caption sidecar to
`/mnt/ssd1/Projects/PeakATail_wd/manuscript/figures/`.

Keep that boundary. The most common way these figures go wrong is a plotter that
recomputes from coverage something the statistical test defined differently.

## 2. Inputs

| Input | Why |
|---|---|
| Indexed BAM (`.bam` + `.bai`) | Coverage comes from region queries; without an index every panel costs a full scan |
| Barcode -> group TSV | Maps each cell barcode to the condition drawn as a track. Must use the barcodes the caller saw |
| Switch results (differential tables) | Which sites to draw AND the usage fractions printed on them |
| Replication table | Restricts to switches that held in an independent sample |
| Reference GTF | The annotation guard (section 3) |

**Set `LC_ALL=C` for every step.** A `tr_TR` locale silently misorders BED/interval data.

### Worked example: mouse spermatogenesis (3 ordered stages, 2 replicate animals)

```
BAM              /mnt/ssd1/Projects/PeakATail_wd/data/benchmark/gse104556/starsolo/Mouse1_scRNAseq/Aligned.sortedByCoord.out.bam
                 (and Mouse2_scRNAseq/… ; both indexed)
barcode -> stage /mnt/ssd1/Projects/PeakATail_wd/results/stage3_spermatogenesis_v2/mouse1/input/stage3_labels.tsv
                 columns: cb, stage  (SPC / RS / ES)
differential     /mnt/ssd1/Projects/PeakATail_wd/results/stage3_spermatogenesis_v2/mouse1/diff/true/differential/fisher_ES_vs_SPC.tsv
                 also fisher_ES_vs_RS.tsv, fisher_RS_vs_SPC.tsv
replication      /mnt/ssd1/Projects/PeakATail_wd/results/stage3_spermatogenesis_v2/summary/switch_replicated_pas_true.tsv
GTF              /mnt/ssd1/Projects/PeakATail_wd/data/references/mouse/Mus_musculus.GRCm38.102.gtf
```

### Worked example: lung-adenocarcinoma cohort (2 cell types, replication across patients)

```
cohort tree      /mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v3
BAMs             /mnt/ssd2/Laugney_Aligned/<SRR>_STAR/<SRR>_Aligned.sortedByCoord.out.bam
                 one GSM spans up to 8 SRR runs; pool them
SRR -> GSM map   /mnt/ssd1/Projects/PeakATail_wd/data/laughney/20250507_laughney_metadata.csv
                 column "Unnamed: 0" = SRR, "Sample.Name" = GSM
barcode -> type  /mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v3/switch/<GSM>/labelled/<GSM>.labels.tsv
                 index = <GSM>_<barcode>, column celltype
differential     /mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v3/switch/<GSM>/true/differential/fisher_<A>_vs_<B>.tsv
replication      /mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v3/replication/primary_noMetBone/pas_K2/replicated.tsv
GTF              /home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.99.gtf
```

Both differential-table families share one schema. The columns that matter:
`gene_id, chrom, start, strand, cluster1, cluster2, qvalue, n_reads_pas_cluster1,
n_reads_pas_cluster2, n_reads_gene_cluster1, n_reads_gene_cluster2, delta_proportion`.

## 3. The two selection rules

Statistics alone will hand you genes that look perfect and are wrong.

**Rule 1 — exactly two replicated sites per gene, across ALL comparisons.**
Only then do "proximal" and "distal" mean anything. A gene with three or more significant
sites cannot have two labelled that way without picking arbitrarily. Mouse set: keeps 722
of 5,361 genes. Apply per *gene*, not per comparison, or a gene with a third site in
another contrast slips through.

**Rule 2 — both sites inside the assigned gene, no other gene in the window.**
The assigned gene's span must contain both sites, and no other annotated gene may overlap
the drawn window including flanks. Leaves 259 of those 722.

### Why rule 2 exists

Two genes passed every statistical filter — replicated in both mice, q down to 5e-33,
large opposite-direction effects — and were drawn before the guard existed:

- **Map3k11**'s "distal site" sits inside **Kcnk7**, a different gene.
- **Pom121** overlaps **Nsun5**.

Both would have shipped as 3'UTR switches that are really artefacts of site-to-gene
assignment in overlapping loci (PeakATail issue #99). What exposed it was a strand
anomaly: only 12% of reads at the Map3k11 window were on its own strand, because an
antisense gene sits there. No amount of statistical significance would have caught this.

Applying the guard also moved the direction balance from 65% to 73% shortening —
overlapping loci were inflating the apparent lengthening class by about a third. The
guard is a result, not just hygiene.

Put it in the **extractor** as an assertion so a bad gene stops the run. Working
implementation at
`/mnt/ssd1/Projects/PeakATail_wd/scripts/manuscript_figures/_switch_track_extract.py`:

```python
def assert_annotation_clean(genes, gtf):
    """Both sites inside the assigned gene; no other gene in the window."""
    for g in genes:
        s = [g["prox"], g["dist"]]
        lo, hi = min(s) - FLANK, max(s) + FLANK
        own = G[G.gene_id == g["gene"]].iloc[0]
        assert own.start <= min(s) and max(s) <= own.end, g["symbol"]
        other = G[(G.chrom == g["chrom"]) & (G.end > lo)
                  & (G.start < hi) & (G.gene_id != g["gene"])]
        assert not len(other), (g["symbol"], f"{len(other)} other gene(s)")
```

## 4. Extraction — match the caller's read acceptance

| Decision | Setting, and why |
|---|---|
| Read filter | `samtools view -F 3844` — the tool's own `CLIP_EXCLUDE_FLAGS`: unmapped, secondary, QC-fail, duplicate, supplementary |
| Strand | Draw only gene-strand reads, because that is how the caller counts. Record both strands and check: 92-99.8% should sit on the gene strand. Near 12% means an antisense neighbour and a rule-2 failure |
| Barcode | Keep a read only if its `CB` tag is in that condition's cell set; strip any suffix so it matches the label file |
| Window | Sites span +/- 1500 bp flank. Do not truncate a wide locus — bin it |
| CIGAR | Advance and add depth on `M = X`; advance without depth on `D N`. Skipping `N` merges introns into coverage |

## 5. Plotting — four decisions that make or break it

| Decision | Do this |
|---|---|
| Normalise | Divide depth by the number of cells in that condition, or a condition with twice the cells looks twice as active |
| Share the scale | One y-max across a gene's conditions, printed on the figure. Per-track autoscale hides the switch, which *is* the change in shape |
| Print tested usage | `n_reads_pas / n_reads_gene` per condition from the differential table. Peak height confounds usage with expression — never measure the percentage off the coverage |
| Bin wide loci | Reduce to ~900 points by max-per-bin. Draws a 445 kb locus honestly at lower resolution rather than truncating |

**Lay the vertical geometry out in inches, not figure fractions.** Fractions do not
compose: three collisions and two box overflows shipped before this changed. Define
`fy(inches_from_top)` and `fh(inches)` against a single `FIG_H` and position through them.

## 6. Build asserts

Each was added after the defect it catches got through:

- **Text overlap** — pairwise bbox test over every text object, *including axes tick
  labels and axis titles*. Figure-level text alone misses the real collisions.
- **Box containment** — a label whose centre is inside a box must fit inside it with margin.
- **Type floor** — nothing below 6 pt.
- **Edge check** — the outer 8 px of the rendered PNG must be blank.
- **Canvas width** — pinned to 7.09 in (180 mm print column).

**Prove each assert fails.** Introduce the defect deliberately, confirm the build stops,
then revert. An assert that has never fired may not work.

## 7. Two atlas-free checks worth running

Curated panels show what a good call looks like, never what a typical one looks like. Run
these on a random unfiltered sample before trusting a call set. Implementations:
`/mnt/ssd1/Projects/PeakATail_wd/scripts/manuscript_figures/_site_step_scan.py` and
`/mnt/ssd1/Projects/PeakATail_wd/scripts/manuscript_figures/fig_switch_qc.py`.

**Site separation.** Distance between the outermost switching sites per gene. A tandem
3'UTR switch happens inside one 3'UTR, so mass beyond ~10 kb is something else. Mouse set:
48% of multi-site genes exceed 10 kb, 9% exceed 100 kb.

**Coverage step at the call.** Mean depth in the 200 bp inside the transcript over the
200 bp outside, strand-aware. A genuine poly(A) site in 3'-tag data drops off a cliff. On
300 random replicated sites: 68.9% drop at least threefold, 24% show no step, 16% are
inverted. That 68.9% uses no annotation atlas at all and lands beside the 0.706 atlas
agreement reported for the default call set — two independent measures converging.

## 8. Starter prompt for Claude Code

Paste into a session opened in the companion project, filling in the five input paths.

```
Build a genome-browser figure of replicated poly(A)-site switches, following
the recipe at
/mnt/ssd1/Projects/PeakATail_wd/scripts/manuscript_figures/SWITCH_TRACK_PLOTS.md

Read these as reference implementations (same machine, absolute paths):
  /mnt/ssd1/Projects/PeakATail_wd/scripts/manuscript_figures/_switch_track_extract.py
  /mnt/ssd1/Projects/PeakATail_wd/scripts/manuscript_figures/fig_switch_tracks.py
  /mnt/ssd1/Projects/PeakATail_wd/scripts/manuscript_figures/_pubstyle.py

My inputs:
  BAM (indexed)      <path>
  barcode -> group   <path>   column for the condition drawn as a track
  switch results     <path>   differential tables with n_reads_pas / n_reads_gene
  replication table  <path>
  reference GTF      <path>

Split into two scripts: an extractor that walks the BAM once and writes
coverage / meta / usage TSVs, and a plotter that reads only those TSVs.

Select genes by two rules, both asserted in the extractor:
  1. exactly two replicated sites for the gene across ALL comparisons,
     so proximal and distal are unambiguous;
  2. both sites inside the assigned gene, with no other annotated gene
     overlapping the drawn window.
Report how many genes survive each rule, and the direction balance
before and after rule 2.

Extraction: samtools view -F 3844; keep only reads on the gene's strand;
keep only barcodes in the condition's cell set; +/-1500 bp flank; handle
CIGAR D and N as gaps. Report the fraction of reads on the gene strand
per gene and flag anything below 80%.

Plot: one row per gene, one track per condition, depth normalised per
cell, ONE shared y-scale per gene printed on the figure, dashed markers
at the two sites. Print the tested usage (n_reads_pas / n_reads_gene)
from the differential table — never measure it off the coverage.
Bin wide loci to ~900 points by max-per-bin; never truncate.

Lay the vertical geometry out in INCHES via fy()/fh(), not figure
fractions. Assert: no text-text overlap (including axes tick labels),
no label spilling its box, nothing under 6 pt, blank outer 8 px of the
PNG, canvas width 7.09 in. Prove each assert fails on a deliberate
regression before trusting it.

Set LC_ALL=C. Write a caption sidecar naming every number's source.
```

Adapt one thing per dataset: **what a track is.** Ordered stages read as a developmental
axis; two cell types read as a contrast; patients or donors are replication and belong
beside the tracks rather than as more of them.

## 9. Reference implementations

| File (absolute) | What it shows |
|---|---|
| `/mnt/ssd1/Projects/PeakATail_wd/scripts/manuscript_figures/_switch_track_extract.py`<br>`/mnt/ssd1/Projects/PeakATail_wd/scripts/manuscript_figures/fig_switch_tracks.py` | Ordered conditions (three stages), two independent replicates as columns. Cleanest template to copy |
| `/mnt/ssd1/Projects/PeakATail_wd/scripts/manuscript_figures/_cohort_track_extract.py`<br>`/mnt/ssd1/Projects/PeakATail_wd/scripts/manuscript_figures/fig_cohort_tracks.py` | Two cell types pooled across patients, per-patient replication drawn beside the tracks. Copy when replication is across samples, and for pooling several BAMs per sample |
| `/mnt/ssd1/Projects/PeakATail_wd/scripts/manuscript_figures/_random_switch_extract.py`<br>`/mnt/ssd1/Projects/PeakATail_wd/scripts/manuscript_figures/fig_random_switches.py` | Random unfiltered sample, with max-per-bin downsampling for wide loci. Copy to audit a call set rather than illustrate it |
| `/mnt/ssd1/Projects/PeakATail_wd/scripts/manuscript_figures/_site_step_scan.py`<br>`/mnt/ssd1/Projects/PeakATail_wd/scripts/manuscript_figures/fig_switch_qc.py` | The two atlas-free checks |
| `/mnt/ssd1/Projects/PeakATail_wd/scripts/manuscript_figures/_pubstyle.py` | Shared palette, type scale, rcParams. Import rather than restating colours |

Rendered figures for reference:
`/mnt/ssd1/Projects/PeakATail_wd/manuscript/figures/fig_switch_tracks.png`,
`fig_cohort_tracks.png`, `fig_random_switches.png`, `fig_switch_qc.png`, each with a
`.caption.md` sidecar naming the source of every number on the canvas.

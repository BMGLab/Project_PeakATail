# Sierra — benchmark tool status

**State: INSTALLED and smoke-tested (2026-08-18)**

## Citation
Patrick R, Humphreys DT, Janbandhu V, Oshlack A, Ho JWK, Harvey RP, Lo KK.
*Sierra: discovery of differential transcript usage from polyA-captured single-cell RNA sequencing data.*
Genome Biology, 2020, 21:167. (First author: Ralph Patrick; VCCRI.)
Repo: https://github.com/VCCRI/Sierra

## Versions installed
- Sierra **0.99.27** (installed from GitHub master via remotes)
- R 4.5.3 (conda-forge r-base inside env)
- regtools **1.0.0** (bioconda; needed to make the junctions BED that FindPeaks requires)
- Env: conda env **bench_sierra** (`/home/biolab/miniconda3/envs/bench_sierra`)

## Install commands that worked (in order)
```bash
CONDA=/home/biolab/miniconda3/bin/conda
$CONDA create  -n bench_sierra --solver=libmamba --override-channels -c conda-forge -c bioconda -y r-base r-remotes
$CONDA install -n bench_sierra --solver=libmamba --override-channels -c conda-forge -c bioconda -y \
    bioconductor-genomicranges bioconductor-genomicalignments bioconductor-rtracklayer \
    bioconductor-rsamtools bioconductor-biostrings bioconductor-bsgenome bioconductor-genomicfeatures \
    bioconductor-summarizedexperiment bioconductor-singlecellexperiment bioconductor-biocparallel \
    bioconductor-genefilter
$CONDA install -n bench_sierra --solver=libmamba --override-channels -c conda-forge -c bioconda -y \
    bioconductor-dexseq bioconductor-gviz bioconductor-biocstyle
$CONDA install -n bench_sierra --solver=libmamba --override-channels -c conda-forge -c bioconda -y \
    r-reshape2 r-plyr r-dplyr r-foreach r-doparallel r-matrix r-mlmetrics r-progress \
    r-ggplot2 r-cowplot r-flock r-magrittr r-r.utils r-data.table r-scales
/home/biolab/miniconda3/envs/bench_sierra/bin/Rscript -e \
    'remotes::install_github("VCCRI/Sierra", upgrade="never", build_vignettes=FALSE)'
$CONDA install -n bench_sierra --solver=libmamba --override-channels -c conda-forge -c bioconda -y regtools
```
Note: Seurat is only in Suggests and was NOT installed (only needed for `NewPeakSeurat`).
To add later: `$CONDA install -n bench_sierra --solver=libmamba --override-channels -c conda-forge -c bioconda -y r-seurat`

## Smoke-test evidence
```
$ Rscript -e 'library(Sierra); packageVersion("Sierra"); exists("FindPeaks")'
Sierra version: 0.99.27
FindPeaks exists: TRUE
```
(Load emits 3 benign warnings: `replacing previous import 'GenomicRanges::setdiff' by 'dplyr::setdiff'` etc.)
Key exports: FindPeaks, CountPeaks, MergePeakCoordinates, AggregatePeakCounts, AnnotatePeaksFromGTF,
ReadPeakCounts, NewPeakSeurat, NewPeakSCE, GetExpressedPeaks, SelectGenePeaks, PeakSeuratFromTransfer.
Also functionally verified `regtools junctions extract` on
`data/testdata/downsampled_aligned.bam` -> 329,972 BED12 junction records.

## Inputs Sierra needs
1. **BAM**: coordinate-sorted, indexed, with cell-barcode and UMI tags. Defaults `CBtag="CB"`, `UMItag="UB"`
   = CellRanger conventions, so the 10x pbmc_10k_v3 `possorted_genome_bam.bam` works as-is (tags configurable
   for Drop-seq etc.).
2. **GTF**: the same annotation used for alignment (for pbmc_10k_v3: CellRanger GRCh38 reference GTF).
3. **Junctions file** (FindPeaks only): either a STAR `SJ.out.tab(.gz)` or a regtools BED12
   (`.bed(.gz)`), auto-detected by filename; read with plain `read.table` — bgzip/tabix NOT required.
   Make it with: `regtools junctions extract -s XS -o junctions.bed possorted_genome_bam.bam`.
   (Strand may come out `?` if the BAM lacks XS attrs; harmless — FindPeaks only uses junction
   coordinates, block sizes, and the count column.)
4. **Barcode whitelist** (CountPeaks only): one barcode per line, matching the CB tag values —
   use CellRanger's `filtered_feature_bc_matrix/barcodes.tsv.gz` (keep the `-1` suffix; a gzipped file is fine).

## Outputs
- **FindPeaks -> peak/PAS coordinate table** (tab-separated, one header line):
  `Gene  Chr  Strand  MaxPosition  Fit.max.pos  Fit.start  Fit.end  mu  sigma  k  exon/intron  exon.pos  LogLik  polyA_ID`
  where `polyA_ID = Gene:Chr:Fit.start-Fit.end:Strand`. Peaks are Gaussian fits (NLS) to read
  coverage; `Fit.max.pos`/`MaxPosition` is the peak summit — the best single-nt proxy for the PAS.
- **CountPeaks -> peak x cell UMI count matrix** in 10x-style MatrixMarket triplet form inside
  `output.dir`: `matrix.mtx.gz`, `sitenames.tsv.gz` (polyA_IDs), `barcodes.tsv.gz`.
  Load back with `ReadPeakCounts(data.dir=...)`.
- Downstream (not needed for PAS benchmarking): AnnotatePeaksFromGTF, DUTest (DEXSeq-based DTU).

## Exact invocation for the 10x pbmc_10k_v3 BAM
```bash
ENV=/home/biolab/miniconda3/envs/bench_sierra
BAM=/path/to/pbmc_10k_v3/possorted_genome_bam.bam
GTF=/path/to/refdata-gex-GRCh38/genes/genes.gtf
WL=/path/to/pbmc_10k_v3/filtered_feature_bc_matrix/barcodes.tsv.gz
OUT=/mnt/ssd1/Projects/PeakATail_wd/results/benchmarks/sierra_pbmc10k

mkdir -p $OUT
$ENV/bin/regtools junctions extract -s XS -o $OUT/junctions.bed $BAM

$ENV/bin/Rscript -e '
library(Sierra)
FindPeaks(output.file   = "'$OUT'/sierra_peaks.txt",
          gtf.file      = "'$GTF'",
          bamfile       = "'$BAM'",
          junctions.file= "'$OUT'/junctions.bed",
          ncores        = 8)
CountPeaks(peak.sites.file = "'$OUT'/sierra_peaks.txt",
           gtf.file        = "'$GTF'",
           bamfile         = "'$BAM'",
           whitelist.file  = "'$WL'",
           output.dir      = "'$OUT'/counts",
           countUMI = TRUE, ncores = 8)'
```
Full FindPeaks signature (defaults): min.jcutoff=50, min.jcutoff.prop=0.05, min.cov.cutoff=500,
min.cov.prop=0.05, min.peak.cutoff=200, min.peak.prop=0.05, fit.method="NLS", ncores=1.
Runtime note (paper/wiki): whole-BAM FindPeaks on a 10k-cell dataset is hours-scale; use ncores
and expect this to be the slow step of the benchmark.

## Dataset / chemistry restrictions
- Designed for **polyA-captured 3'-tag scRNA-seq**: 10x Chromium 3' **v2 and v3 both supported**
  (paper demonstrates on 10x 3'). Peak model assumes read pileups at polyA-proximal ends.
- Operates purely on the **aligned cDNA read (R2) in the BAM** — no access to raw read1/FASTQ needed;
  barcodes/UMIs must already be in BAM tags (CellRanger does this).
- 10x **5' chemistry is not appropriate** (peak model assumes 3' polyA capture).
- Non-10x platforms work if BAM has barcode/UMI tags (set `CBtag`/`UMItag` in CountPeaks).

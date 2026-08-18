# scAPAtrap — install status

**State: installed** (2026-08-18)

## Citation
Wu X, Liu T, Ye C, Ye W, Ji G. *scAPAtrap: identification and quantification of alternative
polyadenylation sites from single-cell RNA-seq data.* **Briefings in Bioinformatics** 2021;
22(4):bbaa273. First author: Xiaohui Wu. GitHub: BMILAB/scAPAtrap (installed HEAD
12ee8f5df34fbe81642fe8564132e3f847f93132, package version 0.2.0, dated 2023-09-10).

## Install commands that worked (in order)
```bash
# 1. conda env with R + all Bioconductor deps + external binaries (one solve, ~4 min)
conda create -n bench_scapatrap -y --override-channels -c conda-forge -c bioconda \
  r-base r-remotes r-reshape2 r-dplyr r-matrix r-magrittr \
  bioconductor-genomicranges bioconductor-iranges bioconductor-rsamtools \
  bioconductor-genomicalignments bioconductor-derfinder bioconductor-bumphunter \
  bioconductor-regioner samtools subread

# 2. scAPAtrap itself from GitHub (pure R, NeedsCompilation: no)
/home/biolab/miniconda3/envs/bench_scapatrap/bin/Rscript -e \
  'remotes::install_github("BMILAB/scAPAtrap", upgrade="never")'

# 3. remaining external tools the pipeline shells out to
conda install -n bench_scapatrap -y --override-channels -c conda-forge -c bioconda umi_tools star
```

## Versions installed
| component | version |
|---|---|
| scAPAtrap (R pkg) | 0.2.0 |
| R | 4.5.3 |
| samtools | 1.24 |
| featureCounts (subread) | v2.1.1 |
| umi_tools | 1.1.6 |
| STAR | 2.5.2b (old solver pick; only used when starting from FASTQ — NOT needed for a CellRanger BAM) |

Env prefix: `/home/biolab/miniconda3/envs/bench_scapatrap`

## Smoke-test evidence
```
$ Rscript -e 'library(scAPAtrap); packageVersion("scAPAtrap")'
scAPAtrap 0.2.0 loaded OK
exports: TRAP.PARAMS, TRAPFILES, convertAPAtrapData, countPeaks, dedupByPos, extractBcAndUb,
findPeaks, findPeaksByStrand, findTails, findTailsByPeaks, findUniqueMap, generateAlignBam,
generateFinalBam, generateRefIndex, generateSAF, generatescExpMa, initScAPAtrap, loadBpCoverages,
reducePeaks, scAPAtrap, separateBamBystrand, setTools, setTrapParams

$ setTools(list(samtools=..., umitools=..., featureCounts=..., star=...), check=TRUE)
setTools check=TRUE passed   # validates all four binaries in the env
```

## Input format
- Coordinate-sorted, indexed 10x BAM **with CB/UB tags** (CellRanger `possorted_genome_bam.bam`
  works directly; `TenX=TRUE` makes umi_tools use `--extract-umi-method=tag --umi-tag=UB --cell-tag=CB`).
- Optional cell-barcode whitelist (CellRanger filtered `barcodes.tsv`, strip the `-1` suffix).
- `trap.params$chrs` must match BAM chromosome naming (CellRanger GRCh38 uses `chr1`-style).
- `trap.params$readlength` must be set to the actual mapped read length (pbmc_10k_v3 R2 = 91 nt;
  the default is 49 — wrong for v3 data, must override).
- Pipeline internally: findUniqueMap (samtools MAPQ-255 filter) -> dedupByPos (umi_tools dedup
  per cell) -> separateBamBystrand -> loadBpCoverages (derfinder) -> findPeaks -> generateSAF ->
  countPeaks (featureCounts + umi_tools count) -> findTails (soft-clipped poly(A) evidence) ->
  generatescExpMa.

## Output format
- `<outputDir>/peaks.saf` — peak/PAS coordinates in SAF format: GeneID, Chr, Start, End, Strand
  (PAS = 3' end of peak; strand-aware).
- `<outputDir>/counts.tsv.gz` — UMI counts per peak per cell (long format: gene(=peakID), cell, count).
- `<outputDir>/scAPAtrapData.rda` — list with `peaks.meta` (data.frame of peak coordinates,
  poly(A)-tail support flag) and `peaks.count` (sparse peak x cell matrix). So: **both PAS/peak
  coordinates AND a peak-by-cell counts matrix**. `convertAPAtrapData` converts to
  Seurat / SingleCellExperiment / movAPA PACdataset.
- Many intermediate BAMs land next to the input (see `TRAPFILES()` after a run); outputDir must
  NOT already exist or the wrapper errors.

## Dataset / chemistry restrictions
- Designed for 3'-enriched scRNA-seq: 10x Chromium **v2 and v3 both supported** (`TenX=TRUE`);
  also usable for CEL-seq/Drop-seq-like data with per-position dedup (`TenX=FALSE`).
- Uses only the genomic read (R2 for 10x); no R1 requirement beyond CB/UB already being in tags.
- `findTails` needs reads that run into the poly(A) tail (soft-clipped A's) — works on 10x v3
  91-nt R2; tail evidence is optional (`tails.search='no'` default) and used only to flag/filter peaks.
- No genome/annotation files needed when starting from a BAM (annotation-free peak calling).

## Exact invocation for the 10x pbmc_10k_v3 CellRanger BAM
```r
# run with: conda run -n bench_scapatrap Rscript run_scapatrap_pbmc10k.R   (LC_ALL=C exported)
library(scAPAtrap)
E <- "/home/biolab/miniconda3/envs/bench_scapatrap/bin"
tools <- setTools(list(samtools      = file.path(E, "samtools"),
                       umitools      = file.path(E, "umi_tools"),
                       featureCounts = file.path(E, "featureCounts"),
                       star          = file.path(E, "STAR")))
tp <- TRAP.PARAMS()
tp$TenX       <- TRUE
tp$readlength <- 91                        # pbmc_10k_v3 R2 length
tp$chrs       <- paste0("chr", c(1:22, "X", "Y"))   # CellRanger GRCh38 naming
tp$thread     <- 12
bc            <- read.delim("barcodes.tsv", header = FALSE)$V1   # filtered barcodes
tp$barcode    <- gsub("-[0-9]+$", "", bc)
scAPAtrap(tools, tp,
          inputBam  = "/path/to/pbmc_10k_v3_possorted_genome_bam.bam",
          outputDir = "scapatrap_pbmc10k",          # must not pre-exist
          logf      = "scapatrap_pbmc10k.log")
# -> scapatrap_pbmc10k/{peaks.saf, counts.tsv.gz, scAPAtrapData.rda}
```

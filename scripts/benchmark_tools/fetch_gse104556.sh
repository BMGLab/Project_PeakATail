#!/usr/bin/env bash
# fetch_gse104556.sh - download benchmark dataset 2: GSE104556 (mouse testis scRNA-seq)
#
# GSE104556 = "A broad single-cell transcriptome view of the male mouse germ line"
#   Lukassen et al. 2018, PMID 30204153. Two biological replicates, ~2550 cells total,
#   10x Genomics Chromium Single Cell 3' v2 chemistry, Illumina HiSeq 2500, mm10.
# SRA study: SRP119327 (BioProject PRJNA413049)
#   SRR6129050  GSM2803334  Mouse1_scRNAseq  PAIRED  199,284,287 read pairs
#   SRR6129051  GSM2803335  Mouse2_scRNAseq  PAIRED  198,383,896 read pairs
#
# IMPORTANT - CellRanger note:
#   GEO/SRA do not provide a BAM for this dataset, only reads (+ author-processed
#   barcodes/genes/matrix on the GEO series page). To benchmark PAS callers that need
#   CB/UB tags, the FASTQs must be re-processed with CellRanger (`cellranger count`)
#   against the 10x mm10 reference (refdata-gex-mm10-2020-A). v2 chemistry:
#   R1 = 26 bp (16 bp cell barcode + 10 bp UMI), R2 = cDNA. fasterq-dump must be run
#   with --include-technical --split-files so the barcode read is kept, and files must
#   be renamed to the CellRanger convention <Sample>_S1_L001_R1_001.fastq.gz etc.
#   The resulting BAM uses mm10/GRCm38 with 'chr'-prefixed names (chr1..chr19, chrX,
#   chrY, chrM) - same 'chr' vs Ensembl naming mismatch as the PBMC dataset: translate
#   before intersecting with PolyASite mouse atlas (bare 1..19, X, Y, MT).
#
# Prereqs: sra-tools >= 3.x (prefetch, fasterq-dump) e.g. via
#   conda create -y -p /mnt/ssd1/Projects/PeakATail_wd/tools/bench_envs/sratools -c bioconda -c conda-forge sra-tools pigz
# Disk: ~110 GB FASTQ + ~100 GB SRA cache; CellRanger run needs ~64 GB RAM.
#
# This script only DOWNLOADS and prepares FASTQs; the cellranger invocation is printed
# at the end (run it manually - it is multi-hour).

set -euo pipefail
export LC_ALL=C   # tr_TR locale on this machine breaks sort

OUTDIR=/mnt/ssd1/Projects/PeakATail_wd/data/benchmark/gse104556
SRRS=(SRR6129050 SRR6129051)
SAMPLES=(Mouse1_scRNAseq Mouse2_scRNAseq)
THREADS=8

mkdir -p "$OUTDIR"/{sra,fastq}
cd "$OUTDIR"

# 1) prefetch the .sra objects (resumable; --max-size raised for ~20 GB runs)
for srr in "${SRRS[@]}"; do
    prefetch --max-size 100G -O "$OUTDIR/sra" "$srr"
done

# 2) dump FASTQs, KEEPING the technical barcode read (R1) - required for 10x v2
for i in "${!SRRS[@]}"; do
    srr=${SRRS[$i]}
    fasterq-dump --split-files --include-technical -e "$THREADS" \
        -O "$OUTDIR/fastq" "$OUTDIR/sra/$srr/$srr.sra"
done

# 3) rename to CellRanger convention and compress.
#    Expected per-run outputs: <SRR>_1.fastq (26 bp: CB+UMI -> R1) and
#    <SRR>_2.fastq (cDNA -> R2). VERIFY read lengths before renaming:
#      head -2 fastq/SRR6129050_1.fastq | awk 'NR==2{print length($0)}'   # expect 26
#    If a run yields 3 files, the extra short one is the sample index (-> I1).
for i in "${!SRRS[@]}"; do
    srr=${SRRS[$i]}; sample=${SAMPLES[$i]}
    mv "fastq/${srr}_1.fastq" "fastq/${sample}_S1_L001_R1_001.fastq"
    mv "fastq/${srr}_2.fastq" "fastq/${sample}_S1_L001_R2_001.fastq"
    pigz -p "$THREADS" "fastq/${sample}_S1_L001_R1_001.fastq" "fastq/${sample}_S1_L001_R2_001.fastq"
done

# 4) (optional) author-processed matrix from GEO, for the cell-barcode whitelist
wget -c -P "$OUTDIR" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104556/suppl/GSE104556_barcodes.tsv.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104556/suppl/GSE104556_genes.tsv.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104556/suppl/GSE104556_matrix.mtx.gz"

cat <<'MSG'
Download done. Next (manual, multi-hour) - produce CB/UB-tagged BAMs with CellRanger:

  # reference (once):
  #   wget -c https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz && tar xzf ...
  cellranger count --id=gse104556_mouse1 \
      --transcriptome=/path/to/refdata-gex-mm10-2020-A \
      --fastqs=/mnt/ssd1/Projects/PeakATail_wd/data/benchmark/gse104556/fastq \
      --sample=Mouse1_scRNAseq --chemistry=SC3Pv2 --localcores=8 --localmem=64
  # repeat with --id=gse104556_mouse2 --sample=Mouse2_scRNAseq
  # BAM appears at <id>/outs/possorted_genome_bam.bam with CB/UB tags.
MSG

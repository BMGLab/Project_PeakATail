# SCAPTURE — install status

**State: INSTALLED, smoke-tested OK** (2026-08-18)

## Citation
Li, Guo-Wei, et al. "SCAPTURE: a deep learning-embedded pipeline that captures polyadenylation
information from 3' tag-based RNA-seq of single cells." **Genome Biology** 22, 221 (**2021**).
(YangLab/SCAPTURE; first author Guo-Wei Li; GPLv3 for academic use.)

## Where things live
- Repo clone: `/mnt/ssd1/Projects/PeakATail_wd/tools/SCAPTURE` (shallow clone of https://github.com/YangLab/SCAPTURE, main; scripts chmod +x done)
- Conda env: `bench_scapture` (`/home/biolab/miniconda3/envs/bench_scapture`, ~4.5 GB)
- Version: SCAPTURE v1.0 (2021/01/25, only released version; git main @ 2026-08-18)

## Install commands that worked
```bash
cd /mnt/ssd1/Projects/PeakATail_wd/tools && git clone --depth 1 https://github.com/YangLab/SCAPTURE.git
/home/biolab/miniconda3/bin/conda tos accept --override-channels \
    -c https://repo.anaconda.com/pkgs/main -c https://repo.anaconda.com/pkgs/r
cd SCAPTURE && CONDA_CHANNEL_PRIORITY=flexible \
    /home/biolab/miniconda3/bin/conda env create -n bench_scapture -f SCAPTURE_env.yml --solver=libmamba
# featureCounts is a documented requirement but is NOT in their yml — added:
CONDA_CHANNEL_PRIORITY=flexible /home/biolab/miniconda3/bin/conda install -n bench_scapture \
    -c bioconda -c conda-forge --solver=libmamba -y 'subread=1.6.4'
chmod +x /mnt/ssd1/Projects/PeakATail_wd/tools/SCAPTURE/scapture*
```
Gotchas hit: (1) machine conda config has `channel_priority: strict`, which makes the fully-pinned
yml unsolvable ("excluded by strict repo priority") — must override with `CONDA_CHANNEL_PRIORITY=flexible`;
(2) Anaconda TOS had to be accepted once for pkgs/main + pkgs/r (done, persists).

## Env contents (relevant pins)
python 3.7.8, R 3.6.3, samtools 1.9, bedtools 2.26.0, stringtie 2.1.4, umi_tools 1.1.0,
featureCounts 1.6.4 (added), tensorflow-gpu 2.0.0 (pip) + cudatoolkit 10.0/cudnn 7.6.5,
homer, dropseq_tools, UCSC utils, GNU parallel.

## Smoke-test evidence
```bash
source /home/biolab/miniconda3/etc/profile.d/conda.sh && conda activate bench_scapture
export PATH=/mnt/ssd1/Projects/PeakATail_wd/tools/SCAPTURE:$PATH
scapture -h          # prints full usage (annotation|PAScall|PASmerge|PASquant), exit 1 on no args
python -c "import tensorflow as tf; print(tf.__version__)"   # -> 2.0.0
featureCounts -v     # -> featureCounts v1.6.4
```
DeepPASS end-to-end check: `CUDA_VISIBLE_DEVICES="" python DeepPASS/Predict.py -m DeepPASS/best_model.h5
-p <3 random 200nt seqs> -o out/` loaded the h5 model and wrote `Predict_Result.txt`
(all 3 random seqs -> class 0/negative, as expected). Ran on CPU.
**GPU caveat**: box has RTX A4000 (Ampere, CC 8.6); env pins CUDA 10.0 which predates Ampere.
Run PAScall with `CUDA_VISIBLE_DEVICES=""` to force CPU for the DeepPASS step (it is fast; fine on CPU).

## Input format
- **BAM**: Cell Ranger possorted BAM with `CB` (cell barcode) and `UB` (UMI) tags. No FASTQ needed.
- **References**: genome .fa (+ .fai), chrom sizes file, and a GTF whose 9th column carries
  `gene_name` and `gene_type` attributes (GENCODE-style; hard-coded perl regex on `gene_type`).
- **Chemistry restrictions**: any 3' tag-based scRNA-seq with CB/UB in BAM; paper used 10x v2 and v3
  PBMCs — both fine. Only R2 (cDNA) matters; pass its length as `-l` (pbmc_10k_v3 R2 = 91 nt).
  DeepPASS models provided for human + mouse only (`--species human|mouse`).
- **Reference-naming caveat for our pbmc_10k_v3 BAM**: BAM uses chr-less Ensembl names ("1","2",...).
  Matching refs already on disk: `/mnt/ssd1/Projects/PeakATail_wd/data/references/humanSTARindex/Homo_sapiens.GRCh38.dna.primary_assembly.fa`
  (+.fai) and `.../humanSTARindex/Homo_sapiens.GRCh38.99.gtf`, chromsize `/mnt/ssd1/Projects/PeakATail_wd/data/references/chrom.sizes.nochr.filt`.
  BUT the Ensembl GTF uses `gene_biotype`, not `gene_type` — fix with
  `sed 's/gene_biotype/gene_type/g' Homo_sapiens.GRCh38.99.gtf > scapture.gtf`
  (or use GENCODE GTF with `sed 's/^chr//'`). This GTF tweak is the one untested step.

## Output format
- **PAScall** (per sample): BED12+3 peak files — `<prefix>.{exonic,intronic,3primeExtended}.peaks.evaluated.bed`;
  cols 1-12 = spliced PAS peak region (BED12), col 13 = # known poly(A)-DB sites supporting the PAS,
  col 14 = DeepPASS prediction ("positive"/"negative"), col 15 = ±100 nt sequence around cleavage site.
  So: **PAS/peak coordinates = BED12, strand-aware; distal end of the peak = cleavage site region.**
- **PASmerge**: `<prefix>.Integrated.bed` (same BED12+3 schema, samples merged).
- **PASquant**: `<prefix>.KeepCell.UMIs.tsv.gz` — **PAS x cell UMI counts matrix**
  (peak IDs are gene-embedded, e.g. GENE-P1/P2 -> usable for APA), plus optional per-celltype bigWigs.

## Exact invocation planned for pbmc_10k_v3
```bash
source /home/biolab/miniconda3/etc/profile.d/conda.sh && conda activate bench_scapture
export PATH=/mnt/ssd1/Projects/PeakATail_wd/tools/SCAPTURE:$PATH   # scapture finds DeepPASS via `which scapture`
export LC_ALL=C
REF=/mnt/ssd1/Projects/PeakATail_wd/data/references
BAM=/mnt/ssd1/Projects/PeakATail_wd/data/benchmark/pbmc_10k_v3/pbmc_10k_v3_possorted_genome_bam.bam
sed 's/gene_biotype/gene_type/g' $REF/humanSTARindex/Homo_sapiens.GRCh38.99.gtf > scapture.gtf
scapture -m annotation -o scapture_anno -g $REF/humanSTARindex/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --gtf scapture.gtf --cs $REF/chrom.sizes.nochr.filt --extend 2000
CUDA_VISIBLE_DEVICES="" scapture -m PAScall -a scapture_anno -o pbmc10k \
    -g $REF/humanSTARindex/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    -b $BAM -l 91 -p 16 --species human
# select confident PAS (single sample, so skip PASmerge):
perl -alne '$,="\t";print @F[0..11] if $F[13] eq "positive";' \
    pbmc10k.exonic.peaks.evaluated.bed pbmc10k.intronic.peaks.evaluated.bed > pbmc10k.PASquant.bed
scapture -m PASquant -b $BAM --pas pbmc10k.PASquant.bed \
    --celllist barcodes.tsv -o pbmc10k.PASquant -p 16   # barcodes from pbmc_10k_v3_filtered_feature_bc_matrix
```
Note (as of writing): the pbmc_10k_v3 BAM in data/benchmark/ was still mid-download by another
process (~6%, wget log live) — verify `samtools quickcheck` before running.

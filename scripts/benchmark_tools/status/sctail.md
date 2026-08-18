# scTail — install & smoke-test status

**State: INSTALLED & SMOKE-TESTED** (2026-08-18)

## Citation
Hou R, Huang Y. *scTail: precise polyadenylation site detection and its alternative
usage analysis from reads 1 preserved 3' scRNA-seq data.* **Genome Biology**, 2025.
doi:10.1186/s13059-025-03710-7 (preprint: bioRxiv 2024.07.05.602174).
Repo: https://github.com/StatBiomed/scTail — docs: https://sctail.readthedocs.io

## Installed version / environment
- scTail **0.1.8** (PyPI) in venv `/mnt/ssd1/Projects/PeakATail_wd/tools/bench_envs/sctail`
  (Python 3.10.12 from /usr/bin/python3.10; env size ~1.4 GB)
- Key deps: torch 2.13.0+cpu (CPU-only wheel, deliberate — CNN inference falls back to
  CPU automatically), numpy 1.26.4, pandas 2.3.3, pysam 0.24.0, pyranges 0.1.4,
  anndata 0.11.4, kipoiseq 0.7.1, attrs 21.4.0 (upstream pins attrs<=21.4.0), scipy 1.15.3,
  scikit-learn 1.7.2, pyfaidx 0.9.0.4. `pip check`: no broken requirements.
- Pretrained CNN PAS classifiers ship inside the wheel:
  `site-packages/scTail/model/{human,mouse}_pretrained_model.pth`

## Install commands that worked (chunked, each well under 8 min)
```bash
/usr/bin/python3.10 -m venv /mnt/ssd1/Projects/PeakATail_wd/tools/bench_envs/sctail
V=/mnt/ssd1/Projects/PeakATail_wd/tools/bench_envs/sctail/bin
$V/pip install --upgrade pip setuptools wheel
$V/pip install torch --index-url https://download.pytorch.org/whl/cpu
$V/pip install "numpy<2.0" "pandas>=2.1.1" "scipy>=1.12.0" "scikit-learn>=1.3.2" "matplotlib>=3.8.1" "tqdm>=4.66.1"
$V/pip install "pysam>=0.22.0" "pyranges>=0.0.129" "anndata>=0.10.3" "pyfaidx>=0.7.2.2" "attrs<=21.4.0"
$V/pip install "kipoiseq>=0.7.1"
$V/pip install scTail
```

## Smoke-test evidence
- `scTail --help` prints: "Welcome to scTail v0.1.8! Command lines available:
  scTail-callPeak / scTail-peakMerge / scTail-count"
- `scTail-callPeak --help`, `scTail-peakMerge --help`, `scTail-count --help` all print
  full optparse usage (verified 2026-08-18).
- `python -c "import scTail, torch"` → "scTail import OK, version 0.1.8",
  torch tensor allocates on CPU.

## Input the tool needs (READ CAREFULLY — this is the benchmark blocker)
scTail's whole premise is calling PAS from **aligned READ 1** ("reads 1 preserved"
3' scRNA-seq): read1 must be **>100 bp long (ideally 150/151 bp)** so that after the
16 bp CB + 10/12 bp UMI + ~32 bp poly(dT)VN there is genomic cDNA immediately
upstream of the PAS. Required preprocessing is **STAR/STARsolo aligning read1 as the
barcode mate**, e.g. (docs, 10x v2):

```
STAR --soloType CB_UMI_Simple --soloBarcodeMate 1 --clip5pNbases 58 0 \
     --soloStrand Reverse --readFilesIn reads_1.fastq.gz reads_2.fastq.gz \
     --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 \
     --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
     --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate \
     --soloCBwhitelist 737K-august-2016.txt
```
(10x v3: UMI 12 bp, whitelist `3M-february-2018.txt`, clip 5' of read1 accordingly.)
Then filter the BAM to reads carrying **GX, CB and UB** tags (pysam), index it.
BAMs >30 GB: split per-cell-barcode chunks with `sinto filterbarcodes`.

Other required inputs: GTF annotation, genome FASTA (+.fai), chromosome-sizes file,
cell-barcode whitelist of cells to keep. `--species human|mouse` **only** (CNN
classifier trained on those two genomes).

### Dataset-chemistry restriction (blocker for pbmc_10k_v3)
- **A CellRanger BAM alone is NOT sufficient**: CellRanger aligns only read2; scTail's
  PAS-calling step needs read1 alignments.
- **Standard 10x runs are incompatible with the PAS-calling step**: pbmc_10k_v3 R1
  FASTQ is 28 bp (CB16+UMI12 only) — nothing to align after clipping. scTail requires
  datasets where R1 was sequenced to 150 bp past the poly(T) (as in the scTail paper's
  own data). So on pbmc_10k_v3 scTail **cannot call PAS de novo**.
- What *is* salvageable for a head-to-head: `scTail-count` only assigns **reads2** to a
  given PAS-cluster BED, so scTail's *quantification* step could be run on
  pbmc_10k_v3 given a PAS BED from elsewhere (e.g. its published PAS atlases or
  another caller) — but that benchmarks counting, not PAS detection. For a fair
  detection benchmark, pick a public "long-R1" 10x dataset (e.g. ones used in the
  scTail paper, GSE deposits listed there).

## Output formats
1. `scTail-callPeak` → per-sample working dir with paraclu intermediates, CNN input/
   prediction files, and **`positive_result.bed`** = PAS clusters passing the CNN filter
   (genomic PAS cluster coordinates, stranded).
2. `scTail-peakMerge` → **`merged_cluster.bed`** (clusters across samples merged,
   `--maxDistance` default 40 bp).
3. `scTail-count` → **`all_cluster.h5ad`** (cell x all-PAS counts, AnnData) and
   **`two_cluster.h5ad`** (genes with two PAS, for APA/BRIE2 analysis).

## Exact invocation for a 10x dataset (given a compliant read1-aligned BAM)
```bash
V=/mnt/ssd1/Projects/PeakATail_wd/tools/bench_envs/sctail/bin
$V/scTail-callPeak -b star_read1_sorted_filtered.bam \
    --gtf gencode.v44.annotation.gtf --cellbarcode barcodes.tsv \
    -f GRCh38.primary_assembly.genome.fa --species human \
    --chromoSize GRCh38.chrom.sizes -o sctail_out --minCount 50 -p 8
echo -e "sctail_out/positive_result.bed" > sctail_out/sample_list.tsv
$V/scTail-peakMerge --sampleList sctail_out/sample_list.tsv -o sctail_out
$V/scTail-count --cellbarcode barcodes.tsv -b star_read2.bam \
    -o sctail_out --PAScluster sctail_out/merged_cluster.bed -p 8
```
Notes: `-d/--device` defaults to GPU card 0, silently falls back to CPU (our env is
CPU-only torch). `LC_ALL=C` for any external sort steps.

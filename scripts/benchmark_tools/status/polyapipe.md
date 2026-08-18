# polyApipe — install status

**State: installed and smoke-tested end-to-end (Python part). R polyApiper NOT installed (optional, see below).**
Date: 2026-08-18

## What it is
polyApipe (MonashBioinformaticsPlatform/polyApipe) — de novo polyA-site/APA caller for 10x scRNA-seq.
Finds peaks of polyA-evidence reads (soft-clipped A-tails) in a CB/UMI-tagged BAM, then counts UMIs
per peak per cell. Rated best de novo sensitivity/accuracy in the 2026 NAR benchmark.

## Citation
No peer-reviewed paper. Cite the poster:
Harrison P, Williams S, Powell D, Albrecht D, Beilharz T.
"Tools for identifying and characterizing alternative polyadenylation in scRNA-Seq."
F1000Research (poster, Oz Single Cell 2019), 2019. doi:10.7490/f1000research.1117076.1
(+ GitHub: https://github.com/MonashBioinformaticsPlatform/polyApipe)

## Install (commands that worked)
```bash
python3 -m venv /mnt/ssd1/Projects/PeakATail_wd/tools/bench_envs/polyapipe
cd /mnt/ssd1/Projects/PeakATail_wd/tools/bench_envs
git clone --depth 1 https://github.com/MonashBioinformaticsPlatform/polyApipe.git polyapipe/src
polyapipe/bin/pip install --upgrade pip
polyapipe/bin/pip install pysam setuptools        # pysam 0.24.0; setuptools needed for distutils shim on py3.13
# `pip install .` of the repo FAILS (setuptools flat-layout auto-discovery error: top-level dirs data/ + polyApiper/).
# Workaround — install the script directly with venv-python shebang:
V=/mnt/ssd1/Projects/PeakATail_wd/tools/bench_envs/polyapipe
sed "1s|.*|#!$V/bin/python3|" $V/src/polyApipe.py > $V/bin/polyApipe.py && chmod +x $V/bin/polyApipe.py
# Runtime helper binaries (umi_tools sdist does not build on py3.13 — its setup.py bootstraps
# setuptools-10.0 which imports removed distutils; used conda instead):
conda create -y -n bench_polyapipe -c conda-forge -c bioconda subread umi_tools
```

## Versions
- polyApipe.py 0.1.0 (`--version`); venv python 3.13.9; pysam 0.24.0
- conda env `bench_polyapipe`: featureCounts v2.1.1 (subread), UMI-tools 1.1.6
- system samtools 1.18 (/usr/local/bin/samtools) — required, found on PATH

## Smoke-test evidence
- `/mnt/ssd1/.../polyapipe/bin/polyApipe.py --help` → full usage text, exit 0
- `--version` → `polyApipe.py 0.1.0`
- Full end-to-end run on bundled demo BAM (repo `data/demo/SRR5259354_demo.bam`) completed
  ("Processed 1200 polyA reads, 1182 included (had CB and UR tags) ... Done!", exit 0),
  producing `demo_run_polyA_peaks.gff`, `demo_run_counts.tab.gz`, `demo_run_polyA.bam`.

## Input format
- One or more coordinate-sorted, indexed BAMs (.bai required) with cell-barcode and UMI tags.
  CellRanger BAMs work as-is. Defaults: `--cell_barcode_tag CB`, `--umi_tag UR` (RAW UMI!).
  For CellRanger output pass `--umi_tag UB` to use corrected UMIs.
- No FASTQs needed; polyA evidence comes from soft-clipped A-tails on the aligned cDNA read.

## Output format
1. `<out>_polyA_peaks.gff` — PAS/peak coordinates. GFF3-ish, one row per peak: a `region_size`
   (default 250 bp) window upstream of the inferred PAS. Peak name `<chrom>_<pos>_<f|r>`;
   for `f` (+ strand) the PAS is the region END, for `r` (- strand) the PAS is the region START.
   Attributes: `peak=`, `peakdepth=` (read support), `misprime=` (True/False internal-priming flag).
   Example: `1  polyAends  polyAends  3199728  3199977  .  -  .  peak="1_3199728_r"; peakdepth="5"; misprime="False";`
2. `<out>_counts.tab.gz` — long-format UMI counts matrix: columns `gene  cell  count`
   where gene = peak ID (umi_tools count output). One file per input BAM if multiple.
3. `<out>_polyA.bam` — the extracted polyA-evidence reads (intermediate, kept).

## Dataset/chemistry restrictions
- Designed for 10x Genomics 3' scRNA-seq; works with v2 and v3 chemistry (only needs CB+UMI tags
  in the BAM and enough read length for A-tail soft-clips). No read1/FASTQ requirement.
  Not for 5' kits / non-polyA protocols. Read-level defaults: minMAPQ 10, minpolyA 5,
  misprime_A_count 8 within misprime_in 10 bp.
- Key tunables for benchmark: `--depth_threshold` (default 1 read/peak — very permissive),
  `--region_size 250`.

## Exact invocation for the 10x pbmc_10k_v3 CellRanger BAM
```bash
export LC_ALL=C
export PATH=/home/biolab/miniconda3/envs/bench_polyapipe/bin:$PATH   # featureCounts + umi_tools
BAM=/path/to/pbmc_10k_v3_possorted_genome_bam.bam                    # must have .bai next to it
mkdir -p /mnt/ssd1/Projects/PeakATail_wd/results/benchmarks/polyapipe
/mnt/ssd1/Projects/PeakATail_wd/tools/bench_envs/polyapipe/bin/polyApipe.py \
    -i "$BAM" \
    -o /mnt/ssd1/Projects/PeakATail_wd/results/benchmarks/polyapipe/pbmc10k \
    --cell_barcode_tag CB --umi_tag UB \
    -t 8
```

## R polyApiper (NOT installed)
Downstream APA analysis of the counts; README (2024) says it is "less well developed, and not
currently under active development" and recommends using the Python outputs directly — which is
all the PAS benchmark needs. If wanted later:
`R -e 'BiocManager::install(c("MonashBioinformaticsPlatform/weitrix")); devtools::install_local("/mnt/ssd1/Projects/PeakATail_wd/tools/bench_envs/polyapipe/src/polyApiper")'`
into ~/R/bench-lib.

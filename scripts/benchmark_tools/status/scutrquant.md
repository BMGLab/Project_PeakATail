# scUTRquant — benchmark tool status

Date: 2026-08-18
State: **partial** (orchestrator installed + smoke-tested; per-rule conda envs and the
358 MB UTRome index are deferred to first `--use-conda` run by the pipeline's own design —
full build not attempted per task instructions).

## Citation
Fansler M.M., Mitschka S., Mayr C. "Quantifying 3'UTR length from scRNA-seq data reveals
changes independent of gene expression." *Nature Communications* 15, 4050 (2024).
doi:10.1038/s41467-024-48254-9. Repo: Mayrlab/scUTRquant, pipeline version v0.5.1
(SQ_VERSION in Snakefile).

## What the tool is
Snakemake pipeline around a patched `kallisto bus` + `bustools`. It does **not call
peaks/PAS de novo** — it quantifies a fixed, prebuilt "UTRome" transcriptome
(GENCODE v39 protein-coding 3' ends augmented with high-confidence cleavage sites called
from the Human Cell Landscape; transcripts truncated to last 500 nt; isoforms whose
cleavage sites are <200 nt apart are merged). PAS positions for benchmarking = 3' ends of
the UTRome isoforms, obtainable from the UTRome GTF or the `rowRanges` of the output SCE.

## Outputs
- `data/sce/utrome_hg38_v1/<dataset>.txs.Rds` — Bioconductor `SingleCellExperiment`:
  sparse counts matrix cells x 3'UTR isoforms; `rowRanges` = GenomicRanges of isoforms
  (genomic coordinates, so isoform 3' end = PAS); `rowData` with isoform annotations.
- `data/sce/utrome_hg38_v1/<dataset>.genes.Rds` — gene-level counts SCE.
- Optional `output_format: ["h5ad"]` writes AnnData instead (no rowRanges).
- QC HTML report per sample under `qc/umi_count/`.
- Intermediates kept: BUS files, `txs.mtx` / `genes.mtx` + barcodes/genes txt under
  `data/kallisto/<target>/<sample_id>/` (counts made with `bustools count --em --genecounts`).
- All output paths are relative to the Snakefile dir `/mnt/ssd1/Projects/PeakATail_wd/tools/scUTRquant/`.

## Chemistry / input restrictions
- 3'-end tag-based scRNA-seq only. `tech:` config is passed to `kallisto bus -x`;
  10xv2 (16 bp CB + 10 bp UMI, whitelist 737K-august-2016.txt) and 10xv3
  (16 bp CB + 12 bp UMI, whitelist 3M-february-2018.txt) are the documented cases;
  10x Multiome GEX supported experimentally (gex_737K-arc-v1.txt).
- FASTQ mode needs full R1 (CB+UMI) plus R2 (cDNA), given as `R1;R2[;R1;R2...]` per sample.
- BAM mode (`file_type: bam` in sample sheet) feeds a CellRanger possorted BAM directly to
  the patched `kallisto bus --bam` (custom build `kallisto=0.46.2sq` from the author's
  personal conda channel `merv`; installed automatically by `--use-conda` via
  envs/kallisto-bustools.yaml, with bustools=0.40.0).
- Strand setting for 10x 3': `strand: "--fr-stranded"`.
- Human target: `utrome_hg38_v1`; mouse: `utrome_mm10_v2`. Custom targets possible via
  txcutr-db and `extdata/targets/targets.yaml`.

## Install commands that worked
```bash
cd /mnt/ssd1/Projects/PeakATail_wd/tools
git clone --depth 1 https://github.com/Mayrlab/scUTRquant.git   # -> tools/scUTRquant, v0.5.1

# env bench_scutrquant (conda solve of snakemake via bioconda timed out twice at ~7.5 min,
# so python+pandas via conda, snakemake via pip — snakemake only needs `conda` on PATH
# at runtime for --use-conda):
/home/biolab/miniconda3/bin/conda create -y -n bench_scutrquant -c conda-forge \
    --override-channels "python=3.11" "pandas>=2,<3"
/home/biolab/miniconda3/envs/bench_scutrquant/bin/pip install "snakemake==7.32.4" "pulp<2.8"
/home/biolab/miniconda3/envs/bench_scutrquant/bin/pip install "setuptools<81"
# (setuptools 84 broke it: Snakefile's min_version() imports pkg_resources, removed in 84;
#  80.10.2 works with a deprecation warning)
```
Versions installed: snakemake 7.32.4 (pip), python 3.11.15, pandas 2.3.3, setuptools 80.10.2.
Failed attempts (honest record): `conda create -n bench_scutrquant -c conda-forge -c bioconda
snakemake-minimal=7.32.4 pandas python=3.11` — killed at 460 s; retry of
`conda install -n bench_scutrquant -c conda-forge -c bioconda snakemake-minimal=7.32.4`
also killed at 440 s (bioconda solve too slow on this box).

## Smoke-test evidence
- `snakemake --version` -> `7.32.4`; `import pandas` -> 2.3.3.
- Dry run against our real pbmc_10k_v3 CellRanger BAM config builds the complete 14-job DAG
  with zero errors:
  `[INFO] scUTRquant v0.5.1 ... Loaded 1 samples ... Job stats: all, bustools_correct,
  bustools_correct_sort, bustools_count_genes, bustools_count_txs, bustools_sort,
  download_10X_whitelists, download_utrome_hg38_v1, generate_gene_merge, generate_tx_merge,
  kallisto_bus, mtxs_to_sce_genes, mtxs_to_sce_txs, report_umis_per_cell -- total 14.
  This was a dry-run (flag -n).`
- (Run with `export PYTHONNOUSERSITE=1` to keep ~/.local site-packages out of the env.)

## pbmc_10k_v3 run requirements (documented, not yet executed)
Input BAM already local: `/mnt/ssd1/Projects/PeakATail_wd/data/benchmark/pbmc_10k_v3/pbmc_10k_v3_possorted_genome_bam.bam` (22 GB, CB/UB tags).
Downloads the pipeline will fetch on first run (sizes verified via HTTP HEAD):
- hg38 UTRome bundle (gtf + kallisto .kdx index + merge tsv):
  https://github.com/Mayrlab/hcl-utrome/releases/download/v1.0.0/utrome.e30.t5.gc39.pas3.f0.9999.w500.tar.gz — **358 MB** compressed.
- tx/gene annotation Rds+csv (Mayrlab/atlas-hs v0.1.0) — ~3.5 MB total.
- 10xv3 barcode whitelist 3M-february-2018.txt.gz — 18.4 MB (all 4 whitelists fetched by the rule).
- Plus 5 rule conda envs (kallisto-bustools from channel `merv`, two r-base=4.3 envs,
  anndata, downloading) — expect several GB and >8 min total; prebuild them in chunks.

Config written (ready to use):
- `/mnt/ssd1/Projects/PeakATail_wd/scripts/benchmark_tools/scutrquant/config.yaml`
  (dataset pbmc_10k_v3, target utrome_hg38_v1, tech 10xv3, --fr-stranded,
   bx_whitelist extdata/bxs/3M-february-2018.txt, min_umis 500, cell_annots null)
- `/mnt/ssd1/Projects/PeakATail_wd/scripts/benchmark_tools/scutrquant/sample_sheet.csv`
  (`pbmc_10k_v3,bam,<absolute BAM path>`)

## Exact invocation for the 10x BAM
```bash
cd /mnt/ssd1/Projects/PeakATail_wd/tools/scUTRquant
export PYTHONNOUSERSITE=1 LC_ALL=C
export PATH=/home/biolab/miniconda3/bin:$PATH   # snakemake needs `conda` for --use-conda

# step 1 (next command; chunkable, rerun until done): build rule envs only
timeout 440 /home/biolab/miniconda3/envs/bench_scutrquant/bin/snakemake \
  --use-conda --conda-frontend conda --conda-create-envs-only -j 4 \
  --configfile /mnt/ssd1/Projects/PeakATail_wd/scripts/benchmark_tools/scutrquant/config.yaml

# step 2: real run (downloads UTRome+whitelists, then kallisto bus --bam on the 22 GB BAM;
# expect ~1-2 h; run in background, NOT under the 8-min limit)
/home/biolab/miniconda3/envs/bench_scutrquant/bin/snakemake \
  --use-conda --conda-frontend conda -j 16 \
  --configfile /mnt/ssd1/Projects/PeakATail_wd/scripts/benchmark_tools/scutrquant/config.yaml
```
Expected outputs: `tools/scUTRquant/data/sce/utrome_hg38_v1/pbmc_10k_v3.txs.Rds` (+ .genes.Rds),
`tools/scUTRquant/qc/umi_count/utrome_hg38_v1/pbmc_10k_v3.umi_count.html`.

## Benchmark note for PeakATail comparison
scUTRquant's "PAS calls" are the fixed UTRome annotation, identical for every dataset; the
dataset-specific signal is which isoforms get nonzero counts. For site-level comparison,
take 3' ends of expressed isoforms from `rowRanges(sce)` (or the UTRome GTF) with
expression thresholds from the txs counts matrix.

---

# ADDENDUM 2026-08-19 — pbmc_10k_v3 run + salvage

State: **success-after-salvage** (quantification fully completed; R/SCE packaging failed on
first run from a host-env module gap, fixed with one pip install and rerun — see below).

## What completed on the first full run (run.log, /usr/bin/time -v)
- STAGE 1 conda-env build: 00:42-02:50 (+03), two attempts (~2 h 08 m total incl. retry).
- STAGE 2 pipeline (`-j 16`): wall **25:52.06**, 312% CPU, peak RSS **4.19 GB**, exit 1.
  - kallisto bus --bam on the 22 GB CellRanger BAM: started 02:51:24; whole
    kallisto+bustools chain (sort, correct, count txs+genes) done by ~03:16 (~25 min).
  - 10 of 14 jobs finished. Outputs present in
    `tools/scUTRquant/data/kallisto/utrome_hg38_v1/pbmc_10k_v3/`:
    txs.mtx (1,082,852 barcodes x 49,410 merged isoforms, 33,845,069 nonzeros, 510 MB),
    genes.mtx, barcodes/genes txt, output.sorted.bus (4.4 GB), run_info.json.
- FAILED at rule `mtxs_to_sce_genes` (and would also have failed `mtxs_to_sce_txs`):
  `ModuleNotFoundError: No module named 'nbformat'` at Snakefile:340.

## Root cause (precise)
The failure was in the **host snakemake env** (`bench_scutrquant`), not the rule's R conda
env: snakemake 7.32.4's `script.py` does `import nbformat` (line 1458) whenever it executes
any `script:` rule. `nbformat` happened to exist only in `~/.local` user site-packages; we
run with `PYTHONNOUSERSITE=1` (correctly, to isolate the env), which hides it — so the env
itself never had it. The pip snakemake install had also skipped it.

## One-line fix (worked)
```bash
PYTHONNOUSERSITE=1 /home/biolab/miniconda3/envs/bench_scutrquant/bin/pip install nbformat
# -> nbformat 5.11.1 into the env's own site-packages
```
Rerun of only the 4 remaining jobs (`--rerun-triggers mtime` so the deleted temp .bus
intermediates did not trigger a kallisto redo), `-j 8`:
log `results/benchmark_tools/pbmc_10k_v3/scutrquant/salvage_rerun.log`.
RERUN RESULT: **SUCCESS** — 4/4 jobs (report_umis_per_cell, mtxs_to_sce_genes,
mtxs_to_sce_txs, all), wall 4:53.13, peak RSS 2.39 GB, exit 0. All pipeline outputs now
exist: `tools/scUTRquant/data/sce/utrome_hg38_v1/pbmc_10k_v3.txs.Rds` (87 MB),
`pbmc_10k_v3.genes.Rds` (79 MB), `qc/umi_count/utrome_hg38_v1/pbmc_10k_v3.umi_count.html`.
So final state is effectively **complete** (partial only in that it needed one manual
env fix mid-run). Cross-validation of the independent extraction below vs the SCE:
SCE is 49,410 isoforms x 12,253 cells (cell count identical to our >=500-UMI filter);
all 40,532 pas.bed isoforms match SCE rownames; PAS positions agree 40,532/40,532 with
3' ends of `unlist(range(rowRanges(sce)))`; SCE nonzero-count isoforms = 40,532 = bed rows.
(SCE keeps UCSC chrom names, e.g. chr7; pas.bed is translated to Ensembl.)

## Benchmark extraction (independent of SCE)
`results/benchmark_tools/pbmc_10k_v3/scutrquant/pas.bed` — BED6 1-bp PAS points,
Ensembl chrom names (chr stripped, chrM->MT), derivation in `#` header. Built by
`scripts/benchmark_tools/scutrquant/extract_pas_bed.py` directly from txs.mtx +
txs.genes.txt + the UTRome GTF (verified: all 49,410 merged-isoform IDs are
representative transcript_ids present in the GTF; PAS = 3' end, +:GTF end / -:GTF start;
spot-check NOC2L minus-strand 3' end chr1:944,203 correct).
- Cells (barcodes with >=500 total UMIs, the run's own min_umis): **12,253**
- Isoforms detected in cells: **40,532** (of 49,410 catalog; 40,696 over all barcodes)
- Total UMIs in cells: 91,173,024
- Companions: `pas_counts.tsv` (per-isoform cells/UMIs incl. all-barcode columns),
  `pas_summary.tsv`.

## MUST-STATE caveat for the comparison figure
scUTRquant is **annotation-based**: its "calls" are the fixed hg38 UTRome catalog
(GENCODE v39 3' ends + HCL cleavage sites, sites <200 nt apart merged) filtered by
detection in this dataset. It cannot discover novel PAS, and precision vs an
annotation-derived atlas is near-tautological — footnote this on any figure using pas.bed.

# Data accession manifest

Raw sequencing data is **not** included in this repository (~460 GB on the
analysis machine). This file records every input actually consumed by the
analyses here, so each result can be traced back to a public accession.

For the wider dataset survey — including candidates that were evaluated and
rejected, with access flags — see [`manuscript/04_datasets.md`](manuscript/04_datasets.md).
That document is the *shortlist*; this one is the *record of what was used*.

---

## 1. Benchmark datasets (head-to-head, seven tools)

### `pbmc_10k_v3` — human PBMC, 10x 3' v3

| field | value |
|---|---|
| Source | 10x Genomics public datasets — `pbmc_10k_v3` |
| URL | https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_v3 |
| Files used | `pbmc_10k_v3_possorted_genome_bam.bam` (44,198,457,238 B) + `.bai`; `pbmc_10k_v3_filtered_feature_bc_matrix.tar.gz` (94,334,700 B) |
| Provenance | CellRanger 3.0.0, reference `refdata-cellranger-GRCh38-3.0.0` (Ensembl 93), STAR 2.5.1b |
| Chromosome naming | bare Ensembl (`1..22, X, Y, MT` + scaffolds), 194 `@SQ` lines |
| Tags | `CB`/`UB` (~98% of reads carry `CB`), plus `CR/UR`, `GX/GN`, `RE`, `xf` |
| Read length | R2 = 91 bp |
| Cells | 11,769 filtered barcodes |
| Local path | `data/benchmark/pbmc_10k_v3/` (git-ignored) |

### `GSE104556` — Lukassen adult mouse testis, 10x 3'

| field | value |
|---|---|
| Accession | GEO **GSE104556** |
| SRA runs used | **SRR6129050**, **SRR6129051** (mouse1, mouse2) |
| Retrieval | `scripts/benchmark_tools/fetch_gse104556.sh` (`prefetch` + `fasterq-dump`) |
| Realignment | STARsolo, whitelist `737K-august-2016.txt` (10x v2) |
| Local path | `data/benchmark/gse104556/` (git-ignored) |
| Role | scAPA's flagship dataset; spermatogenesis 3'UTR shortening as positive control |

### PacBio Kinnex PBMC — long-read truth set

| field | value |
|---|---|
| Source | PacBio public cloud, Kinnex single-cell RNA datasets |
| URL | https://downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/ |
| Sets used | `primary_10x3p`, `secondary_gemx3p` (10x 3' v3.1 and GEM-X v4) |
| Local path | `data/benchmark/kinnex_pbmc_10x3p/` (git-ignored) |
| Caveat | Long-read only — no matched Illumina library from the same GEM well, and a different donor from the 10x public PBMCs. Used as an **atlas-style, site-level** truth, not cell-matched truth. |

## 2. Application dataset

### `GSE123904` — Laughney et al. 2020 human LUAD

| field | value |
|---|---|
| Accession | GEO **GSE123904** |
| Tissue | Human lung adenocarcinoma: primary, normal lung, brain/bone/adrenal metastases |
| Chemistry | 10x 3' (v2 era, HiSeq 2500) |
| Scale | ~5k cells/sample across 40+ samples, 17 patients |
| Local path | `data/laughney/` and `/mnt/ssd2/Laugney_Aligned/` (git-ignored); sample sheet tracked at `data/laughney/sample_sheet_laughney.csv` |
| Sweep outputs | `RERUN_2026-08_fixed` experiment tree |

## 3. Reference / annotation resources

| resource | version / file | notes |
|---|---|---|
| PolyASite 2.0 | `polyasite2.GRCh38.96.rep_sites.bed6` | representative sites; atlas-level PAS reference for precision/recall |
| Ensembl GTF (human) | `Homo_sapiens.GRCh38.99.gtf` | bare-Ensembl names, matches the pbmc_10k_v3 BAM directly |
| Gene-end BED | `data/references/gene_end.bed` | tracked in this repo |
| 10x barcode whitelist | `737K-august-2016.txt` | 10x 3' v2 |

### Derived shared references (regenerable)

`score_tool.py` regenerates and asserts these; they are the shared denominators
that make cross-tool scores comparable:

- `detected_genes.txt` — 14,851 unique gene IDs from column 5 of `runs/grid/lg_annotate/annotatedpas.bed`
- `genebodies.merged.bed` — `gene_end.bed` sorted + `bedtools merge` (null-shuffle inclusion regions)
- `pas2.in_detected_genes.bed` — PolyASite 2.0 representative sites restricted to strand-matched gene bodies of those 14,851 genes (`bedtools intersect -s -u`); **285,220 sites**, the recall denominator

## 4. What is deliberately absent

- Raw FASTQ/BAM for every dataset above — all are retrievable from the public accessions listed.
- Intermediate pipeline outputs (`results/`, 126 GB) and the sweep trees (`archive/`, 46 GB).
- Controlled-access datasets (EGA) — none were used in the analyses in this repository.

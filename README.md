# PeakATail — working directory (BioLab)

Analysis working directory for the **PeakATail** single-cell APA / 3'UTR isoform tool
(CLI name: `ema`; code lives in `tools/PeakATail`, upstream: https://github.com/BMGLab/PeakATail,
docs: https://bmglab.github.io/PeakATail/). Goal: a Q1 methods+application manuscript —
validation and benchmarking of PeakATail, then cell-type-specific APA biology in
lung cancer (Laughney GSE123904), tumor-associated macrophages, and PBMC scRNA-seq.

Reorganized on 2026-08-12 — see `REORGANIZATION.md` for the full old→new path map.

## Layout

```
data/         inputs (git-ignored): Laughney h5ads, STAR index, references
  laughney/     Laughney-specific inputs (barcodes, samplesheets, rds, test BAM)
  references/   reference files (gene_end.bed — formerly ./test.bed)
scripts/      analysis code (tracked)
  laughney/     Laughney notebooks & Rmds (scanpy QC, LuCA subsetting, downstream)
  pipeline/     Nextflow pipeline: main.nf, main_ek.nf, config, envs, pat_down.R
  misc/         utilities (createBW.sh, gtf2bed.py, random.chr.sh, steps.tmp.sh protocol log)
results/      analysis outputs (git-ignored)
  emaout/, ema_merge/, star/, utr_sampled/   April-2025 pipeline outputs
  yktest/       YK's complete April-2025 ema runs (laughney/ 5.9G, macs/ 1.2G + bigwigs)
  figures/      all figures (scanpy PDFs, UMAP panels, PAS-gene featureplots)
  report_data/  copy of the July-2025 partial run staged for a report
  laughney_rerun_2026-08_fixed -> /mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed
manuscript/   manuscript planning & drafts (tracked) — start at manuscript/README.md
tools/        PeakATail checkout (branch: biolab-manuscript = develop @ 3ce8cdc)
archive/      moved-aside items (git-ignored): old nextflow scratch (work/ 32G),
              amirtest/ (14G), empty placeholder files, rotated logs
```

## Key external locations

- **Fresh cohort rerun (2026-08):** `/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed/`
  (17-sample Laughney cohort, strategy grid + PolyASite benchmarks, reannotate sweep;
  per-celltype switch tasks failed on a `null/pasbed.bed` bug — needs fixed main.nf + `-resume`).
- **Validation runs:** `/mnt/ssd2/Laugney_Aligned/peakatail_experiments/VALIDATION_de26d16/` (11/11 PASS audit).
- Laughney h5ad store: `/mnt/ssd0/20250728_laughney_h5ad_files_ek/`.

## Running the Nextflow pipeline

```bash
nextflow run scripts/pipeline/main.nf -w <scratch-dir> -resume   # from the repo root
```

Known-good ownership quirks: `scripts/figures/`, `report_data/` (contents copied into
`results/`), and root-owned `data/null`, `data/work` cannot be moved/removed without sudo —
see REORGANIZATION.md.

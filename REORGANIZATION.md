# Repository reorganization — 2026-08-12

Everything was **moved, never deleted** (git `mv` for tracked paths). Branch: `reorg-manuscript`.
Old notebooks/scripts that reference old relative paths (e.g. `emaout/`, `../figures`) may need
path updates — check this table for the new location.

## Old → new path map

| Old path | New path | Notes |
|---|---|---|
| `202505_macs_analysis.R` | `scripts/misc/202505_macs_analysis.R` | macrophage Seurat script |
| `20250721_laughney_macs_cell_barcodes.csv` | `data/laughney/…` | |
| `20250724_laughney_pasgenes_featureplots.pdf` | `results/figures/…` | |
| `20250724_luca_subset_laughney.rds` | `data/laughney/…` | 189M Seurat object |
| `laughney_3utr_selected.bam(.bai)` | `data/laughney/…` | |
| `output.png` | `results/figures/laughney_umap_metadata_panels.png` | renamed (was generic) |
| `sample_sheet.csv` | `scripts/pipeline/sample_sheet.csv` | ≠ `data/sample_sheet.csv` (different file, untouched) |
| `sample_sheet_laughney.csv` | `data/laughney/…` | rows point at /mnt/lun2,lun3 BAMs |
| `test.bed` | `data/references/gene_end.bed` | renamed — md5-identical to the pipeline's gene_end.bed |
| `main.nf`, `main_ek.nf`, `nextflow.config`, `envs/`, `bamfiles.yaml` | `scripts/pipeline/…` | `params.samplesheet`/`params.rscript` updated accordingly |
| `scripts/pat_down.R` | `scripts/pipeline/pat_down.R` | referenced by main.nf `params.rscript` |
| `scripts/20250506_laughney_scanpy_ek.ipynb` | `scripts/laughney/…` | |
| `scripts/20250723_laughney_downstream.rmd` | `scripts/laughney/…` | reads `amirtest/laugh/emaout` → now `archive/amirtest/laugh/emaout` |
| `scripts/202507_23_laughney_luca_subset_ek.ipynb` | `scripts/laughney/…` | |
| `scripts/20250724_laughney_subset_luca.rmd` | `scripts/laughney/…` | |
| `scripts/202511_laughney_manipulations.ipynb` | `scripts/laughney/…` | |
| `scripts/20250728_laughney_scdownstrem_samplesheet(_paths).csv` | `data/laughney/…` | inputs, not code |
| `scripts/createBW.sh`, `gtf2bed.py`, `random.chr.sh`, `environment.yml`, `steps.tmp.sh` | `scripts/misc/…` | steps.tmp.sh is the command-protocol log |
| `emaout-old/` | `results/emaout-old/` | empty stubs of an aborted run |
| `yktest/` | `results/yktest/` | complete April ema runs: laughney/ (5.9G), macs/ (1.2G) |
| `scripts/figures/` | **copied** to `results/figures/` | source dir owner-locked (ebrukocakaya), original left in place |
| `report_data/` | **copied** to `results/report_data/` | source dir owner-locked (amiramiritabat), original left in place |
| `try.bed`, `regions_3utr_py.bed`, `requirements.txt`, `environment.yml` (root), `scripts/laughney.ipynb`, `scripts/Untitled-1.ipynb` | `archive/empty-placeholders/` | all 0-byte |
| `amirtest/` | `archive/amirtest/` | 14G superseded sandbox (2 tracked files moved via git) |
| `work/` | `archive/work/` | 32G Nextflow scratch — delete to reclaim space |
| `.nextflow/` | `archive/nextflow-cache/` | |
| `.nextflow.log*` (10 files) | `archive/nextflow-logs/` | |
| `emaout/negbed.bed`, `emaout/negmatrix.mtx` | — | were already deleted from disk; deletion committed |

New: `results/laughney_rerun_2026-08_fixed` → symlink to
`/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed`.

## `data/` grouped into subfolders

`data/` was a flat dump of ~50 items; it is now:

| Old path | New path | Notes |
|---|---|---|
| `data/humanSTARindex` (symlink), `data/chrom.sizes.nochr.filt`, `data/hg38.chrom.sizes` | `data/references/` | joins `gene_end.bed`; `main.nf` `params.genomeDir/gtf/chromSizes` updated |
| `data/*.fastq.gz` (20 symlinks → /mnt/lun2 macrophage + A549 co-culture reads) | `data/macrophage/` | `scripts/pipeline/sample_sheet.csv` read1/read2 paths updated |
| `data/2025*_laughney*.h5ad`, `20250507_laughney_metadata.csv`, `PATIENT_LUNG_ADENOCARCINOMA_ANNOTATED.h5` | `data/laughney/` | 17G of Laughney inputs |
| `data/test.bam`, `KI270373.1.bam`, `KI270744.1.bam(.bai)`, `downsampled_aligned.bam(.bai)` | `data/testdata/` | small dev/test slices |
| `data/.nextflow`, `data/.nextflow.log*` | `archive/nextflow-logs/data-run/` | logs moved; `.nextflow` dir is root-owned and stayed |

Absolute paths inside old notebooks that pointed at `…/data/<file>.h5ad` now need
`…/data/laughney/<file>.h5ad`. Also fixed: `main.nf`/`main_ek.nf` hardcoded
`scripts/gtf2bed.py` → `scripts/misc/gtf2bed.py`.

**Pre-existing breakage (not caused by this reorg):** both `.nf` files source
`tools/PeakATail/.emaenv/bin/activate`, a virtualenv that no longer exists — it was part of
the discarded "yk from biolab" commit. This legacy pipeline needs a real env (or should be
retired in favor of the ssd2 Nextflow pipeline) before it can run again.

## Items that need sudo (couldn't be moved as `biolab`)

Directory renames need write permission on the directory itself; these are owned by others:

- `scripts/figures/` (ebrukocakaya) — contents already copied to `results/figures/`; safe to `sudo rm -rf` after verifying.
- `report_data/` (amiramiritabat) — contents copied to `results/report_data/`; same.
- `data/null/` (root, 16K), `data/work/` (root, 4K) — stray root-owned dirs from an old run; `sudo rm -rf` when convenient.
- `data/.nextflow/`, `data/Laugney_scdownstream/` (6.2M), `data/ge_apa_merged_mudata/` (425M) — root-owned, so they stayed at the `data/` top level. Move with sudo: `Laugney_scdownstream` → `data/laughney/`, `ge_apa_merged_mudata` → `data/` is fine, `.nextflow` → `archive/`.

One-time permanent fix for the shared-repo permission problem (new git objects currently
require a workaround because ~90 `.git` paths belong to ebrukocakaya):

```bash
sudo chown -R biolab:biolab /mnt/ssd1/Projects/PeakATail_wd/.git
sudo chmod -R a+rwX /mnt/ssd1/Projects/PeakATail_wd/.git
```

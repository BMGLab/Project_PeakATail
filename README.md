# Project_PeakATail

Analysis code, figure source data and environment specifications for the
**PeakATail** manuscript — validation and benchmarking of single-cell poly(A)
site detection, and cell-type-resolved alternative polyadenylation (APA)
analyses in lung adenocarcinoma and mouse testis.

The **tool itself** lives in a separate repository:
[BMGLab/PeakATail](https://github.com/BMGLab/PeakATail) (CLI name `ema`,
docs at <https://bmglab.github.io/PeakATail/>). This repository is the
companion analysis record for the paper.

## Layout

```
manuscript/          plans, reports and the figure index
  figures/             final figures as PNG + vector PDF
  github/              issue/PR text staged for the tool repo
source_data/         the tables each manuscript figure is plotted from
scripts/
  manuscript_figures/  one script per figure
  benchmark_tools/     seven-tool head-to-head harness + per-tool status files
  pipeline/            Nextflow pipeline (Laughney application analyses)
  laughney/            notebooks and Rmds for the LUAD cohort
  misc/                utilities
envs/                conda environment exports + TOOL_VERSIONS.md
data/                only small tracked references (sample sheets, gene_end.bed)
DATA_ACCESSIONS.md   every input traced to a public accession
```

Raw data is **not** in this repository (~460 GB on the analysis machine).
Every input is traceable through [`DATA_ACCESSIONS.md`](DATA_ACCESSIONS.md).

## Reproducing a figure

Each figure has one script in `scripts/manuscript_figures/` and one or more
tables in `source_data/`. The tables are the plotted values, so figures can be
regenerated without re-running the upstream pipeline:

```bash
python scripts/manuscript_figures/benchmark_headtohead.py
```

Note that the figure scripts contain absolute paths to the analysis machine
(`/mnt/ssd1`, `/mnt/ssd2`, `/mnt/ssd0`). Re-running them end-to-end requires
either the raw data in those locations or path edits; `source_data/` exists so
that the plotted numbers are inspectable without either.

## Benchmark

Seven single-cell APA callers were run on identical input BAMs per dataset:
PeakATail, polyApipe, scAPAtrap, SCAPTURE, scTail, scUTRquant and Sierra.
Versions, pins and environment files are in
[`envs/TOOL_VERSIONS.md`](envs/TOOL_VERSIONS.md); the input contract and
scoring protocol are in `scripts/benchmark_tools/README.md`; per-tool install
and smoke-test transcripts are in `scripts/benchmark_tools/status/`.

Results are written up in `manuscript/09_headtohead_results.md` and
`manuscript/07_curated_benchmark_report.md`. **These include negative results
for PeakATail** — they are reported as found.

## History

This repository's history was rewritten once, before publication, to remove a
451 MB intermediate matrix (`amirtest/macro/emaout/negmatrix.mtx`) that
exceeded GitHub's file-size limit. All 41 commits are otherwise preserved;
commit SHAs differ from the internal working repository.

## Before release

Items to settle before making this repository public and minting a DOI:

- [ ] Replace the placeholder author list in [`CITATION.cff`](CITATION.cff) with the final manuscript author list and ORCIDs.
- [ ] Add a `LICENSE` file (CITATION.cff currently declares MIT).
- [ ] Push branch `biolab-manuscript` of BMGLab/PeakATail so the pinned benchmark commit `18678ef` resolves publicly — it is currently local-only, which blocks exact reproduction.
- [ ] Confirm whether the three post-benchmark read-module fixes (up to `2e7fc0d`) change any reported number; the current `source_data/` predates them.
- [ ] Flip the repository to public, then enable the Zenodo webhook and cut a release to mint the DOI.

## License

MIT (see `CITATION.cff`; a `LICENSE` file is still to be added).

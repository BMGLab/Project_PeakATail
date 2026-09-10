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

## History and syncing

This repository's history was rewritten once, before publication, to remove a
451 MB intermediate matrix (`amirtest/macro/emaout/negmatrix.mtx`) that
exceeded GitHub's file-size limit. All 41 commits are otherwise preserved;
commit SHAs differ from the internal working repository.

Because of that rewrite this repository cannot be updated with a plain `git
pull` from the internal working directory. Use
[`maintenance/sync_from_working_dir.sh`](maintenance/sync_from_working_dir.sh),
which re-runs the identical filter (deterministic, so unchanged commits keep
their SHAs) and replays the packaging commits on top.

## Before this repository is made public

It currently contains the **unsubmitted manuscript draft** and internal working documents that
came across in the sync — `manuscript/00_draft/` (MAIN.md, SUPPLEMENTARY.md), and notably
`PI_DECISIONS.md` and `REVIEW_SIMULATION.md`, the latter being a *simulated* referee report
written in-project, not a real review. Decide deliberately which of these should be visible
before flipping the repository to public; nothing here is public yet.

The Availability statement in the manuscript promises this repository carries the benchmark
scorer, the replication filter, the long-read truth-set build scripts, the label-confirmation
and universe policy files, the figure scripts with their audit tables and captions, the
pre-registration documents, and the run manifests. As of this commit all of those are present:

| promised | where |
|---|---|
| benchmark scorer | `scripts/benchmark_tools/score_tool.py` |
| replication filter | `scripts/` (`replication_filter.py`) |
| Kinnex truth + decoy point BEDs | `truth_sets/kinnex/` (gzipped) |
| Kinnex build script | `scripts/benchmark_tools/` |
| pre-registration documents | `manuscript/24_`, `26_`, `27_` |
| figure scripts and captions | `scripts/manuscript_figures/`, `manuscript/figures/*.caption.md` |
| figure source data | `source_data/` |
| numbers audit, competitor run table | `audit/` |
| run manifests | `run_manifests/` |

## Before release

Items to settle before making this repository public and minting a DOI:

- [ ] Replace the placeholder author list in [`CITATION.cff`](CITATION.cff) with the final manuscript author list and ORCIDs.
- [ ] Add a `LICENSE` file (CITATION.cff currently declares MIT).
- [ ] Push branch `biolab-manuscript` of BMGLab/PeakATail so the pinned benchmark commit `18678ef` resolves publicly — it is currently local-only, which blocks exact reproduction.
- [ ] Confirm whether the three post-benchmark read-module fixes (up to `2e7fc0d`) change any reported number; the current `source_data/` predates them.
- [ ] Flip the repository to public, then enable the Zenodo webhook and cut a release to mint the DOI.

## License

MIT (see `CITATION.cff`; a `LICENSE` file is still to be added).

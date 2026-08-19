# Benchmarked tool versions

Every tool in the head-to-head benchmark was installed and smoke-tested on the
same machine on 2026-08-18 and consumed the **same input BAM** per dataset (see
`scripts/benchmark_tools/README.md` for the input contract). Per-tool install
transcripts, smoke-test evidence and exact invocations live in
`scripts/benchmark_tools/status/<tool>.md`.

| Tool | Version | Source / pin | Environment file |
|---|---|---|---|
| **PeakATail** (CLI `ema`) | 0.2.0 | github.com/BMGLab/PeakATail, commit `18678ef` (branch `biolab-manuscript`) | installed from source |
| polyApipe | 0.1.0 | polyApipe.py; featureCounts 2.1.1 (subread), UMI-tools 1.1.6 | `bench_polyapipe.yml` |
| scAPAtrap | 0.2.0 | BMILAB/scAPAtrap, commit `12ee8f5df34fbe81642fe8564132e3f847f93132` (2023-09-10) | `bench_scapatrap.yml` |
| SCAPTURE | v1.0 (2021-01-25, only release) | git main @ 2026-08-18; TensorFlow 2.0.0, featureCounts 1.6.4 | `bench_scapture.yml` |
| scTail | 0.1.8 | pip; PyTorch | (venv, see `status/sctail.md`) |
| scUTRquant | v0.5.1 (`SQ_VERSION`) | Mayrlab/scUTRquant, `--depth 1` clone @ 2026-08-18 | `bench_scutrquant.yml` |
| Sierra | 0.99.27 | GitHub master via `remotes`, R lib `~/R/bench-lib` | `bench_sierra.yml` |
| STAR / STARsolo | see `bench_star.yml` | used to produce placeholder CB/UB tags for umi_tools-based arms | `bench_star.yml` |

## Reproducing an environment

```bash
conda env create -f envs/bench_sierra.yml   # etc.
```

The `.yml` files are full `conda env export --no-builds` dumps, so they pin
package versions but not builds. They were exported from the live environments
on 2026-08-19, after all benchmark runs completed.

## Pipeline

The Nextflow pipeline used for the Laughney application analyses is in
`scripts/pipeline/` (`main.nf`, `nextflow.config`, `envs/star_env.yaml`).

## Caveat on the PeakATail pin

The benchmark ran against commit `18678ef` of the `biolab-manuscript` branch.
That branch also carries three later read-module fixes (up to `2e7fc0d`) that
were made *after* the benchmark runs and are therefore **not** reflected in the
numbers in `source_data/`. See README "Before release".

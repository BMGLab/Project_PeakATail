# PeakATail (reference tool) — status

Date: 2026-08-18. Tool: `/mnt/ssd1/Projects/PeakATail_wd/tools/PeakATail`,
branch `biolab-manuscript`, commit `18678ef`, `ema --version` -> **0.2.0**.
(NB `/home/biolab/Projects/PeakATail_wd` is the same directory — same
device+inode — so tracebacks may show either path.)

## State found, and repairs performed (all into the tool's own `.venv`)

The bundled `.venv` (Python **3.10.12**) was stale relative to the
biolab-manuscript branch:

1. Missing declared deps (all in `pyproject.toml [project] dependencies`).
   Installed 2026-08-18, commands that worked:
   ```
   ./.venv/bin/pip install 'sortedcontainers>=2.4'
   ./.venv/bin/pip install 'click>=8.1' 'rich>=13.7' 'questionary>=2.0' 'tqdm>=4.66'
   ./.venv/bin/pip install 'PyYAML>=6.0' 'psutil>=5.9' 'statsmodels>=0.14' 'pyarrow>=15.0'
   ```
2. **Stale console script**: `.venv/bin/ema` still points at the OLD egg-info
   entry (`ema = ema.main:main`, see `ema.egg-info/entry_points.txt`); any
   invocation — even `ema --help` — falls into the legacy pipeline body and
   crashes with `could not open alignment file 'None'`. The branch declares
   `ema = ema.cli:main` (click group) in pyproject.
3. Refresh attempt **failed honestly**:
   `./.venv/bin/pip install --no-deps -e .` →
   `ERROR: Package 'peakatail' requires a different Python: 3.10.12 not in '>=3.11'`.

## Working invocation (the workaround used by run_peakatail.sh)

The click CLI itself imports and runs fine on 3.10:

```
./.venv/bin/python -c 'from ema.cli import main; main()' <subcommand> [args...]
```

Smoke-test evidence (2026-08-18):

- `... main()' --help` prints the click group with subcommands
  `collapse merge parse-gtf reannotate run switch wizard`.
- `... main()' --version` → `version 0.2.0`.
- `... main()' run --help` prints the full option set (`--bam-dir --gtf
  --output --threads --barcode-tag --cb-len --seq-len --ignore-chro
  --atlas --atlas-mode --peak-strategy ...`).

`PEAKATAIL_NO_TIMESTAMP=1` makes `--output` deterministic
(`ema/cli/common.py:resolve_output_dir`).

## Proper fix (next command, NOT run — likely >8 min, chunk it)

```
cd /mnt/ssd1/Projects/PeakATail_wd/tools/PeakATail
/usr/bin/python3.11 -m venv .venv311
./.venv311/bin/pip install 'pysam>=0.22' 'numpy>=1.26' 'pandas>=2.1' 'scipy>=1.11'   # batch 1
# batch 2: anndata scanpy igraph leidenalg louvain pybedtools statsmodels pyarrow ...
./.venv311/bin/pip install -e . --no-deps
```
Then point run_peakatail.sh's `VENV_PY` at `.venv311/bin/python`.

## Input / output formats

- Input: one or more CB/UB-tagged BAMs (`--bam-dir` single-BAM convenience or
  `--bam-files`/YAML datasets) + a GTF whose chrom naming matches the BAM.
- PAS output: `<out>/pasbed.bed` — BED6 peak **intervals**
  (`chrom start end id score strand`), chrom naming follows the input BAM
  (bare Ensembl for pbmc_10k_v3). Also `annotatedpas.bed` (+gene columns),
  `annotated_matrix.mtx` (cell x PAS counts), stage dirs `01_peak_calling` …
  `08_differential`. Verified on the Laughney RERUN_2026-08_fixed runs.

## Exact invocation for the 10x BAM

See `../run_peakatail.sh` (the reference wrapper; preflight refuses the
still-downloading BAM). Core command:

```
/usr/bin/time -v -o runtime_mem.txt \
  .venv/bin/python -c 'from ema.cli import main; main()' run \
    --bam-dir  data/benchmark/pbmc_10k_v3/pbmc_10k_v3_possorted_genome_bam.bam \
    --gtf      /home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.99.gtf \
    --output   results/benchmark_tools/pbmc_10k_v3/peakatail/run \
    --threads 16 --barcode-tag CB --cb-len 16 --seq-len 91 \
    --ignore-chro MT --plot-engine none --no-progress
```
(no `--atlas`: snapping onto the atlas we later score against would be circular)

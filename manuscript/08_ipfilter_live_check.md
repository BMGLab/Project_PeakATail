# Internal-priming filter: live wiring check (issues #69/#71)

Date: 2026-08-19. Tool: `tools/PeakATail`, branch `biolab-manuscript`, commit
`18678ef`, `ema` 0.2.0, invoked via the py3.10 venv workaround
(`.venv/bin/python -c 'from ema.cli import main; main()'`; see
`scripts/benchmark_tools/status/peakatail.md`).

## Verdict (citable)

> PeakATail's internal-priming filter is live and off by default: a run with no
> `--ip-filter` flag executes no internal-priming code and its PAS BEDs are
> byte-identical to an `--ip-filter --ip-filter-mode annotate` run's;
> `annotate` mode keeps every PAS and records a `True`/`False`
> `internal_priming` flag as the 10th column of the per-dataset
> `annotatedpas.bed` (here 138/3,664 PAS flagged, 3.8%); `filter` mode drops
> exactly the flagged PAS (3,664 -> 3,526 strand-level; unified `pasbed.bed`
> 1,862 -> 1,773, a strict subset, -89 rows, 0 added).

Per issue #69's question: the default is genuinely ip-off, confirmed both in
code (`ema/main.py`: `ip_filter = bool(getattr(args, "ip_filter", False))`;
`_apply_pas_filters` returns `None` without importing filter code when off)
and behaviorally (the no-flag run's log has zero `peak_filters` lines, writes
no `peak_filters_stats.json`, and leaves the `internal_priming` column empty).

## Setup

- Input BAM: `data/testdata/downsampled_aligned.bam` (1.25 GB, bare-Ensembl
  chroms, CB tag 12 bp, UB 8 bp, reads 101 bp), **RG-stripped** to
  `results/ipfilter_check/input_noRG.bam` first — see "Environment repairs".
- GTF: `/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.99.gtf`
- Genome FASTA (for `--ip-filter`):
  `/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.dna.primary_assembly.fa`
  (`.fai` present)
- Launcher: `results/ipfilter_check/run.sh` (nohup, DONE.ok/FAILED.err
  markers; three runs in parallel, 8 threads each).

Common flags for all three runs:

```
PEAKATAIL_NO_TIMESTAMP=1 .venv/bin/python -c 'from ema.cli import main; main()' run \
  --bam-dir results/ipfilter_check/input_noRG.bam \
  --gtf /home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.99.gtf \
  --output results/ipfilter_check/<name> \
  --threads 8 --barcode-tag CB --cb-len 12 --seq-len 101 \
  --ignore-chro MT --plot-engine none --no-progress
```

| run | extra flags | wall time | max RSS |
|---|---|---|---|
| `off` | (none) | 10:08 | 0.65 GB |
| `annotate` | `--ip-filter --genome-fasta <fa> --ip-filter-mode annotate` | 9:36 | 0.63 GB |
| `filter` | `--ip-filter --genome-fasta <fa> --ip-filter-mode filter` | 9:23 | 0.64 GB |

All three completed the full 8-stage pipeline (01_peak_calling ...
08_differential).

## Observations

### 1. Default = off (issue #69)

- `off` log: **zero** `peak_filters` mentions; no
  `04_pas_gene_assignment/peak_filters_stats.json`.
- `off/03_gtf_annotation/default/annotatedpas.bed` has the 10-column D9
  layout but the `internal_priming` column is `""` for all 3,664 rows
  (column present, never populated — the documented "stage never ran"
  encoding).
- `diff off/pasbed.bed annotate/pasbed.bed` -> identical (1,862 rows each);
  root `annotatedpas.bed` (legacy 9-col gene-annotation format) also
  byte-identical; per-dataset `01_peak_calling/default/pasbed.bed`
  identical across **all three** runs (3,664 rows — peak calling is
  deterministic and unaffected by the filter, which runs downstream of it).

### 2. `annotate` mode annotates, drops nothing (issue #71 wiring)

- `annotate/04_pas_gene_assignment/peak_filters_stats.json`:
  pos total=1,731 passed=1,731 flagged=60 (3.5%); neg total=1,933
  passed=1,933 flagged=78 (4.0%); removed=0 both strands.
- `annotate/03_gtf_annotation/default/annotatedpas.bed` column 10
  (`internal_priming`): 3,526 `False` / 138 `True` — 60+78=138 matches the
  stats JSON exactly.

### 3. `filter` mode drops exactly the flagged PAS

- stats JSON: pos 1,731 -> 1,671 (60 removed), neg 1,933 -> 1,855
  (78 removed) — removed == flagged on both strands.
- Root `pasbed.bed`: 1,862 (off) -> 1,773 (filter). `comm` on sorted rows:
  89 rows only-in-off, **0** rows only-in-filter — a strict subset (the 138
  strand-level PAS collapse to 89 rows after PAS unification).
- Root legacy `annotatedpas.bed`: 2,029 -> 1,938 gene-assignment rows.
- Note: the per-dataset `annotatedpas.bed` still lists all 3,664 PAS with
  their flags even in filter mode (it is built from the pre-filter
  per-dataset pasbed) — the drop is realized in the unified/root artifacts
  and everything downstream of gene assignment.

### 4. Ledger caveat

No `provenance/pas_ledger.tsv` was written by any of the three runs: on this
branch the single-BAM path does not emit the provenance ledgers (they are
written on the multi-sample/atlas path). The ip flag's queryable surfaces in
a single-BAM run are the per-dataset `annotatedpas.bed` 10th column and
`peak_filters_stats.json`.

## Environment repairs performed (2026-08-19, would otherwise invalidate the check)

1. **`pyfaidx` was missing from the tool venv.** Crucially,
   `ema/experimental/internal_priming.py` treats a missing pyfaidx as a
   *silent no-op* (logs an error, copies input to output, flags nothing) — an
   ip run without it would falsely look like "filter does nothing". Installed
   with `./.venv/bin/pip install 'pyfaidx>=0.7'` (got 0.9.0.4) before any run.
2. **Mixed-length RG ids crash the pipeline.** The raw test BAM carries three
   read groups (`Mono_day0`, `A549trd_M2mac_day8`, `untrd_mac_day7`);
   `ema/countmatrix/read.py` builds barcodes as `f"{rg}_{cb}"` and
   `ema/countmatrix/cb_encode.py::encode_cb_batch` requires equal-length
   strings -> `ValueError: cannot reshape array` in every run. Worked around
   by stripping `@RG`/`RG:Z:` into `results/ipfilter_check/input_noRG.bam`
   (all reads then share the uniform `default_` prefix). Candidate upstream
   bug report; unrelated to the ip-filter wiring itself.

## Artifacts

Everything under `/mnt/ssd1/Projects/PeakATail_wd/results/ipfilter_check/`:
`run.sh`, `{off,annotate,filter}/` run dirs, `{off,annotate,filter}.{log,time}`,
`input_noRG.bam(.bai)`, `DONE.ok`.

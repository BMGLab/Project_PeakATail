# source_data — the numbers behind each figure

One or more tables per manuscript figure, copied verbatim from the figure
generation run. These are the plotted values: they let a reader check every
number in the paper without access to the raw data or the pipeline.

Naming follows the figure it belongs to — e.g. `benchmark_headtohead.tsv`
backs `manuscript/figures/benchmark_headtohead.pdf`, produced by
`scripts/manuscript_figures/benchmark_headtohead.py`.

`FIGURE_SOURCE_README.md` is the generation-time README carried over from the
figure output directory; it documents how each table was produced.

Notable tables:

- `benchmark_consolidated.tsv` — all tools, all datasets, all scoring arms in one table
- `benchmark_headtohead.tsv` — the head-to-head summary figure
- `fdr_calibration_null_pvalues.tsv` — full permutation-null p-value distribution behind the FDR-calibration result
- `spermatogenesis_control_*.tsv` — panels of the mouse testis positive control

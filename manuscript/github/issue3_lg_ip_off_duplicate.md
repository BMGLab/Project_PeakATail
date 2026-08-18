In `RERUN_2026-08_fixed/runs/grid/`, the arm meant to test "internal-priming filter off" produced output **byte-identical** to the default arm:

```
md5sum lg_annotate/pasbed.bed lg_ip_off/pasbed.bed
20cb455024723dfa4256c1329f262386  (both)
```

Their `benchmark_vs_polyasite_v3.json` files are numerically identical too. So the grid's IP contrast is really only two-sided: `lg_ip_filter` (20,609 PAS, −8.9%) vs default (22,633).

Two possible causes worth distinguishing:
1. The sweep config generated identical CLI invocations for the two arms (config bug), or
2. The tool's default already runs with IP filtering off, making an "off" arm inherently a duplicate (in which case the grid design just needs relabeling).

Related docs problem: `docs/cli/run.md` still carries **contradictory rows** for `--ip-filter`/`--genome-fasta` — the D9 annotate-mode table *and* stale "currently no-op" rows. HANDOFF says D9 wiring landed (`ema/experimental/internal_priming.py`). Before the manuscript Methods claims integrated IP handling, someone should verify the live behavior with a small BAM and delete whichever doc rows are wrong.

@TRextabat — can you check which of (1)/(2) it was in `gen_configs.py` / `peakcall_grid.tsv`?

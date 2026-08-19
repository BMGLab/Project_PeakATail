Running `ema run` on a **stock CellRanger BAM** (10x public `pbmc_10k_v3`, `--barcode-tag CB --cb-len 16`) produces a **completely empty run**: 0 peaks, 0 barcodes, 0-byte pos/neg/pas beds — after 93 minutes of streaming the 44 GB BAM — then crashes much later in `find_close` with an opaque `pandas.errors.EmptyDataError: No columns to parse from file`.

## Cause

`ema/countmatrix/read.py`:

```python
cb = read.get_tag(barcode)
if len(cb) != barcode_len:
    return 0, 0, 0, 0, 0   # silent drop
```

CellRanger tags corrected barcodes with a GEM-group suffix: `CB:Z:TTGCATTGTTTCGTTT-1` → **18 characters**. With `cb_len=16`, the length check silently discards **100% of reads**. STARsolo BAMs (no suffix) are unaffected — which is why all internal runs worked. The README's claim of CellRanger compatibility is currently false in practice.

## Fix applied on our side (branch `biolab-manuscript`, commit `ad59cdd`)

Strip a trailing `-<digits>` when the raw tag is over-length by exactly the suffix:

```python
if len(cb) != barcode_len:
    dash = cb.rfind("-")
    if dash == barcode_len and cb[dash + 1:].isdigit():
        cb = cb[:dash]
```

Happy to turn this into a PR onto `develop` with a regression test (synthetic BAM with suffixed tags).

## Two follow-ups this exposes

1. **Fail fast**: if the BAM stream ends with zero usable CB reads, abort with a clear message ("0 reads carried a parseable <tag> barcode — check --barcode-tag/--cb-len") instead of writing empty artifacts and failing later in an unrelated stage.
2. **Multi-GEM-group inputs**: aggregated CellRanger BAMs carry `-1`, `-2`, … suffixes that distinguish samples; plain stripping merges them. The right general fix keeps the suffix as part of the sample identity (like the existing `rg_cb` composite) — worth deciding before the PR.

@TRextabat

**Update — second CellRanger-input bug in the same code path** (commit `6df8eed` locally): unmapped reads (kept, barcoded, in CellRanger BAMs) have `reference_end = None`, crashing `read_check` with a `TypeError` two hours into the run. Both bugs share a root cause: the read loop was only ever exercised on STARsolo BAMs pre-filtered to mapped 3′UTR reads. A regression test with a stock CellRanger BAM slice would catch this class.

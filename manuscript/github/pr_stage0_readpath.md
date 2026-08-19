## What this fixes

Makes stock **CellRanger BAMs** work (they currently produce empty or corrupted runs) and repairs a regression introduced by an earlier draft of these fixes:

1. GEM-group CB suffix (`AAAC...-1`, 18 chars) failed the `cb_len` check → **100% of reads silently dropped** (93 min to an empty run).
2. Unmapped kept reads (`reference_end=None`) → TypeError hours in.
3. RG IDs containing underscores corrupted the `<sample>_<CB>` composite (split on first `_`) → all barcodes collapsed into one column. The earlier draft *discarded* such RGs — itself a data-corrupting regression for `ema merge` output (two merged samples would share columns). Replaced with a **bijective sanitiser**; test proves two underscore-bearing RGs stay distinct.

Branch is green: 0 new failures vs develop (mock drift fixed).

## Honest scope limit (pre-existing, documented, not a regression)

Multi-sample *merged* CellRanger input still isn't fully handled by `DatasetManager.prepare`; single-BAM CellRanger runs (the benchmark case) are correct. Tracked for Stage 1.

## Tests

New CellRanger-shaped unit tests + RG bijection tests; suite parity with develop. Adversarially re-verified (SOUND, with the scope limit above surfaced by the verifier).

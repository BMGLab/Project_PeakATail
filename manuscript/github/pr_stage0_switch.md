## What this fixes

Two `switch length` correctness defects that together invalidate shipped PDUI:

1. **Strand inversion, residual path**: the `(gene,"_gene_",0,1,1)` fallback sentinel used by `--isoform-agg per_isoform --utr-unmatched gene` ignores strand — **1667/1667 fallback rows inverted** in shipped output (per_gene itself was fixed in #66). Also `switch length` silently fell back to a plus-strand convention when coordinates were missing; it now RAISES, and gains `--pasbed`.
2. **PDUI computed on TF-IDF weights**: `build_count_dfs` read `.X`, which clustering overwrites in place. 100% of shipped proximal/distal "read counts" are non-integer; on informative cells r = 0.47 vs true counts, 10.6% of ΔPDUI signs flip. Now reads a raw `--counts-layer` (default `counts`, stashed before normalisation) and RAISES on non-integral matrices unless `--allow-non-count-matrix`.

Plus the permanent guard: any `pdui_classic.tsv` row must satisfy strand `'+'` ⇒ proximal < distal, `'-'` ⇒ proximal > distal — the check that produced the 2299/2299 inversion evidence.

## Impact on the blocked SWITCH_CELLTYPE rerun (#67)

**Hold the rerun until this merges** — otherwise it returns inverted PDUI computed on TF-IDF weights. Unblock path: backfill `layers['counts']` from `06_preprocessing/*/preprocessed.h5ad` (verified byte-compatible), then only B2+combine+switch re-execute.

## Tests

The 3 xfail(strict) tests now pass for real; net +15 vs develop, failures a strict subset of develop's 7. Docs updated (breaking change called out). Adversarially re-verified (SOUND).

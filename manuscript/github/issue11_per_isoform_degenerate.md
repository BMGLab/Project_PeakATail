# Issue: `switch length --isoform-agg per_isoform` emits degenerate proximal==distal pairs; the strand guard misclassifies them

**Labels:** bug, switch-length. Found and reproduced on the final-benchmark testis runs (code `4efeb12`, both mice).

## Symptom
`ema switch length --strategy classic --isoform-agg per_isoform --utr-unmatched drop` exits rc=1 on both
GSE104556 mice: mouse1 18,116 / 13,709,930 rows flagged (7 genes), mouse2 15,004 / 15,148,xxx (8 genes).

## Diagnosis (verifier-reproduced)
All flagged rows have `proximal_pas_id == distal_pas_id` with equal starts — 14 (gene, transcript) units
broadcast to all 1,294 cells; every read-bearing flagged row has `proximal_reads == distal_reads` and
`pdui` exactly 0.5. There are **zero** true strand inversions (0 different-strand rows). The per-isoform
aggregation selects the same PAS as both endpoints when a transcript has one usable PAS.

The guard at `ema/switch_test/runner.py:~397` uses `~(proximal_start < distal_start)`, so equality lands
in the "inverted / strand" bucket and the error message blames strand selection — the wrong cause.

## Proposed fix
- In per-isoform aggregation, drop (or mark `single_pas`) transcripts whose proximal and distal endpoints
  resolve to the same PAS instead of emitting a degenerate pair.
- Split the guard into `inverted` (proximal_start > distal_start, a real bug signal) and `degenerate`
  (equal ids/starts), with messages naming each.
- Add a regression test: transcript with a single PAS → no degenerate row, run completes.

`--isoform-agg per_gene` passes the same guard cleanly on both mice (strand convention verified on all
~10.4k genes per mouse, 0 violations). Repro commands and counts:
`results/stage3_spermatogenesis_final/REPORT_PROVISIONAL.md` §5 on biolab.

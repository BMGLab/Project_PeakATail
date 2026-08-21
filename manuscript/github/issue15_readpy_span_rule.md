# Issue: `read_check`'s span rule silently discards 13.7 % of valid-CB reads (96 % of them spliced) before both the clip detector and the count matrix — and fabricates the 3' end of ~10 % of the rest

**Labels:** bug, accuracy, quantification. Measured by `results/algo_headroom/A2_clip_sensitivity/`
and `A1_end_pileup/` on the frozen v2 tree `9dfdefb`; the genome-wide census and the
truth:decoy quality of the recovered reads were checked by the independent verifier
(`results/algo_headroom/VERIFY/`, verdict FIXED).

Related: **#98** (per-isoform degenerate pairs), **#99** (PAS→gene assignment) — both are downstream
of the matrix this rule feeds. The code already carries `# TODO: this block will be fixed`.

## 1. Symptom

`ema/countmatrix/read.py:101-109`:

```python
span = read_end - read_start
if span > seq_len:
    return 0, 0, 0, 0, 0
elif span < seq_len:
    read_end = read_start + seq_len
```

`read_end` is `pysam`'s `reference_end`, so `span` is the alignment's **reference** span, not the
read's length. Any read whose *reference* footprint exceeds `--seq-len` — i.e. **any spliced
alignment** — is dropped outright, before clip detection, before peak calling and before the count
matrix. Any read shorter than `--seq-len` has its 3' end **overwritten** with `read_start + seq_len`.

## 2. Measured diagnosis

* The rule discards **13.74 % of all valid-CB reads genome-wide, 96.0 % of them spliced**
  (`A2 tables/T2b`, `T3d`). On the PBMC chr19+21 slice the figure is 24.19 %, and on chr19:1–20 Mb it
  reaches 23.84 % — the discard rate tracks gene density, so it is worst exactly where 3' UTRs are.
* **88.8 M reads never enter the count matrix** for this reason on PBMC 10k v3.
* For clip detection specifically the rule costs **56,335 qualifying clip reads (+1.76 % genome-wide;
  +4.64 % on the slice)**. Small — but it is the **only** relaxation we tested that is positive on
  **both** truths at **every** operating point (+0.7 % to +1.4 % atlas, +1.9 % to +2.1 % Kinnex
  recall at matched precision; projected ~+0.4 % genome-wide after the 2.6× slice haircut).
* **The recovered reads are as good as the ones we already use.** The increment's truth:decoy ratio
  against the Kinnex long-read truth is **3.0 : 1**, against 3.67 : 1 for the existing evidence
  (`A2 tables/T7b`). Every *other* relaxation we tested came in at 0.83–0.89 : 1, i.e. more than half
  internal priming — which is why lowering `--polya-min-clip` is rejected and this is not.
* At *verified* PAS, the span rule alone removes **19,362 reads = 23.6 % of all CB reads there,
  99.6 % of them spliced** (`A2 tables/T12b`).
* The `elif span < seq_len` branch **fabricates the end coordinate of 8.2–11.7 % of the reads it
  keeps** (`A1` stage 0). Their reported 3' end is `start + 91` — a constant offset from the
  alignment start, i.e. a shifted copy of coverage. This is one reason the end-position channel
  carries no information beyond local depth after depth matching (0.95–1.01×; `VERIFY` §2.1).
* **0 hard-clipped alignments in 663,383,342 records** — nothing is hidden by hard clipping, so the
  soft-clip path is the whole story.

## 3. Proposed change

- [ ] Replace the reference-span test with a **query-length** test: reject on
      `read.query_length > seq_len` (or `infer_read_length()`), not on `reference_end - reference_start`.
      A spliced alignment of a 91 nt read is a 91 nt read.
- [ ] Do not extend short alignments to a fabricated `start + seq_len` on the **clip** path; the clip
      site is defined by the soft-clip boundary, which is known exactly.
- [ ] If the extension is load-bearing for the coverage path, make it explicit and flagged
      (`--extend-short-reads/--no-extend-short-reads`) and record the fabricated-end count in the run
      log, so downstream readers know which coordinates are inferred.
- [ ] Behind a flag on `peakAtail-prime` (`--span-rule {legacy,query-length}`, default `legacy` until
      the 24 §3.1 criterion is met), because it changes **both** the candidate set and the count
      matrix.
- [ ] Regression tests: a spliced alignment whose reference span exceeds `seq_len` but whose query
      length does not must be accepted; the accepted-read count on the committed fixture must be
      asserted.
- [ ] While in this file: `read_check` applies **no `-F 3844`**, so secondary alignments of
      multimappers are counted once per alignment on the coverage path (25 % of records in some
      chr19 windows; `A1` incidental findings). Decide and document whether that is intended —
      `pas_support.tsv` already carries both `clip_reads` and `clip_reads_f3844`, so the two
      conventions coexist in one output today.

## 4. How to reproduce

    export LC_ALL=C
    cd /mnt/ssd1/Projects/PeakATail_wd
    sed -n '95,115p' tools/pa-polya-run-9dfdefb3/ema/countmatrix/read.py
    cat results/algo_headroom/A2_clip_sensitivity/tables/T2b* T7b* T12b*
    sed -n '1,40p' results/algo_headroom/A2_clip_sensitivity/VERIFIER_ADDENDUM.md

## 5. What it is expected to buy

**For detection: little — projected ~+0.4 % recall at matched precision**, and A2's slice figure
already had a 2.6× haircut applied to get there. It is on the list only because it is the single
relaxation that is positive on both truths at every operating point and whose added evidence is the
same quality as the evidence we already trust.

**For quantification: potentially a lot, and it is unmeasured.** 88.8 M reads and 13.7 % of all
valid-CB reads never reach the matrix, non-randomly — spliced reads, so gene-body- and
gene-density-dependent. Every count, every UMI, every switch test and both open downstream issues
(#98, #99) sit on top of that matrix. **Nobody has measured what changes when those reads are
included**; that measurement should be part of this issue, not an afterthought.

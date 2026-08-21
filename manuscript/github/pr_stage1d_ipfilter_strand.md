# PR: fix(ip-filter) — test the downstream side of the cleavage site on the `-` strand

**Branch:** `fix/ip-filter-strand` (3 commits: `3142270` fix + tests, `e79005b` docs + changelog)
**Base:** `feat/polya-evidence` at `4efeb12` — **merge after that branch** (it is the code the final Stage-2 benchmark ran); then `develop`.
**Fixes:** Stage-1d issue §1 (`manuscript/github/issue10_stage1d_ipfilter_memory.md`; found by the final-gate verifier, `manuscript/15_final_gate.md` §6).
**Scope:** `ema/experimental/internal_priming.py` only (+ tests + docs). No CLI/YAML/default/stats-key change. `+`-strand output byte-identical.

## What was wrong

`filter_internal_priming()` fetched one forward-coordinate window `[pos-10, pos+30)` for **both** strands
(`pos` = BED `end` on `+`, BED `start` on `-`) and scanned it for `AAAAAA` / A-fraction ≥ 0.7 on `+`,
`TTTTTT` / T-fraction ≥ 0.7 on `-`. Oligo-dT internal priming is caused by a genome-encoded A-stretch
**downstream of the cleavage site in transcript orientation**. On `-`, transcript-downstream is
genomically *upstream* (lower coordinates), so the shipped window tested 30 nt upstream / 10 nt downstream
of the cleavage site — mostly the wrong side:

```text
shipped (both strands)                         genomic:  [pos-10 ............ pos+30)
  + strand  5'---UUUUUUUUU|cleavage>AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA   correct: 10 up / 30 down
  - strand  3'   ...down...<cleavage|UUUUUUUUUUUUUUUUUUUUUUUUUUUUUU   wrong:   30 up / 10 down (transcript)

fixed (`-` mirrored)                            genomic:  [pos-30 ............ pos+10)
  - strand  TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT<cleavage|UUUUUUUUU---5'   10 up / 30 down, reverse-complemented
```

## The fix

* `ip_window(pos, strand, left, right)` → genomic `[start, end)`: `[pos-left, pos+right)` on `+`,
  `[pos-right, pos+left)` on `-` (start clamped at 0; the contig slice truncates the end).
* `call_internal_priming(window_seq, strand, a_stretch, a_fraction)` → `(flag, tested_seq)`: pure;
  upper-cases, reverse-complements on `-`, then applies the unchanged rule (`'A'*a_stretch in seq` or
  A-fraction ≥ `a_fraction`) to the transcript-oriented sequence. Equivalent to the old T-scan on the
  forward strand, now on the correct window.
* `check_internal_priming(contig, pos, strand, ...)` composes the two for a `str` or `pyfaidx` record;
  `filter_internal_priming()` calls it per row. Flags, `stats` keys, modes (`annotate` / `filter`),
  `--ip-window-left/right`, `--ip-a-stretch`, `--ip-a-fraction` and their defaults are unchanged.
* `+` strand is byte-identical to 4efeb12 (randomised check, 20,000 cases incl. `a_stretch ∈ {3,6}`,
  `a_fraction ∈ {0.5,0.7}`, `left ∈ {0,5,10}`, `right ∈ {10,30}`, contig-edge positions; and 0 changed
  flags on 444,000 real `+` sites below).

## Tests (`tests/test_internal_priming_strand.py`, 19 new)

Synthetic **C/G-only** genome (no A/T run or A/T-rich window can occur by accident) with implanted signals:

| case | assertion |
|---|---|
| (a) `+`, 6-A run downstream of BED `end` | flagged; same run upstream → not flagged; `+` == shipped rule over a sweep of run positions |
| (b) `-`, 6-T run genomically **upstream** of BED `start` (transcript-downstream), at `[pos-25, pos-19)` | flagged — and the re-implemented 4efeb12 rule is asserted to **miss** it (regression) |
| (b) `-`, 6-T run genomically downstream at `[pos+15, pos+21)` | not flagged — and the 4efeb12 rule is asserted to **flag** it (wrong side) |
| (b) strand symmetry | a `-` site on `g` == the `+` site at the mirrored coordinate on `revcomp(g)` |
| (c) fraction rule, both strands | 40-nt window with 28 A (0.70, no 6-run) flagged; 27 A (0.675) not |
| (d) contig edges | `-` at `start=5` (window clamped `[0,15)`), `+` run ending at the contig end, position beyond the contig (empty window → `False`), contig missing from FASTA (row kept, `False`) |
| (e) end-to-end | tiny FASTA + BED6 with 4 sites (one per case) through `filter_internal_priming()` in `annotate` and `filter` mode, and through `peak_filters.apply_filters()` (the seam `ema.main` calls per strand BED) |

Full suite (`python -m pytest -q tests`, pyfaidx 0.9): baseline 4efeb12 = 2 failed / 1215 passed / 9 skipped /
4 xfailed / 1 xpassed; this branch = **2 failed / 1234 passed** / 9 skipped / 4 xfailed / 1 xpassed. The 2 failures are
the pre-existing environment failures in `tests/test_pyproject_install.py` (entry point not installed in this venv).

## Measured impact on the final Stage-2 run (no caller re-run)

Both rules applied to the emitted `pasbed.bed` of the final run (code 4efeb12) with the tool's own pyfaidx path
(`results/benchmark_tools/ipfix_impact/ip_strand_impact.py`; per-site flags in `*.flags.tsv`, counts in
`*.summary.tsv`). `pos` convention identical to the caller's; chromosome names are Ensembl in both BED and FASTA
(0 rows skipped). Sanity: the old rule re-applied to the IP-arm survivors flags **0** sites on both strands, and
applied to the no-IP arm (same caller, same seed) it reproduces the IP arm's precision default **exactly**
(44,394 sites, P@100 0.7167) — so the no-IP arm is the pre-filter candidate set and the post-hoc numbers are exact.

| dataset / set | strand | n | old flagged | new flagged | both | newly flagged | newly un-flagged |
|---|---|---:|---:|---:|---:|---:|---:|
| PBMC IP-arm survivors | + | 167,923 | 0 | 0 | 0 | 0 | 0 |
| PBMC IP-arm survivors | − | 160,464 | 0 | 5,599 (3.49 %) | 0 | 5,599 | 0 |
| PBMC pre-filter candidates (no-IP arm) | + | 203,152 | 35,229 (17.3 %) | 35,229 | 35,229 | 0 | 0 |
| PBMC pre-filter candidates (no-IP arm) | − | 199,708 | 39,244 (19.6 %) | 33,626 (16.8 %) | 28,027 | 5,599 | **11,217** |
| testis mouse1 IP-arm survivors | + | 36,774 | 0 | 0 | 0 | 0 | 0 |
| testis mouse1 IP-arm survivors | − | 35,487 | 0 | 1,505 (4.24 %) | 0 | 1,505 | 0 |
| testis mouse2 IP-arm survivors | + | 36,412 | 0 | 0 | 0 | 0 | 0 |
| testis mouse2 IP-arm survivors | − | 35,946 | 0 | 1,302 (3.62 %) | 0 | 1,302 | 0 |

Within the PBMC tier-1 ≥2-molecule subset (the pre-registered precision default): `−` candidates 30,484,
old flagged 9,603, new flagged 7,472, newly flagged 386, newly un-flagged 2,517.

**Pre-registered precision default, P@100 (PolyASite 2.0 rep sites + protein-coding TES, strand-matched, 100 bp;
`scripts/benchmark_tools/score_tool.py` with the `stage2_final_launch.sh` arguments; outputs in
`results/benchmark_tools/ipfix_impact/score_*.tsv`):**

| set | n scored | P@100 | R@100 (detected-gene atlas) | F1@100 (detected) |
|---|---:|---:|---:|---:|
| PBMC, as run (old rule) | 44,394 | 0.7167 | 0.1707 | 0.2757 |
| PBMC, newly-flagged 386 removed | 44,008 | **0.7184** | 0.1698 | 0.2747 |
| PBMC, full corrected filter on the candidates (386 removed, 2,517 restored) ≈ a re-run | 46,524 | **0.7062** | 0.1754 | 0.2811 |
| testis mouse1, as run | 25,991 | 0.7414 | — | — |
| testis mouse1, newly-flagged 297 removed | 25,695 | 0.7458 | 0.2013 | 0.3170 |
| testis mouse2, as run | 26,164 | 0.7554 | — | — |
| testis mouse2, newly-flagged 226 removed | 25,938 | 0.7585 | 0.2042 | 0.3218 |

(The verifier's estimate was 0.7167 → ~0.7182 for the first PBMC row.) The sites the old rule removed for the
wrong reason (transcript-upstream T-runs) are, as a group, slightly below the average precision of the set, so a
re-run with the corrected filter lands at ≈0.706 rather than 0.717 on PBMC while recovering ~2,100 sites and
+0.5 pp recall. The mouse "re-run equivalent" cannot be derived post hoc (no no-IP final mouse arm); the mouse
rows only remove the newly flagged sites.

**Benchmark conclusions are unchanged**: the pre-registered gate (P@100 ≥ 0.50 on PBMC and both testis mice)
passes in every row above, and the `+` strand — half of every call set — is untouched. What *does* change is
the per-site `internal_priming` flag / drop set on `−`: ~3.5–4.2 % of surviving `−` sites should have been
flagged and ~28 % of the removed `−` sites (11,217 of 39,244 on PBMC) should not have been. **Any use of the
per-site flag from a pre-fix run must re-run the filter** — `ema reannotate --genome-fasta <fa>` on the run
directory does it without re-calling peaks — rather than trust the run's flag.

## Files

```
 CHANGELOG.md                          |  29 ++++
 docs/cli/run.md                       |  24 ++-
 docs/strategies/peak-calling.md       |   3 +-
 ema/experimental/internal_priming.py  | 109 ++++++++----
 tests/test_internal_priming_strand.py | 305 ++++++++++++++++++++++++++++++++++
 5 files changed, 439 insertions(+), 31 deletions(-)
```

_Commits: `3142270` fix + tests; `e79005b` docs + changelog; `9f1a23a` docs table repair (verifier). Based on `feat/polya-evidence`; merge after the Stage 1b+1c PR._

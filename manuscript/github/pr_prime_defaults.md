# defaults: the internal-priming veto runs when you type nothing — with a byte-for-byte v2 no-op guarantee

**Base: `develop` (this branch already merges `develop`'s current tip, `5073d63`).** Branch:
`peakAtail-prime`, 31 commits + 1 merge on top of `9dfdefb` (= v2, the code behind the current
manuscript numbers).

This is a **defaults-and-correctness PR, not an accuracy PR.** Its central claim is a no-op
guarantee: **the branch's default output is byte-identical to the v2 manuscript arm** (md5-equal
`pas.bed` / `pas_tier1.bed` / `pas_tier2.bed` on all four benchmark libraries — PBMC 10k v3, both
GSE104556 testis mice, PBMC 4k), so merging it cannot move a single published number. The
pre-registered accuracy criterion (`manuscript/24` §3.1, default vs default) therefore **FAILS with
every delta exactly 0.000000 — an identity, not a regression**. Amendment 3 of that document records
why: the criterion was written against a v2 default (`--ip-filter-mode filter`) that never existed;
every benchmark arm passed it explicitly. What this PR changes is what a user gets by **typing
nothing**.

---

## 1. The no-op guarantee, proven on the exact commit being pushed

Re-run on this PR's HEAD (post-merge), against fresh reference runs of the frozen v2 worktree
(`pa-polya-run-9dfdefb3`), `scripts/prime/identity_check.py` three-tier protocol:

| run pair | data files | `pas.bed` md5 |
|---|---|---|
| v2 code vs prime + `V2_COMPAT_FLAGS`, PBMC chr19+21 slice | **50/50 byte-identical** | `227720c57cfab4baf3e021fa2313fc68` (both) |
| v2 code vs prime + `V2_COMPAT_FLAGS`, GSE104556 mouse1 chr18+19 slice | **50/50 byte-identical** | `7677c8894b59a431784ad0cad99560fd` (both) |
| **prime at bare defaults** (only `--genome-fasta`) vs the v2 *manuscript* arm (v2 + `--ip-filter --genome-fasta --ip-filter-mode filter`), PBMC slice | all four derived BEDs md5-equal | `bb65b0e0e6cef78c4dc0e3300a3ed6e3` (both); counts `pas 15,925 / tier1 8,524 / tier2 7,401 / tier1>=2mol 2,883` |

In the compat pairs the run journal is identical after timestamp/worker-order normalisation, and
`run_config.json` / `run_manifest.json` differ **only by 29 new keys** (the recorded prime options
and the resolved internal-priming block) plus the documented `cleavage_offset` respelling (`0` →
`none`, the same parsed no-op) — zero keys lost, zero values drifted.

The guarantee is enforced in the suite, not just claimed here:

* `tests/test_prime_v2_compat_golden.py` — the real caller over a committed fixture BAM, every
  output byte SHA-256-checked against goldens generated from the frozen v2 worktree; shown to bite
  when a behavioural change lands without a flag.
* `tests/test_prime_compat_flags.py` — ties the documented compat command line
  (`V2_COMPAT_FLAGS`, 14 flag/value pairs), the library-level pin and the CHANGELOG prose to each
  other in both directions, and diffs RunConfig **defaults** (not just field names) against the
  frozen v2 tree, so a silently moved default cannot hide again.
* `tests/test_prime_ip_default_mode.py` — drives the real filter seam with `args` seeded exactly as
  a no-flag `ema run` seeds it and asserts a flagged site is gone from the BED on disk.

There is **no `--compat` flag** (an earlier help text claimed one; corrected in `06cd227`). The v2
incantation is the 14-pair `V2_COMPAT_FLAGS` command line, and restoring v2's *call set* alone needs
only `--ip-filter-default off`.

## 2. What changes for a user who types nothing

One behavioural default: **`--ip-filter-default auto`** (new flag, default `auto`) plus
`--ip-filter-mode` defaulting to the sentinel `auto`, which resolves to `filter`. When a readable
`--genome-fasta` is present, the internal-priming veto that every published benchmark arm passed
explicitly now runs by default; with no FASTA the caller warns loudly, naming the cost, and emits
the v2 flagless output byte-for-byte (verified: md5 `227720c5...`, same as v2).

Default vs the old flagless behaviour, full PBMC 10k v3 BAM, ≥2-molecule arms
(`results/prime_bench/fourway_headline.txt` §2; matched-budget controls
`results/prime_bench/adversarial/veto_matched_pbmc10k.tsv`):

* **Headline, never quoted alone:** atlas P@100 0.5907 → 0.7062 (**+19.5 %**) at R_det@100
  0.1911 → 0.1754 (**−8.2 %**) — but those two arms differ in call count (62,110 vs 46,524), so the
  like-for-like statement is the matched-budget one below.
* **At matched call count (both 46,524), prime's default dominates the flagless configuration on
  all five read-outs at once:**

  | read-out | v2 flagless @ 46,524 | prime default @ 46,524 |
  |---|---:|---:|
  | atlas P@100 | 0.6640 | **0.7062** |
  | R_det@100 | 0.1675 | **0.1754** |
  | Kinnex ≥5-UMI P@25 | 0.6782 | **0.7647** |
  | Kinnex ≥20-UMI P@25 | 0.5228 | **0.5584** |
  | Kinnex IP-decoy@25 (lower better) | 0.2810 | **0.1300** |

* **The one row prime does not win, stated rather than hidden:** forced to the flagless
  configuration's own 62,110-call budget, prime must reach into its 1-molecule tier and lands
  marginally behind (P 0.5818 / R 0.1870 vs 0.5907 / 0.1911). That row is **87 % decided by the
  tie-break at the cut** (121,041 sites tied at 1 molecule, 15,586 kept) and is not a like-for-like
  comparison — recorded with its tie-exposure columns in the same TSV.

The measurement behind making the veto the default at all: it is worth +7.7 % to +12.8 % relative
recall at matched atlas precision — bigger than every detector change tested put together, because
it is the only stage that reads genomic sequence
(`results/algo_headroom/VERIFY/tables/v8_ipveto_value.tsv`).

## 3. BREAKING CHANGE — read this before merging

**A flagless v2 user who supplies `--genome-fasta` silently loses 17.09 % of their calls** — the
internally-primed ones (PBMC full set: 402,765 → 333,920, −68,845). Nothing about their command
line changes; the output does. That is the intended behaviour and it is the definition of a
breaking default change, so it must ride a version bump and a CHANGELOG entry that says exactly
this. The CHANGELOG entry is in the branch ("BREAKING — the flagless call set changes"); the
version itself is untouched (`pyproject.toml` stays 0.2.0, the entries are filed under
*Unreleased*), so the bump lands with the release tag at whatever number the maintainer chooses —
this PR's position is only that it must not ship inside a patch release. Escapes, in increasing
strictness:

* `--ip-filter-mode annotate` — the veto still runs but drops nothing (flags are recorded);
* `--ip-filter-default off` (or `--no-ip-filter`) — v2's call set;
* the full 14-pair `V2_COMPAT_FLAGS` line — v2, byte-for-byte, every file.

## 4. Also on by default — all measured call-set-neutral

| change | flag | what it does | why it is safe |
|---|---|---|---|
| feature sidecar | `--pas-features on` | appends 25 documented columns (sequence context, IP covariates, hexamers, local candidate context, `clip_offset_mean`) to `pas_support.tsv` | adds, drops and moves no PAS; `off` is v2's sidecar byte-for-byte; costs 32 µs / 0.49 kB per PAS |
| per-site IP flag | (rides the sidecar: `ip_tool_flag`, `ip_tool_afrac`, …) | the per-site internal-priming flag #95 §3 asked for, in both modes | column-append only |
| exact clip-rate QC | `--clip-rate-sampling pass` | counts every accepted read and qualifying clip during the peak-calling pass; exact rate per (contig, strand) | zero extra I/O; replaces the head-of-chr1 sample that reads **4.2× high** (2.2565 % vs true 0.5364 % on PBMC) — an over-estimate that would mask a destroyed poly(A) channel, the one condition the QC exists to catch. `head` = v2 including its log line; `strided` kept as the measured-unreliable obvious design |
| inferred cleavage | `--emit-inferred-cleavage on` | appends an `inferred_cleavage` column; the reported coordinate never moves | `pasbed.bed` byte-identical; the measured per-library offset is zero (−0.334 bp PBMC / +0.214 bp mouse1), which is why `--cleavage-offset` defaults `none` |

## 5. Ships OFF, with the measured negative results attached

Per the pre-registration's flag policy: where the measurement says an idea does not pay, the flag
ships off and the numbers ship with it.

* **`--read-geometry {fixed,keep,true}`, default `fixed` (= v2).** `true` stops discarding spliced
  reads (24.2 % of valid-CB reads on the PBMC slice, 98 % spliced) and buys **+31.2 % raw
  count-matrix mass** (PBMC slice; +23.5 % mouse1; ~+15.9 % genome-wide) — but it **fails the
  pre-registered detection criterion**: ΔR_det +0.0031 on PBMC (< +0.010 required) and ΔP@100
  **−0.0292 on mouse1** (5.9× past the −0.005 allowance). The recovered reads carry real but less
  precise evidence: admitting them slides along the curve instead of lifting it. The flag exists
  for quantification use; making it default needs a quantification criterion the pre-registration
  does not have.
* **`--pas-score {none,calibrated,select}`, default `none`.** The calibrated per-site score
  transfers across species (fitted on mouse1 alone, up-and-right on held-out PBMC and mouse2 at
  matched n: PBMC P@100 0.7062 → 0.7538), but **it partly learned the atlas's own annotation
  prior**: against Kinnex long reads at matched call count t5 P@25 falls 0.7647 → 0.7050, and at
  matched atlas precision 0.7647 → 0.6624 (`manuscript/28` §2 records it as losing ~15 % of
  long-read-verified sites); on separating long-read termini from decoys it does not beat the
  caller's own inverted `ip_tool_afrac` (AUC 0.7185 vs 0.7790), while an otherwise identical
  Kinnex-trained model does (0.8041). No fixed probability threshold clears the criterion on all
  three datasets (the ranking transfers; the calibration does not). Full transcript:
  `results/prime/TASK_D_pas_score.md`.
* Also off, each with its recorded reason: `--read-exclude-flags` (a clean slide along the curve),
  `--cleavage-offset auto` (measured no-op), `--pas-gene-rescue` (no rescue is free),
  `--dynamic-threshold-clamp` (§6 — off *is* the fix's honesty: on cannot change a run that
  does not crash, and off stays v2).

## 6. Bug fix that rides along: the `--dynamic-threshold` IndexError, guarded

The 2026-08 parameter sweep proved (`results/paramsweep/VERDICTS.md` §5) that
`--dynamic-threshold --lambda-fold-change 2.0` — the parameter reference's own "find more PAS"
setting — **aborts the caller with a bare `IndexError`** at
`l_end = data_array[-current_threshold]`: nothing bounds the dynamically-raised threshold by the
live window's length. The sweep's verifier reproduced the abort on the mouse slice too
(species-independent); it is a v2 defect, reachable only with `--dynamic-threshold` (off by
default), so no published number is affected. This PR adds **`--dynamic-threshold-clamp {off,on}`,
default `off` = v2 to the character** (same expression, same abort — the cardinal rule wins over
the fix); `on` bounds the index and logs a census of how often it fired.
`tests/test_prime_dynamic_threshold_clamp.py` reproduces the crash on a 14-read synthetic BAM
through both the monolithic and pipeline loops (verified to raise the identical IndexError at
`peackcalling.py:528` on the frozen v2 worktree), proves `on` completes those identical runs, and
proves `on` == `off` byte-for-byte on any dynamic-threshold run that does not crash. Issue draft
with the full diagnosis accompanies this PR.

## 7. Flag reference

| flag | values | default | v2 value |
|---|---|---|---|
| `--ip-filter-default` | `off`, `auto` | **`auto`** | `off` |
| `--ip-filter-mode` | `auto`, `annotate`, `filter` | **`auto` (→ `filter`)** | `annotate` |
| `--no-ip-filter` | flag | off | (n/a; `--ip-filter` absent in v2) |
| `--pas-features` | `off`, `on` | **`on`** | `off` |
| `--emit-inferred-cleavage` | `off`, `on` | **`on`** | `off` |
| `--clip-rate-sampling` | `head`, `strided`, `pass` | **`pass`** | `head` |
| `--cleavage-offset` | `none`, `auto`, signed int | `none` | `0` (same no-op) |
| `--read-geometry` | `fixed`, `keep`, `true` | `fixed` | `fixed` |
| `--read-exclude-flags` | int (SAM mask) | `0` | `0` |
| `--pas-score` | `none`, `calibrated`, `select` | `none` | (absent) |
| `--pas-score-model` / `--pas-score-min` | name/path, float | `prime1` / `-1` | (absent) |
| `--pas-gene-rescue` / `--pas-gene-rescue-min-mol` | `off`,`inside` / int | `off` / `0` | (absent) |
| `--dynamic-threshold-clamp` | `off`, `on` | `off` | `off` (the IndexError) |

Canonical machine-readable copy: `ema.cli.config_schema.V2_COMPAT_FLAGS`;
`tests/test_prime_compat_flags.py` keeps it, the CHANGELOG prose and the library pin in lockstep.

## 8. Tests

Full suite on this HEAD: **1,570 passed**, 9 skipped, 4 xfailed, 1 xpassed, plus exactly the 2
known `tests/test_pyproject_install.py` environment failures (branch-cut baseline was 1,277
passed; every addition is a new test, nothing pre-existing moved). The three-path agreement suite
(monolithic / pipeline / tiled) covers every new per-read seam, and the four compat/golden test
files in §1 and §6 are the no-op guarantee in executable form.

## 9. Suggested review order (~31 commits, grouped)

1. **The contract first:** `c7752d4` (PRIME_PLAN: the five changes, their flags, the measurement
   behind each default), `f6e6111` (v2 byte-identity golden on the fixture, before any behaviour
   change).
2. **The one behavioural default:** `f42c985` (TASK E bundle: `--ip-filter-default auto`, exact
   clip-rate, cleavage-offset column, gene seam measurement), then `f1f8174` (auto turns the veto
   ON, it does not make it DROP), then **`6954082`** — the most important commit in the branch to
   review: the shipped default was a measured no-op (veto ran, mode `annotate`, dropped 0), and
   this makes `--ip-filter-mode` default resolve to `filter`, with the run log and `run_config.json`
   recording the resolution.
3. **The guarantee machinery:** `a71f1dc` (the compat command line is checkable, not prose),
   `06cd227` (no `--compat` flag exists; docs corrected).
4. **Call-set-neutral additions:** `325140a` → `f09b020` (pas-features: schema, merge-seam fixes,
   memory, cost), `e7259ec`/`8f9dc3e` (cleavage-offset column and its measured no-op).
5. **The negative results, implemented then shelved:** `ecf488b` → `9e98dc1` (read-geometry),
   `f4b7820` → `a877acf` (pas-score; `f95e7c8` is the transfer measurement, `58bbc44` the author's
   own two-defect correction), `6edac79` (PRIME_PLAN rev 2 — the verifier's veto-stays-a-hard-gate
   ruling).
6. **The crash guard:** `1f80f9b` (§6; includes the wiring an interrupted session left unreachable
   and the discard of its unimplemented `--warn-inert-options` stub); `5e18e6f` (the §3 BREAKING
   statement, added to the CHANGELOG itself).
7. **The merge:** `9ff958a` (develop's `7b50a5f` switch-diff NameError fix and `5073d63` logging
   teardown fix; both verified off the caller output path).

## 10. Issue cross-references

* **#95** — this PR completes §3's first ask: a **per-site internal-priming flag in both modes**
  (`ip_tool_flag` + covariates in `pas_support.tsv`, on by default). The default flip itself is
  built on the measurement recorded in #95's addendum (the shipped rule flags 17.09 % genome-wide;
  applying it takes the flagless arm from P 0.5907 to 0.7062 and decoy 0.3122 to 0.1300 — arms of
  different call count; §2's matched-budget table is the like-for-like statement) — and on
  that addendum's own withdrawal: `--ip-rule kinnex` stays dead. §1 (minus-strand window) was #96;
  §2 (memory/CPU) was #97; §3's other bookkeeping items (zero-count PAS drop, unified col-5,
  seq-len log line) remain open.
* **#99** — §2 (clip-rate warning samples only the head of chr1) is completed by
  `--clip-rate-sampling pass`, and the 2× wrong docstring constant is corrected (1.152 % →
  0.5730 % measured). §1 (gene assignment in overlapping loci) is measured on this branch (the
  lncRNA seam, CHANGELOG §4 of the TASK E entry) but not changed: no rescue was free.
* **#94 / #98** — not addressed here (switch-test statistics and per_isoform degenerate pairs;
  both remain open). The merged develop tip includes `7b50a5f`, which fixes a *different* switch
  defect (the within_utr/between_utr NameError) — adjacent, not a resolution of either issue.

---

🤖 Generated with [Claude Code](https://claude.com/claude-code)

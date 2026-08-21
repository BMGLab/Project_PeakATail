# `scripts/prime/` — dev harness for the `peakAtail-prime` branch

Everything here supports one job: making the **shipped / v2 / prime** comparison in
`manuscript/24_prime_preregistration.md` reproducible. Read 24 before changing a default;
read `tools/pa-prime/PRIME_PLAN.md` for what the branch is actually changing.

**Always `export LC_ALL=C`.** The box's `tr_TR` locale silently corrupts GNU `sort` and
`bedtools` ordering. Every script here does it; anything you write must too.

## Scripts

| script | usage |
|---|---|
| `common.sh` | sourced by the others: paths, the launcher's reference-argument sets, `assert_code`, `guard_load`. Not executable on its own. |
| `make_slice.sh [out.bam]` | verifies (or builds) the PBMC chr19+21 dev slice. Reuses the perf team's slice after checking with `samtools idxstats` that exactly contigs 19 and 21 carry reads — it does not assume. |
| `run_slice.sh --label L [--code DIR] [--threads N] [--bam B] [--gtf G] [--seq-len N] [-- extra ema flags]` | runs the caller on the slice into `/mnt/ssd0/emaout/peakatail_benchmark/prime/<label>/`. |
| `score_pas.sh <pas.bed> <label> [--species human\|mouse] [--outdir D] [--slice-contigs 19,21] [--no-kinnex]` | scores with `score_tool.py` + the launcher's reference arguments; prints n / P@10 / P@25 / P@100 / R_det@100 / F1 / null, plus Kinnex t5 & t20 P@25 and the IP-decoy rate at 25 bp. |
| `identity_check.py REF_DIR NEW_DIR [-o report.tsv] [--strict-volatile]` | the compatibility-mode check: byte-identity of two caller output trees. |

## Two things these scripts refuse to do

**Run against the wrong code.** Five worktrees of this package exist (`PeakATail`, `pa-prime`,
`pa-perf`, `pa-polya`, and the frozen `pa-polya-run-*` snapshots). `assert_code` imports `ema`,
`ema.strategies.clip_seeded` and `ema.countmatrix.polya` and asserts all three resolve inside the
requested tree before a run starts. Picking up the wrong worktree has invalidated results on this
project before.

This is not hypothetical, and it fails **silently** without the assertion: `ema` is importable with
**no `PYTHONPATH` at all**, resolving to
`/home/biolab/Projects/PeakATail_wd/tools/PeakATail/ema/__init__.py` — the `develop` worktree,
through an editable install and the bind mount (`/home/biolab/Projects/PeakATail_wd` and
`/mnt/ssd1/Projects/PeakATail_wd` are the same directory). A forgotten `PYTHONPATH` therefore runs
**develop**, not prime, and produces plausible output. Demonstrated: `assert_code` on a tree with no
`ema/` fails with `AssertionError: WRONG TREE: ema -> /home/biolab/.../tools/PeakATail/ema/__init__.py`.
Every run writes the resolved path to `<label>/modpath.txt`; check it before believing a number.

**Compete with the Stage-3 chain.** `guard_load` refuses to start if load1 is above 60 or less than
20 GB is available, and `run_slice.sh` refuses `--threads > 8`. The replication chain under
`/mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v3/` owns 24 cores and is never touched.

## What `run_slice.sh` fixes, and why

Flags are copied verbatim from `scripts/benchmark_tools/stage2_final_launch.sh`
(`--barcode-tag CB --cb-len 16 --seq-len 91 --ignore-chro MT --plot-engine none --no-progress
--peak-strategy clip_seeded`), so a slice run differs from a manuscript arm only in the input BAM
and whatever extra flags you pass. `PEAKATAIL_NO_TIMESTAMP=1` is exported so the output leaf is
exactly `run/` — without it the caller appends `_<ts>` and no two runs can be compared by path.
`PYTHONDONTWRITEBYTECODE=1` is exported so a run against a frozen read-only snapshot never writes
`.pyc` into it. The scorer-ready `pas.bed` is derived with the launcher's own awk +
`normalize_chroms.py --to ensembl`.

## What `identity_check.py` treats as volatile — and what it refuses to

Three tiers, because "byte-identical" is only meaningful once you say which bytes:

* **DATA** — everything not named below. Must be byte-identical. A difference **FAILS**.
* **CONFIG** (`run_config.json`, `run_manifest.json`, `*.peak_jobs.json`) — real content
  (resolved settings, artifact hashes, the job list) wrapped in paths and clock values that can
  never match. Compared after substituting the run roots and blanking ISO timestamps, run ids,
  sha256 artifact hashes and the `wall_s` / `maxrss_mb` telemetry. A residual difference **FAILS**:
  it means a setting drifted.
* **LOG** (`peakatail_*.log`, `*.log`) — the run journal, full of elapsed times, per-worker RSS and
  nondeterministic worker completion order. Normalised (paths, timestamps, `elapsed: N.Ns`,
  `N s`/`N GB RSS`) and compared as a **sorted line set**; reported with the residual line count,
  **never fatal**. Files whose names embed the run timestamp are paired by normalised name.

`.h5ad` files are compared byte-wise first and, only if the bytes differ, dataset-by-dataset with
`h5py` (HDF5 containers are not guaranteed byte-stable). Exit status is 0 only when every file is
SAME / SAME_NORM / SAME_H5 / SAME_LOG / LOG_DIFF and neither side has an unmatched file.

Files whose names embed a run timestamp are paired by their normalised name. If two files in the
*same* tree normalise to one name a dict would silently drop one of them and the check would pass
with a file uncompared, so that case is reported as **NAME_COLLISION** and is fatal (added in the
second validation pass; no real run has yet produced one).

`--strict-volatile` disables all three normalisations and demands raw byte equality everywhere.

## Validation of the harness (2026-08-21, before any prime behaviour change)

Two runs of the **same** slice, one importing the frozen v2 snapshot `tools/pa-polya-run-9dfdefb3`
and one importing the unmodified `peakAtail-prime` branch:

```
slice_v2_ref          6:35.35 wall, 1,163,768 kB peak RSS, pas 18,865 | tier1 10,886 | tier2 7,979 | tier1>=2mol 3,643
slice_prime_compat    6:23.33 wall, 1,163,748 kB peak RSS, pas 18,865 | tier1 10,886 | tier2 7,979 | tier1>=2mol 3,643
identity_check.py --> VERDICT: IDENTICAL (54 files; 50 SAME, 3 SAME_NORM, 1 SAME_LOG), exit 0
```

Negative control on the same pair, with three deliberate corruptions of a copy — one base
coordinate moved in `pasbed.bed`, `seqlen` 91→92 in `run_config.json`, `pas_gene.tsv` deleted:

```
DIFF        pasbed.bed          849b7073 != 3f235ad1
DIFF_NORM   run_config.json     (the 91->92 survives normalisation, as it must)
DIFF_NORM   run_manifest.json   (it embeds the resolved config and artifact hashes)
ONLY_IN_REF pas_gene.tsv
VERDICT: NOT IDENTICAL, exit 1
```

Transcripts: `results/prime/identity_slice_v2_vs_prime.txt`,
`results/prime/identity_negative_control.txt`, per-file report
`results/prime/identity_slice_v2_vs_prime.tsv`.

## Second validation pass (2026-08-21, independent re-run)

Everything above was re-run from scratch rather than taken on trust, and three checks were added.
Full table: `results/prime/harness_validation_pass2.tsv`; branch state:
`results/prime/branch_baseline.tsv`.

* the full test suite was re-run on the branch: **1,277 passed**, 9 skipped, 4 xfailed, 1 xpassed,
  the 2 known `tests/test_pyproject_install.py` environment failures;
* the frozen v2 snapshot `tools/pa-polya-run-9dfdefb3` was confirmed clean at `9dfdefb`, and its
  `ema/` tree confirmed **identical** to `tools/pa-prime/ema/` (`diff -r`), so the positive control
  really is "same source, two trees, two output directories";
* the positive control was re-run (54 files, IDENTICAL, exit 0) and the four derived scorer-ready
  BEDs re-hashed (equal);
* a **second, independent negative control** (different corruptions from the first): one character
  changed inside `pas_support.tsv` *at unchanged file length* -> `DIFF`; an extra output file ->
  `ONLY_IN_NEW`; `polya_min_clip` 6->5 and `seqlen` 91->92 inside `run_config.json` -> `DIFF_NORM`
  (a settings drift must survive normalisation); exit 1
  (`results/prime/identity_negative_control_pass2.txt`);
* the `NAME_COLLISION` guard was added and exercised synthetically;
* `score_pas.sh` was re-run on the v2 IP arm and reproduced pass 1 exactly (P@10 0.5727,
  P@100 0.7392, R_det@100 0.2080, F1 0.3246, Kinnex t20 0.6011, decoy 0.0933);
* the branch gained a fixture-level compat test,
  `tools/pa-prime/tests/test_prime_v2_compat_golden.py`, which hashes every byte the caller writes
  for two arms over `tests/fixtures/cellranger_pbmc_tiny.bam` against goldens generated from the
  frozen v2 worktree. It was shown to bite: deleting the `span > seq_len` rule in `read.py`
  (Change 1 landing without a flag) fails both arms.

### The pair of slice runs quoted above was re-made from scratch in the second pass

Both arms were run again in one session — nothing was inherited — so the byte-identity claim rests on
runs whose code provenance was stamped and checked at launch:

```
pass2_v2ref       code tools/pa-polya-run-9dfdefb3 @ 9dfdefb (frozen v2)   6:48.78 wall, 1,162,864 kB
pass2_primehead   code tools/pa-prime @ 6edac79 (branch tip, 3 commits)    6:42.84 wall, 1,165,352 kB
both              pas 18,865 | tier1 10,886 | tier2 7,979 | tier1>=2mol 3,643
identity_check.py --> VERDICT: IDENTICAL (54 files; 50 SAME, 3 SAME_NORM, 1 SAME_LOG), exit 0
```

Transcript: `results/prime/identity_pass2_v2_vs_primehead.txt`, per-file report `...tsv`.

**Determinism, checked rather than assumed.** The two v2-code slice runs made an hour apart in
different sessions (`slice_v2_ref` and `pass2_v2ref`) are themselves byte-identical (54 files,
exit 0), and their derived `pas.bed` carries the same md5 `227720c57cfab4baf3e021fa2313fc68`. That is
what licenses `manuscript/24` §3.1's statement that prime-vs-v2 deltas are not noisy estimates and
need no run-to-run tolerance.

**Harness guards, exercised on their failure paths.** `run_slice.sh --threads 9` refuses (exit 2);
`assert_code` on a tree without `ema/` fails with `WRONG TREE`; `make_slice.sh` re-verified the
perf team's slice by `samtools idxstats` (19: 41,608,052 reads, 21: 9,290,404, nothing else) and
re-confirmed its md5.

*Caveat on the negative control:* it compares a **copied** tree, so the run root recorded inside the
log no longer equals the directory being compared and the journal reports a path difference. On real
runs (the only way the checker is used) the recorded root and the compared root are the same.

## Slice scoring, and why slice numbers are not manuscript numbers

`--slice-contigs 19,21` restricts every reference — atlas, TES, chrom.sizes, gene bodies, detected
atlas, Kinnex truth and decoys — to those contigs, so recall has a slice-sized denominator. The
restricted references come out at 25,675 atlas rep sites and 14,664 detected-gene atlas sites, which
matches `results/algo_headroom/A1_end_pileup/` exactly. **The slice flatters the caller** (A1 §1:
+4% precision, +19% relative recall vs genome-wide) and no headline claim may rest on it. Every
manuscript number comes from a full-BAM run on PBMC 10k v3 and both testis mice.

## Cross-check of the scorer against an independent implementation

`results/algo_headroom/A1_end_pileup/` scored the v2 default arm on this slice with its own bedtools
pipeline, independently of `score_tool.py`. Running the same arm through this harness
(`run_slice.sh --label slice_v2_ref_ipfilt ... -- --ip-filter --genome-fasta <hg38> --ip-filter-mode
filter`, then `score_pas.sh ... --slice-contigs 19,21`) reproduces it:

| | A1 (independent) | this harness | Δ |
|---|---:|---:|---:|
| P@10 | 0.5727 | **0.5727** | 0.0000 |
| P@100 | 0.7385 | 0.7392 | +0.0007 |
| R_det@100 | 0.2085 | 0.2080 | −0.0005 |
| F1_det@100 | 0.3252 | 0.3246 | −0.0006 |
| Kinnex t20 P@25 | 0.5990 | 0.6011 | +0.0021 |
| Kinnex decoy@25 | 0.0929 | 0.0933 | +0.0004 |
| n | 2,895 | 2,883 | −12 |

The 12-call difference is expected and is **not** a harness error: A1 restricted a **genome-wide**
v2 run to contigs 19+21, whereas `run_slice.sh` runs the caller on a **2-contig BAM**, so the
caller's global per-cell filters (`min_read 1500`, `min_pas_per_cell`, gene detection) see 7 % of the
library and retain a slightly different cell set. **A slice run is not a genome-wide run restricted
to the slice** — compare slice runs with slice runs, and never quote a slice number as a
genome-wide one. Full table: `results/prime/harness_validation.tsv`.

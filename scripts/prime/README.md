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

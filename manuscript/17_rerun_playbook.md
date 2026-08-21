# Re-run playbook — after Amir merges the Stage 1b/1c + Stage 1d PRs (PI directive 2026-08-21)

Goal: pull the merged PeakATail version, freeze it, and redo every analysis whose numbers depend on
the caller, **against the unchanged pre-registered gates** (13 §1–3, 14, 15, 16). Nothing below
re-opens a gate or a definition; only the code commit changes. Do not start until the IP-filter
(Stage 1d-A) and memory/CPU (Stage 1d-B) PRs are merged into `develop`; run order is fixed.

## 0. Pull and freeze
```bash
cd /mnt/ssd1/Projects/PeakATail_wd
git -C tools/PeakATail fetch origin && git -C tools/PeakATail log --oneline origin/develop -5   # confirm the merges
export SRC=tools/PeakATail SRC_COMMIT=$(git -C tools/PeakATail rev-parse origin/develop) OUT_TAG=final_v2
# the launcher creates the detached worktree tools/pa-polya-run-<sha8> itself; record the sha in 06_roadmap.md
```
Run the tool's test suite from the new worktree first (`PYTHONPATH=tools/pa-polya-run-<sha8> .venv/bin/python -m pytest -q tests`;
expect only the 2 known `test_pyproject_install` env failures).

## 1. Stage-2 benchmark (4 arms, ~4 h wall, ≥300 GB free RAM unless 1d-B landed)
```bash
LC_ALL=C nohup bash scripts/benchmark_tools/stage2_final_launch.sh > results/benchmark_tools/stage2_final_logs/launch_v2_$(date +%Y%m%d_%H%M).log 2>&1 &
```
Outputs: `results/benchmark_tools/pbmc_10k_v3/peakatail_clipseeded_final_v2{,_ipfilt}`, `.../gse104556/peakatail_clipseeded_final_v2/mouse{1,2}`.
Then the gate verifier (same checks as 15: scorer re-run byte-identical, own bedtools pipeline, provenance,
tier tag, IP-flag re-application with the **corrected** strand rule → expect 0 newly flagged) and a
`15_final_gate_v2.md`. Expected movement (from the post-hoc estimate in `github/pr_stage1d_ipfilter_strand.md`):
PBMC default P@100 0.717 → ~0.706 with ~+2,100 sites; mice ≈ +0.4 pp; gates unchanged. Any other movement
must be explained (1d-B is required to be output-identical).

## 2. Kinnex trusted-novel (minutes)
```bash
LC_ALL=C tools/PeakATail/.venv/bin/python scripts/reliability/trusted_novel_pas.py call ... --pas-bed results/benchmark_tools/pbmc_10k_v3/peakatail_clipseeded_final_v2_ipfilt/pas_PRESPEC_precision_default.bed ...   # same args as results/reliability/trusted_novel_final_pbmc/01_call.sh
```
Same pre-registered definition; report `16_trusted_novel_kinnex_v2.md`. If a v2 definition is to be tested
it must be pre-registered *before* this step and validated on a held-out truth (see 16 §"What the paper says").

## 3. Stage 3 Laughney (cohort mode; ~10 h peak calling + ~8 h switch queue + replication)
Re-run `scripts/stage3/` with the new snapshot: cohort_run (`launch_cohort.sh` with the new PYTHONPATH),
label confirmation is **unchanged** (labels are GEX-derived; reuse `labels/confirmed_labels.tsv`), universe
rebuilt from the new sidecars (`stage3_build_universe.py`), switch queue (`run_cohort_switch_all.sh`),
replication (`run_replication.sh`, patient unit, MetBone excluded in primary). Verifier scripts live in
`/mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v2/tmp/verify_switch/`.

## 4. Spermatogenesis control + switch calibration on the new caller
Re-run the mouse1/mouse2 `switch length` SPC>RS>ES control and the 20-permutation calibration
(`results/fdr_calibration_v2/` drivers) on the v2 mouse clusters.h5ad with Fisher cells-mode, `--marker-top-n 0`.

## 5. Figures and text
Regenerate `final_benchmark`, `trusted_novel_funnel`, `fdr_calibration_v2` and the Stage-3 figures from the v2
outputs (scripts read the score TSVs; point them at `*_final_v2*`), then replace the v1 numbers in
`01_outline_and_journals.md` — every number must again trace to a verified `*_v2.md`.

## What does NOT need re-running
Competitor runs (Sierra, polyApipe, SCAPTURE, scAPAtrap, scUTRquant), the Kinnex truth sets, the curated
atlas references, the label confirmation, and the calibration *design*. The v1 results stay in the repo as
the record of the pre-fix caller.

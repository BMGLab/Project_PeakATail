## Motivation

For the manuscript submission, reviewers must be able to (a) see the test suite pass automatically on every change, and (b) cite the software from a machine-readable record. Today `.github/workflows/` only deploys docs — the pytest suite runs nowhere automatically — there is no `CITATION.cff`, and the README's command table lists 8 of the 12 actual `ema` commands. This PR fixes all three. Closes #__RELENG_ISSUE__ — the remaining release items there (PyPI, bioconda, Zenodo DOI) are release-cutting steps that need maintainer credentials, to be filed/handled separately at tag time.

## What changed

- **`.github/workflows/ci.yml`** (new): runs the pytest suite on pushes and PRs targeting `develop` (plus `workflow_dispatch`), matrix over Python **3.11 / 3.12**, `pip install -e '.[test]'`, pip caching via `actions/setup-python`, `bedtools` installed so the bedtools-guarded tests actually run instead of skipping. Per-test 300 s timeout comes from the existing `[tool.pytest.ini_options]` via `pytest-timeout`; `timeout-minutes: 45` bounds each job.
- **`pyproject.toml`**: adds a minimal `[test]` extra (`pytest`, `pytest-timeout`) so CI doesn't pull the lint/type-check toolchain; the existing `dev` extra is unchanged.
- **`CITATION.cff`** (new): CFF 1.2.0, title *PeakATail*, version 0.2.0, MIT, repo URL. Authors are "PeakATail developers (BMGLab)" as a collective plus Amir Amiri Tabat — an in-file note marks the author list as a placeholder to be finalized with the manuscript; no ORCIDs were invented. Validates clean with `cffconvert --validate` (schema 1.2.0).
- **`README.md`**: the command table now lists **all 12 commands** (adds `ema reannotate`, `ema switch trend`, `ema switch combine`, `ema collapse`; cross-checked against the Click registry in `ema/cli/__init__.py` + `ema/cli/switch.py` and `docs/cli/index.md`), and the four `switch` doc links are fixed from `cli/switch/<cmd>/` to the real `cli/switch-<cmd>/` pages. Commands without a dedicated docs page yet (`trend`, `combine`, `collapse`) link to the CLI reference index.

No runtime code is touched — the diff is CI config, packaging metadata, and docs.

## Test evidence

Full suite run locally against this branch exactly as CI will run it (module path asserted to resolve to this branch before collection):

```
MODULE UNDER TEST: .../PeakATail-relwork/ema/__init__.py   (asserted)
5 failed, 925 passed, 9 skipped, 1 xfailed, 230 warnings in 90.78s
```


**The 5 failures are pre-existing on `origin/develop` (18678ef) and untouched by this branch** — this PR changes no runtime code or tests, and each failure is a test/code drift already present at HEAD (e.g. `test_downstream_parallel` passes `per_dataset_dir=`, a kwarg `run_one_dataset_downstream()` no longer accepts):

- `tests/test_cli_switch.py::test_isoform_agg_per_gene_dispatches_to_per_gene_branch` (`_FakeAdata` stub lacks `.var`)
- `tests/test_downstream_parallel.py::TestRunOneDatasetDownstream::{test_raises_on_empty_sub_indices, test_full_pipeline_calls_with_mocks}` (stale `per_dataset_dir=` kwarg)
- `tests/test_ip_annot_filter_d6.py::{test_ip_filter_removes_internally_primed_peak_keeps_clean_one, test_ip_filter_default_mode_annotates_instead_of_dropping}` (`KeyError: 'filtered'`)

The first CI run on `develop` will therefore be red on exactly these 5 until they're fixed — that visibility is the point of adding CI; fixing them is a small follow-up for whoever owns the D6/downstream test contracts. (The local `xfail` is `test_entry_point_runs`, an artifact of a stale `ema 0.1a1` console script in the shared venv; a fresh CI `pip install -e '.[test]'` wires `ema.cli:main` and it passes there.)

## Follow-ups (separate, need maintainer credentials)

Tag-triggered PyPI release workflow → bioconda recipe → Zenodo DOI on `v0.3.0`.

---
*Opened by Claude (AI assistant) for the manuscript effort; the gh session on this machine authenticates as @yasinkaymaz. Contents (workflow, CFF, README table) were reviewed against the live CLI registry before posting.*

🤖 Generated with [Claude Code](https://claude.com/claude-code)

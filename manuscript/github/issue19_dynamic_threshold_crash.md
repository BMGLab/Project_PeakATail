# Issue: `--dynamic-threshold` at its documented setting crashes the caller with a bare `IndexError` — nothing bounds the look-back index by the live window

**Labels:** bug, crash, peak-calling. Found by the 2026-08 parameter sweep
(`results/paramsweep/VERDICTS.md` §5) and reproduced by its independent accuracy verifier
(`results/paramsweep/VERIFY/rerun_identity.tsv`), against the frozen v2 tree `9dfdefb`.

Related: **#95** (the sweep was run to close out the parameter reference's tuning map), **#99**
(same class: what a flag's help text says vs what the code path does).

## Symptom

`--dynamic-threshold --lambda-fold-change 2.0` — the parameter reference's own first entry under
"find more PAS" — aborts the run:

```
IndexError: list index out of range
```

raised from inside a spawned chromosome worker (`chrom_parallel.py`), with no message naming a
flag, a contig or a remedy. The sweep lost a whole arm to it before anyone read the source.

## Exact location (v2 `9dfdefb` line numbers)

`ema/countmatrix/peackcalling.py:528`:

```python
l_end = data_array[-current_threshold]
```

`data_array` is the live window of read 3' ends; with `--dynamic-threshold` the threshold is
recomputed per read (`peackcalling.py:518–525`):

```python
current_threshold = max(floor_threshold, int(local_lambda * lambda_fold_change))
```

and **nothing bounds it by `len(data_array)`**. When the local read density pushes the threshold
past the number of ends still in the window, the negative index runs off the front of the list.
The 3-stage pipeline's finder carries the identical expression at `peak_pipeline.py:313` (its
threshold recompute is `peak_pipeline.py:304–310`).

Reachability: `--lambda-fold-change` is read at exactly two places in the tree, both inside
`if dynamic_threshold:`, and `--dynamic-threshold` is off by default — so **no published PeakATail
number is affected**; the defect is real for exactly the users the tuning map sends there.

## Reproduction

* **Real data (sweep provenance, `VERDICTS.md` §5):** the PBMC chr19+21 dev slice with
  `--peak-strategy clip_seeded --dynamic-threshold --lambda-fold-change 2.0` aborts as above. The
  sweep's table initially said the same arm "completes on mouse"; the verifier re-ran it and
  showed that claim FALSE — the GSE104556 mouse1 chr18+19 slice crashes at the same line
  (`VERIFY/rerun_identity.tsv`), so the defect is **species-independent**, not a property of one
  library's depth.
* **Minimal synthetic (14 reads):** one forward-strand cluster whose starts sit ~80 bp before
  their ends. With `default_threshold = floor_threshold = 3` the signal fires at read 3; at read
  11 the background deque tops the estimator's 10-read minimum and the recomputed threshold jumps
  past the live window's 10 ends → `data_array[-N]` with `N > 10` → `IndexError`. The BAM builder
  and both loop paths (monolithic and pipeline) are committed as
  `tests/test_prime_dynamic_threshold_clamp.py` on the `peakAtail-prime` branch; run against the
  frozen v2 tree the same fixture raises the identical `IndexError` at `peackcalling.py:528`.
* `--dynamic-threshold --lambda-fold-change 1.2` runs (PBMC 15,925 → 15,670 PAS — *fewer*, and
  the candidate ceiling gets worse), which is the sweep's separate FAIL verdict on the lead
  itself: unusable at the documented setting, worth +0.00014 recall at the setting that runs.

## Proposed fix — contained in the `peakAtail-prime` PR

`--dynamic-threshold-clamp {off,on}`, default `off`:

* `off` (default) is v2 to the character — same expression, same abort — because that branch's
  cardinal rule is byte-for-byte v2 reproducibility and a fix that changes what a run emits needs
  its own gate first.
* `on` bounds the look-back index to the live window (returning the oldest end still in it) and
  logs a census of how often it fired. The single copy of the rule
  (`ema.countmatrix.dynamic_threshold.resolve_l_end`) carries an exhaustively-tested identity
  guarantee: for every in-range threshold both settings return the same element, so `on` can only
  change a run that would otherwise have aborted.

If a maintainer prefers the clamp unconditional on `develop` (there is no byte-identity contract
there), the guard module drops in as-is; the flag is the branch's constraint, not the fix's.

## Two documentation defects the same sweep proved — filed here so they travel together

(`VERDICTS.md` §8.3; both byte-identity results, both species:)

1. **`--pas-gap` does nothing on a single-BAM run**, despite its schema help text ("minimum gap
   between PAS within a peak"). It is consumed only by `merge_pas_beds` in the multi-dataset
   unified path (`main.py:1516`); `--pas-gap 25` and `200` are byte-identical to the baseline.
   The help text should say it is a *multi-dataset merge* parameter, or the flag should be moved
   to that command's surface.
2. **`--min-cells` and `--min-pas-per-cell` filter the AnnData, not the call set.** Both act in
   `preprocessing()`, after `pasbed.bed` is written; `--min-cells 1` and `--min-pas-per-cell 10`
   leave `pasbed.bed` byte-identical. Any reading of them as call-set filters — including in a
   parameter sweep — is wrong, and the docs should say which outputs they touch.

## Evidence

`results/paramsweep/VERDICTS.md` §5 (the crash), §4 and §6–8 (the inert-flag byte-identity table),
`results/paramsweep/VERIFY/rerun_identity.tsv` (the mouse reproduction the sweep itself missed);
fix and regression tests: `peakAtail-prime` commit `1f80f9b`.

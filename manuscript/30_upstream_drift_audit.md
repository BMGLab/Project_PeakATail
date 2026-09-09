# Upstream drift audit: what changed in the tool after the frozen commits

Run 2026-09-09, after `origin/develop` moved 59 commits past the manuscript's v2 commit
`9dfdefb`. Four independent audit lenses over MAIN.md and SUPPLEMENTARY.md against the live
repo, each finding adversarially verified. 34 corrections survived verification, 14 were
refuted; 18 non-overlapping edits applied (the other 16 are subsumed by a longer correction
covering the same passage).

## The decision: STAY FROZEN at 9dfdefb. Do not re-run.

No merged upstream fix invalidates any number the paper reports. The decisive fact, verified
at both refs rather than argued:

`pasbed.bed` is written by `find_close` at `pasbed.saveas(mergebed)` **before**
`pasbed.closest(...)` and before the tier filter — at 9dfdefb (`find_close.py` L91-97) and
at origin/develop (L301 saveas, L310 closest, L371 tier filter). The benchmark's scored call
sets are built by awk over that file (`scripts/benchmark_tools/stage2_final_launch.sh:46-56`).
**The annotation stage therefore cannot touch the scored call set.** Every atlas-agreement
precision, every call count (46,524 / 26,255 / 26,526 / 20,672), the Kinnex concordances,
the cross-donor concordances and the matched-budget curves are provably immune to PR #104.

Per merged PR:

| PR | What it changed | Effect on a reported number |
|---|---|---|
| #100 | peakAtail-prime defaults; 0.2.0 -> 0.3.0 | none — its own identity check proves prime's bare default is md5-identical to the v2 manuscript arm on `pas.bed`/`pas_tier1.bed`/`pas_tier2.bed` across all four libraries |
| #102 | `--dynamic-threshold` bounds (#101) | none — the flag is used by no reported run |
| #104 | PAS->gene assignment in overlapping loci (#99) | none on scored call sets, per the ordering above; would change gene-level *labels*, which is why no ranked gene list is published |
| #105 | `--marker-top-n` default 200 -> 0; Fisher denominator from the full matrix (#94, partial) | none — every reported run already passed `--marker-top-n 0` explicitly |
| #106 | `per_isoform` degenerate pairs (#98) | none — no reported result uses that arm |
| #107 | `reannotate` atomic write | none |

Re-running would additionally invalidate the paper's calibration evidence, which rests on 20
label-permutation nulls per mouse computed at `9dfdefb`, and would break the pre-registration
that names that commit.

## What the audit corrected

- **The count and the tense were both wrong.** MAIN.md said "Two defects are open publicly"
  and named #99 and #98. Both are CLOSED/COMPLETED (4 and 3 September 2026). They are still
  present in *this paper's results*, because every fix merged after `9dfdefb` — verified with
  `git merge-base --is-ancestor` for 71e3ade, 3f77717, b25843c and c978e81, all of which are
  ancestors of develop and none of 9dfdefb. So the caveats stay; only the tense changes.
- **#99 has two halves, and the draft conflated them.** PR #104 fixed the gene assignment;
  PR #100 fixed the separate clip-rate head-sampling estimator. The #104 fix also carries a
  residual the maintainer records: 96 of 465 relabels still fall back to terminus proximity.
  That residual, not an open issue, is what now justifies withholding the ranked gene list.
- **#94 is the one cited issue still open**, and the draft did not say so. PR #105 closed 2
  of its 6 checkboxes and the maintainer deliberately left it open. It bears on **no** result
  here: every reported run used `--marker-top-n 0`.
- **Defaults on develop have moved** and the paper described them in the present tense:
  `marker_top_n` 200 -> 0, `ip_filter_mode` `annotate` -> `auto`, `cleavage_offset` 0 -> none.
  These are now pinned to the frozen commit rather than stated as current behaviour.
- **prime is no longer "an unmerged branch"** — PR #100 merged it.
- **0.3.0 is bumped in metadata but NOT released.** `git tag` shows only `0.1.a1`/`0.1.a2`;
  the develop CHANGELOG still says "convert these `Unreleased` headings at tag time". The
  corrected text says "version-bumped to 0.3.0 on develop", never "released as 0.3.0".

## Two live defects the manuscript does not cite

Both are ours, opened during this work, and both postdate the frozen commits:

- **PR #109** (open) — `--isoform-agg within_utr` aborted on every invocation in every
  released version. Relevant to the Fig 5 3'UTR-scoped caveat.
- **Issue #110** (open) — `between_utr` output leaves `gene_id`, `chrom`, `start`, `end`,
  `strand` empty.

## Method note

Applying the audit mechanically was itself unsafe: several verifier `corrected_proposal`
fields contained reviewer commentary ("EDIT 1 — FROM: … TO: …") rather than replacement
prose, and a first pass pasted that into the draft. The applier now rejects any replacement
matching a commentary pattern and falls back to the original proposal, with a hard gate that
re-scans both files afterwards. Numeric tokens were diffed before and after: the only ones
removed are the issue numbers 98 and 99 moving into new phrasing.

# Companion repository: state, blockers, and the runbook to public + DOI

**Repository:** https://github.com/BMGLab/Project_PeakATail
**Assessed:** 2026-09-02. Every status below was checked mechanically; the
command that produced each finding is shown so it can be re-run.

**Nothing in this document has been executed.** No commit, no push, no `gh`
write, no visibility change. All GitHub calls made while writing it were
read-only (`gh api repos/...`, `git ls-remote`).

---

## 1. Where the repository actually stands

```
gh api repos/BMGLab/Project_PeakATail
  -> visibility  PRIVATE
     license     null            (no LICENSE file at all)
     defaultBranch  main
     createdAt   2026-08-19T19:41:07Z
     pushedAt    2026-08-19T19:43:10Z     <-- untouched for two weeks
```

43 commits on `main`: 41 filtered upstream commits plus two packaging commits,
`081bdb0` ("Publication packaging: source data, environments, accession
manifest") and `bc0eae1` ("Add sync script…").

The last upstream commit it carries is `ca3200a`, the filtered image of the
working directory's `17de20f` (2026-08-19 22:18). **The working directory has
moved 84 commits since then**, up to `f5aa77b` (2026-09-02 23:21) — the entire
manuscript draft, the final figure roster, the v2 re-runs, the second-donor
test, the prime benchmark and Stage 3 final. The public repository is two weeks
and a whole manuscript behind.

The history is a `git-filter-repo` rewrite and **cannot be updated by a plain
pull**. `maintenance/sync_from_working_dir.sh` is the only supported path.
`git-filter-repo` 2.47.0 is installed at `/home/biolab/.local/bin/git-filter-repo`.

---

## 2. Blocker status — the README's "Before release" list

The repository's own README lists five items. Two are now **resolved by
evidence**, one is resolved **only if the PI accepts a proposed file**, and two
remain genuinely open.

### Blocker 1 — placeholder authors in `CITATION.cff` · **OPEN — file drafted, PI decision required**

Current file names four authors, led by the placeholder entity
`"PeakATail manuscript authors (BMGLab)"`, and includes **Ebru Koçakaya**, who
is not on the manuscript author list the PI set on 2026-09-02.

A validated replacement is staged at **`submission/CITATION.cff.companion`**
(parses as YAML; passes the CFF 1.2.0 JSON schema with zero errors).

**This needs a PI decision before it is used.** Software authorship and paper
authorship may legitimately differ, and the git record is lopsided:

```
git shortlog -sne --all
  111  ebrukocakaya <kocakayaebru@gmail.com>
   10  yasinkaymaz  <yasinkaymaz@gmail.com>
    2  Yasin Kaymaz <yasinkaymaz@gmail.com>
    1  Sercan Ozturk <n.sercan.ozturk@gmail.com>
    1  ebrukocakaya <sevtapduman@gmail.com>
```

The staged file follows the README's own instruction ("match the manuscript
author list") and therefore lists A.A.T. and Y.K. only. **Removing a named
person from a citation record is not a decision an assistant may make**, so it
is raised here rather than applied. Whatever is decided, this file is what
Zenodo reads when it mints the DOI — it becomes the authorship record of the
archived software.

*Who:* PI decides; maintainer commits.

### Blocker 2 — no `LICENSE` file · **OPEN — PI legal decision**

Confirmed: `gh api` reports `license: null`; there is no `LICENSE` in the tree.
README says "MIT (see `CITATION.cff`; a `LICENSE` file is still to be added)"
and `CITATION.cff` declares `license: MIT` — **a licence claim with no licence
text**, which is the weakest possible position and, separately, does not satisfy
Genome Biology's OSI-compliance requirement.

See **`submission/LICENSE_DECISION.md`** for the candidates, verified comparator
evidence, a recommendation and copy-ready texts
(`submission/LICENSE.CANDIDATE-*.txt`). That brief also documents a second,
independent defect it found: the **tool** repo publishes MIT on `develop` and
GPL-3.0 on `main` simultaneously.

*Who:* PI (and, if Ege University asserts ownership of staff software, the
university's technology-transfer/legal office).

### Blocker 3 — benchmark pin `18678ef` on the unpushed `biolab-manuscript` branch · **RESOLVED**

The premise no longer holds. `18678ef` is publicly reachable and has been for
some time:

```
cd tools/PeakATail
git branch -r --contains 18678ef        # 19 remote branches
  -> origin/develop  (the default branch), origin/peakAtail-prime,
     origin/feat/polya-evidence, origin/fix/ip-filter-strand, and 15 others

git merge-base --is-ancestor 18678ef origin/develop   -> YES
git merge-base --is-ancestor 18678ef 4efeb125          -> YES   (v1)
git merge-base --is-ancestor 18678ef 9dfdefb           -> YES   (v2)
```

`18678ef` is an **ancestor of both frozen manuscript commits**, and both of
those are themselves on the public default branch:

```
git ls-remote --heads origin
  refs/heads/develop          5073d63e…      (origin/develop; contains v1 and v2)
  refs/heads/peakAtail-prime  5e18e6f0…      (pushed; PR #100 head)
```

Nothing needs to be pushed. `biolab-manuscript` (local tip `6697b0d`) is still
unpushed, but that is now irrelevant — everything the benchmark pins is
reachable without it. **Do not push `biolab-manuscript` to fix this**; it would
add a stale branch for no benefit.

*Action:* delete the blocker from the README; correct `envs/TOOL_VERSIONS.md`,
which still says the pin lives on `biolab-manuscript` (see §3).

### Blocker 4 — "do the three post-benchmark read-module fixes (up to `2e7fc0d`) change any reported number?" · **RESOLVED — no**

`envs/TOOL_VERSIONS.md` warns that the benchmark ran before three read-module
fixes and that `source_data/` therefore predates them. **That is false.** The
three fixes exist on `origin/develop` under different SHAs (they were re-landed
through PRs, so the local `biolab-manuscript` SHAs are not the ones that
shipped), and each is patch-identical to its local twin and an ancestor of both
frozen manuscript commits:

| local SHA (unreachable) | shipped SHA on `origin/develop` | `git patch-id` | ancestor of v1 `4efeb125` | ancestor of v2 `9dfdefb` |
|---|---|---|---|---|
| `ad59cdd` strip CellRanger GEM-group suffix from CB tags | `9de09a5` | IDENTICAL | YES | YES |
| `6df8eed` skip unmapped/CIGAR-less reads in `read_check` | `19790c8` | IDENTICAL | YES | YES |
| `2e7fc0d` underscore-bearing RG tags corrupt `sample_cb` | `eb13529` | IDENTICAL | YES | YES |

Reproduce:

```bash
cd tools/PeakATail
for pair in "ad59cdd 9de09a5" "6df8eed 19790c8" "2e7fc0d eb13529"; do set -- $pair
  a=$(git show $1 | git patch-id --stable | cut -d' ' -f1)
  b=$(git show $2 | git patch-id --stable | cut -d' ' -f1)
  echo "$1 vs $2: $([ "$a" = "$b" ] && echo IDENTICAL || echo DIFFERENT)"
  git merge-base --is-ancestor $2 9dfdefb && echo "  in v2: YES" || echo "  in v2: NO"
done
```

So every number in the manuscript was computed on code that **already contains**
all three fixes. Worth noting for its own sake: the local SHAs `ad59cdd`,
`6df8eed` and `2e7fc0d` are on **no branch at all**, local or remote
(`git branch -a --contains 2e7fc0d` is empty). They survive only in this working
copy's object store, so the README's citation of `2e7fc0d` is a reference to a
commit no one else can resolve — another reason to delete the item rather than
carry it forward.

*Action:* delete the blocker from the README and delete the "Caveat on the
PeakATail pin" section from `envs/TOOL_VERSIONS.md`.

### Blocker 5 — flip public, enable Zenodo, cut a release · **OPEN — gated on 1 and 2**

Verified: Zenodo's GitHub integration lists **only public repositories**, so the
flip must precede the DOI. Steps in §5.

*Who:* PI / a BMGLab organisation admin.

---

## 3. What `source_data/`, `envs/` and the metadata files need

The sync script rebuilds history from the working directory, but the working
directory **git-ignores `results/`** (`.gitignore` line 4), so figure source
tables are not carried by the sync. They must be refreshed by hand — this is
the script's own closing note, and it is mechanically necessary:

```
git ls-files results/figures/manuscript/ | wc -l   ->  0
```

### `source_data/` — 70 files missing, 3 stale, 0 orphaned

Compared `results/figures/manuscript/*.tsv` (104 files, 14 MB) against
`source_data/*.tsv` (34 files):

- **70 new tables absent from the repo.** The entire final figure roster:
  `fig1_overview*` (5), `fig2_accuracy*` (3), `fig3_tradeoff*` (3),
  `fig4_calibration*` (6), `fig5_spermatogenesis*` (8), `fig6_cohort*` (8),
  `figS1_datasets*` (2), `figS2_peakqc*` (7), `figS3_nulldesign*` (2),
  `figS4_seconddonor*` (3), `figS7_novelfunnel*` (5), `figS8_compute*` (3),
  `figS9_calibration_extended*` (6), `figS10_clustering*` (2),
  `figS11_gatehistory*` (2), `figS12_versions*` (4), plus
  `kinnex_truth_validation.tsv`.
- **3 tables changed** and must be overwritten: `benchmark_consolidated.tsv`,
  `benchmark_headtohead.tsv`, `pbmc_novelty.tsv`.
- **31 unchanged.** **0 orphaned** — every table in the repo is still produced,
  including those backing the 9 retired figures.
- `source_data/FIGURE_SOURCE_README.md` is a copy of
  `results/figures/manuscript/README.md` and is now stale (md5
  `5f273c2b…` vs `b95cc05e…`). Re-copy it.
- `source_data/README.md` (hand-written) still describes the **retired** figure
  names as the examples — `benchmark_headtohead.tsv` backing
  `manuscript/figures/benchmark_headtohead.pdf`. Rewrite its "Notable tables"
  section against the final `fig1…fig6` / `figS1…figS12` roster.

A plain `cp` of all 104 tables plus the README is the whole refresh:

```bash
cp "$WD"/results/figures/manuscript/*.tsv                source_data/
cp "$WD"/results/figures/manuscript/README.md            source_data/FIGURE_SOURCE_README.md
```

### `envs/` — re-export mechanically possible, and two factual corrections needed

All six exported environments still exist on the machine
(`conda env list`: `bench_polyapipe`, `bench_scapatrap`, `bench_scapture`,
`bench_scutrquant`, `bench_sierra`, `bench_star`), so re-export is a loop:

```bash
for t in polyapipe scapatrap scapture scutrquant sierra star; do
  conda env export -n "bench_$t" --no-builds > "envs/bench_$t.yml"
done
```

**Re-export only if a tool environment actually changed since 2026-08-19.** The
existing exports are the ones the benchmark ran under; re-exporting a drifted
environment would replace a true record with a false one. If nothing changed,
leave the six `.yml` files alone — they are correct as they stand.

`envs/TOOL_VERSIONS.md` needs two edits that are **not** optional:

1. The PeakATail row reads
   `commit 18678ef (branch biolab-manuscript)`. Replace with the manuscript's
   own frozen commits — `4efeb1252e6d7c7d51b252b443b5be5947ae000f` (v1) and
   `9dfdefb3eb353b0817ef79c4eb9ace6d6c8aab53` (v2, the manuscript's record) —
   and state that both are on the public default branch `develop`.
2. Delete the "Caveat on the PeakATail pin" section outright (§2, blocker 4).
   As written it tells readers the published numbers omit three bug fixes they
   in fact include — the single most damaging wrong sentence in the repository.

There is no `envs/` entry for PeakATail itself ("installed from source"). The
tool repo carries a `uv.lock` (410 KB) at each frozen commit, which *is* the
exact dependency pin; consider naming it in the table instead of "installed
from source".

### `DATA_ACCESSIONS.md` — three stale numbers

It is otherwise excellent (it is the source that will fill the availability
statement's dataset paragraph), but three values predate the v2 re-run and
disagree with the manuscript:

| field | `DATA_ACCESSIONS.md` | manuscript (MAIN.md Methods) |
|---|---|---|
| detected genes, `pbmc_10k_v3` | 14,851 | **14,949** |
| recall denominator `pas2.in_detected_genes.bed` | 285,220 sites | **285,136** sites |
| GSE123904 scale | "17 patients" | **17 libraries from 14 patients** |

The manuscript's values are the audited ones (NUMBERS_LEDGER.tsv). Fix the
manifest, not the paper. The third is a substantive correction, not a typo: the
cohort has 17 libraries from 14 patients, three of whom contributed two
libraries.

### One gap the manifest cannot fill by itself

The Kinnex row records the **source URL**
(`https://downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/`, sets
`primary_10x3p` and `secondary_gemx3p`) but **no accession** — PacBio
distributes these from its own cloud, not from SRA/GEO. MAIN.md still carries an
open `[[CITE]]` for a Kinnex data citation "with accession/identifier"
(reviewer point M8). If no accession exists, the honest resolution is to cite
the URL, the two set names and the access date, and say so — **not** to invent
an accession. See `availability_statement.md`, which is written that way.

---

## 4. Pre-publication audit (already run — results below)

Run against the 339 tracked files that would actually be published (i.e.
excluding the three filtered paths).

- **Secrets: clean.** No token, key or credential pattern
  (`ghp_`/`gho_`/`github_pat_`/`AKIA`/`BEGIN … PRIVATE KEY`/`api_key=`/
  `password=`/`secret=`) in tracked non-filtered content.
- **E-mail in file content: clean.** The only real address is
  `yasin.kaymaz@ege.edu.tr` (×2); every other hit is an R idiom
  (`seu@meta.data`, `laugh@meta.data`).
- **E-mail in git author metadata: four personal Gmail addresses become
  public with the history** — `kocakayaebru@gmail.com`,
  `n.sercan.ozturk@gmail.com`, `sevtapduman@gmail.com`,
  `yasinkaymaz@gmail.com`. This is normal for an open repository and is not a
  blocker, but the three people who are not the PI should be told before the
  flip. If any of them objects, the fix must happen *during* the filter step
  (`git filter-repo --mailmap`), not afterwards — the history is immutable once
  it is public and archived.
- **File sizes: safe after filtering.** Exactly one tracked file exceeds
  GitHub's 100 MB hard limit — `archive/amirtest/macro/emaout/negmatrix.mtx`
  (451.33 MB) — and it is the file the filter exists to remove. Nothing else
  exceeds even the 50 MB warning threshold. Largest survivor:
  `archive/amirtest/macro/emaout/negbed.bed` at 25.24 MB, also filtered; then
  `scripts/laughney/202511_laughney_manipulations.ipynb` at 2.50 MB.
- **Absolute machine paths:** 141 tracked files reference `/mnt/ssd0|1|2`.
  Already disclosed in the README ("figure scripts contain absolute paths…
  `source_data/` exists so that the plotted numbers are inspectable without
  either"). No action needed; do not strip them, the disclosure is more honest.
- **New root file since the last sync:** `PeakATail_Parameters_Reference.pdf`
  (100 KB, added 2026-08-23). It will become public. Confirm with A.A.T. that
  he is content to publish it.

---

## 5. The runbook

`WD=/mnt/ssd1/Projects/PeakATail_wd`. Phases are ordered by dependency; do not
reorder. Everything through Phase 3 is reversible; Phase 4 onward is not.

### Phase 0 — decisions that gate everything (PI)

| # | Decision | Input |
|---|---|---|
| 0.1 | Licence for **both** repositories | `submission/LICENSE_DECISION.md` |
| 0.2 | Copyright holder line | `LICENSE_DECISION.md` §6 |
| 0.3 | Companion-repo authorship: does the analysis record credit Ebru Koçakaya (and Sercan Özturk)? | §2 blocker 1 |
| 0.4 | Version tag for the archived release, e.g. `v1.0.0-manuscript` | — |
| 0.5 | Are the three non-PI contributors content for their commit e-mails to be public? | §4 |

Nothing below can be finished without 0.1–0.4.

### Phase 1 — sync the content (maintainer, on this machine)

1. **Commit the working directory first.** The sync clones the WD repo, so only
   committed, tracked content travels. At the time of writing, uncommitted:
   `.vscode/settings.json`, `manuscript/figures/fig4_calibration.pdf`,
   `scripts/manuscript_figures/_molsweep.py`, `_pubstyle.py`, `submission/`.
   Decide for each whether it should be public. (`submission/` — this directory —
   contains PI-facing decision briefs; it is probably *not* wanted in a public
   repo. Consider `.gitignore`-ing it or keeping it unstaged.)

2. **Ensure temp space.** The script clones the full 229 MB `.git` including the
   451 MB blob, then filters it. Point `TMPDIR` at a disk with ≥2 GB free.

3. **Clone the publication repo fresh** (do not reuse a stale checkout):

   ```bash
   git clone git@github.com:BMGLab/Project_PeakATail.git /path/to/pub
   cd /path/to/pub
   git log --oneline -1        # expect bc0eae1
   ```

4. **Run the sync:**

   ```bash
   maintenance/sync_from_working_dir.sh "$WD"
   ```

   It clones the WD, re-filters with
   `--invert-paths --path archive/amirtest --path amirtest --path tools/PeakATail`,
   fetches the filtered tip of `reorg-manuscript`, and rebases the two packaging
   commits onto it. It never pushes.

5. **Verify determinism before trusting the rebase.** filter-repo is
   deterministic only for a fixed version and argument set; if SHAs shifted, the
   merge-base would be wrong and the rebase would replay the wrong range. Check
   that the previous sync point survived:

   ```bash
   git merge-base --is-ancestor ca3200a HEAD && echo "OK: history is a superset"
   git log --oneline | wc -l     # expect exactly 127
   git log --oneline -3          # expect bc0eae1-equivalent on top
   ```

   **127 is exact, not an estimate.** `reorg-manuscript` carries 125 commits;
   the previous sync pruned 0 commits as empty (41 in, 41 out); and none of the
   84 new commits touches *only* filtered paths, so none becomes empty either
   (checked commit by commit). 125 + 2 packaging = 127. A different number means
   the filter behaved differently from last time — investigate before pushing.

   If `ca3200a` is *not* an ancestor, stop — the filter has drifted and the
   history must be rebuilt deliberately rather than rebased.

   Conflict risk is low: `README.md` and `.gitignore` are the only files both
   sides touch, and neither changed upstream in the 84 new commits (223 files
   added, 14 modified).

6. **Refresh `source_data/` and, only if an environment drifted, `envs/`** — §3.

### Phase 2 — publication metadata (maintainer, one commit)

7. Add `LICENSE` — copy the chosen `submission/LICENSE.CANDIDATE-*.txt` and
   substitute the copyright line per `LICENSE_DECISION.md` §6.
8. Replace `CITATION.cff` with `submission/CITATION.cff.companion`, adjusted for
   decisions 0.1/0.3. Re-validate:
   `python submission/validate_cff.py CITATION.cff submission/cff-schema-1.2.0.json`
9. `envs/TOOL_VERSIONS.md` — the two corrections in §3.
10. `DATA_ACCESSIONS.md` — the three stale numbers in §3.
11. `README.md` — replace "MIT (see `CITATION.cff`; a `LICENSE` file is still to
    be added)" with the real statement; delete resolved blockers 3 and 4; keep
    1, 2 and 5 only until they are done, then delete the section.
12. `source_data/README.md` — rewrite against the final figure roster.
13. **Do not add a `.zenodo.json`.** Verified against Zenodo's documentation:
    *"If both files are present in your repository, only the `.zenodo.json`
    metadata will be used for the GitHub release archiving. The `CITATION.cff`
    metadata will be ignored by Zenodo."* Neither repo has one today. Adding one
    would silently override the CITATION.cff you just curated.

### Phase 3 — review, then push (maintainer)

14. `git log --oneline | head -20`, `git status`, `git diff --stat ca3200a..HEAD`.
15. Re-run the §4 audit on the final tree.
16. `git push --force-with-lease origin main` — the push **must** be a force,
    because the rebase rewrote the two packaging commits onto a new base. This
    is expected and safe for a repo with no other collaborators; confirm nobody
    else has a clone with local work first.

### Phase 4 — go public (PI / BMGLab org admin) — **irreversible in practice**

17. Settings → General → Danger Zone → Change visibility → Public. Verify:
    `gh api repos/BMGLab/Project_PeakATail --jq '.visibility, .license.spdx_id'`
    should return `PUBLIC` and the chosen SPDX id.
18. Confirm the tool repo is consistent — `LICENSE_DECISION.md` §7 item 2
    (the GPL-3.0 file still served on `BMGLab/PeakATail` `main`) should be fixed
    in the same session, or a referee will find two licences on the tool repo
    the same day the companion goes public.

### Phase 5 — Zenodo DOI (PI)

19. Sign in to Zenodo with the GitHub account that administers BMGLab.
20. Profile menu → **GitHub** → **Sync now**. Only public repositories appear;
    if `Project_PeakATail` is missing, the org has not granted Zenodo access —
    an organisation owner must approve the Zenodo OAuth app for BMGLab.
21. Toggle `BMGLab/Project_PeakATail` **on**. This installs the release webhook.
    Refresh to confirm it is listed as enabled.
22. **Only then** cut the release — the webhook does not archive releases that
    predate it. Tag per decision 0.4, with a release title and notes naming the
    manuscript and the two frozen tool commits.
23. Zenodo ingests the tarball and mints two DOIs: a **version DOI** for this
    release and a **concept DOI** that always resolves to the latest version.
    **Cite the version DOI in the manuscript** — it is the one that pins what a
    reader gets. Record the concept DOI in the README badge.
24. Check the deposition page: the author list, title and licence should have
    come from `CITATION.cff`. Fix any field on Zenodo before publishing the
    record; after publishing, metadata is editable but files are not.

### Phase 6 — back-fill (maintainer)

25. Replace `[[PLACEHOLDER: Zenodo DOI pending — repository must be public
    before a DOI can be minted]]` in `manuscript/00_draft/MAIN.md`
    (Declarations → Availability of data and materials) with the version DOI.
    STATUS.md item A3 closes here.
26. Uncomment and fill `doi:` in **both** `CITATION.cff` files; fill both
    `version:` placeholders with the tag. Re-validate both.
27. Add the Zenodo DOI badge to both READMEs.
28. Consider giving the **tool** repo the same treatment — Genome Biology asks
    that "the version used in the manuscript" be archived in a DOI-assigning
    repository, and that is the tool, not only the analysis code. A tag at
    `9dfdefb` archived to Zenodo would satisfy it directly.

---

## 6. Critical path, and what it means for the manuscript

```
0.1 licence ─┐
0.3 authors ─┼─> Phase 1 sync ─> Phase 2 metadata ─> Phase 3 push
0.4 tag     ─┘                                            │
                                                          v
                                            Phase 4 public ─> Phase 5 DOI ─> Phase 6 back-fill
```

Only the PI can start it (0.1, 0.3, 0.4) and only the PI can finish it (Phases 4
and 5 need organisation-admin and Zenodo rights). Phases 1–3 are perhaps a
half-day of mechanical work once the decisions exist.

This is **STATUS.md blocker A3**, one of the three items standing between the
draft and submission, and it is the only one of the three whose cost is
administrative rather than analytical: A1 (Kinnex de-duplication) and A2
(Table T7 / S5 / S6) need computation, A3 needs decisions. It is therefore the
one that can be finished first, and it should be — a referee sent to a private
repository at review time is a predictable, avoidable complaint, and reviewer
point M9vi already anticipates it.

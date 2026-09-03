# Repository packaging and the availability statement

Prepared 2026-09-02/03 for the PeakATail Genome Biology (Method track)
submission. Everything listed below is a **draft for the PI and the
maintainer**; nothing was applied to either repository and nothing in
`manuscript/` was touched.

> **Scope.** This index covers only the repository-packaging deliverables. Other
> files in `submission/` — `references.bib`, `CITATION_MAP.tsv`,
> `CITATION_GAPS.md`, `_build_citation_map.py`, `cover_letter.md`,
> `SUBMISSION_CHECKLIST.md` — come from sibling submission tasks and are not
> described here. `CITATION_GAPS.md` is worth reading alongside
> `availability_statement.md`: the manuscript's open `[[CITE]]` for the Kinnex
> data citation is resolved there in prose form.

> **Probably do not publish this directory.** It contains PI-facing decision
> briefs. If the working directory is synced to the public companion repo, add
> `submission/` to `.gitignore` first — see `COMPANION_REPO_PLAN.md` §5 step 1.

## Files from this task

| File | What it is |
|---|---|
| `LICENSE_DECISION.md` | The licence brief: candidates, what each implies, what the six benchmark comparators actually use (queried live), the GPL-dependency constraint, a recommendation, and the exact copyright lines. **No LICENSE file was placed** — that is the PI's legal decision. |
| `LICENSE.CANDIDATE-MIT.txt` | Canonical MIT text, verbatim from `gh api licenses/mit`. Substitute `[year] [fullname]`. |
| `LICENSE.CANDIDATE-BSD-3-Clause.txt` | Canonical BSD-3-Clause text. Substitute `[year], [fullname]`. |
| `LICENSE.CANDIDATE-GPL-3.0.txt` | Canonical GPL-3.0 text. Copied verbatim, no substitution. |
| `CITATION.cff` | Drop-in for the **tool** repo (BMGLab/PeakATail). Validated. |
| `CITATION.cff.companion` | Drop-in for the **companion** repo (BMGLab/Project_PeakATail). Validated. Carries an open authorship question for the PI. |
| `validate_cff.py` + `cff-schema-1.2.0.json` | Re-run the validation yourself: `python validate_cff.py CITATION.cff cff-schema-1.2.0.json` |
| `COMPANION_REPO_PLAN.md` | Current status of all five "Before release" blockers, the `source_data/`/`envs/` refresh delta, a pre-publication audit, and the ordered runbook to public + Zenodo DOI. |
| `availability_statement.md` | Drop-in prose for the manuscript's availability declaration, with a provenance table for every accession and URL. |

## The three findings that change the plan

1. **Two blockers are already resolved.** The pinned benchmark commit `18678ef`
   is publicly reachable on `origin/develop` (and is an ancestor of both frozen
   manuscript commits), and the three read-module fixes the repo warns are
   *missing* from the published numbers are in fact **present** in both v1 and
   v2 under re-landed SHAs. `COMPANION_REPO_PLAN.md` §2 has the proofs.

2. **The tool repo is not unlicensed — it is doubly licensed.**
   `BMGLab/PeakATail` publishes MIT on `develop` and GPL-3.0 on `main` at the
   same time, both public. This must be fixed whatever the PI decides.
   `LICENSE_DECISION.md` §0.

3. **The genuinely unlicensed repo is the companion one**, which is also the one
   that must go public before Zenodo can mint a DOI. It is on the critical path
   for STATUS.md blocker A3.

## Suggested order

1. PI reads `LICENSE_DECISION.md` → decides licence + copyright holder.
2. PI resolves the companion-repo authorship question
   (`COMPANION_REPO_PLAN.md` §2, blocker 1).
3. Maintainer runs `COMPANION_REPO_PLAN.md` §5 Phases 1–3.
4. PI runs Phases 4–5 (public, then Zenodo).
5. Maintainer back-fills the DOIs (Phase 6) and pastes
   `availability_statement.md` Part A into MAIN.md.

## Verification state

- Both `CITATION.cff` files: parse as YAML, pass the CFF 1.2.0 JSON schema with
  zero errors. The `doi:` key is commented out on purpose — CFF's DOI pattern
  rejects a placeholder string, proven by a negative control.
- Comparator licences, repo visibility and commit reachability: read-only
  `gh api` / `git ls-remote` / `git merge-base`, 2026-09-02.
- No git commit, no push, no `gh` write, no visibility change was made.

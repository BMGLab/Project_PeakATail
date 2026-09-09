# Licence decision brief — PeakATail (tool repo) and Project_PeakATail (companion repo)

Prepared 2026-09-02 for the PI. **No LICENSE file has been placed by this brief.**
Choosing a licence is a legal decision for the PI and, if Ege University asserts
institutional ownership of software written by its staff, for the university's
technology-transfer/legal office. This document supplies the facts, the
verified comparator evidence, a recommendation with its reasoning, and
copy-ready licence texts. It is not legal advice.

Candidate texts staged next to this file, fetched verbatim from the GitHub
Licenses API (`gh api licenses/<id>`) on 2026-09-02:

| file | md5 | bytes |
|---|---|---|
| `LICENSE.CANDIDATE-MIT.txt` | `fc17e78fb80a4fb0d94717082c01f844` | 1070 |
| `LICENSE.CANDIDATE-BSD-3-Clause.txt` | `bcd178f572e975c5dd719efbd8b00c9b` | 1500 |
| `LICENSE.CANDIDATE-GPL-3.0.txt` | `5b4473596678d62d9d83096273422c8c` | 35150 |

---

## 0. Correction to the brief: the tool repo is not unlicensed

The task that produced this document stated that the tool repo has no licence.
**That is no longer true, and the real problem is different and worse.** Verified
2026-09-02:

```
gh api repos/BMGLab/PeakATail  -> visibility PUBLIC, license.spdx_id = "MIT",
                                  defaultBranchRef = develop
git cat-file -p origin/develop:LICENSE | head -3   -> "MIT License / Copyright (c) 2026 BMGLab"
git cat-file -p origin/main:LICENSE    | head -3   -> "GNU GENERAL PUBLIC LICENSE / Version 3"
```

`LICENSE` and `CITATION.cff` are both **tracked** on `develop`. MIT arrived in
commit `d146590` (2026-05-11, "chore(license): standardise on MIT"), whose message
records that it was resolving an *earlier* inconsistency — `pyproject.toml` said
MIT while the on-disk LICENSE was GPL-3.0-only.

### The conflict you actually have

That commit fixed `develop` and never touched `main`. Both branches are public.

| branch | LICENSE text | position |
|---|---|---|
| `develop` (default branch, and the manuscript's code line) | MIT | 322 commits ahead of the merge base |
| `main` | GPL-3.0 | 10 commits ahead of the merge base, 322 behind |

So **github.com/BMGLab/PeakATail publishes two contradictory licences at the same
time.** GitHub's UI and API report only the default branch's (MIT), which makes
the GPL file on `main` easy to miss and easy for a reviewer or a downstream user
to find. A referee who checks `main` — the branch name most people try first —
sees GPL-3.0 on a paper whose availability statement will say something else.

This must be resolved regardless of which licence is chosen. Two contradictory
licence grants on one public repository is the one outcome that is worse than
having no licence at all: a recipient can plausibly claim they relied on
whichever branch they read.

Frozen manuscript commits `4efeb125` (v1) and `9dfdefb` (v2) are both on the
`develop` line, i.e. **under the MIT file**, not the GPL one. Anyone reproducing
the paper is working from MIT-licensed code today.

### The companion repo is the genuinely unlicensed one

```
gh api repos/BMGLab/Project_PeakATail -> visibility PRIVATE, license = null
```

Its `README.md` says "License: MIT (see `CITATION.cff`; a `LICENSE` file is still
to be added)" and its `CITATION.cff` declares `license: MIT`. That is a licence
*claim* with no licence *text* — the weakest position available: no operative
grant, and a README sentence that a court would read as an intention rather than
a licence. It is also the repo that must go public and be archived for a DOI, so
it is on the critical path. See `COMPANION_REPO_PLAN.md`.

---

## 1. What Genome Biology requires

From *Being open: our policy on source code* (Genome Biology editorial,
doi:10.1186/s13059-016-1040-y), verified 2026-09-02:

- source code must be "deposited in a recognized public source code repository";
- it must be "released under a license complying with an Open Source Definition,
  as defined by the Open Source Initiative" — **the journal accepts all
  OSI-compliant licences**, copyleft and permissive alike;
- for Method/Software articles, authors must "archive the version used in the
  manuscript in a DOI-assigning repository, such as Figshare or Zenodo";
- the licence and access route must be stated in the manuscript.

**All three candidates below are OSI-approved.** The journal does not constrain
the choice; it only requires that a choice be made, be OSI-compliant, and be
stated. The constraint is therefore entirely about what the lab wants.

---

## 2. The three realistic candidates

SPDX status verified against the current SPDX licence list (`spdx/license-list-data`,
list version `b885de7`) on 2026-09-02.

### MIT — SPDX `MIT` (current, OSI-approved)

- 168 words. Permits use, copying, modification, merging, publication,
  distribution, sublicensing and sale, with no conditions beyond preserving the
  copyright notice and the licence text.
- **Commercial redistribution: freely allowed, including in closed-source
  products.** A company may take PeakATail, modify it, embed it in a proprietary
  single-cell analysis platform and never publish the changes.
- No patent grant (explicit or, arguably, implied). If patenting any part of the
  method is contemplated, MIT is silent on patents — Apache-2.0 is the permissive
  licence that grants them explicitly. Raise this only if a patent is in play.
- No trademark/endorsement clause.

### BSD-3-Clause — SPDX `BSD-3-Clause` (current, OSI-approved)

- 223 words, three numbered conditions. Clauses 1 and 2 say the notice must be
  retained in source redistributions and reproduced in the documentation of
  binary redistributions — MIT achieves the same in one sentence. **Clause 3 is
  the real difference: it forbids using the name of the copyright holder or its
  contributors to endorse or promote derived products without prior written
  permission.**
- Commercial redistribution: identical to MIT — freely allowed, closed source
  permitted.
- The only substantive reason to prefer it over MIT is clause 3: it stops a
  company from marketing "powered by PeakATail — Ege University" without asking.
  If institutional-name protection matters to the university, this is the cheap
  way to get it; otherwise it is MIT with extra words.

### GPL-3.0 — SPDX `GPL-3.0-only` or `GPL-3.0-or-later`

- **Note: the bare SPDX id `GPL-3.0` is DEPRECATED.** Use `GPL-3.0-only` (this
  version only) or `GPL-3.0-or-later` (this or any later version). Writing
  `GPL-3.0` in `pyproject.toml` or `CITATION.cff` is a defect in new metadata.
  `-or-later` is the FSF's own recommendation and the more common choice.
- Strong copyleft: anyone who distributes the software or a **derivative work**
  must do so under GPL-3.0 and must provide corresponding source. Includes an
  explicit patent grant and an anti-tivoisation clause.
- **Commercial redistribution: allowed, but only under GPL.** Selling it is fine;
  embedding it in a closed-source product is not. This is the practical
  difference: it forecloses adoption by commercial pipeline vendors and by
  companies whose legal policy bans copyleft in shipped products, and it
  complicates inclusion in permissively-licensed aggregations.
- Using the tool to produce results, or calling its CLI from a pipeline, is
  **not** distribution of a derivative work — GPL does not restrict academic use
  or publication of results by anyone, under any of these options.

---

## 3. What comparable tools in this benchmark actually use

Queried live via `gh api repos/<owner>/<repo>` on 2026-09-02. Repository
identities taken from `scripts/benchmark_tools/status/*.md` and
`envs/TOOL_VERSIONS.md` in the companion repo — not guessed.

| Tool | Repository | Licence (as GitHub reports it) | Class |
|---|---|---|---|
| polyApipe | MonashBioinformaticsPlatform/polyApipe | **LGPL-2.1** | weak copyleft |
| Sierra | VCCRI/Sierra | **GPL-3.0** | strong copyleft |
| scAPAtrap | BMILAB/scAPAtrap | **AGPL-3.0** | strongest copyleft (network clause) |
| SCAPTURE | YangLab/SCAPTURE | **"Other"** — see below | dual / non-free |
| scTail | StatBiomed/scTail | **Apache-2.0** | permissive (+ patent grant) |
| scUTRquant | Mayrlab/scUTRquant | **GPL-3.0** | strong copyleft |

**SCAPTURE is not OSI-compliant as written.** Its LICENSE file (fetched verbatim)
reads:

> Copyright ©2021 Shanghai Institute of Nutrition and Health. All Rights Reserved.
> Licensed GPLv3 for open source use or contact YangLab (yanglab@picb.ac.cn) for commercial use.
> Permission to use, copy, modify, and distribute this software and its documentation for
> educational, research, and not-for-profit purposes, without fee and without a signed licensing
> agreement, is hereby granted […]

That is a dual-licence with a not-for-profit restriction in the permission grant,
which is why GitHub classifies it `NOASSERTION`/"Other" rather than GPL-3.0. It
would **not** satisfy Genome Biology's OSI-compliance requirement as the primary
grant. Useful as a cautionary example: do not write a bespoke licence.

**Reading of the panel.** Five of six comparators are copyleft; only scTail
(Apache-2.0) is permissive, and it is the smallest and newest project in the set
(5 stars) and the one tool that could not even be run here. The field's centre of
gravity is GPL-3.0. Choosing MIT makes PeakATail the *most* permissive tool in
its own comparison table — an outlier, but an outlier in the direction that costs
the project nothing and maximises who can adopt it. Choosing GPL-3.0-or-later
would put it squarely in the majority of its own benchmark.

For context, the resource the benchmark scores against — PolyASite 2.0 — is a
data atlas, not code; its terms govern the atlas, not this tool's licence.

---

## 4. The one technical constraint: GPL runtime dependencies

Verified from the installed environment
(`tools/PeakATail/.venv/lib/python3.*/site-packages/*.dist-info/METADATA`), not
from memory:

| dependency (declared in `pyproject.toml`) | installed version | licence |
|---|---|---|
| `leidenalg` | 0.10.2 | **GPLv3+** |
| `louvain` | 0.8.2 | **GPLv3+** |
| `igraph` | 0.11.9 | **GPL** |
| `pysam` | 0.23.0 | MIT |
| `pybedtools` | 0.12.0 | MIT |
| `scanpy` | 1.11.1 | BSD-3-Clause |
| `anndata` | 0.11.4 | BSD-3-Clause |
| `statsmodels` | — | BSD |

`pyproject.toml` declares `leidenalg`, `louvain` and `igraph` as **hard runtime
dependencies**. They are not imported directly anywhere in `ema/`; they are
reached indirectly through `sc.tl.leiden(..., flavor='igraph')` in
`ema/clustering/strategies/leiden_libsize.py:86` and
`ema/clustering/strategies/leiden_tfidf.py:229`.

Why this matters, stated without overclaiming: the FSF's position is that a
Python program which imports a GPL library forms a single combined work, so
distributing an MIT-licensed package that *requires* GPL libraries is a
recognised grey area. It is also extremely common — scanpy itself is BSD-3-Clause
and treats `leidenalg` as a dependency — and to our knowledge has never been
litigated for pip dependencies. Two clean ways out:

1. **Move `leidenalg`, `louvain` and `igraph` from `dependencies` to an optional
   extra** (e.g. `[project.optional-dependencies] clustering = [...]`), the way
   scanpy does, and let the clustering strategies raise a clear
   "install `peakatail[clustering]`" error if they are missing. The base install
   is then MIT-with-permissive-deps and the question disappears.
2. Licence the tool GPL-3.0-or-later, which makes the question moot by construction.

**Option 1 is unusually attractive here** because the manuscript *dropped the
clustering claim* — clustering is no longer a load-bearing part of the paper, so
demoting its dependencies costs the paper nothing. This is a code change, not a
licence change, and it is the PI's/maintainer's call; it is flagged here only
because it is the single fact that could change the licence answer.

---

## 5. Recommendation

**Adopt MIT for both repositories, and pair it with the dependency change in
§4 option 1.**

Reasoning, strongest first:

1. **It ratifies reality instead of creating a third record.** The tool repo's
   default branch, its `pyproject.toml` (`license = "MIT"`, PEP 639 syntax), its
   README badge, its `CITATION.cff`, GitHub's own licence detection, and the
   companion repo's README and `CITATION.cff` all already say MIT. Confirming MIT
   costs one commit on `main`. Switching to GPL means editing six declarations
   across two repositories and re-checking every downstream mention — more work,
   more chance of a fresh inconsistency, and a public repo that said MIT for four
   months before changing its mind.
2. **It matches the paper's own thesis.** This is a Method paper whose value is
   that other people measure with the tool and benchmark against it. Precision-first
   PAS calling wants to be embeddable — in nf-core and Galaxy wrappers, in Bioconda,
   and in commercial single-cell analysis stacks. MIT removes every barrier to that;
   GPL keeps commercial vendors out, which is the population most able to give the
   method reach.
3. **It aligns with the stack it sits on.** scanpy and anndata are BSD-3-Clause;
   pysam and pybedtools are MIT. The permissive scientific-Python convention is
   the one PeakATail already lives inside.
4. **It is OSI-approved**, so Genome Biology is satisfied either way — the journal
   is not a tiebreaker.

**Choose differently if** either of these is true, and both are PI calls this
brief cannot make:

- *The lab wants to prevent a company from taking the code closed.* Then
  **GPL-3.0-or-later**. It is the licence four of the six comparators use, it
  resolves §4 by construction, and it costs adoption by commercial pipelines.
  Do not choose bare `GPL-3.0` — it is a deprecated SPDX identifier.
- *Ege University wants explicit protection of its name against a redistributor
  implying endorsement.* Then **BSD-3-Clause** — MIT's terms plus that one clause.

**AGPL-3.0 is not recommended** even though scAPAtrap uses it: its network clause
targets hosted services, PeakATail is a CLI, and AGPL is the licence most often
blocked outright by institutional and corporate policy.

**Both repositories should carry the same licence.** A permissive tool with a
copyleft analysis companion (or the reverse) invites exactly the confusion §0
describes.

---

## 6. Exact copyright lines

`YYYY` is **2026** (first publication year; the existing MIT file on `develop`
already says 2026). The **copyright holder is a legal decision**: for
university-employed authors in Türkiye the holder may be the institution, the
authors personally, or both. The existing file names "BMGLab", which is a lab
name and not a legal entity — that is worth correcting whichever licence is
chosen. Realistic options, in decreasing institutional formality:

- `Ege University`
- `Ege University and the PeakATail authors`
- `Amir Amiri Tabat and Yasin Kaymaz`

### MIT — line 3 of `LICENSE.CANDIDATE-MIT.txt`

The staged file carries the canonical placeholders `[year]` and `[fullname]`.
Substitute in place:

```
Copyright (c) 2026 <holder>
```

e.g. `Copyright (c) 2026 Ege University`. One command:

```bash
sed -i 's/\[year\] \[fullname\]/2026 Ege University/' LICENSE
```

### BSD-3-Clause — line 3 of `LICENSE.CANDIDATE-BSD-3-Clause.txt`

Note the comma, which MIT does not have:

```
Copyright (c) 2026, <holder>
```

```bash
sed -i 's/\[year\], \[fullname\]/2026, Ege University/' LICENSE
```

### GPL-3.0 — no substitution in the LICENSE file

The `<year>` / `<name of author>` tokens in `LICENSE.CANDIDATE-GPL-3.0.txt`
appear only in the appendix "How to Apply These Terms to Your New Programs"
(lines 635 and 655), which is **template text that stays verbatim**. The GPL
LICENSE file is copied byte-for-byte. The copyright notice goes instead into a
header on each source file and into the README:

```
PeakATail — single-cell poly(A) site detection and APA analysis
Copyright (C) 2026 Ege University

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.
```

---

## 7. Once the choice is made — every place it must be stated

Do all of these in one commit per repository, so no branch or file ever
disagrees with another.

**Tool repo (github.com/BMGLab/PeakATail):**

1. `LICENSE` on `develop` — replace only if the choice is not MIT.
2. **`LICENSE` on `main` — mandatory in every scenario.** Reconcile `main` with
   `develop` so the GPL-3.0 file stops being publicly served. Either merge/replace
   the file, or, better, retire `main` as a stale branch. This is the §0 defect
   and it is independent of the licence choice.
3. `pyproject.toml` — `license = "<SPDX id>"` (PEP 639 expression syntax; the file
   already carries a comment warning not to add a legacy `License ::` classifier
   alongside it, which remains correct).
4. `CITATION.cff` — `license:` field. Use `submission/CITATION.cff`, which is
   validated and already set to `MIT`; change the field if the choice differs.
5. `README.md` — licence badge and the "License" section (both currently say MIT).
6. `docs/` — grep for licence mentions before publishing.

**Companion repo (github.com/BMGLab/Project_PeakATail):**

7. Add the `LICENSE` file — it currently has none at all.
8. `CITATION.cff` — `license:` field (currently `MIT`, alongside placeholder authors).
9. `README.md` — replace "a `LICENSE` file is still to be added" with the real statement.

**Manuscript:**

10. The availability statement must name the licence explicitly — Genome Biology
    requires it stated in the manuscript, and the current draft does not name one.
    `submission/availability_statement.md` carries a marked slot for it.

**Verification after the change** — this must come back with one answer, not two:

```bash
cd tools/PeakATail
for b in origin/main origin/develop; do
  printf '%-18s ' "$b"; git cat-file -p "$b:LICENSE" 2>/dev/null | head -1
done
gh api repos/BMGLab/PeakATail        --jq '.license.spdx_id'
gh api repos/BMGLab/Project_PeakATail --jq '.license.spdx_id'
```

---

## 8. What this brief deliberately did not do

- **No `LICENSE` file was placed** in either repository, and no existing licence
  file was edited. The candidates are staged beside this document only.
- **No commit, push, or `gh` write** of any kind. All GitHub calls were read-only
  (`gh api repos/...`, `gh api licenses/...`, `git ls-remote`).
- **No copyright holder was chosen** and no licence text was paraphrased — the
  three staged files are byte-verbatim from the GitHub Licenses API, with md5s
  recorded above so they can be re-verified.

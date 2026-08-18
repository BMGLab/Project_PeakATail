Documentation contradictions that must be resolved before the manuscript Methods section is written (each one is a claim a reviewer can check against the docs site):

1. **UMI handling** — `docs/concepts/data-flow.md` states "UMI tags are not used; PeakATail counts raw read 3′ ends", while the README and docs index require `UB:Z` tags in the input BAM. Which is true? If UMIs are genuinely unused, the Methods must not imply UMI deduplication; if they are used, data-flow.md is wrong.
2. **`--ip-filter` rows** — `docs/cli/run.md` contains both the D9 annotate-mode table and stale "currently no-op" rows for the same flags (see the duplicate-arm issue for the experimental side of this).
3. **README command table** — 8 of 12 commands listed; `reannotate`, `switch trend`, `switch combine`, `collapse` are documented in `docs/cli/` but absent from the README.
4. **Citation placeholder** — README citation block lists a placeholder author; align with the intended author list and add `CITATION.cff` (tracked in the release-engineering issue).

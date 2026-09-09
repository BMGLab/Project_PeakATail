# .verify/ — audit trail for the adversarial verification pass (2026-09-03)

Evidence and tooling for the independent re-check of `submission/`. Nothing here is
part of the submission; it exists so any claim in this package can be re-tested.

| file | what it is |
|---|---|
| `check_bib.py` | parses `references.bib` and resolves every DOI against `api.crossref.org`; prints a per-entry title / first-author / author-count / year comparison |
| `crossref_records.json` | the raw Crossref record fetched for each of the 33 DOI-bearing entries |
| `check_epmc.py` | second, independent source: resolves every DOI against the Europe PMC REST API and compares journal / volume / pages / year / PMID |
| `epmc_records.json` | the raw Europe PMC records |
| `fix_authors.py` | the script that expanded 10 abbreviated author lists to the publisher-deposited full given names; only rewrites when the surname sequence matches exactly and Crossref strictly adds information |
| `canon_mit.txt`, `canon_bsd-3-clause.txt`, `canon_gpl-3.0.txt` | canonical licence texts from the GitHub licenses API; the staged `LICENSE.CANDIDATE-*.txt` files are byte-identical to these |
| `gnu_gpl3.txt` | the GPL-3.0 text from gnu.org, used as a second source (identical to the above after whitespace normalisation) |
| `cff-schema-fresh.json` | CFF 1.2.0 schema fetched fresh; byte-identical to the staged `cff-schema-1.2.0.json` |
| `references.bib.bak`, `CITATION_MAP.tsv.bak` | pre-edit copies, for diffing |
| `all_dois.txt`, `bib_dois.txt` | every DOI asserted anywhere in `submission/`, and the subset in the .bib |

Re-run the two reference checks with:

    export LC_ALL=C
    cd /mnt/ssd1/Projects/PeakATail_wd/submission
    python3 .verify/check_bib.py     # Crossref
    python3 .verify/check_epmc.py    # Europe PMC

Expected: every entry `OK` except `herrmann2020polyasite` (flagged `CHECK:YEAR` — Crossref's
`issued` date is the 2019 online-first date; the NAR 48(D1) database issue is 2020, which Europe
PMC confirms) and the two GEO `@misc` entries, which carry no DOI by design and are verified
against NCBI E-utilities `esummary` instead.

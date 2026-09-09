<!--
SUBMISSION_CHECKLIST.md — Genome Biology (Method track) requirements mapped onto the actual state
of this project. Compiled 2026-09-02 from: the live Genome Biology submission guidelines (see §0 for
exactly what was read and how), manuscript/00_draft/{MAIN.md,SUPPLEMENTARY.md,STATUS.md,
PI_DECISIONS.md,NUMBERS_LEDGER.tsv}, manuscript/figures/FIGURES_MANIFEST.md, manuscript/github/
PUBLISHED.md, submission/LICENSE_DECISION.md, submission/CITATION_GAPS.md and
submission/COMPANION_REPO_PLAN.md. Every count in a "current state" paragraph was measured from the
files on disk (method in §3), not copied from a summary; where a measurement disagrees with a
summary, both are shown.

This file does not modify the manuscript. Items needing a text edit are listed as actions for
whoever holds MAIN.md / SUPPLEMENTARY.md (locked during the legend-sync pass).
-->

# Genome Biology, Method track — submission checklist for PeakATail

**Status vocabulary:** `DONE` · `PI ACTION` (a decision or text only the PI/authors can supply) ·
`BLOCKED-BY` (needs work or an external event first, named).

**Verification vocabulary:** each requirement carries how it was established.
`[LIVE]` read from the journal's current guidelines pages on 2026-09-02.
`[PRACTICE]` corroborated from Genome Biology papers as published.
`[SECOND-HAND]` verified by another workstream in this repo, not re-read here.
`[INFERRED]` general publishing practice — **not** confirmed against a journal page; treat as a
prompt to check, not as a requirement.

---

## 0. How the journal's requirements were verified (read this before trusting §1–§18)

Genome Biology's guidelines have moved from `genomebiology.biomedcentral.com` to
`link.springer.com/journal/13059`. Every `biomedcentral.com` guideline URL now 301-redirects there,
and those Springer pages answer a direct fetch with a cookie-authorization redirect
(`idp.springer.com/authorize`) instead of content. `web.archive.org` and `pmc.ncbi.nlm.nih.gov` are
both unreachable from this environment (domain block; CAPTCHA).

The pages below were therefore read through a text-extraction proxy (`r.jina.ai`) of the live
Springer URLs on **2026-09-02**. The content is the live page's, but it passed through one
intermediary — so for anything a submission would fail on, **re-read the page in a browser before
acting**:

- `link.springer.com/journal/13059/submission-guidelines` — article types, cover letter, figures,
  tables, additional files, peer-review model.
- `link.springer.com/journal/13059/submission-guidelines/methodology` — the Method-track page.
- `link.springer.com/journal/13059/submission-guidelines/research` — the Research page (quoted here
  only as the alternative track, for the abstract question in §2).
- `link.springer.com/brands/bmc/editorial-policies` — data/materials availability, preprints,
  competing interests, authorship, duplicate submission.

Corroboration from published practice: **scTail** (Genome Biology, 2025; DOI 10.1186/s13059-025-03710-7)
carries the JATS subject label **"Methodology"** and a single-paragraph, unstructured abstract of
~110 words; **SCAPTURE** (Genome Biology, 2021; DOI 10.1186/s13059-021-02437-5) likewise has an
unstructured abstract. Both were read via the Europe PMC REST API.

Not verified anywhere, flagged `[INFERRED]` below: manuscript file formats accepted by the
submission system, ORCID requirements at submission, reporting-standard checklists, the AI/LLM
disclosure policy text, and the suggested/opposed-reviewer fields.

---

## 1. Article type and section structure

**Requirement `[LIVE]`.** The Method track's guidelines page is titled **"Methodology"**:
*"Genome Biology publishes outstanding new methods that will be of utility to a wide audience"*,
and such articles *"should describe novel methods that are shown to be a clear advance over
existing state-of-the-art methods in a side-by-side demonstration, where possible"*, with
benchmarking against a dataset where ground truth is known recommended. Section order:
**Title page → Abstract → Keywords → Background → Results → Discussion → Conclusions → Methods →
Abbreviations → Declarations → References → Figures, tables and additional files.**
`[PRACTICE]` published papers of this type carry the label "Method".

**Current state.** `MAIN.md` runs Title/authors → Abstract → Keywords → **Key points** → Background →
Results → Discussion → Conclusions → Methods → Declarations → Figure legends → Supplementary
information. The side-by-side demonstration the track asks for is the six-tool benchmark on shared
BAMs, one scorer, shuffled nulls and matched call counts.

| # | Item | Status | What is needed |
|---|---|---|---|
| 1.1 | Method/Methodology is the right track | `DONE` | Fits the type's own wording; the letter argues it. |
| 1.2 | Section order as listed | `PI ACTION` | Only two deviations: no **Abbreviations** section (§1.3), and a **Key points** block between Keywords and Background — GB does not list it. Keep it only if the handling editor accepts it; otherwise fold the load-bearing bullets into Background/Discussion. Text edit, owner of MAIN.md. |
| 1.3 | Abbreviations section | `PI ACTION` | Absent. The draft uses PAS, APA, UMI, CB, FDR, PDUI, AMI, ARI, IP, R_det, P@100. Add a short list, or confirm with the editor that it may be omitted. |
| 1.4 | Declarations placed after Methods | `DONE` | Present and correctly positioned. |
| 1.5 | References section | `BLOCKED-BY` §7.2 | The file has no reference list yet; citations are still 47 `[[CITE: …]]` markers, though the bibliography behind them is built and verified (`submission/references.bib`). |

---

## 2. Abstract — the one hard limit this draft currently breaks

**Requirement `[LIVE]`, Methodology page, quoted:** *"The Abstract should not exceed 100 words.
Please minimize the use of abbreviations and do not cite references in the abstract. The abstract
should be unstructured."*
**Requirement `[LIVE]`, Research page, quoted:** *"The Abstract should not exceed 250 words"*, with
three sections **Background / Results / Conclusions**.
`[PRACTICE]` both Genome Biology method papers checked (scTail 2025, SCAPTURE 2021) have short
unstructured abstracts — consistent with the 100-word rule, not with a 250-word structured one.

**INDEPENDENTLY RE-VERIFIED 2026-09-03 by two routes that do not use the r.jina.ai proxy.**
(a) A web search of the official Methodology guidelines page returns the rule verbatim: *"The
Abstract should not exceed 100 words. Please minimize the use of abbreviations and do not cite
references in the abstract. The abstract should be unstructured."* (b) First-party journal output,
counted here from the Europe PMC records: **scTail (Genome Biology 2025) = 97 words, unstructured;
SCAPTURE (Genome Biology 2021) = 85 words, unstructured.** Both sit under 100 and neither carries
Background/Results/Conclusions headings. A direct fetch of link.springer.com still fails (JS
"Client Challenge"), and genomebiology.biomedcentral.com now 301-redirects to Springer — so the
proxy caveat below was honest, but the rule itself no longer rests on the proxy.

**Current state (RE-MEASURED 2026-09-03; supersedes the "272 / 285" figures previously printed
here, which do not reproduce).** The draft's abstract is **structured** (Background / Results /
Conclusions). Counting whitespace tokens after stripping bold markers and HTML comments:
**276 words** for the abstract proper — which reproduces `STATUS.md`'s 276 exactly — or **281** if
em-dash-joined tokens are split, and **291** if the trailing `**Keywords:**` line is counted with
it. Under every convention it is **≈2.8× the Methodology limit and ~10% over the Research limit.**
The draft's own note assumed "limit ≤ 250", i.e. the Research rule.

| # | Item | Status | What is needed |
|---|---|---|---|
| 2.1 | Abstract format for the chosen track | `PI ACTION` — **decide first, it changes the rewrite** | If submitting as **Method/Methodology**: rewrite to ≤100 words, unstructured. If the editor agrees the paper is better handled as **Research**: keep the structure and cut ~26+ words. Everything else in §2 depends on this call. |
| 2.2 | Abstract length | `BLOCKED-BY` 2.1 | 276 → 100 (Methodology) or → ≤250 (Research). At 100 words only ~3 claims survive: the evidence model and the pre-registered default's standing, the calibrated switch test, and cross-patient replication with a null. The four failed/negative results cannot all fit; keep the failed pre-registered gate, which is the one referees called the paper's distinguishing feature. |
| 2.3 | No references cited in abstract | `DONE` | None cited. |
| 2.4 | Minimize abbreviations | `PI ACTION` | Current abstract uses APA, PAS, R_det, PBMC. A 100-word rewrite should carry at most APA and PAS. |
| 2.5 | Keywords | `DONE` | Eight keywords present (three to ten is the usual BMC range `[INFERRED]`). |
| 2.6 | Key points block | `PI ACTION` | Not a listed section for this track (§1.2). If it stays, it is not covered by any verified limit; if it goes, the negatives it carries must be visible elsewhere. |

---

## 3. Manuscript length

**Requirement `[LIVE]`.** Neither the Methodology page nor the Research page states a main-text word
limit, a figure limit, a table limit or a reference limit. The **only** numeric text limit found is
the abstract's.

**Current state (measured here, `LC_ALL=C`, HTML comments stripped).** `MAIN.md` = **12,970 words**
(13,361 with comments). By section: title block 52 · abstract + keywords 293 · key points 350 ·
Background 924 · Results 4,415 · Discussion 690 · Conclusions 134 · Methods 3,236 · Declarations 347 ·
figure legends 1,963 · supplementary-legend section 566. Main-text prose only
(Background → Methods) = **9,399 words**. `SUPPLEMENTARY.md` = 3,052 words.

*Method:* whitespace tokens after stripping `<!-- … -->` comments, split at `#`-headings —
`python3 -c "import re;t=open('MAIN.md').read();print(len(re.sub(r'<!--.*?-->','',t,flags=re.S).split()))"`
under `LC_ALL=C`. The abstract was counted separately with bold markers removed (§2).

| # | Item | Status | What is needed |
|---|---|---|---|
| 3.1 | Length within journal limits | `DONE` (no limit found) | Reproduces `STATUS.md`'s table exactly. **Note:** the "11,923 words" figure quoted in the task brief does not reproduce under any scope measured here — the auditable number is 12,970 total / 9,399 main-text prose. Use those. |
| 3.2 | Figure legends length | `PI ACTION` | 1,963 words for six legends is long by any standard. Every caveat in them is load-bearing (they carry the "travelling caveats"), so this is a judgement call, not a defect — but expect an editor to ask. |

---

## 4. Figures

**Requirement `[LIVE]`.** Formats: *"EPS, PDF, Microsoft Word, PowerPoint, TIFF, JPEG, PNG, BMP,
CDX"*. Resolution: *"approximately 300 dpi at the final size"*. Size: *"Individual figure files
should not exceed 10 MB"*.

**Current state (RE-MEASURED 2026-09-03; the "~21 MB" figure previously printed here does NOT
reproduce and has been corrected).** 6 main figures and 10 supplementary figures are built, each as
a vector PDF (fonttype-42, no Type 3, `pdffonts`-checked) plus a 600 dpi PNG. Exact byte counts:

| set | files | PNG total | note |
|---|---|---|---|
| 6 submitted main figures | fig1–fig6 | **4.73 MB** | excludes `fig1_overview_detailed.png` |
| `fig1_overview_detailed.png` | 1 | 1.66 MB | auxiliary variant — confirm whether it is submitted |
| 10 supplementary figures | S1–S4, S7–S12 | **12.38 MB** | largest is `figS12_versions.png` at **1.98 MB** |
| all 16 main+supp PNGs | 16 | **18.77 MB** | still under 20 MB, but with little headroom |
| all 10 supplementary **PDFs** | 10 | **0.41 MB** | the format to bundle from |

Largest single figure file of any kind is 1.98 MB, well inside the 10 MB per-figure cap. Two
supplementary slots, **S5** and **S6**, are declared "in preparation" and nothing in the draft
cites them.

| # | Item | Status | What is needed |
|---|---|---|---|
| 4.1 | Accepted format | `DONE` | PDF (vector) is on the accepted list; PNG is the fallback. |
| 4.2 | Resolution ≥ ~300 dpi | `DONE` | 600 dpi PNGs; PDFs are vector. |
| 4.3 | ≤10 MB per figure | `DONE` | Largest single file 1.98 MB. |
| 4.4 | One file per figure, numbered in citation order | `DONE` | Frozen numbering in `FIGURES_MANIFEST.md`; Figs 1–6 and S1–S12 with S5/S6 as the gaps. |
| 4.5 | S5 and S6 | `BLOCKED-BY` A1/A2 (`PI_DECISIONS.md`) | S5 needs the de-duplicated long-read truth; S6 needs the shuffled-position null for the re-anchored window. **Decide by freeze:** build both, or drop them and renumber S7→S5 … S12→S10 in one pass. Do not submit with two "in preparation" slots — the simulated referee (M10) treats that as structurally incomplete. |
| 4.6 | Figure legends in the manuscript, not on the images | `DONE` | The 2026-09-02 surgery pass moved legends off the images into script-written sidecars; MAIN.md carries them. |
| 4.7 | Colour accessibility | `DONE` | Okabe–Ito palette throughout (house rule, verified in the manifest). |

---

## 5. Tables

**Requirement `[LIVE]`.** *"Tables should be numbered and cited in the text in sequence using Arabic
numerals."* Large tables belong in additional files `[INFERRED]`.

| # | Item | Status | What is needed |
|---|---|---|---|
| 5.1 | Main-text tables | `DONE` | There are none; nothing to number. |
| 5.2 | Supplementary Table T7 (competitor versions, parameters, run commands, run-log caveats) | `BLOCKED-BY` collation (A2) | Referenced from Methods and does not exist. Content is in the verified run registry and run logs; version numbers must be read from the logs, never reconstructed. Every competitor claim in Figs 2/3 is vouched for by this table. |
| 5.3 | Table naming convention | `PI ACTION` | "T7" is a project-internal label. Renumber supplementary tables to the journal's convention (Additional file *n* / Table S1…) at assembly. |

---

## 6. Additional files (supplementary information)

**Requirement `[LIVE]`.** *"Maximum file size for additional files is 20 MB each."*
Naming/citation convention (Additional file 1, cited in text, listed at the end of the manuscript):
`[INFERRED]` — the detail page could not be opened; check it in a browser.

| # | Item | Status | What is needed |
|---|---|---|---|
| 6.1 | Additional file 1 exists and is cited | `DONE` | `SUPPLEMENTARY.md` is titled "Additional file 1: Supplementary figure legends" and is cited from MAIN.md's Supplementary information section. |
| 6.2 | ≤20 MB per additional file | `PI ACTION` at assembly | **Corrected 2026-09-03: this was a false alarm.** The 10 supplementary PNGs total **12.38 MB**, not ~21 MB, so they are already inside the 20 MB cap; all 16 main+supp PNGs come to 18.77 MB, also inside it but with little headroom. The advice is still worth taking on quality grounds rather than size: bundle the supplement from the **PDFs** (0.41 MB for all ten, vector) rather than from 600 dpi rasters. |
| 6.3 | Supplementary figures S5/S6 | `BLOCKED-BY` — see 4.5 | Same decision. |
| 6.4 | Supplementary table T7 | `BLOCKED-BY` — see 5.2 | |
| 6.5 | **Supplementary figures ship as separate files, not embedded** | `DONE` (2026-09-03) | `submission/4_BUILD/output/PeakATail_supplementary.docx` was built from `SUPPLEMENTARY.md` and carries the **legends only — zero embedded images**, verified by re-opening the file (`submission/4_BUILD/verify_docx.py`). The ten supplementary figures (S1–S4, S7–S12) upload as their own files beside it, exactly as Figs 1–6 do for the main manuscript; bundle them from the vector PDFs per 6.2. Same rule for the main manuscript file: `PeakATail_manuscript_submission.docx` contains **0** embedded images. Build and regeneration commands: `submission/BUILD_NOTES.md`. |

---

## 7. References and citations

**Requirement `[INFERRED]`** (numbered reference list in the journal's style; the style page was not
reachable). `[LIVE]` on what may be cited: *"Only articles, clinical trial registration records and
abstracts that have been published or are in press, or are available through public e-print/preprint
servers, may be cited."*

**Current state (measured here, and cross-read against the bibliography workstream's outputs in this
directory).** MAIN.md still carries **47 `[[CITE: …]]` markers** over 26 phrasings (23 distinct
targets after merging three double-phrased pairs) and has **no reference list in the file yet**.
The bibliography itself, however, is largely built: `submission/references.bib` holds 35 entries,
each verified against a live record (Europe PMC / Crossref / NCBI E-utilities) with a DOI and an
access date; `submission/CITATION_MAP.tsv` maps every token to its entry; `submission/CITATION_GAPS.md`
reports **21 of 23 targets and 45 of 47 token instances resolved**.

| # | Item | Status | What is needed |
|---|---|---|---|
| 7.1 | Bibliography resolved | `DONE` for 21/23 targets | See `CITATION_GAPS.md`. Nothing was resolved by resemblance; two targets are deliberately left open rather than guessed. |
| 7.2 | Reference list inserted into the manuscript | `BLOCKED-BY` the MAIN.md lock | The `[[CITE]]` markers must be replaced by numbered citations and a reference list generated from `references.bib`. Text edit for the MAIN.md owner after the legend-sync pass. |
| 7.3 | G1 — spermatogenesis 3′UTR-shortening gene panel | `PI ACTION` — **blocking** | The 16-gene panel is a hard-coded list in `scripts/manuscript_figures/spermatogenesis_control_steps/step7_genes.py` with no recorded source. A literature-derived panel with no literature cannot be cited, and the paper reports its non-reproduction as a negative result. Either the PI supplies the sources, or the claim is rewritten so it does not rest on an uncitable panel. |
| 7.4 | G3 — Kinnex x3p / GEM-X dataset citation and accession | `BLOCKED-BY` A1 | Must come from the truth-build record. The referee simulation (M8) calls its absence a reproducibility defect: *"publicly available PacBio Kinnex long-read PBMC data" is not a citation a reader can follow.* |
| 7.5 | Never fabricate a reference | `DONE` (policy, and observed) | Every `[[CITE]]` target is named and unambiguous; the two unresolved ones stay unresolved. |

---

## 8. Data availability

**Requirement `[LIVE]`, BMC editorial policies.** Every manuscript needs an **"Availability of data
and materials"** section stating where the data supporting the conclusions are, with
community-recognised repositories used where they exist.

**Current state.** The section names: 10x Genomics public `pbmc_10k_v3` and `pbmc4k`; GEO
**GSE104556** (mouse testis); GEO **GSE123904** (lung adenocarcinoma); PolyASite 2.0 (GRCh38/GRCm38
representative sites); and the PacBio Kinnex truth/decoy BEDs distributed with the companion
repository.

| # | Item | Status | What is needed |
|---|---|---|---|
| 8.1 | Availability section present | `DONE` | With accessions for both GEO datasets. |
| 8.2 | Kinnex source identified | `BLOCKED-BY` A1 | See 7.4 — the accession and download route are still missing. |
| 8.3 | Third-party data are available | `DONE` | All datasets are previously published or publicly distributed. |
| 8.4 | Derived data (truth/decoy BEDs, audit TSVs) | `BLOCKED-BY` A3 | They live in the companion repository, which is still private. |

---

## 9. Code availability, licence and archived version — the Method track's own gate

**Requirement `[SECOND-HAND]`** (verified by `submission/LICENSE_DECISION.md` against the Genome
Biology editorial *Being open: our policy on source code*, doi:10.1186/s13059-016-1040-y — not
re-read here): source code must be deposited in a recognised public repository; released under a
licence complying with the **Open Source Initiative**'s Open Source Definition (any OSI licence,
copyleft or permissive); for Method/Software articles **the version used in the manuscript must be
archived in a DOI-assigning repository** such as Zenodo or Figshare; and the licence and access
route must be **stated in the manuscript**.
`[LIVE]` the Methodology page repeats: *"Source code for tools described in Method or Software
articles should be deposited in a public repository, with an OSI-compliant license."*

| # | Item | Status | What is needed |
|---|---|---|---|
| 9.1 | Tool code public | `DONE` | github.com/BMGLab/PeakATail, public. |
| 9.2 | Manuscript versions pinned | `DONE` | `4efeb125` (v1) and `9dfdefb3eb353b0817ef79c4eb9ace6d6c8aab53` (v2, the record) named in the Availability section. |
| 9.3 | OSI-compliant licence, unambiguous | `PI ACTION` — **do not submit before this is fixed** | The repository currently publishes **two contradictory licences at once**: MIT on `develop` (the default branch and the manuscript's code line) and GPL-3.0 on `main`. Full evidence, options and a recommendation (adopt MIT on both, one commit on `main`) are in `submission/LICENSE_DECISION.md`. |
| 9.4 | Companion repo licensed | `PI ACTION` | github.com/BMGLab/Project_PeakATail has **no licence at all** (verified in the licence brief). |
| 9.5 | Licence named in the manuscript | `PI ACTION` | The Availability section names no licence; the journal requires it stated. Text edit for the MAIN.md owner once 9.3 is decided (`LICENSE_DECISION.md` §"Manuscript" lists this as item 10 and reserves a slot for the wording — coordinate with that workstream rather than writing a second version). |
| 9.6 | Companion repository public | `BLOCKED-BY` A3 (PI/maintainer) | Still `visibility: PRIVATE`, untouched since 2026-08-19; the runbook to public + DOI, with its blockers, is `submission/COMPANION_REPO_PLAN.md`. Referees expect working access at review time (referee simulation M9vi). The repository's history is a `git-filter-repo` rewrite: `maintenance/sync_from_working_dir.sh` is the only supported sync path — a plain pull will not update it. |
| 9.7 | Archived version with a DOI | `BLOCKED-BY` 9.6 | Zenodo cannot mint a DOI for a private repository. The manuscript placeholder and the cover letter both carry this slot. |
| 9.8 | Open-development record | `DONE` | PRs #92, #93, #96, #97, #100 and issues #94, #95, #98, #99, #101 (`manuscript/github/PUBLISHED.md`), cited in the Availability section and in the cover letter. |

---

## 10. Authors, ORCID and authorship

**Requirement `[LIVE]`/`[INFERRED]`.** Author list, affiliations and a designated corresponding
author with an email are required `[LIVE, cover-letter/submission page]`; ORCID for the submitting/
corresponding author is standard Springer Nature practice `[INFERRED]` — the ORCID text found on the
GB guidelines page concerns *reviewer* identification, so treat authors' ORCID as strongly expected
but unverified. Authorship criteria (substantial contribution, approval of the submitted version)
are `[LIVE]` in the BMC editorial policies.

| # | Item | Status | What is needed |
|---|---|---|---|
| 10.1 | Author list and affiliations | `DONE` | Amir Amiri Tabat (Computer Engineering, Faculty of Computer and Informatics Sciences, Ege University) and Yasin Kaymaz (Bioengineering, Faculty of Engineering, Ege University), set verbatim by the PI 2026-09-02. |
| 10.2 | Corresponding author + email + ORCID | `DONE` | Y.K., yasin.kaymaz@ege.edu.tr, ORCID 0000-0002-9725-7536. |
| 10.3 | A.A.T. ORCID | `PI ACTION` | Placeholder in MAIN.md; supply or register one, or confirm it is omitted. |
| 10.4 | All authors approve submission | `PI ACTION` | Asserted in the cover letter; must actually be obtained from A.A.T. before sending. |

---

## 11. Declarations (all seven subheadings are mandatory)

**Requirement `[LIVE]`.** Required subheadings, in this order: **Ethics approval and consent to
participate · Consent for publication · Availability of data and materials · Competing interests ·
Funding · Authors' contributions · Acknowledgements** (plus optional "Authors' information").
*"If any of the sections are not relevant to your manuscript, please include the heading and write
'Not applicable' for that section."*

| # | Item | Status | What is needed |
|---|---|---|---|
| 11.1 | All seven headings present | `DONE` | All seven are in MAIN.md. |
| 11.2 | Heading order | `PI ACTION` (trivial) | The draft puts Availability of data and materials first; the journal lists it third. Reorder at assembly. |
| 11.3 | Ethics approval and consent | `DONE` | "Not applicable — previously published, publicly available datasets only; no new human or animal data were collected." Correct form for public-data-only work `[LIVE: the "Not applicable" instruction]`. |
| 11.4 | Consent for publication | `DONE` | "Not applicable". |
| 11.5 | Competing interests | `PI ACTION` | Placeholder. Required in the manuscript **and** declared in the cover letter — the two must match. |
| 11.6 | Funding | `PI ACTION` | Placeholder. If none, "Not applicable"/"No funding was received" — never invent a grant number. |
| 11.7 | Authors' contributions | `PI ACTION` | Placeholder; a CRediT split is drafted in-comment (A.A.T. — software, methodology, investigation; Y.K. — conceptualization, supervision, methodology, writing). PI to confirm or amend. |
| 11.8 | Acknowledgements | `PI ACTION` | Placeholder; "Not applicable" is acceptable. |

---

## 12. Reporting standards, ethics and research-integrity declarations

| # | Item | Status | What is needed |
|---|---|---|---|
| 12.1 | Reporting checklist (ARRIVE/CONSORT/etc.) | `DONE — not applicable` `[INFERRED]` | No new experiments, no animals handled, no human subjects; the mouse data are a published GEO dataset. |
| 12.2 | Ethics for public human data | `DONE` | Covered by 11.3. |
| 12.3 | AI/LLM disclosure | `PI ACTION` `[INFERRED — verify]` | Springer Nature's stated position (from search results, **not** from a page that would open here: the AI-policy URL 404'd) is that LLMs cannot be authors and that LLM use must be documented in Methods, with AI-assisted copy-editing exempt. Given how this manuscript and its analysis records were produced, the PI should read the current policy and decide what, if anything, is disclosed. Do not guess the wording. |
| 12.4 | Pre-registration disclosure | `DONE` | Unusually strong: gates registered before the final run, git-verified timestamps, disclosed post-hoc sweep, amendments — all in Methods and Fig S11. |

---

## 13. Preprint policy

**Requirement `[LIVE]`, BMC editorial policies.** Posting on a preprint server *"does not constitute
previous publication"*; authors *"should disclose details of preprint posting, including DOI and
licensing terms, upon submission"*.

| # | Item | Status | What is needed |
|---|---|---|---|
| 13.1 | Preprint decision | `PI ACTION` | No preprint is referred to anywhere in the draft or the project records read here, but that is an absence of evidence, not a confirmation — the PI should confirm. If one is or will be posted (e.g. bioRxiv), disclose its DOI and licence at submission, and make sure it does not contradict the cover letter's "not published elsewhere" sentence (it does not, under this policy, but the disclosure is still required). |

---

## 14. Peer review, suggested and opposed reviewers

**Requirement `[LIVE]`.** *"Genome Biology operates a single-anonymous peer review system, where
reviewers are aware of the names and affiliations of the authors, but reviewer reports provided to
authors are anonymous."* Whether the system asks for suggested/opposed reviewers, and how many, is
`[INFERRED]` — not verified.

| # | Item | Status | What is needed |
|---|---|---|---|
| 14.1 | Suggested reviewers | `PI ACTION` | Placeholder in the cover letter. Give 3–5 with institutional emails, none a collaborator or co-author within the usual exclusion window; single-anonymous review means they will see the author names. |
| 14.2 | Opposed reviewers | `PI ACTION` | Placeholder in the cover letter. Name them with a one-line, non-polemical reason. |
| 14.3 | Anticipated referee concerns | `DONE` | `REVIEW_SIMULATION.md` is a full simulated referee report (12 majors, 14 minors) with `PI_DECISIONS.md` A1–A3, B1–B7, C1–C10 as the response queue. The cover letter pre-empts the four most likely concerns (atlas-metric gameability, recall deficit, the ten-permutation floor, donor-mismatched long-read truth) in the authors' own words. |

---

## 15. Cover letter

**Requirement `[LIVE]`, quoted.** The cover letter must contain: *"An explanation of why your
manuscript should be published in Genome Biology; An explanation of any issues relating to journal
policies; A declaration of any potential competing interests; Confirmation that all authors have
approved the manuscript for submission; Confirmation that the content of the manuscript has not been
published, or submitted for publication elsewhere."*

| # | Item | Status | What is needed |
|---|---|---|---|
| 15.1 | Drafted | `DONE` | `submission/cover_letter.md` (~740 words including letterhead; typesets to about one page). |
| 15.2 | Why Genome Biology | `DONE` | Two grounds: the side-by-side benchmark the track asks for, and a reusable evaluation discipline. |
| 15.3 | Journal-policy issues | `PI ACTION` | Once §9 is settled, consider one sentence naming the licence and the DOI. The letter currently states the code record and carries the Zenodo placeholder. |
| 15.4 | Competing interests declared | `PI ACTION` | Placeholder in the letter; must match 11.5 word for word. |
| 15.5 | All authors approve | `DONE` (stated) / `PI ACTION` (obtain) | See 10.4. |
| 15.6 | Not published or under consideration elsewhere | `DONE` | Stated. Verify no preprint or parallel submission contradicts it (§13). |
| 15.7 | Suggested / opposed reviewers | `PI ACTION` | The two placeholders of §14. |

---

## 16. File manifest for the submission

`[INFERRED]` — assembled from what the guidelines describe rather than from a page listing required
files; confirm against the submission system's own upload list.

| # | File | Status | Note |
|---|---|---|---|
| 16.1 | Cover letter | `DONE` | `submission/cover_letter.md` → convert; fill three placeholders. |
| 16.2 | Main manuscript (title page → declarations → references → figure legends) | `PI ACTION` | **Conversion pass DONE 2026-09-03:** `submission/4_BUILD/output/PeakATail_manuscript_submission.docx` — A4, Times New Roman 12 pt, 1.5 spacing, continuous line numbers, "Page X of Y" footer, plain section numbers (1 Background … 5 Methods), title page with superscript affiliation markers and the corresponding-author ORCID, legends collected at the end, **0 embedded images**. Regenerate with `bash submission/4_BUILD/build.sh` (see `submission/BUILD_NOTES.md`). Still `PI ACTION` for content, not format: the Abstract rewrite of §2, the §1 order/abbreviation fixes, and the reference list from `submission/references.bib` per §7.2 — the file still carries 48 visible `[[CITE]]` tokens and no reference list. |
| 16.3 | Figures 1–6 | `DONE` | Six vector PDFs, one per figure, in `manuscript/figures/`. |
| 16.4 | Additional file 1 — supplementary figure legends + supplementary figures | `BLOCKED-BY` §4.5, §5.2, §6.2 | S5/S6 decision, T7 collation, and the 20 MB bundling rule. |
| 16.5 | Additional file — Supplementary Table T7 | `BLOCKED-BY` §5.2 | |
| 16.6 | ORCIDs for both authors | `PI ACTION` | §10.3. |
| 16.7 | Zenodo DOI for the archived companion snapshot | `BLOCKED-BY` §9.6–9.7 | Referenced from both the manuscript and the cover letter. |
| 16.8 | Additional file 1 as a word-processor file | `DONE` (format) | `submission/4_BUILD/output/PeakATail_supplementary.docx`, built from `SUPPLEMENTARY.md`; legends only, 0 embedded images (§6.5). Content blockers unchanged: S5/S6 (§4.5) and T7 (§5.2). |
| 16.9 | Co-author reading copy — **not a submission file** | `DONE` | `submission/4_BUILD/output/PeakATail_manuscript_reading_copy.docx`: identical text with Figs 1–6 embedded at their first citation point and each legend beneath. For internal circulation only; never upload it. |

---

## 17. What actually blocks submission, in order

1. **Abstract rewrite (§2).** The only verified journal limit the draft breaks. ≤100 words and
   unstructured for the Method/Methodology track — a substantive rewrite, not a trim.
2. **Citations (§7).** The bibliography is built and verified (21/23 targets), but the manuscript
   still carries 47 `[[CITE]]` markers and no reference list, and two targets are open: the
   spermatogenesis gene panel with no recorded literature (G1, a PI call) and the Kinnex dataset
   accession (G3, tied to A1).
3. **Licence resolution (§9.3–9.5).** Two contradictory licences on one public repository, plus an
   unlicensed companion repository and no licence named in the manuscript. `LICENSE_DECISION.md`
   has the evidence and a recommendation; the decision is the PI's.
4. **Companion repository public + Zenodo DOI (§9.6–9.7).** Gates the archived-version requirement
   and referee access (`PI_DECISIONS.md` A3).
5. **S5/S6 and Table T7 (§4.5, §5.2).** Build or drop-and-renumber; T7 needs collation only.
6. **The five PI placeholders (§10.3, §11.5–11.8).** ORCID, competing interests, funding, CRediT,
   acknowledgements.
7. **De-duplicated, panel-wide long-read axis (`PI_DECISIONS.md` A1).** Not a journal requirement —
   the referee simulation's top request, and the change most likely to strengthen the headline claim.
   A scheduling decision, not a submission gate, but it also unblocks §7.4 (the Kinnex citation) and
   §4.5 (Fig S5).

---

## 18. Verified-vs-inferred summary

**Verified on the live guidelines pages (2026-09-02, via the proxy route of §0):** article type name
and description; section order; abstract limits for both Methodology (≤100, unstructured) and
Research (≤250, structured); the seven declaration subheadings and the "Not applicable" instruction;
figure formats, ~300 dpi and the 10 MB per-figure limit; the 20 MB per-additional-file limit; table
numbering; single-anonymous peer review; the five cover-letter requirements; the citable-sources
rule; the OSI-licence sentence for Method/Software articles. **Verified in BMC's editorial
policies:** the required availability section, code access statement, DOI-archived software,
preprints, competing interests, authorship criteria, no-duplicate-submission.

**Second-hand:** the Genome Biology source-code editorial (doi:10.1186/s13059-016-1040-y) as quoted
in `submission/LICENSE_DECISION.md`.

**Inferred, not verified — check before relying on any of it:** accepted manuscript file formats;
author ORCID requirements; reference style; additional-file naming and in-text citation convention;
the keyword count range; reporting-standard checklists; the AI/LLM disclosure policy; and the
suggested/opposed-reviewer fields in the submission system.

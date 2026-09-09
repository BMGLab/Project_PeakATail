# BUILD_NOTES.md — how the PeakATail submission Word documents are built

Built 2026-09-03 from the verified draft at `manuscript/00_draft/`.
Everything in `submission/*.docx` is **derived**: nothing here is authored, and the build
never writes into `manuscript/`. Re-run the build after any edit to the draft or the figures.

---

## 1. What each document is for

| File | Audience | Contains | Figures |
|---|---|---|---|
| `PeakATail_manuscript_submission.docx` | **The journal.** The official first-submission copy uploaded as the manuscript file. | Title page → Abstract → Key points → 1 Background → 2 Results → 3 Discussion → 4 Conclusions → 5 Methods → Declarations → Figure legends → Supplementary information | **None embedded.** Figures upload as six separate files (`manuscript/figures/fig{1..6}_*.pdf`, vector; 600 dpi PNG is the fallback). Legends are collected at the end, per journal convention. |
| `PeakATail_manuscript_reading_copy.docx` | **Co-authors.** Circulation copy, comfortable to read start to finish. | Same text, same order, minus the collected legend section | **Six embedded**, each at its first citation point with its legend immediately beneath. |
| `PeakATail_supplementary.docx` | **The journal**, as Additional file 1. | `SUPPLEMENTARY.md`: the S1–S12 legends and the Supplementary tables section. | **None embedded.** Supplementary figures ship as separate files too — see §6 below and `SUBMISSION_CHECKLIST.md` §6.5. |

Neither manuscript file contains a table: the draft has none in the main text
(`SUBMISSION_CHECKLIST.md` §5.1), so "tables inline" in the reading copy is a no-op.
The only table-shaped item is Supplementary Table T7, which is still in preparation.

## 2. Regenerating — the one command

```bash
export LC_ALL=C OMP_NUM_THREADS=1
bash /mnt/ssd1/Projects/PeakATail_wd/submission/4_BUILD/build.sh
```

The script is idempotent, prints every step, and ends with a 38-check verification pass
that exits non-zero on any failure. Tooling used: `/usr/bin/pandoc` 2.9.2.1 and system
`python3` with `python-docx` 1.2.0 and `lxml`.

## 3. What the build does, step by step

Files live in `submission/4_BUILD/`.

**0. Pin the style base.** `pandoc --print-default-data-file reference.docx > reference_base.docx`
— pandoc's own default reference document, so every style its docx writer expects survives.

**1. One reference.docx per profile.** `make_reference.py {submission,reading,supplementary}`
rewrites `docDefaults`, the heading/body styles and the page setup, drops pandoc's blue
heading colour, and appends the three custom styles the build markdown uses
(`TitlePageLine`, `FigureImage`, `FigureLegend`).

| profile | page | margins | body | line spacing | H1 / H2 | legend |
|---|---|---|---|---|---|---|
| submission | A4 | 1.00 in | Times New Roman 12 pt | 1.5 | 14 / 12 pt bold | — |
| reading | A4 | 0.79 in (2 cm) | Times New Roman 11 pt | 1.15 | 14 / 12 pt bold | 10 pt, indented |
| supplementary | A4 | 0.87 in (2.2 cm) | Times New Roman 11 pt | 1.15 | 14 / 12 pt bold | — |

**2. Preprocess the Markdown.** `prepare_md.py` writes `main_submission.md`,
`main_reading.md` and `supplementary.md`. It is assembly only — no claim, number, citation
or author is altered. What it does:

* Drops HTML comments. One is stripped from `MAIN.md` (the author-block note recording the
  PI's verbatim instruction) and none from `SUPPLEMENTARY.md`. That comment carries the
  sixth `[[PLACEHOLDER]]` occurrence, so the .docx shows **11** placeholder occurrences
  where the Markdown shows 12 — all **6 distinct** placeholders are still present and
  visible, because each is also mirrored in Declarations → *Outstanding items*.
* Turns the H1 title into a `Title`-styled block; puts the author line, each affiliation
  line and the corresponding-author line in their own `TitlePageLine` paragraph, so the
  two affiliations do not fold into one paragraph on a Markdown soft break, and so
  pandoc's `^1^` / `^2,\*^` become real superscript runs.
* Promotes `## Abstract` and `## Key points` to H1, so H1 = major section and H2 =
  subsection throughout.
* Writes plain section numbers into the heading text — **pandoc 2.9's docx writer ignores
  `--number-sections`**, so the numbers have to be literal. Numbered: Background 1,
  Results 2 (2.1–2.7), Discussion 3, Conclusions 4, Methods 5 (5.1–5.10). Unnumbered:
  Abstract, Key points, Declarations, Figure legends, Supplementary information — the
  front and back matter journals never number. To number Declarations too, add
  `"Declarations"` to `NUMBERED` in `prepare_md.py`.
* Inserts a page break after the title page, and (submission copy only) before
  *Figure legends*.
* Reading copy only: places each of Figs 1–6 after the body paragraph that **first cites
  it**, followed by that figure's legend taken verbatim from the manuscript's own
  `# Figure legends` section, then drops that collected section because every legend is
  now inline. All six anchors fall in Results:

  | Fig | first cited in | anchor paragraph opens |
  |---|---|---|
  | 1 | 2.1 | "PeakATail calls poly(A) sites (PAS) from the one element…" |
  | 2 | 2.1 | "The internal-priming filter is measured on its own output…" (cites Fig 2a) |
  | 3 | 2.3 | "Tightening the matching window to 10 bp…" |
  | 4 | 2.4 | "Before reporting a single switch we asked whether `peakatail switch diff`…" |
  | 5 | 2.5 | "We next asked whether the calibrated test recovers a known biological program…" |
  | 6 | 2.6 | "We then applied the same discipline where labels are noisier…" |

**3. pandoc → docx.** Verbatim:

```bash
PANDOC=/usr/bin/pandoc
WD=/mnt/ssd1/Projects/PeakATail_wd
FMT='markdown-implicit_figures-smart'

"$PANDOC" -f "$FMT" -t docx \
  --reference-doc="$WD/submission/4_BUILD/reference_submission.docx" \
  --resource-path="$WD/manuscript/figures:$WD" \
  -o "$WD/submission/4_BUILD/output/PeakATail_manuscript_submission.docx" \
  "$WD/submission/4_BUILD/main_submission.md"

"$PANDOC" -f "$FMT" -t docx \
  --reference-doc="$WD/submission/4_BUILD/reference_reading.docx" \
  --resource-path="$WD/manuscript/figures:$WD" \
  -o "$WD/submission/4_BUILD/output/PeakATail_manuscript_reading_copy.docx" \
  "$WD/submission/4_BUILD/main_reading.md"

"$PANDOC" -f "$FMT" -t docx \
  --reference-doc="$WD/submission/4_BUILD/reference_supplementary.docx" \
  --resource-path="$WD/manuscript/figures:$WD" \
  -o "$WD/submission/4_BUILD/output/PeakATail_supplementary.docx" \
  "$WD/submission/4_BUILD/supplementary.md"
```

Both reader extensions are deliberate and **must not be dropped**:

* `-implicit_figures` — a paragraph holding only an image would otherwise become a pandoc
  figure with an empty caption; we want a plain inline image in a `FigureImage` paragraph.
* `-smart` — **required for correctness.** With `smart` on (pandoc's default for
  `markdown`), pandoc rewrites `--` as an en dash and straight quotes as curly ones. That
  silently corrupted the bare command-line flag `--peak-workers` in Results into
  `–peak-workers`. With `smart` off the verified text passes through character for
  character; the cost is that apostrophes and quotation marks stay straight, which is the
  right trade for a file whose numbers and flags are audited.

**4. Postprocess.** `postprocess_docx.py` adds what pandoc 2.9 cannot emit:

* **Continuous line numbers** (`w:lnNumType countBy="1" restart="continuous" distance="360"`)
  on the submission copy only — Genome Biology asks for them; the reading copy is easier to
  read without. Note the omitted `w:start`: Word defaults it to 1, but LibreOffice *adds*
  `w:start` to the count, so `w:start="1"` makes LibreOffice label the first line "2".
* A centred **"Page X of Y" footer** (`PAGE` / `NUMPAGES` fields) on all three, with
  `w:suppressLineNumbers` on the footer paragraph so it stays out of the line count.
* `w:pgNumType start="1"` and the document properties (title, authors).

**5. Verify.** `verify_docx.py` re-opens each .docx with python-docx and checks section
order, image counts and sizes, token counts, superscripts, ORCID, line numbering, the
footer fields, and that **every source paragraph is present in the output** (no-loss).

## 4. Measured state of the current build

| | submission | reading copy | supplementary |
|---|---|---|---|
| file size | 49.2 kB | 4.15 MB | 20.5 kB |
| paragraphs | 153 | 157 | 28 |
| words (as counted in the .docx) | 14,883 | 14,881 | 3,402 |
| embedded images | **0** | **6** | **0** |
| pages (LibreOffice render, A4) | 38 | 27 | 5 |
| `[[CITE]]` occurrences | 48 | 48 | 0 |
| `[[PLACEHOLDER]]` occurrences | 11 (6 distinct) | 11 (6 distinct) | 0 |
| continuous line numbers | yes | no | no |
| page-number footer | yes | yes | yes |

The reading copy's 2-word deficit against the submission copy is exactly the dropped
"Figure legends" heading; every legend itself is present inline.

Embedded image geometry (600 dpi PNGs scaled **by width**, never by pixel count):
all six are 6.40 in wide (text column 6.69 in) at 6.40 × 7.13, 5.82, 5.96, 3.00, 6.36 and
5.91 in — every one inside the 10.12 in usable page height, so no figure is clipped or
forced to rotate. Effective resolution at that width is ~665 dpi.

The six embedded parts are **byte-identical** to `manuscript/figures/fig{1..6}_*.png`
(md5-checked): pandoc stores the PNGs unchanged and only records a display size, so the
reading copy carries the full 600 dpi rasters. A LibreOffice *PDF export* of the reading
copy downsamples them to 300 dpi JPEG — that is an export setting of the render used for
the visual check, not a property of the .docx.

## 5. Verification result (all 38 checks pass)

```
python3 submission/4_BUILD/verify_docx.py /mnt/ssd1/Projects/PeakATail_wd
```

Confirms, per document: heading sequence; **0 / 6 / 0** embedded image parts and the same
number of drawing references; image widths ≤ text column and heights ≤ usable page height;
`[[CITE]]` and `[[PLACEHOLDER]]` counts equal to the build Markdown; continuous line
numbering present only where intended; `PAGE`+`NUMPAGES` in the footer part; the author
line once with superscript runs exactly `['1', '2,*']`; both affiliations as separate
paragraphs each with a superscript marker; `ORCID: 0000-0002-9725-7536` and
`yasin.kaymaz@ege.edu.tr` present; and every source paragraph found in the document text
(148 / 152 / 28 blocks checked, 0 missing).

## 6. Files that ship alongside these documents

Both the main figures and the supplementary figures upload **as separate files** — neither
Word document embeds them (the reading copy is a circulation convenience, not a submission
file).

* Main figures: `manuscript/figures/fig{1..6}_*.pdf` (vector) — PNG at 600 dpi is the fallback.
  `fig1_overview_detailed.*` is the archival data-rich variant and is **not** submitted.
* Supplementary figures: `manuscript/figures/figS{1,2,3,4,7,8,9,10,11,12}_*.pdf` — ten files,
  0.41 MB as PDFs, bundled with `PeakATail_supplementary.docx` as Additional file 1.
  S5 and S6 are named pending slots with real blockers and are cited nowhere.

## 7. Known open items these documents inherit from the draft

The build is faithful assembly; it does not fix anything the draft still owes. Carried
through verbatim and visible in the .docx:

* 48 `[[CITE]]` tokens and **no reference list.** Resolving them is typesetting, and the
  draft's own *Outstanding items* says each token "is left exactly where it stands"; the
  verified bibliography is in `submission/references.bib` (35 entries) with
  `submission/CITATION_MAP.tsv`. Nothing was invented here.
* 6 `[[PLACEHOLDER]]` items (A.A.T. ORCID, Zenodo DOI, competing interests, funding, CRediT
  contributions, acknowledgements), each in place *and* mirrored in *Outstanding items*.
* The Abstract still exceeds the journal's limit (`SUBMISSION_CHECKLIST.md` §2) — a
  rewrite, not a formatting job.
* Supplementary Figs S5 / S6 and Supplementary Table T7 remain "in preparation".

## 8. If a document needs tuning

| Want | Change |
|---|---|
| Double spacing instead of 1.5 | `PROFILES["submission"]["line"] = 480` in `make_reference.py` |
| 11 pt body | `PROFILES["submission"]["body_pt"] = 11` |
| Number Declarations too | add `"Declarations"` to `NUMBERED` in `prepare_md.py` |
| Bigger / smaller inline figures | `READING_IMG_WIDTH_IN` in `prepare_md.py` (keep ≤ 6.69 in) |
| Line numbers on the reading copy | drop `--no-line-numbers` from its `postprocess_docx.py` call in `build.sh` |
| Figure placed elsewhere | it follows the first paragraph matching `\bFig N`; move the citation or special-case it in `prepare_md.py` |

Then re-run `build.sh` — it rebuilds and re-verifies everything.

## 9. Scratch

PDF renders used for the visual check live in the git-ignored
`logs/docxbuild_20260903/` (LibreOffice profile + PDFs); delete freely. Nothing was
written to `/mnt/ssd2` or to `results/benchmark_tools/*final*`, no figure data or
manuscript claim was touched, and nothing was committed.

# Typography spec for the three Word documents

Set by the PI on 2026-09-05 and applied identically to all three outputs. Every value
below is asserted against the built files, not against the build configuration: the
checker resolves each run's *effective* size through its paragraph style chain to
`docDefaults`, so a run that merely inherits the right size still counts as correct.

| Element | Spec | Where it is set |
|---|---|---|
| Main text | 11 pt | `docDefaults` + `Normal`, `BodyText`, `FirstParagraph`, `BlockText` in `make_reference.py` |
| Headings (H1-H3) | bold 12 pt | heading styles in `make_reference.py`; H3 keeps italic as the only remaining level cue |
| Figure and table legends | 10 pt | `FigureLegend` style (reading copy) and a run-level pass in `postprocess_docx.py` (submission copy, where legends sit in one collected section under a body style) |
| Table body text | 9 pt | run-level pass in `postprocess_docx.py` |
| Body paragraphs | justified both edges | `jc="both"` on the four body styles |
| Headings, title, title-page lines | left aligned | `jc="left"` set explicitly so they do not inherit justification from `Normal` |

Typeface is Times New Roman throughout, Courier New for code; A4. Line spacing and page
margin are the only things that still differ by profile: the submission copy is double
spaced (line 480/240) with a 25.4 mm margin for a reviewer to mark up, the reading and
supplementary copies are 1.15 spaced with tighter margins.

## Two things this spec does not currently reach

**Tables: none exist yet.** The 9 pt rule is implemented and verified to run, but it
formatted 0 paragraphs, because the manuscript contains no table. The only one planned,
Supplementary Table T7 (competitor tool versions, parameters, run commands and run-log
caveats), is marked *in preparation* in `SUPPLEMENTARY.md`. The rule will bind the first
time a table is added; nothing further is needed then.

**Two defects were found by reading the built files back, not by the build log.** Both are
fixed, and both are the reason this verification resolves properties rather than trusting
the builder:

1. The legend rule matched on paragraph text opening `Fig N |`, which is also exactly how
   every supplementary *heading* opens. Twelve supplementary headings were being shrunk
   from 12 pt to 10 pt. Paragraph style now wins over the text pattern.
2. Five body paragraphs in the main text carry pandoc's bare `Normal` style and so were
   left ragged-right while the other 86 were justified. `Normal` is now justified at the
   root, with headings and title styles opted out explicitly.

## Reproducing the check

    bash submission/4_BUILD/build.sh
    python3 submission/4_BUILD/check_typography.py submission/4_BUILD/output/PeakATail_manuscript_submission.docx \
                                    submission/4_BUILD/output/PeakATail_manuscript_reading_copy.docx \
                                    submission/4_BUILD/output/PeakATail_supplementary.docx

Last verified 2026-09-05: all three documents pass; body 11 pt, headings bold 12 pt,
18 legend paragraphs at 10 pt in each main copy, 91/91 and 14/14 body paragraphs justified.

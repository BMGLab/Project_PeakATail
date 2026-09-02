# Publication design directives — PI, 2026-09-02 (binding for the main-figure design pass)

Status: RECORDED, awaiting execution. The footer→legend surgery must finish first (same files).
The design-pass verifier must check off every item below explicitly.

1. **Publication designs, not analysis reports** (PI, verbatim intent). One message per panel;
   panels serving the record move to supplements/legends; generous whitespace; unified visual
   system (palette, markers, type hierarchy) across all six mains; minimum print type size;
   overlap of text/lines/plots is a hard failure (bbox checks + visual review at print size).
2. **Subtitles off the plots** — every figure- and panel-level descriptive subtitle line moves to
   the legend sidecar (no-loss rule). Applies to all 16 figures.
3. **Fig 1 simplified for a non-technical reader**: 3–4 plain-language panels (the data and the
   poly(A)-tail evidence; clustering tail-bearing ends into a site; the internal-priming trap and
   the genome-sequence filter; the trusted output feeding calibrated cell-type comparisons).
   At most 1–2 headline-number badges. No trade curve (see item 5). Data-rich version stays
   reachable behind a style switch. Technical verifier: simplification must not create a false
   claim (internal priming is the danger spot).
4. **Development history out of the main figures** (PI: end users don't care how we improved the
   tool): Fig 2 loses "shipped (pre-fix)" and the arrowed dev path — only the current tool's
   user-choosable operating points (precision default; ≥1-molecule sensitivity arm) plus
   competitors and matched-N. Fig 3d compute shows the CURRENT version vs competitors only —
   no v1→v2 arrow. History lives in S8/S12 with one Results sentence pointing there.
   Deliberate exception, PI-informed: Fig 4's calibration finding stays (a scientific result,
   reframed as "which test configurations control FDR", not as self-repair).
5. **Redundancy across mains eliminated** (PI: Fig 1e ≈ Fig 3a): the trade curve appears in
   EXACTLY ONE place — either Fig 3a stays its home, or it folds into Fig 2's P/R plane and
   Fig 3 sharpens to robustness (resolution/replication) — design pass decides, lean = fold into
   Fig 2. Pairwise redundancy test across all main-figure panels.
6. **Sentence-case labels**: axis labels, panel titles and prose tick labels start with a capital
   letter. Canonical identifiers exempt (polyApipe, scAPAtrap, scUTRquant, gene symbols, CLI
   flags, file names keep their casing). Grep-level check in verification.

Downstream after the pass: legend re-sync into manuscript/00_draft (MAIN.md/SUPPLEMENTARY.md),
"(Fig X)" already gone from subtitles there, then the two Word documents (submission copy;
reading copy with mains + tables embedded).

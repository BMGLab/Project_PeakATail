#!/usr/bin/env python3
"""
Regression tests for two figure defects the per-script guards cannot see
(design judge, 2026-09-03; findings 2 and 3).

FINDING 3 -- MATHTEXT UNDER THE 6 pt FLOOR
    Every figure script asserts that no Text artist is smaller than
    `_pubstyle.TYPE["annotation_min"]` (6 pt).  That guard reads
    `Text.get_fontsize()`, which is the size of the STRING, not of the glyphs
    matplotlib actually draws inside it.  A mathtext subscript or superscript
    is drawn at 0.7x the surrounding size, so

        ax.set_xlabel("detected-gene recall R$_{det}$ @100 bp")   # 7.5 pt

    passes the guard while printing "det" at 5.25 pt, and a default log axis
    (LogFormatterSciNotation, "$\\mathdefault{10^{-3}}$") prints its exponent at
    0.7 x 6.5 = 4.55 pt.  Both shipped in the 2026-09-02 renders.  This test
    reads the glyph runs out of the PDFs themselves, so it sees what the guards
    cannot.

FINDING 2 -- A LEGEND DRAWN TWICE
    `fig.legend(...)` already registers the legend on the figure (it appends to
    `fig.legends`).  A following `fig.add_artist(lg)` registers the SAME artist
    a second time in `fig.artists`, and `Figure.draw` walks both lists, so the
    legend is drawn twice.  In the raster the two draws land on identical
    pixels and nothing looks wrong; in the vector PDF every text and marker of
    that legend is emitted twice.  Figs 2 and 3 shipped that way.  This test is
    a source check, because the defect is invisible in the raster the eye
    reviews.

Run:
    export LC_ALL=C OMP_NUM_THREADS=1
    python3 -m pytest scripts/manuscript_figures/tests/test_figure_typography.py -v
"""
import re
import zlib
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
SCRIPTS = HERE.parent
FIGDIR = SCRIPTS.parents[1] / "manuscript" / "figures"

TYPE_FLOOR_PT = 6.0

# The 16 figures of the submission scheme (FIGURES_MANIFEST.md).
SUBMITTED = [
    "fig1_overview", "fig2_accuracy", "fig3_tradeoff", "fig4_calibration",
    "fig5_spermatogenesis", "fig6_cohort",
    "figS1_datasets", "figS2_peakqc", "figS3_nulldesign", "figS4_seconddonor",
    "figS7_novelfunnel", "figS8_compute", "figS9_calibration_extended",
    "figS10_clustering", "figS11_gatehistory", "figS12_versions",
]

# `fig1_overview_detailed` is the archival data-rich render behind FIG1_STYLE=detailed.
# It is NOT a submitted figure and figures/README.md keeps it "deliberately NOT held to
# the design laws of that pass -- its callout density is the reason it is kept", so its
# 5.0-5.8 pt type is a property of the record, not a defect.  Excluded on purpose.
NOT_SUBMITTED = ["fig1_overview_detailed"]


def glyph_run_sizes(pdf: Path):
    """Effective point size of every text-showing operator in a matplotlib PDF.

    Walks each content stream tracking the q/Q graphics stack, the CTM scale set
    by `cm`, and the font size set by `/F<n> <size> Tf`; a mathtext sub/superscript
    shows up here as its own smaller Tf size.  Returns {size: count}.
    """
    raw = pdf.read_bytes()
    sizes = {}
    tok = re.compile(rb"(?:/F\d+\s+([-+0-9.]+)\s+Tf)|"
                     rb"(?:((?:[-+0-9.]+\s+){6})cm)|"
                     rb"(\bq\b)|(\bQ\b)|(\bTJ\b)|(\bTj\b)")
    for m in re.finditer(rb"stream\r?\n(.*?)endstream", raw, re.S):
        try:
            d = zlib.decompress(m.group(1))
        except Exception:
            continue
        if b"Tf" not in d:
            continue
        cur, stack, scale = None, [], 1.0
        for t in tok.finditer(d):
            if t.group(1) is not None:
                cur = float(t.group(1))
            elif t.group(2) is not None:
                a, b = (float(v) for v in t.group(2).split()[:2])
                scale = (a * a + b * b) ** 0.5
            elif t.group(3):
                stack.append(scale)
            elif t.group(4):
                scale = stack.pop() if stack else 1.0
            elif cur:
                k = round(cur * scale, 3)
                sizes[k] = sizes.get(k, 0) + 1
    return sizes


@pytest.mark.parametrize("stem", SUBMITTED)
def test_no_glyph_below_the_type_floor(stem):
    """No glyph in a submitted figure may print below 6 pt -- mathtext included."""
    pdf = FIGDIR / f"{stem}.pdf"
    assert pdf.exists(), f"{pdf} is missing; regenerate the figure"
    sizes = glyph_run_sizes(pdf)
    assert sizes, f"{stem}: no text found -- the PDF parser needs updating, not the figure"
    small = {s: n for s, n in sizes.items() if s < TYPE_FLOOR_PT - 1e-6}
    assert not small, (
        f"{stem}: {sum(small.values())} glyph run(s) below {TYPE_FLOOR_PT} pt: {small}. "
        "A mathtext sub/superscript prints at 0.7x its label size, so either drop the "
        "mathtext (R_det, not R$_{det}$) or label the log axis with _pubstyle.plain_log()."
    )


@pytest.mark.parametrize("stem", NOT_SUBMITTED)
def test_archival_render_is_exempt_but_still_parsed(stem):
    """The archival render is exempt from the floor; only assert it still has text."""
    pdf = FIGDIR / f"{stem}.pdf"
    if not pdf.exists():
        pytest.skip(f"{pdf} not built in this worktree")
    assert glyph_run_sizes(pdf), f"{stem}: no text found"


@pytest.mark.parametrize("stem", SUBMITTED)
def test_no_legend_is_registered_twice(stem):
    """`fig.add_artist(lg)` on a legend `fig.legend()` already returned draws it twice."""
    src = (SCRIPTS / f"{stem}.py").read_text(encoding="utf-8")
    named = set(re.findall(r"^\s*(\w+)\s*=\s*fig\.legend\(", src, flags=re.M))
    readded = set(re.findall(r"^\s*fig\.add_artist\(\s*(\w+)\s*\)", src, flags=re.M))
    both = named & readded
    assert not both, (
        f"{stem}: legend(s) {sorted(both)} are registered by fig.legend() and then again "
        "by fig.add_artist(), so they are drawn twice -- identical pixels in the PNG, "
        "duplicated text and markers in the vector PDF. Drop the add_artist() call."
    )

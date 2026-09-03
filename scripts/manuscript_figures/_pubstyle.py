"""Shared publication style for all manuscript figures (DESIGN_DIRECTIVES.md item 1).

Import as:  from _pubstyle import PAL, TOOL_STYLE, TYPE, apply_rc, sentence_case
Every main-figure script MUST use these constants so the six figures read as one system.
"""
import matplotlib

# Okabe-Ito derived, colourblind-safe (validated 2026-08-19; six-checks pass)
PAL = {
    "peakatail": "#0072B2",   # our tool, always
    "good":      "#009E73",   # positive/confirmed
    "bad":       "#D55E00",   # negative/artefact/warning
    "alt":       "#CC79A7",
    "accent":    "#E69F00",
    "light":     "#56B4E9",
    "neutral":   "#999999",   # nulls, reference bands — never a data series
    "ink":       "#1A1A1A",
    "muted":     "#5A6570",
    "grid":      "#D9DEE3",
}
# One marker/colour identity per tool across every figure (canonical casing preserved)
TOOL_STYLE = {
    "PeakATail":  dict(color="#0072B2", marker="D"),
    "polyApipe":  dict(color="#D55E00", marker="o"),
    "scAPAtrap":  dict(color="#E69F00", marker="o"),
    "SCAPTURE":   dict(color="#009E73", marker="^"),
    "Sierra":     dict(color="#CC79A7", marker="s"),
    "scUTRquant": dict(color="#999999", marker="o", mfc="white"),  # catalog-based: hollow
}
# Type scale (points, at print size — canvases are laid out at final width, <= 7.09 in preferred)
TYPE = {"panel_letter": 10, "panel_title": 8.5, "axis_label": 7.5, "tick": 6.5,
        "annotation": 6.5, "annotation_min": 6.0}

def apply_rc():
    matplotlib.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "DejaVu Sans",
        "font.size": TYPE["axis_label"], "axes.titlesize": TYPE["panel_title"],
        "axes.labelsize": TYPE["axis_label"], "xtick.labelsize": TYPE["tick"],
        "ytick.labelsize": TYPE["tick"], "legend.fontsize": TYPE["annotation"],
        "savefig.dpi": 600, "axes.linewidth": 0.7,
    })

_CANONICAL = ("PeakATail", "polyApipe", "scAPAtrap", "SCAPTURE", "scUTRquant", "Sierra",
              "PolyASite", "Kinnex", "CellRanger", "STARsolo", "GRCh38", "GRCm38",
              "R_det", "P@", "F1", "PBMC", "UMI", "BAM", "GEX", "APA", "PAS", "FDR",
              "pbmc4k", "pbmc_10k_v3", "clip_seeded", "lambda_gradient")

def sentence_case(label: str) -> str:
    """Capitalize the first alphabetic character unless the label starts with a
    canonical identifier whose casing must be preserved (DESIGN_DIRECTIVES.md item 6)."""
    s = label.lstrip()
    for c in _CANONICAL:
        if s.startswith(c):
            return label
    for i, ch in enumerate(label):
        if ch.isalpha():
            return label[:i] + ch.upper() + label[i + 1:]
    return label


def plain_log(ax, axis="y", compact=False):
    """Label a log axis in plain text instead of matplotlib's 10^n mathtext.

    Mathtext superscripts and subscripts are drawn at 0.7x the surrounding font
    size, so a 6.5 pt tick label carries a 4.55 pt exponent -- below the 6 pt
    print floor of DESIGN_DIRECTIVES.md item 1.  The per-script type guards
    measure Text artists, not the glyph runs inside a mathtext string, so they
    cannot see it (design-judge finding 3, 2026-09-03).  This keeps every tick
    label at the tick size.

    Ticks themselves are untouched (the LogLocator still chooses them): only the
    label text changes.  `compact` switches thousands to the 10k / 1M form the
    call-count axes already use; minor ticks stay unlabelled, as the default log
    formatter leaves them on every axis this is applied to.
    """
    from matplotlib.ticker import FuncFormatter, NullFormatter

    def fmt(v, _pos=None):
        if v <= 0:
            return ""
        if compact and v >= 1e6:
            return f"{v / 1e6:g}M"
        if compact and v >= 1e3:
            return f"{v / 1e3:g}k"
        if v >= 1:
            return f"{v:g}"
        return f"{v:.10f}".rstrip("0")

    a = ax.yaxis if axis == "y" else ax.xaxis
    a.set_major_formatter(FuncFormatter(fmt))
    a.set_minor_formatter(NullFormatter())
    return ax

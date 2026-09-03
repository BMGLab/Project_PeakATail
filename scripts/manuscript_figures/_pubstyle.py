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

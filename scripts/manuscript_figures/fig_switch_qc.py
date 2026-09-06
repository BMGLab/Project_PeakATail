"""Two atlas-free checks on PeakATail's replicated switch calls.  Not a manuscript figure.

(a) How far apart are the sites a gene is said to switch between?  A tandem 3'UTR
    switch happens inside one 3'UTR, so spans far beyond ~10 kb are unlikely to be
    that, and point at sites being assigned to a gene they do not terminate.
(b) Does a called site sit at a real transcript 3' end?  In 3'-tag data a genuine
    poly(A) site shows a sharp coverage drop: depth in the 200 bp inside the
    transcript over depth in the 200 bp outside.  This uses no atlas at all, so it
    is independent of the precision figure the paper quotes.

Build: python3 /tmp/step_scan.py 300   (writes site_step_scan.tsv)
       python3 scripts/manuscript_figures/fig_switch_qc.py
"""
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _pubstyle import PAL, TYPE, apply_rc

apply_rc()

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
OUTDIR = WD / "manuscript/figures"
TSVDIR = WD / "results/figures/manuscript"
STEM = "fig_switch_qc"
INK, MUTED, GRID = PAL["ink"], PAL["muted"], PAL["grid"]
BLUE, GOOD, BAD = PAL["peakatail"], PAL["good"], PAL["bad"]

d = pd.read_csv(WD / "results/stage3_spermatogenesis_v2/summary/switch_replicated_pas_true.tsv",
                sep="\t")
g = pd.DataFrame({"span": d.groupby("gene_id").start_m1.agg(lambda s: s.max() - s.min()),
                  "n": d.groupby("gene_id").start_m1.nunique()})
g = g[g.n > 1]
S = pd.read_csv(TSVDIR / "site_step_scan.tsv", sep="\t")
S = S[S.inside >= 1.0]

FIG_W, FIG_H = 7.09, 3.40
fig = plt.figure(figsize=(FIG_W, FIG_H))
texts = []


def T(x, y, s, **kw):
    kw.setdefault("transform", fig.transFigure)
    kw.setdefault("color", INK)
    t = fig.text(x, y, s, **kw)
    texts.append(t)
    return t


def style(ax):
    for k, sp in ax.spines.items():
        sp.set_visible(k in ("left", "bottom"))
        sp.set_color(GRID)
        sp.set_linewidth(0.7)
    ax.tick_params(labelsize=TYPE["tick"], colors=MUTED, length=2.5)
    ax.yaxis.grid(True, color=GRID, lw=0.5)
    ax.set_axisbelow(True)


T(0.030, 0.940, "a", fontsize=TYPE["panel_letter"], fontweight="bold", va="baseline")
T(0.058, 0.940, "How far apart are the switching sites?",
  fontsize=TYPE["panel_title"], va="baseline")
axA = fig.add_axes([0.075, 0.300, 0.385, 0.545])
bins = np.logspace(np.log10(30), np.log10(4e6), 34)
axA.hist(g.span.clip(lower=31), bins=bins, color=BLUE, alpha=0.85, lw=0)
axA.set_xscale("log")
axA.axvline(10000, color=BAD, lw=1.1, ls=(0, (3, 1.6)), zorder=4)
frac = (g.span > 10000).mean()
axA.annotate(f"{frac:.0%} of genes sit\nbeyond 10 kb", xy=(10000, axA.get_ylim()[1] * 0.94),
             xytext=(26000, axA.get_ylim()[1] * 0.92), fontsize=6.3, color=BAD,
             va="top", ha="left", linespacing=1.35)
axA.set_xlabel("Distance between the outermost switching sites in a gene (bp)",
               fontsize=TYPE["axis_label"])
axA.set_ylabel("Genes", fontsize=TYPE["axis_label"])
style(axA)
T(0.075, 0.150, f"{len(g):,} genes with more than one replicated site.\n"
                "A tandem 3′UTR switch happens inside one 3′UTR, so\n"
                "the mass beyond 10 kb is unlikely to be one.",
  fontsize=TYPE["annotation"], color=MUTED, va="top", linespacing=1.45)

T(0.525, 0.940, "b", fontsize=TYPE["panel_letter"], fontweight="bold", va="baseline")
T(0.553, 0.940, "Does a call sit at a real transcript end?",
  fontsize=TYPE["panel_title"], va="baseline")
axB = fig.add_axes([0.570, 0.300, 0.385, 0.545])
st = S.step.clip(lower=0.03, upper=3e4)
axB.hist(st, bins=np.logspace(np.log10(0.03), np.log10(3e4), 34), color=GOOD, alpha=0.85, lw=0)
axB.set_xscale("log")
axB.axvline(3, color=INK, lw=1.1, ls=(0, (3, 1.6)), zorder=4)
axB.axvline(1, color=BAD, lw=1.0, ls=(0, (1.6, 1.6)), zorder=4)
clear, inv = (S.step >= 3).mean(), (S.step < 0.8).mean()
axB.annotate(f"{clear:.0%} drop sharply", xy=(3, axB.get_ylim()[1] * 0.95),
             xytext=(7, axB.get_ylim()[1] * 0.93), fontsize=6.3, color=INK, va="top")
axB.annotate(f"{inv:.0%} inverted:\nmore coverage\noutside than in",
             xy=(1, axB.get_ylim()[1] * 0.55), xytext=(0.038, axB.get_ylim()[1] * 0.55),
             fontsize=6.3, color=BAD, va="top", ha="left", linespacing=1.35)
axB.set_xlabel("Read depth inside the transcript ÷ outside, at the call",
               fontsize=TYPE["axis_label"])
axB.set_ylabel("Sites", fontsize=TYPE["axis_label"])
style(axB)
T(0.570, 0.150, f"{len(S)} replicated sites sampled at random. No atlas is\n"
                f"used here, yet {clear:.1%} land at a sharp 3′ end, against the\n"
                "0.706 atlas agreement the paper reports.",
  fontsize=TYPE["annotation"], color=MUTED, va="top", linespacing=1.45)

fig.canvas.draw()
r = fig.canvas.get_renderer()
allt = list(texts) + [t for ax in fig.axes for t in ax.texts if t.get_text().strip()]
bad = []
for i in range(len(allt)):
    bi = allt[i].get_window_extent(renderer=r)
    for j in range(i + 1, len(allt)):
        bj = allt[j].get_window_extent(renderer=r)
        if bi.x0 < bj.x1 - 1 and bj.x0 < bi.x1 - 1 and bi.y0 < bj.y1 - 1 and bj.y0 < bi.y1 - 1:
            bad.append((allt[i].get_text()[:22], allt[j].get_text()[:22]))
assert not bad, f"text overlap: {bad[:4]}"

for ext, kw in (("png", dict(dpi=600)), ("pdf", {})):
    fig.savefig(OUTDIR / f"{STEM}.{ext}", **kw)
    print("wrote", OUTDIR / f"{STEM}.{ext}")
from PIL import Image
im = np.array(Image.open(OUTDIR / f"{STEM}.png").convert("L"))
e = 8
ink = int((im[:e, :] < 200).sum() + (im[-e:, :] < 200).sum()
          + (im[:, :e] < 200).sum() + (im[:, -e:] < 200).sum())
print(f"edge check: {im.shape[1]}x{im.shape[0]} px, ink in the outer {e} px = {ink}")
assert ink == 0, "content is clipped at the canvas edge"

print(f"span>10kb {frac:.1%} of {len(g):,} genes | step>=3 {clear:.1%}, inverted {inv:.1%} "
      f"of {len(S)} sites")

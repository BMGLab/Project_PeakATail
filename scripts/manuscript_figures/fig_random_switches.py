"""A random, unfiltered look at PeakATail's replicated switch calls.

Not a manuscript figure.  The curated browser panel selects for legibility; this
one is a uniform random sample of 12 genes (seed 20260906) from the 5,361 genes in
the replicated spermatogenesis set, with NOTHING filtered out: genes with three or
more sites, sites hundreds of kilobases apart, small effects and overlapping loci
all stay in.  It is here to answer whether the calls look real in general.

Every site drawn replicated in both mice by construction, so the mouse-2 check is
already inside the input; what is on trial here is whether the coverage supports
the call.

Build: python3 scripts/manuscript_figures/_random_switch_extract.py
       python3 scripts/manuscript_figures/fig_random_switches.py
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
STEM = "fig_random_switches"

STAGES = ["SPC", "RS", "ES"]
STAGE_COL = {"SPC": PAL["light"], "RS": PAL["peakatail"], "ES": PAL["accent"]}
INK, MUTED, GRID = PAL["ink"], PAL["muted"], PAL["grid"]
UP_C, DN_C = PAL["good"], PAL["bad"]
NBIN = 900          # per-panel drawing resolution; wide loci are binned, max per bin

cov = pd.read_csv(TSVDIR / "random_switch_coverage.tsv.gz", sep="\t")
meta = pd.read_csv(TSVDIR / "random_switch_meta.tsv", sep="\t")
sites = pd.read_csv(TSVDIR / "random_switch_sites.tsv", sep="\t")
meta = meta.sort_values("span").reset_index(drop=True)
GENES = list(meta.symbol)
assert len(GENES) == 12, GENES

FIG_W, FIG_H = 7.09, 8.80
fig = plt.figure(figsize=(FIG_W, FIG_H))
texts = []


def T(x, y, s, **kw):
    kw.setdefault("transform", fig.transFigure)
    kw.setdefault("color", INK)
    t = fig.text(x, y, s, **kw)
    texts.append(t)
    return t


def fy(i):
    return 1.0 - i / FIG_H


def fh(i):
    return i / FIG_H


TOP, PH, PGAP = 0.78, 1.10, 0.16
TRACK_H, STRIP_H = 0.255, 0.13
CX, CW = [0.055, 0.545], 0.400

T(0.030, fy(0.22), "A random, unfiltered sample of replicated switch calls",
  fontsize=9.5, fontweight="bold", va="baseline")
T(0.030, fy(0.42), "12 genes drawn uniformly at random (seed 20260906) from the 5,361 genes in the "
                   "replicated set. Nothing was filtered:",
  fontsize=TYPE["annotation"], color=MUTED, va="baseline")
T(0.030, fy(0.575), "genes with three or more sites, sites far apart, and small effects are all kept. "
                    "Mouse 1 coverage; every site shown replicated in both mice.",
  fontsize=TYPE["annotation"], color=MUTED, va="baseline")

for gi, sym in enumerate(GENES):
    col, row = gi % 2, gi // 2
    x0 = CX[col]
    g0 = TOP + row * (PH + PGAP)
    m = meta[meta.symbol == sym].iloc[0]
    lo, hi = int(m.lo), int(m.hi)
    st_ = sites[sites.symbol == sym]

    T(x0, fy(g0 + 0.10), sym, fontsize=7.2, fontweight="bold", style="italic", va="baseline")
    T(x0 + CW, fy(g0 + 0.10),
      f"chr{m.chrom} ({m.strand})  {int(m.n_sites)} sites  "
      f"{m.span/1000:.1f} kb  max |Δ| {m.max_abs_dprop:.2f}",
      fontsize=6.0, color=MUTED, va="baseline", ha="right")

    vecs = {}
    for s_ in STAGES:
        c = cov[(cov.symbol == sym) & (cov.stage == s_)]
        v = np.zeros(hi - lo)
        if len(c):
            v[c.pos.values - lo] = c.depth.values
        v = v / max(int(m[f"n_cells_{s_}"]), 1)
        if len(v) > NBIN:                       # browser-style binning, max per bin
            k = len(v) // NBIN
            v = v[:k * NBIN].reshape(NBIN, k).max(axis=1)
        vecs[s_] = v
    ymax = max(v.max() for v in vecs.values()) or 1.0
    n = len(next(iter(vecs.values())))
    xs = np.linspace(lo, hi, n)

    for si, s_ in enumerate(STAGES):
        ax = fig.add_axes([x0, fy(g0 + 0.18 + (si + 1) * TRACK_H), CW, fh(TRACK_H * 0.86)])
        ax.fill_between(xs, 0, vecs[s_], color=STAGE_COL[s_], lw=0, alpha=0.92, zorder=3)
        for _, r in st_.iterrows():
            ax.axvline(r.pos, color=UP_C if r.sign > 0 else DN_C, lw=0.7,
                       ls=(0, (2.2, 1.5)), zorder=4)
        ax.set_xlim(lo, hi)
        ax.set_ylim(0, ymax * 1.16)
        ax.set_xticks([])
        ax.set_yticks([])
        for k, sp in ax.spines.items():
            sp.set_visible(k == "left")
        ax.spines["left"].set_color(GRID)
        ax.spines["left"].set_linewidth(0.7)
        ax.text(0.004, 0.80, s_, transform=ax.transAxes, fontsize=6.0, fontweight="bold",
                color=STAGE_COL[s_], ha="left", va="center")
        if si == 0:
            ax.text(0.996, 0.80, f"0–{ymax:.2g}/cell", transform=ax.transAxes,
                    fontsize=6.0, color=MUTED, ha="right", va="center")

    axs = fig.add_axes([x0, fy(g0 + 0.18 + 3 * TRACK_H + STRIP_H), CW, fh(STRIP_H * 0.7)])
    axs.set_xlim(lo, hi)
    axs.set_ylim(0, 1)
    axs.axis("off")
    axs.plot([lo, hi], [0.72, 0.72], color=GRID, lw=0.8, zorder=1)
    for _, r in st_.iterrows():
        axs.plot([r.pos], [0.72], marker="v", ms=3.4,
                 color=UP_C if r.sign > 0 else DN_C, zorder=3, clip_on=False)

BOT = TOP + 6 * (PH + PGAP) - PGAP + 0.20
T(0.030, fy(BOT), "Ticks mark every replicated site. Green, usage rises with differentiation; "
                  "orange, usage falls. Tracks share one scale within a gene",
  fontsize=TYPE["annotation"], color=MUTED, va="baseline")
T(0.030, fy(BOT + 0.155), "(printed at its top right) and are binned to the panel width, so a wide "
                          "locus is drawn at lower resolution, not truncated.",
  fontsize=TYPE["annotation"], color=MUTED, va="baseline")

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
assert not bad, f"text overlap: {bad[:5]}"
small = [t.get_text()[:22] for t in allt if t.get_fontsize() < TYPE["annotation_min"]]
assert not small, f"below the {TYPE['annotation_min']} pt floor: {small}"

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

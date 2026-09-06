"""Contact sheet: every tandem 3'UTR switch, so all 130 can be eyeballed at once.

Not a manuscript figure. One cell per gene: mouse-1 coverage for the three
spermatogenic stages on a shared scale, with the proximal and distal sites marked.
Cells are ordered so the ones that hold up in the reads come first.

A gene is marked PASS when both its sites sit at a sharp 3' end (depth inside the
transcript over depth outside >= 3) AND the coverage maximum actually moves between
the two sites in the expected direction across stages. Those two checks are read
evidence, independent of the statistical test that called the switch.

Build: python3 scripts/manuscript_figures/_tandem_sheet_extract.py
       python3 scripts/manuscript_figures/fig_tandem_sheet.py
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
STEM = "fig_tandem_sheet"

STAGES = ["SPC", "RS", "ES"]
SCOL = {"SPC": PAL["light"], "RS": PAL["peakatail"], "ES": PAL["accent"]}
INK, MUTED, GRID = PAL["ink"], PAL["muted"], PAL["grid"]
PROX_C, DIST_C, BAD = PAL["good"], PAL["bad"], PAL["bad"]

cov = pd.read_csv(TSVDIR / "tandem_sheet_coverage.tsv.gz", sep="\t")
M = pd.read_csv(TSVDIR / "tandem_sheet_meta.tsv", sep="\t")
M["clean_ends"] = (M.step_prox >= 3) & (M.step_dist >= 3)
M["pass"] = M.clean_ends & M.peak_shift
M = M.sort_values(["pass", "clean_ends", "max_abs_dprop"], ascending=[False, False, False])
M = M.reset_index(drop=True)

NCOL, NROW = 4, 7
PER = NCOL * NROW
PAGES = int(np.ceil(len(M) / PER))
CELL_W, CELL_H = 0.222, 0.108        # figure fractions
X0, Y0 = 0.045, 0.905
GAPX, GAPY = 0.018, 0.020
TRACK = 0.026

PASS_N = int(M["pass"].sum())
print(f"{len(M)} genes, {PASS_N} pass both read checks; {PAGES} pages of {PER}")

for page in range(PAGES):
    sub = M.iloc[page * PER:(page + 1) * PER]
    fig = plt.figure(figsize=(7.09, 9.10))
    texts = []

    def T(x, y, s, **kw):
        kw.setdefault("transform", fig.transFigure)
        kw.setdefault("color", INK)
        t = fig.text(x, y, s, **kw)
        texts.append(t)
        return t

    T(0.045, 0.968, f"Tandem 3′UTR switches, all {len(M)} of them  ·  page {page+1} of {PAGES}",
      fontsize=10, fontweight="bold", va="baseline")
    T(0.045, 0.950,
      f"Mouse 1 coverage, three stages on a shared scale per gene. Ordered best-first: "
      f"{PASS_N} of {len(M)} pass both read checks.",
      fontsize=TYPE["annotation"], color=MUTED, va="baseline")
    T(0.045, 0.934,
      "PASS = both sites at a sharp 3′ end and the coverage maximum moves between them. "
      "Green tick proximal, orange tick distal.",
      fontsize=TYPE["annotation"], color=MUTED, va="baseline")

    for k, (_, g) in enumerate(sub.iterrows()):
        c, r = k % NCOL, k // NCOL
        x = X0 + c * (CELL_W + GAPX)
        ytop = Y0 - r * (CELL_H + GAPY)
        lo, hi = int(g.lo), int(g.hi)

        vecs = {}
        cg = cov[cov.name == g["name"]]
        for s_ in STAGES:
            v = np.zeros(hi - lo)
            cs = cg[cg.stage == s_]
            if len(cs):
                v[cs.pos.values - lo] = cs.depth.values
            vecs[s_] = v / max(int(g[f"n_cells_{s_}"]), 1)
        ymax = max(v.max() for v in vecs.values()) or 1.0

        ok = bool(g["pass"])
        T(x, ytop + 0.010, g["name"], fontsize=6.6, fontweight="bold", style="italic",
          va="baseline", color=INK if ok else MUTED)
        T(x + CELL_W, ytop + 0.010,
          ("PASS" if ok else ("ends ok" if g.clean_ends else "weak")),
          fontsize=5.6, va="baseline", ha="right",
          color=PROX_C if ok else (MUTED if g.clean_ends else BAD),
          fontweight="bold" if ok else "normal")

        for si, s_ in enumerate(STAGES):
            ax = fig.add_axes([x, ytop - (si + 1) * TRACK, CELL_W, TRACK * 0.86])
            ax.fill_between(np.arange(lo, hi), 0, vecs[s_], color=SCOL[s_], lw=0, alpha=0.9)
            for pos, col in ((int(g.prox), PROX_C), (int(g.dist), DIST_C)):
                ax.axvline(pos, color=col, lw=0.7, ls=(0, (2, 1.4)), zorder=4)
            ax.set_xlim(lo, hi)
            ax.set_ylim(0, ymax * 1.15)
            ax.set_xticks([]); ax.set_yticks([])
            for kk, sp in ax.spines.items():
                sp.set_visible(kk == "left")
            ax.spines["left"].set_color(GRID); ax.spines["left"].set_linewidth(0.6)
            ax.text(-0.012, 0.5, s_, transform=ax.transAxes, fontsize=5.2,
                    color=SCOL[s_], ha="right", va="center", fontweight="bold")

        axs = fig.add_axes([x, ytop - 3 * TRACK - 0.012, CELL_W, 0.010])
        axs.set_xlim(lo, hi); axs.set_ylim(0, 1); axs.axis("off")
        for pos, col in ((int(g.prox), PROX_C), (int(g.dist), DIST_C)):
            axs.plot([pos], [0.75], marker="^", ms=2.8, color=col, clip_on=False)
        axs.text(lo, 0.05, f"{g.span/1000:.1f} kb", fontsize=5.2, color=MUTED,
                 ha="left", va="center")
        axs.text(hi, 0.05, ("short" if g.direction == "shortening" else "long")
                 + f" |Δ|{g.max_abs_dprop:.2f}", fontsize=5.2, color=MUTED,
                 ha="right", va="center")

    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    allt = list(texts) + [t for ax in fig.axes for t in ax.texts if t.get_text().strip()]
    bad = []
    for i in range(len(allt)):
        bi = allt[i].get_window_extent(renderer=rend)
        for j in range(i + 1, len(allt)):
            bj = allt[j].get_window_extent(renderer=rend)
            if bi.x0 < bj.x1 - 1 and bj.x0 < bi.x1 - 1 and bi.y0 < bj.y1 - 1 and bj.y0 < bi.y1 - 1:
                bad.append((allt[i].get_text()[:16], allt[j].get_text()[:16]))
    assert not bad, f"page {page+1} text overlap: {bad[:4]}"
    out = OUTDIR / f"{STEM}_p{page+1}.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print("wrote", out)

M.to_csv(TSVDIR / "tandem_sheet_scored.tsv", sep="\t", index=False)
print("wrote", TSVDIR / "tandem_sheet_scored.tsv")

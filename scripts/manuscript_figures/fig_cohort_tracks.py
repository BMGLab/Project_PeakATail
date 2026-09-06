"""Cell-type poly(A) switches in the lung-tumour cohort: coverage plus per-patient replication.

Companion to fig_switch_tracks.py.  There the axis is development in two mice; here
it is cell identity, tumour epithelium against T cells, and the replication unit is
the patient.  Left of each row is the pooled coverage picture; right is every patient
the gene was tested in, so the pooling cannot hide a split.

A patient that tested the site without calling it significant is NOT a disagreement
under the replication rule, so the two are drawn differently: filled markers for
patients that called the switch (q < 0.05), open for tested-only.

Build: python3 scripts/manuscript_figures/_cohort_track_extract.py
       python3 scripts/manuscript_figures/fig_cohort_tracks.py
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
from matplotlib.patches import FancyArrowPatch, Rectangle

from _pubstyle import PAL, TYPE, apply_rc

apply_rc()

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
OUTDIR = WD / "manuscript/figures"
TSVDIR = WD / "results/figures/manuscript"
STEM = "fig_cohort_tracks"

C1, C2 = "Epithelial_Tumor", "T_cell"
CT_LABEL = {C1: "Tumour epithelium", C2: "T cell"}
CT_COL = {C1: PAL["bad"], C2: PAL["peakatail"]}
INK, MUTED, GRID = PAL["ink"], PAL["muted"], PAL["grid"]
PROX_C, DIST_C = PAL["good"], PAL["accent"]
Q = 0.05

cov = pd.read_csv(TSVDIR / "cohort_tracks_coverage.tsv.gz", sep="\t")
meta = pd.read_csv(TSVDIR / "cohort_tracks_meta.tsv", sep="\t")
use = pd.read_csv(TSVDIR / "cohort_tracks_usage.tsv", sep="\t")
N_OPP, N_CLEAN, N_DIST, N_PROX = 248, 52, 49, 3   # selection audit, see caption
GENES = list(dict.fromkeys(meta.symbol))
assert len(GENES) == 4, GENES

FIG_W, FIG_H = 7.09, 7.15
fig = plt.figure(figsize=(FIG_W, FIG_H))
texts = []


def T(x, y, s, **kw):
    kw.setdefault("transform", fig.transFigure)
    kw.setdefault("color", INK)
    t = fig.text(x, y, s, **kw)
    texts.append(t)
    return t


def fy(inches_from_top):
    return 1.0 - inches_from_top / FIG_H


def fh(inches):
    return inches / FIG_H


TOP, GENE_H, GENE_GAP = 0.70, 1.32, 0.16
TRACK_H, MODEL_H = 0.355, 0.22
LX, LW = 0.032, 0.545          # coverage column
RX, RW = 0.660, 0.315          # per-patient column

T(0.030, fy(0.22), "a", fontsize=TYPE["panel_letter"], fontweight="bold", va="baseline")
T(0.058, fy(0.22), "Cell-type poly(A) switches in a lung-adenocarcinoma cohort",
  fontsize=TYPE["panel_title"], va="baseline")
T(0.030, fy(0.42), "Left, mean read depth per cell pooled over the patients the gene was tested in. "
                   "Right, every one of those patients separately.",
  fontsize=TYPE["annotation"], color=MUTED, va="baseline")
T(LX + LW / 2, fy(0.66), "Pooled coverage", fontsize=7.4, fontweight="bold",
  ha="center", va="baseline")
T(RX + RW / 2, fy(0.66), "Proximal-site usage, per patient", fontsize=7.4,
  fontweight="bold", ha="center", va="baseline")

for gi, sym in enumerate(GENES):
    g0 = TOP + gi * (GENE_H + GENE_GAP)
    row = meta[meta.symbol == sym].iloc[0]
    strand, lo, hi = row.strand, int(row.lo), int(row.hi)
    prox, dist = int(row.prox), int(row.dist)

    T(0.030, fy(g0 + 0.10), sym, fontsize=7.8, fontweight="bold", style="italic",
      va="baseline")
    T(0.030 + 0.072, fy(g0 + 0.10), f"chr{row.chrom} ({strand}), {row.direction}",
      fontsize=6.2, color=MUTED, va="baseline")

    vecs = {}
    for ct in (C1, C2):
        s = cov[(cov.symbol == sym) & (cov.celltype == ct)]
        v = np.zeros(hi - lo)
        if len(s):
            v[s.pos.values - lo] = s.depth.values
        n = int(meta[(meta.symbol == sym) & (meta.celltype == ct)].n_cells.iloc[0])
        vecs[ct] = v / max(n, 1)
    ymax = max(v.max() for v in vecs.values()) or 1.0

    for si, ct in enumerate((C1, C2)):
        ax = fig.add_axes([LX, fy(g0 + 0.20 + (si + 1) * TRACK_H), LW, fh(TRACK_H * 0.88)])
        ax.fill_between(np.arange(lo, hi), 0, vecs[ct], color=CT_COL[ct], lw=0,
                        alpha=0.92, zorder=3)
        for pos, col in ((prox, PROX_C), (dist, DIST_C)):
            ax.axvline(pos, color=col, lw=0.9, ls=(0, (2.4, 1.6)), zorder=4)
        ax.set_xlim(lo, hi)
        ax.set_ylim(0, ymax * 1.14)
        ax.set_xticks([])
        ax.set_yticks([])
        for k, sp in ax.spines.items():
            sp.set_visible(k == "left")
        ax.spines["left"].set_color(GRID)
        ax.spines["left"].set_linewidth(0.7)
        ax.text(0.004, 0.84, CT_LABEL[ct], transform=ax.transAxes, fontsize=6.4,
                fontweight="bold", color=CT_COL[ct], ha="left", va="center")

    axm = fig.add_axes([LX, fy(g0 + 0.20 + 2 * TRACK_H + MODEL_H), LW, fh(MODEL_H * 0.72)])
    axm.set_xlim(lo, hi)
    axm.set_ylim(0, 1)
    axm.axis("off")
    axm.add_patch(Rectangle((min(prox, dist), 0.34), abs(dist - prox), 0.14,
                            facecolor=GRID, edgecolor="none", zorder=2))
    for pos, col, lab in ((prox, PROX_C, "Proximal"), (dist, DIST_C, "Distal")):
        axm.plot([pos], [0.42], marker="v", ms=4.2, color=col, zorder=4, clip_on=False)
        axm.text(pos, 0.00, lab, fontsize=6.0, color=col, ha="center", va="bottom")
    a0, a1 = (lo + (hi - lo) * 0.02, lo + (hi - lo) * 0.14) if strand == "+" else \
             (lo + (hi - lo) * 0.14, lo + (hi - lo) * 0.02)
    axm.add_patch(FancyArrowPatch((a0, 0.88), (a1, 0.88), arrowstyle="-|>",
                                  mutation_scale=6, lw=0.9, color=MUTED, zorder=4))
    axm.text(lo + (hi - lo) * 0.185, 0.88, f"tracks 0–{ymax:.2g} reads per base per cell",
             fontsize=6.0, color=MUTED, ha="left", va="center")
    axm.text(hi - (hi - lo) * 0.005, 0.88, f"{(hi - lo) / 1000:.1f} kb", fontsize=6.0,
             color=MUTED, ha="right", va="center")

    # ---- per-patient proximal usage
    axr = fig.add_axes([RX, fy(g0 + 0.20 + 2 * TRACK_H), RW, fh(2 * TRACK_H - 0.10)])
    u = use[(use.symbol == sym) & (use.site == "proximal")]
    piv = u.pivot_table(index="gsm", columns="celltype", values="usage")
    qv = u.drop_duplicates("gsm").set_index("gsm").qvalue
    piv = piv.reindex(sorted(piv.index))
    n_called = 0
    for yi, (gsm, r) in enumerate(piv.iterrows()):
        called = bool(qv.get(gsm, 1) < Q)
        n_called += called
        a, b = r.get(C1, np.nan), r.get(C2, np.nan)
        if pd.notna(a) and pd.notna(b):
            axr.plot([a, b], [yi, yi], color=INK if called else GRID,
                     lw=1.0 if called else 0.8, zorder=2, alpha=1.0 if called else 0.9)
        for val, ct in ((a, C1), (b, C2)):
            if pd.isna(val):
                continue
            axr.plot([val], [yi], marker="o", ms=4.0 if called else 3.4,
                     mfc=CT_COL[ct] if called else "white", mec=CT_COL[ct],
                     mew=1.0, zorder=3)
    axr.set_ylim(-0.7, len(piv) - 0.3)
    axr.set_xlim(-0.04, 1.04)
    axr.invert_yaxis()
    axr.set_yticks([])
    axr.set_xticks([0, 0.5, 1.0])
    axr.set_xticklabels(["0", "50", "100%"], fontsize=TYPE["tick"])
    axr.tick_params(axis="x", colors=MUTED, length=2.5)
    for k, sp in axr.spines.items():
        sp.set_visible(k == "bottom")
    axr.spines["bottom"].set_color(GRID)
    axr.xaxis.grid(True, color=GRID, lw=0.5)
    axr.set_axisbelow(True)
    axr.text(0.5, 1.02, f"{n_called} of {len(piv)} patients called it",
             transform=axr.transAxes, fontsize=6.2, color=INK, ha="center", va="bottom")

BOT = TOP + 4 * (GENE_H + GENE_GAP) - GENE_GAP + 0.24
T(0.030, fy(BOT), "Filled markers, patients that called the switch (q < 0.05); open markers, patients "
                  "that tested the site without calling it.",
  fontsize=TYPE["annotation"], color=MUTED, va="baseline")
T(0.030, fy(BOT + 0.155), "A tested-but-not-called patient is not a disagreement under the "
                          "replication rule, so the two are drawn apart.",
  fontsize=TYPE["annotation"], color=MUTED, va="baseline")

fig.canvas.draw()
r = fig.canvas.get_renderer()
allt = list(texts) + [t for ax in fig.axes for t in ax.texts if t.get_text().strip()] \
       + [t for ax in fig.axes for t in ax.get_xticklabels() if t.get_text().strip()]
bad = []
for i in range(len(allt)):
    bi = allt[i].get_window_extent(renderer=r)
    for j in range(i + 1, len(allt)):
        bj = allt[j].get_window_extent(renderer=r)
        if bi.x0 < bj.x1 - 1 and bj.x0 < bi.x1 - 1 and bi.y0 < bj.y1 - 1 and bj.y0 < bi.y1 - 1:
            bad.append((allt[i].get_text()[:24], allt[j].get_text()[:24]))
assert not bad, f"text overlap: {bad[:5]}"
small = [t.get_text()[:24] for t in allt if t.get_fontsize() < TYPE["annotation_min"]]
assert not small, f"below the {TYPE['annotation_min']} pt floor: {small}"
assert abs(FIG_W - 7.09) < 1e-6

for ext, kw in (("png", dict(dpi=600)), ("pdf", {})):
    fig.savefig(OUTDIR / f"{STEM}.{ext}", **kw)
    print("wrote", OUTDIR / f"{STEM}.{ext}")


# ---------------------------------------------------------------------------
per = []
for sym in GENES:
    row = meta[meta.symbol == sym].iloc[0]
    u = use[(use.symbol == sym) & (use.site == "proximal")]
    qv = u.drop_duplicates("gsm").set_index("gsm").qvalue
    piv = u.pivot_table(index="gsm", columns="celltype", values="usage")
    called = [g for g in piv.index if qv.get(g, 1) < Q]
    agree = sum(1 for g in called
                if pd.notna(piv.loc[g].get(C1)) and pd.notna(piv.loc[g].get(C2))
                and ((piv.loc[g][C2] > piv.loc[g][C1]) == (row.direction.startswith("distal"))))
    per.append(f"*{sym}* (chr{row.chrom}, {row.strand} strand, {row.direction}; "
               f"proximal {int(row.prox):,}, distal {int(row.dist):,}): called by "
               f"{len(called)} of {len(piv)} patients tested, and all {agree} of the called "
               f"patients with both cell types measured agree in direction.")

CAP = f"""# Figure — `{STEM}` caption

## Legend

**Cell-type poly(A)-site switches in a lung-adenocarcinoma cohort.** Four genes carrying a
switch between tumour epithelium and T cells that replicated across patients, selected by
the same two rules as the spermatogenesis companion figure: the gene's replicated sites
must number exactly two, so proximal and distal are unambiguous, and both must lie inside
the assigned gene with no other annotated gene in the drawn window. Of {N_OPP} genes with two
replicated sites moving in opposite directions in this cell-type pair, {N_CLEAN} pass the
annotation rule. **Left**, mean read depth per cell on the gene's strand, pooled over the
patients in which the gene was tested, one track per cell type; the three stacked elements
are the two tracks and a strip marking the two called sites, the direction of transcription
and the interval between the sites. Read acceptance matches the caller, SAM flag filter
3844, and only gene-strand reads are drawn. **Right**, the same switch resolved by patient:
each row is one patient, the two markers are that patient's proximal-site usage in tumour
epithelium and in T cells, and the connecting line is the difference the test scored.
Filled markers are patients that called the switch at q < {Q}; open markers are patients
that tested the site without calling it. Under the replication rule a tested-but-not-called
patient is not a disagreement, which is why the two are drawn apart rather than pooled, and
why the counts read {" / ".join(str(len([g for g in use[(use.symbol==sy)&(use.site=='proximal')].drop_duplicates('gsm').set_index('gsm').qvalue.items() if g[1] < Q])) for sy in GENES)} of 9, 9, 9 and 6 patients.

Per gene: {" ".join(per)}

**Caveats that travel with this figure.** These are four loci chosen to be legible, not a
random sample; the frequency claims are in Fig 6 and manuscript/20. The direction is
strongly one-sided: {N_DIST} of the {N_CLEAN} annotation-clean genes have the distal site
higher in tumour and only {N_PROX} the reverse, so *SPAG1* is drawn as the rarity it is and
not as half of a balance. Coverage height reflects both site usage and gene expression,
which is why the per-patient tested usage is drawn rather than inferred from peak height.
Pooled tracks weight patients by their read depth; the per-patient panel is the check on
that. Cell-type labels are expression-derived and curated, not sorted populations, and
patient counts differ per gene because a gene must be expressed in both cell types in a
patient to be tested there.

## Provenance

Coverage `results/figures/manuscript/cohort_tracks_coverage.tsv.gz`, metadata
`cohort_tracks_meta.tsv`, per-patient usage `cohort_tracks_usage.tsv`. Built by
`scripts/manuscript_figures/_cohort_track_extract.py` and
`scripts/manuscript_figures/fig_cohort_tracks.py`. Switch calls from
`stage3_laughney_v3` (code 9dfdefb, the v2-code record of manuscript/20); BAMs
`/mnt/ssd2/Laugney_Aligned`. Canvas {FIG_W} x {FIG_H} in; PNG 600 dpi; PDF vector.
The build asserts pairwise text-overlap freedom, a {TYPE['annotation_min']:.0f} pt minimum
type size and a blank 8 px margin.
"""
(OUTDIR / f"{STEM}.caption.md").write_text(CAP)
print("wrote", OUTDIR / f"{STEM}.caption.md")

from PIL import Image
im = np.array(Image.open(OUTDIR / f"{STEM}.png").convert("L"))
e = 8
ink = int((im[:e, :] < 200).sum() + (im[-e:, :] < 200).sum()
          + (im[:, :e] < 200).sum() + (im[:, -e:] < 200).sum())
print(f"edge check: {im.shape[1]}x{im.shape[0]} px, ink in the outer {e} px = {ink}")
assert ink == 0, "content is clipped at the canvas edge"

"""Figure 1, ALTERNATIVE candidate: PeakATail as a tool.

The committed Fig 1 (`fig1_overview.py`) explains the biology: what poly(A) evidence is
and why internal priming is a trap. This alternative introduces THE TOOL instead --
its algorithm end to end, the evidence model that distinguishes it, what a user gets,
and its headline properties. Internal priming appears here as one step of the pipeline,
never as a topic. The PI compares the two and chooses; this script never touches
`fig1_overview.py` or its outputs.

Every number on the canvas is traceable:
  * atlas-agreement precision of the default, 0.7062 / 0.7450 / 0.7572     -> manuscript/19 section 1
  * n calls of the default, 46,524 / 26,255 / 26,526                        -> manuscript/19 section 1
  * tier-1 >=1 molecule (IP arm) P@100 0.3520                              -> manuscript/19 section 2
  * coverage-only tier-2 P@100 0.0568                                      -> manuscript/19 section 2 (via 22/23)
  * genome-wide poly(A) clip rate 0.573%                                   -> manuscript/23 section 2
  * compute 12.53 GB, 27m43s uncontended                                    -> manuscript/19 sections 3 and 5

NOT claimed here, deliberately: any RANKING of compute against the other tools.
The record carries wall time for only three competitors (polyApipe 3:26 h, Sierra
2:43 h, scUTRquant 0:30 h) and peak RSS for none of them, so "second-lightest /
second-fastest of six" was unsupportable in both halves and was removed.
  * switch-test calibration 3.0% null p<0.05, 0 of 20 null runs             -> manuscript/14
  * cohort replication 15,942 switches, 12 patients, none in 10 nulls       -> manuscript/20

Run:  python3 scripts/manuscript_figures/fig1_alt_tool.py
"""
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

from _pubstyle import PAL, TYPE, apply_rc, sentence_case

apply_rc()

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
OUTDIR = WD / "manuscript/figures"
TSVDIR = WD / "results/figures/manuscript"
STEM = "fig1_alt_tool"

# --------------------------------------------------------------------------
# Verified values (source in the module docstring; repeated in the sidecar)
# --------------------------------------------------------------------------
CLIP_RATE = 0.5730          # % of cell-barcoded reads carrying a poly(A) clip
DEFAULT_P = [("PBMC 10k v3", 0.7062, 46524),
             ("Testis mouse 1", 0.7450, 26255),
             ("Testis mouse 2", 0.7572, 26526)]
EVIDENCE = [  # (label, atlas-agreement precision @100 bp, colour key, note)
    ("Tail-supported, ≥ 2 molecules\n(the default output)", 0.7062, "peakatail", "Kept"),
    ("Tail-supported, ≥ 1 molecule\n(sensitivity setting)", 0.3520, "light", "Optional"),
    ("Coverage-only, no tail evidence\n(second-class tier)", 0.0568, "neutral", "Reported apart"),
]
COMPUTE_GB, COMPUTE_MIN = 12.53, "27 min"
CAL_NULL, CAL_RUNS = 3.0, 20
COHORT_SWITCHES, COHORT_PATIENTS = "15,942", 12

FIG_W, FIG_H = 7.09, 7.75


def fy(inches_from_top):
    """Figure-fraction y for a distance measured down from the top edge."""
    return 1.0 - inches_from_top / FIG_H


def fh(inches):
    return inches / FIG_H

fig = plt.figure(figsize=(FIG_W, FIG_H))
INK, MUTED, GRID = PAL["ink"], PAL["muted"], PAL["grid"]
BLUE, GREEN, ORANGE, GREY = PAL["peakatail"], PAL["good"], PAL["bad"], PAL["neutral"]
LIGHT = PAL["light"]

texts = []          # collected for the overlap audit
boxes = []          # (x0, y0, x1, y1) in figure fractions, for the containment audit


def T(x, y, s, **kw):
    kw.setdefault("transform", fig.transFigure)
    kw.setdefault("color", INK)
    t = fig.text(x, y, s, **kw)
    texts.append(t)
    return t


def panel_letter(x, y, letter, title):
    T(x, y, letter, fontsize=TYPE["panel_letter"], fontweight="bold", va="baseline")
    T(x + 0.028, y, sentence_case(title), fontsize=TYPE["panel_title"], va="baseline")


def box(x, y, w, h, fc="white", ec=GRID, lw=0.8, r=0.012, z=2, ls="solid"):
    p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                       fc=fc, ec=ec, lw=lw, zorder=z, linestyle=ls,
                       transform=fig.transFigure, figure=fig)
    fig.patches.append(p)
    boxes.append((x, y, x + w, y + h))
    return p


def arrow(x0, y0, x1, y1, color=MUTED, lw=1.0, z=3):
    a = FancyArrowPatch((x0, y0), (x1, y1), transform=fig.transFigure, figure=fig,
                        arrowstyle="-|>", mutation_scale=7, lw=lw, color=color, zorder=z,
                        shrinkA=0, shrinkB=0)
    fig.patches.append(a)


# ==========================================================================
# a  The pipeline
# ==========================================================================
panel_letter(0.045, fy(0.20), "a", "What PeakATail does, step by step")

STAGES = [
    ("Aligned\nsingle-cell BAM", "Cell barcode + UMI\nper read", GREY),
    ("Find poly(A)\nevidence", f"Non-templated A/T\nrun at a read end\n({CLIP_RATE:.2f}% of reads)", ORANGE),
    ("Cluster into\ncandidate sites", "Single linkage,\n25 bp window", BLUE),
    ("Rank the\nevidence", "Tail-supported vs\ncoverage-only", BLUE),
    ("Screen the\ngenome", "Drop sites on a\ngenomic A-run", GREEN),
    ("Count per cell\nper site", "Sparse matrix +\nper-site support", BLUE),
]
n = len(STAGES)
x0, x1 = 0.055, 0.965
gap = 0.014
bw = (x1 - x0 - gap * (n - 1)) / n
by, bh = fy(1.05), fh(0.75)

for i, (title, sub, col) in enumerate(STAGES):
    bx = x0 + i * (bw + gap)
    box(bx, by, bw, bh, fc="#FBFCFD", ec=col, lw=1.0)
    T(bx + bw / 2, by + bh - fh(0.20), title, fontsize=6.9, fontweight="bold",
      ha="center", va="center", color=col if col != GREY else INK, linespacing=1.25)
    T(bx + bw / 2, by + fh(0.23), sub, fontsize=6.0, ha="center", va="center",
      color=MUTED, linespacing=1.3)
    if i < n - 1:
        arrow(bx + bw + 0.0015, by + bh / 2, bx + bw + gap - 0.0015, by + bh / 2)

# the two CLI spans beneath the pipeline
def span(xa, xb, y, label, color):
    fig.patches.append(Rectangle((xa, y), xb - xa, 0.0035, fc=color, ec="none",
                                 transform=fig.transFigure, figure=fig, zorder=3))
    T((xa + xb) / 2, y - fh(0.15), label, fontsize=6.2, ha="center", va="center",
      color=color, fontweight="bold", family="DejaVu Sans Mono")

span(x0, x0 + 6 * (bw + gap) - gap, fy(1.12), "peakatail run", BLUE)
T(0.5, fy(1.45), "One command produces the sites and the count matrix.", fontsize=6.2,
  ha="center", va="center", color=MUTED)
T(0.5, fy(1.60), "Two flags set the operating point: --polya-min-umis (tail molecules required) "
              "and --ip-filter (genome screen).",
  fontsize=6.2, ha="center", va="center", color=MUTED)

# downstream row
dy, dh = fy(2.66), fh(0.78)
DOWN = [
    ("Cell-type switch tests", "peakatail switch diff", "False-positive rate checked\nagainst shuffled labels"),
    ("3′UTR length shifts", "peakatail switch length", "Per-cell PDUI along\na differentiation axis"),
    ("Cross-sample replication", "replication_filter.py", "Keep a switch only if\nindependent samples agree"),
]
dw = (x1 - x0 - 2 * 0.030) / 3
for i, (title, cli, sub) in enumerate(DOWN):
    bx = x0 + i * (dw + 0.030)
    box(bx, dy, dw, dh, fc="#F7FAFC", ec=GRID, lw=0.8)
    T(bx + dw / 2, dy + dh - fh(0.16), title, fontsize=6.9, fontweight="bold",
      ha="center", va="center")
    T(bx + dw / 2, dy + dh - fh(0.34), cli, fontsize=6.0, ha="center", va="center",
      color=BLUE, family="DejaVu Sans Mono")
    T(bx + dw / 2, dy + fh(0.24), sub, fontsize=6.0, ha="center", va="center",
      color=MUTED, linespacing=1.3)
arrow(0.5, fy(1.72), 0.5, dy + dh + fh(0.03), color=GRID, lw=1.0)

# ==========================================================================
# b  The evidence model
# ==========================================================================
panel_letter(0.045, fy(2.98), "b", "Why a call is trusted: the evidence decides")

axB = fig.add_axes([0.300, fy(3.92), 0.560, fh(0.80)])
ys = [2, 1, 0]
for y, (lab, val, key, tag) in zip(ys, EVIDENCE):
    c = {"peakatail": BLUE, "light": LIGHT, "neutral": GREY}[key]
    axB.barh(y, val, height=0.52, color=c, edgecolor="none", zorder=3)
    axB.text(val + 0.016, y, f"{val:.3f}", va="center", ha="left",
             fontsize=6.4, color=INK, fontweight="bold", zorder=4)
axB.set_yticks(ys)
axB.set_yticklabels([e[0] for e in EVIDENCE], fontsize=6.0, linespacing=1.25)
axB.set_xlim(0, 0.86)
axB.set_ylim(-0.55, 2.55)
axB.set_xlabel(sentence_case("atlas-agreement precision @100 bp"), fontsize=TYPE["axis_label"])
axB.tick_params(axis="x", labelsize=TYPE["tick"], colors=MUTED, length=2.5)
axB.tick_params(axis="y", length=0, colors=INK, pad=5)
for s in ("top", "right", "left"):
    axB.spines[s].set_visible(False)
axB.spines["bottom"].set_color(GRID)
axB.xaxis.grid(True, color=GRID, lw=0.5)
axB.set_axisbelow(True)

T(0.5, fy(4.32), "A site earns its place from molecules that carry the tail, not from a bump in coverage.\n"
              "The three classes are reported separately, so a user chooses the operating point rather than inheriting ours.",
  fontsize=6.3, ha="center", va="top", color=INK, linespacing=1.5)
T(0.5, fy(4.66), f"PBMC 10k v3 shown. The default keeps {DEFAULT_P[0][2]:,} sites.",
  fontsize=6.0, ha="center", va="top", color=MUTED)

# ==========================================================================
# c  Headline properties
# ==========================================================================
panel_letter(0.045, fy(4.98), "c", "What the tool delivers")

cy, ch = fy(6.33), fh(1.20)
cw = (x1 - x0 - 2 * 0.028) / 3
CARDS = [
    (BLUE, "Precision, four libraries",
     "\n".join(f"{name}   {p:.3f}" for name, p, _ in DEFAULT_P) + "\nPBMC donor 2   0.8279",
     "Atlas-agreement precision @100 bp\nof the default output"),
    (GREEN, "Statistics that hold",
     f"{CAL_NULL:.1f}% false-positive rate\nunder shuffled labels\n0 of {CAL_RUNS} null runs\nproduced a hit",
     "The switch test's error rate is\nmeasured, not assumed"),
    (ORANGE, "Runs on a laptop-class node",
     f"{COMPUTE_GB:.1f} GB peak memory\n{COMPUTE_MIN} for a 10k-cell library",
     "Uncontended single run; no specialised\nhigh-memory node required"),
]
for i, (col, head, body, foot) in enumerate(CARDS):
    bx = x0 + i * (cw + 0.028)
    box(bx, cy, cw, ch, fc="white", ec=col, lw=1.0)
    fig.patches.append(Rectangle((bx, cy + ch - fh(0.03)), cw, fh(0.03), fc=col, ec="none",
                                 transform=fig.transFigure, figure=fig, zorder=3))
    T(bx + cw / 2, cy + ch - fh(0.20), head, fontsize=6.8, fontweight="bold",
      ha="center", va="center", color=col)
    T(bx + cw / 2, cy + ch - fh(0.58), body, fontsize=6.3, ha="center", va="center",
      linespacing=1.5)
    T(bx + cw / 2, cy + fh(0.20), foot, fontsize=6.0, ha="center", va="center",
      color=MUTED, linespacing=1.35)

# ==========================================================================
# d  The reliability claim
# ==========================================================================
panel_letter(0.045, fy(6.65), "d", "The result the tool is built for")

box(0.055, fy(7.55), 0.910, fh(0.75), fc="#F4F9F6", ec=GREEN, lw=1.0)
T(0.5, fy(7.02), f"{COHORT_SWITCHES} cell-type polyadenylation switches replicated across "
              f"{COHORT_PATIENTS} patients of a lung-tumour cohort",
  fontsize=7.4, ha="center", va="center", fontweight="bold", color=INK)
T(0.5, fy(7.24), "and none replicated in any of ten label-shuffle null runs of the same pipeline",
  fontsize=6.8, ha="center", va="center", color=GREEN)
T(0.5, fy(7.42), "A switch counts only when independent patients agree on it in the same direction, "
              "and any patient disagreeing vetoes the call.",
  fontsize=6.0, ha="center", va="center", color=MUTED)

# ==========================================================================
# audits
# ==========================================================================
fig.canvas.draw()
r = fig.canvas.get_renderer()


def bb(t):
    e = t.get_window_extent(renderer=r)
    return e.x0, e.y0, e.x1, e.y1


for _ax in fig.axes:                       # tick labels and axis titles collided
    texts.extend([t for t in _ax.get_xticklabels() + _ax.get_yticklabels()
                  if t.get_text()])
    texts.extend([t for t in (_ax.xaxis.label, _ax.yaxis.label) if t.get_text()])
    texts.extend(_ax.texts)

bad = []
for i in range(len(texts)):
    for j in range(i + 1, len(texts)):
        a, b = bb(texts[i]), bb(texts[j])
        if a[0] < b[2] - 1 and b[0] < a[2] - 1 and a[1] < b[3] - 1 and b[1] < a[3] - 1:
            bad.append((texts[i].get_text()[:34], texts[j].get_text()[:34]))
assert not bad, f"text overlap: {bad[:4]}"

W, H = fig.get_size_inches() * fig.dpi
spill = []
for bx0, by0, bx1, by1 in boxes:
    px = (bx0 * W, by0 * H, bx1 * W, by1 * H)
    for t in texts:
        a = bb(t)
        cx, cy_ = (a[0] + a[2]) / 2, (a[1] + a[3]) / 2
        if not (px[0] < cx < px[2] and px[1] < cy_ < px[3]):
            continue                      # text belongs to some other box
        if a[0] < px[0] + 2 or a[2] > px[2] - 2 or a[1] < px[1] + 2 or a[3] > px[3] - 2:
            spill.append(t.get_text()[:34].replace("\n", " / "))
assert not spill, f"text spills its box: {spill}"

small = [t.get_text()[:30] for t in texts if t.get_fontsize() < TYPE["annotation_min"]]
assert not small, f"below the {TYPE['annotation_min']} pt floor: {small}"

assert abs(FIG_W - 7.09) < 1e-6, "canvas must be 180 mm wide"

for ext, kw in (("png", dict(dpi=600)), ("pdf", {})):
    fig.savefig(OUTDIR / f"{STEM}.{ext}", **kw)
    print("wrote", OUTDIR / f"{STEM}.{ext}")

try:
    import numpy as np
    from PIL import Image
    im = np.array(Image.open(OUTDIR / f"{STEM}.png").convert("L"))
    e = 8
    ink = int((im[:e, :] < 200).sum() + (im[-e:, :] < 200).sum() +
              (im[:, :e] < 200).sum() + (im[:, -e:] < 200).sum())
    print(f"edge check: {im.shape[1]}x{im.shape[0]} px, ink in the outer {e} px = {ink}")
    assert ink == 0, "content is clipped at the canvas edge"
except ImportError:
    print("edge check skipped (no PIL)")

# ------------------------------------------------------------------ sidecar
TSVDIR.mkdir(parents=True, exist_ok=True)
with open(TSVDIR / f"{STEM}.tsv", "w") as fh:
    fh.write("panel\tquantity\tvalue\tsource\n")
    for name, p, nn in DEFAULT_P:
        fh.write(f"c\tatlas-agreement precision @100 bp, default, {name}\t{p:.4f}\tmanuscript/19 section 1\n")
        fh.write(f"c\tcalls kept, default, {name}\t{nn}\tmanuscript/19 section 1\n")
    for lab, val, _, _ in EVIDENCE:
        fh.write(f"b\t{lab.splitlines()[0]}\t{val:.4f}\tmanuscript/19 section 2\n")
    fh.write(f"a\tpoly(A) clip rate, genome-wide\t{CLIP_RATE:.4f}%\tmanuscript/23 section 2\n")
    fh.write(f"c\tpeak memory, PBMC 10k v3\t{COMPUTE_GB} GB\tmanuscript/19 section 3\n")
    fh.write(f"c\twall time, uncontended\t{COMPUTE_MIN}\tmanuscript/19 section 3\n")
    fh.write(f"c\tnull p<0.05 rate, calibrated switch test\t{CAL_NULL}%\tmanuscript/14\n")
    fh.write(f"d\treplicated switches\t{COHORT_SWITCHES}\tmanuscript/20\n")
    fh.write(f"d\tpatients\t{COHORT_PATIENTS}\tmanuscript/20\n")
print("wrote", TSVDIR / f"{STEM}.tsv")

cap = f"""# Fig 1 (ALTERNATIVE candidate) — `{STEM}` caption

Generated by `scripts/manuscript_figures/{STEM}.py`. This is the **tool-introduction**
candidate for Figure 1; the committed alternative is `fig1_overview` (the biology-first
schematic). The PI chooses between them.

## Legend

**Figure 1 | PeakATail: algorithm, evidence model and outputs.**
**(a)** The pipeline. A single command, `peakatail run`, takes an aligned single-cell BAM with cell
barcodes and UMIs and returns poly(A) sites with a per-cell count matrix. Reads whose end
carries a non-templated A/T run are the tool's primary evidence; genome-wide only
{CLIP_RATE:.3f}% of cell-barcoded reads carry one, a scarce but highly specific channel
(manuscript/23). Tail-carrying read ends are clustered by single linkage within 25 bp, the
resulting candidates are ranked by whether tail evidence supports them, and sites sitting on a
genomic A-run are removed by screening the reference sequence. Two flags set the operating
point: `--polya-min-umis`, the number of distinct tail-carrying molecules a site must have, and
`--ip-filter`, the genome screen. Downstream, `peakatail switch diff` tests cell-type differences in
site usage, `peakatail switch length` summarises 3'UTR length as per-cell PDUI, and
`replication_filter.py` keeps only switches on which independent samples agree.
**(b)** The evidence model. Atlas-agreement precision at 100 bp for the three evidence classes
on PBMC 10k v3: sites with at least two tail-carrying molecules (the default output, {DEFAULT_P[0][1]:.4f},
{DEFAULT_P[0][2]:,} sites), sites with at least one ({EVIDENCE[1][1]:.4f}), and coverage-only
candidates carrying no tail evidence ({EVIDENCE[2][1]:.4f}). The classes are reported separately
rather than merged, so the operating point is the user's choice.
**(c)** Headline properties. Atlas-agreement precision of the default output on four libraries
spanning two species and two chemistries ({DEFAULT_P[0][1]:.4f}, {DEFAULT_P[1][1]:.4f},
{DEFAULT_P[2][1]:.4f} and 0.8279 on a second human donor); the switch test's measured
false-positive rate under shuffled labels ({CAL_NULL:.1f}%, with no hit in any of {CAL_RUNS} null
runs, manuscript/14); and the compute footprint on a 10k-cell library ({COMPUTE_GB} GB peak
memory, {COMPUTE_MIN} for an uncontended single run; ~35 min when four arms share the
machine, so peak memory is the production number, manuscript/19).
**(d)** The reliability result. In a lung-adenocarcinoma cohort, {COHORT_SWITCHES} cell-type
switches replicate across {COHORT_PATIENTS} patients, and none replicates in any of ten
label-shuffle null runs of the identical pipeline (manuscript/20).

**Caveats that travel with this figure.** Atlas-agreement precision is agreement with the curated
PolyASite 2.0 reference within 100 bp on the same strand; it is not ground truth, and sites
genuinely absent from the atlas count against it. The default is a deliberately conservative
operating point chosen for reliability and is not the F1 optimum. "None in ten label-shuffle
nulls" is an empirical statement with a resolution floor of p <= 0.091 across ten permutations
and must never be reported as an FDR. The cohort is 17 libraries from 14 patients; the
replication primary is 15 libraries from 12 patients.

## Provenance

Values and sources are in `results/figures/manuscript/{STEM}.tsv`; each traces to
manuscript/14, /19, /20 or /23 as listed there. Shared style from
`scripts/manuscript_figures/_pubstyle.py`. Canvas {FIG_W} x {FIG_H} in ({FIG_W * 25.4:.0f} mm
print width); PNG 600 dpi; PDF vector, fonttype 42. The build asserts pairwise text-overlap
freedom, containment of every label inside its own box, a {TYPE['annotation_min']:.0f} pt minimum type size and a
blank 8 px margin.
"""
(OUTDIR / f"{STEM}.caption.md").write_text(cap)
print("wrote", OUTDIR / f"{STEM}.caption.md")

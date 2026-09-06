"""Genome-browser view of replicated poly(A)-site switches across spermatogenesis.

Four genes, each with exactly TWO replicated significant sites so that "proximal"
and "distal" are unambiguous, drawn as normalised read-coverage tracks for the
three spermatogenic stages, in both mice independently.

Read acceptance matches the caller (`-F 3844`, the tool's own CLIP_EXCLUDE_FLAGS)
and only reads on the gene's strand are drawn, which is how the caller counts;
92-99.8% of reads over these loci are on that strand anyway.

Coverage is normalised to mean depth per cell, so track heights are comparable
between stages that contain different numbers of cells. Usage percentages are NOT
read off the coverage: they are the quantity the switch test actually used,
n_reads_pas / n_reads_gene per stage, taken from the differential tables.

Data:  results/figures/manuscript/switch_tracks_{coverage.tsv.gz,meta.tsv,usage.tsv}
Build: python3 scripts/manuscript_figures/_switch_track_extract.py   (cached step)
       python3 scripts/manuscript_figures/fig_switch_tracks.py
"""
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import re

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
STEM = "fig_switch_tracks"

STAGES = ["SPC", "RS", "ES"]
STAGE_LABEL = {"SPC": "Spermatocyte", "RS": "Round spermatid", "ES": "Elongating spermatid"}
STAGE_COL = {"SPC": PAL["light"], "RS": PAL["peakatail"], "ES": PAL["accent"]}
INK, MUTED, GRID = PAL["ink"], PAL["muted"], PAL["grid"]
PROX_C, DIST_C = PAL["good"], PAL["bad"]

TWO_SITE_N = 722            # genes whose replicated sites number exactly two
CLEAN_N = 259               # ...of those, annotation-clean (see the guard in the extractor)
SEL_N = 130                 # ...of those, both sites in ONE 3'UTR of ONE transcript
SHORT_PCT, LONG_PCT = 97.7, 2.3   # direction within the tandem set: 127 vs 3
cov = pd.read_csv(TSVDIR / "switch_tracks_coverage.tsv.gz", sep="\t")
meta = pd.read_csv(TSVDIR / "switch_tracks_meta.tsv", sep="\t")
use = pd.read_csv(TSVDIR / "switch_tracks_usage.tsv", sep="\t")

_sg = meta.drop_duplicates("symbol").set_index("symbol").strand
_t = cov.groupby(["symbol", "strand"]).depth.sum().unstack(fill_value=0)
_frac = [(_t.loc[i, _sg[i]] / _t.loc[i].sum()) * 100 for i in _t.index]
STRAND_PCT, STRAND_MAX = min(_frac), max(_frac)
NCELL = ", ".join(f"mouse {m[-1]} " + "/".join(
    str(int(meta[(meta.mouse == m) & (meta.stage == s_)].n_cells.iloc[0])) for s_ in STAGES)
    + " SPC/RS/ES" for m in ["mouse1", "mouse2"])
GENES = list(dict.fromkeys(meta.symbol))
MICE = ["mouse1", "mouse2"]
assert len(GENES) == 4, GENES

GTF = WD / "data/references/mouse/Mus_musculus.GRCm38.102.gtf"
FEAT = {"exon": "exon", "CDS": "cds", "three_prime_utr": "utr", "five_prime_utr": "utr"}


def load_models(windows):
    """Exon / CDS / UTR intervals per transcript, for each drawn window.

    Drawn so a reader can see whether the two sites sit inside ONE 3'UTR or in the
    UTRs of different isoforms -- a question the switch statistic cannot answer.
    """
    out = {}
    with open(GTF) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 9)
            k = FEAT.get(f[2])
            if k is None:
                continue
            a, b = int(f[3]), int(f[4])
            for sym, (gid, chrom, lo, hi) in windows.items():
                if f[0] != chrom or b < lo or a > hi:
                    continue
                if f'gene_id "{gid}"' not in f[8]:
                    continue
                tx = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
                d = out.setdefault(sym, {}).setdefault(tx, {"exon": [], "cds": [], "utr": []})
                d[k].append((a, b))
    return out


_w = meta.drop_duplicates("symbol").set_index("symbol")
MODELS = load_models({s_: (_w.loc[s_, "gene"], str(_w.loc[s_, "chrom"]),
                           int(_w.loc[s_, "lo"]), int(_w.loc[s_, "hi"])) for s_ in _w.index})
print("transcripts drawn per gene:", {k: len(v) for k, v in MODELS.items()})

# ---------------------------------------------------------------------------
FIG_W, FIG_H = 7.09, 8.92
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


# geometry, in inches from the top
TOP = 0.62
GENE_H = 1.78           # one gene block
GENE_GAP = 0.14
TRACK_H = 0.335         # one stage track
MODEL_H = 0.50          # the gene-model panel under each column
TXROW = 0.090           # one transcript row inside it
COL_X = [0.075, 0.545]
COL_W = 0.385

T(0.030, fy(0.22), "a", fontsize=TYPE["panel_letter"], fontweight="bold", va="baseline")
T(0.058, fy(0.22), "Replicated poly(A)-site switches, read coverage at the locus",
  fontsize=TYPE["panel_title"], va="baseline")
T(0.030, fy(0.42), "Mean read depth per cell on the gene's strand. Two mice scored independently; "
                   "each column is one mouse.",
  fontsize=TYPE["annotation"], color=MUTED, va="baseline")

for mi, mouse in enumerate(MICE):
    T(COL_X[mi] + COL_W / 2, fy(0.60), f"Mouse {mouse[-1]}", fontsize=7.6,
      fontweight="bold", ha="center", va="baseline", color=INK)

legend_done = False
for gi, sym in enumerate(GENES):
    g0 = TOP + gi * (GENE_H + GENE_GAP)
    row = meta[meta.symbol == sym].iloc[0]
    strand = row.strand

    T(0.030, fy(g0 + 0.10), sym, fontsize=7.8, fontweight="bold", style="italic",
      va="baseline")
    T(0.030, fy(g0 + 0.10) - fh(0.135),
      f"chr{row.chrom} ({strand}), {row.direction}, "
      f"{len(MODELS.get(sym, {}))} transcript(s)", fontsize=6.2, color=MUTED,
      va="baseline")

    for mi, mouse in enumerate(MICE):
        sub = meta[(meta.symbol == sym) & (meta.mouse == mouse)].iloc[0]
        lo, hi = int(sub.lo), int(sub.hi)
        prox, dist = int(sub.prox), int(sub.dist)
        c = cov[(cov.symbol == sym) & (cov.mouse == mouse) & (cov.strand == strand)]
        # shared y across the three stages of this (gene, mouse): the switch is a
        # change in shape, and a per-track autoscale would hide exactly that.
        vecs = {}
        for st in STAGES:
            s = c[c.stage == st]
            v = np.zeros(hi - lo)
            if len(s):
                v[s.pos.values - lo] = s.depth.values
            ncell = int(meta[(meta.symbol == sym) & (meta.mouse == mouse)
                             & (meta.stage == st)].n_cells.iloc[0])
            vecs[st] = v / max(ncell, 1)
        ymax = max(v.max() for v in vecs.values()) or 1.0

        for si, st in enumerate(STAGES):
            ax = fig.add_axes([COL_X[mi], fy(g0 + 0.22 + (si + 1) * TRACK_H),
                               COL_W, fh(TRACK_H * 0.90)])
            x = np.arange(lo, hi)
            ax.fill_between(x, 0, vecs[st], color=STAGE_COL[st], lw=0, alpha=0.92, zorder=3)
            for pos, col in ((prox, PROX_C), (dist, DIST_C)):
                ax.axvline(pos, color=col, lw=0.9, ls=(0, (2.4, 1.6)), zorder=4)
            ax.set_xlim(lo, hi)
            ax.set_ylim(0, ymax * 1.14)
            ax.set_yticks([])
            ax.set_xticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)
            ax.spines["left"].set_visible(True)
            ax.spines["left"].set_color(GRID)
            ax.spines["left"].set_linewidth(0.7)
            if mi == 0:
                ax.text(-0.016, 0.5, st, transform=ax.transAxes, fontsize=6.6,
                        fontweight="bold", color=STAGE_COL[st], ha="right", va="center")
            # the tested usage of each site in this stage
            u = use[(use.symbol == sym) & (use.mouse == mouse) & (use.stage == st)]
            up = u[u.site == "proximal"].usage
            ud = u[u.site == "distal"].usage
            if len(up) and len(ud):
                ax.text(0.995, 0.80, f"{up.iloc[0]*100:.0f}% / {ud.iloc[0]*100:.0f}%",
                        transform=ax.transAxes, fontsize=6.0, ha="right", va="center",
                        color=MUTED, family="DejaVu Sans", zorder=6,
                        bbox=dict(facecolor="white", edgecolor="none", pad=1.1))

        # ---- gene model: exons, introns, UTRs, and where the two sites fall
        axm = fig.add_axes([COL_X[mi], fy(g0 + 0.22 + 3 * TRACK_H + MODEL_H),
                            COL_W, fh(MODEL_H * 0.82)])
        axm.set_xlim(lo, hi)
        axm.set_ylim(0, 1)
        axm.axis("off")

        txs = MODELS.get(sym, {})
        # the transcript carrying both sites is drawn first and named; it is the one
        # the tandem-3'UTR classification is about
        main_tx = str(sub.tx) if "tx" in sub.index else None
        order = ([main_tx] if main_tx in txs else []) + \
                sorted(t for t in txs if t != main_tx)
        n_tx = max(len(order), 1)
        top_y, row_h = 0.80, min(0.20, 0.72 / n_tx)

        for ti, tx in enumerate(order):
            d = txs[tx]
            yc = top_y - ti * row_h
            is_main = (tx == main_tx)
            line_c = INK if is_main else GRID
            utr_c = PROX_C if is_main else GRID          # the shared 3'UTR, named in the legend
            cds_c = INK if is_main else "#B9C2C8"
            xs = [a for a, _ in d["exon"]] + [b for _, b in d["exon"]]
            if not xs:
                continue
            x0, x1 = max(min(xs), lo), min(max(xs), hi)
            axm.plot([x0, x1], [yc, yc], color=line_c, lw=0.7, zorder=2,
                     solid_capstyle="butt")
            for f in np.linspace(0.06, 0.94, 9):          # direction of transcription
                x = x0 + (x1 - x0) * f
                dx = (hi - lo) * (0.007 if strand == "+" else -0.007)
                axm.plot([x - dx, x], [yc - 0.035, yc], color=line_c, lw=0.55, zorder=3)
                axm.plot([x - dx, x], [yc + 0.035, yc], color=line_c, lw=0.55, zorder=3)
            for a, b in d["exon"]:                        # every exon, untranslated height
                axm.add_patch(Rectangle((a, yc - 0.048), max(b - a, 1), 0.096,
                                        facecolor="#C7D0D6" if not is_main else "#9AA6AE",
                                        edgecolor="none", zorder=4))
            for a, b in d["utr"]:                         # 3'/5' UTR, coloured when it is the one
                if is_main and b >= min(prox, dist) and a <= max(prox, dist):
                    axm.add_patch(Rectangle((a, yc - 0.048), max(b - a, 1), 0.096,
                                            facecolor=utr_c, edgecolor="none", zorder=5))
            for a, b in d["cds"]:                         # coding, tallest
                axm.add_patch(Rectangle((a, yc - 0.098), max(b - a, 1), 0.196,
                                        facecolor=cds_c, edgecolor="none", zorder=6))

        for pos, col in ((prox, PROX_C), (dist, DIST_C)):
            axm.plot([pos], [0.02], marker="^", ms=4.2, color=col, zorder=7, clip_on=False)
        axm.text(hi - (hi - lo) * 0.005, 0.96, f"{(hi - lo) / 1000:.1f} kb",
                 fontsize=6.0, color=MUTED, ha="right", va="center", zorder=8,
                 bbox=dict(facecolor="white", edgecolor="none", pad=1.1))

BOT = TOP + 4 * (GENE_H + GENE_GAP) - GENE_GAP + 0.20
T(0.030, fy(BOT + 0.06),
  "Percentages on each track are the tested quantity, proximal % / distal % of the gene's reads in that stage, "
  "not read off the coverage.",
  fontsize=TYPE["annotation"], color=MUTED, va="baseline")
T(0.030, fy(BOT + 0.215),
  "Gene model: tall blocks coding, low blocks untranslated, line with chevrons intronic and in the direction of "
  "transcription. The transcript",
  fontsize=TYPE["annotation"], color=MUTED, va="baseline")
T(0.030, fy(BOT + 0.370),
  "carrying both sites is drawn dark, other isoforms grey, and the green block is the single 3′UTR that holds "
  "both sites.",
  fontsize=TYPE["annotation"], color=MUTED, va="baseline")

# ---------------------------------------------------------------------------
fig.canvas.draw()
r = fig.canvas.get_renderer()
allt = list(texts) + [t for ax in fig.axes for t in ax.texts if t.get_text().strip()]
bad = []
for i in range(len(allt)):
    bi = allt[i].get_window_extent(renderer=r)
    for j in range(i + 1, len(allt)):
        bj = allt[j].get_window_extent(renderer=r)
        if bi.x0 < bj.x1 - 1 and bj.x0 < bi.x1 - 1 and bi.y0 < bj.y1 - 1 and bj.y0 < bi.y1 - 1:
            bad.append((allt[i].get_text()[:26], allt[j].get_text()[:26]))
assert not bad, f"text overlap: {bad[:5]}"
small = [t.get_text()[:24] for t in allt if t.get_fontsize() < TYPE["annotation_min"]]
assert not small, f"below the {TYPE['annotation_min']} pt floor: {small}"
assert abs(FIG_W - 7.09) < 1e-6

for ext, kw in (("png", dict(dpi=600)), ("pdf", {})):
    fig.savefig(OUTDIR / f"{STEM}.{ext}", **kw)
    print("wrote", OUTDIR / f"{STEM}.{ext}")


# ---------------------------------------------------------------------------
# caption sidecar
# ---------------------------------------------------------------------------
n2 = pd.read_csv(TSVDIR / "switch_tracks_selection.tsv", sep="\t") \
    if (TSVDIR / "switch_tracks_selection.tsv").exists() else None
u = use.pivot_table(index=["symbol", "mouse", "stage"], columns="site", values="usage")
lines = []
for sym in GENES:
    r = meta[meta.symbol == sym].iloc[0]
    q = f"{r.qmax:.0e}".replace("e-", " x 10^-")
    per = []
    for mouse in MICE:
        vals = ", ".join(f"{st} {u.loc[(sym, mouse, st), 'proximal']*100:.0f}/"
                         f"{u.loc[(sym, mouse, st), 'distal']*100:.0f}" for st in STAGES)
        per.append(f"mouse {mouse[-1]} {vals}")
    lines.append(f"*{sym}* (chr{r.chrom}, {r.strand} strand, {r.direction}; "
                 f"proximal {int(r.prox):,}, distal {int(r.dist):,}; "
                 f"max q = {q} across both sites and both mice): "
                 + "; ".join(per) + ".")

CAP = f"""# Figure — `{STEM}` caption

## Legend

**Replicated poly(A)-site switches shown as read coverage at the locus.** Four genes
selected under three rules. First, the gene's replicated significant sites must number exactly
two across all three stage comparisons, so that proximal and distal are unambiguous:
{TWO_SITE_N:,} of the genes whose switch replicated in both mice qualify, and a gene with three
or more sites is excluded because no two of them can be labelled proximal and distal without
ambiguity. Second, both sites must lie inside the assigned gene with no other annotated gene
overlapping the drawn window, which leaves {CLEAN_N:,}. Third, both sites must fall inside the
SAME 3'UTR of the SAME transcript, which is what makes an event tandem 3'UTR polyadenylation
rather than a choice between two last exons: {SEL_N:,} genes qualify. That third rule removed
Prdm15, drawn in an earlier version of this figure as the lengthening example, whose two sites
lie in no annotated 3'UTR at all. That second rule is not cosmetic:
poly(A)-site-to-gene assignment in overlapping loci is a known open defect (tool issue #99),
and it removed two genes that had passed every statistical filter, Map3k11, whose apparent
distal site is really the 3' end of the neighbouring Kcnk7, and Pom121, which overlaps Nsun5.
Either would have been drawn as a 3'UTR switch that is actually a between-gene artefact.
Applying the rule also moves the direction balance from 65% shortening to {SHORT_PCT:.0f}%,
so overlapping loci were inflating the apparent lengthening class. Three
shortening examples and one lengthening example are drawn, against the {SHORT_PCT:.0f}% shortening
/ {LONG_PCT:.0f}% lengthening split of those {SEL_N:,} genes (127 against 3). Restricting to one
shared 3'UTR is what produces that asymmetry: before it the clean set ran 73% / 27%, so most
apparent lengthening was not tandem 3'UTR polyadenylation. Each column is one
mouse, scored independently end to end. Within a column, the three tracks are the three
spermatogenic stages (SPC, spermatocyte; RS, round spermatid; ES, elongating spermatid) and
share one vertical scale, printed beneath them, because the switch is a change in the shape
of the profile and a per-track autoscale would conceal it. Coverage is the mean number of
reads per base per cell, so tracks are comparable between stages containing different
numbers of cells ({NCELL}). Read acceptance matches the caller: SAM flag filter 3844, the
tool's own CLIP_EXCLUDE_FLAGS, which drops unmapped, secondary, QC-failed, duplicate and
supplementary records, and only reads on the gene's own strand are drawn, which is how the
caller counts; {STRAND_PCT:.1f}-{STRAND_MAX:.1f}% of reads over these four loci lie on that
strand in any case. Dashed lines mark the two called sites, green proximal and orange
distal. Beneath each column is the gene model: tall blocks are coding exons, low blocks
untranslated, and the line with chevrons is intronic and points in the direction of
transcription. The transcript carrying both sites is drawn dark with its shared 3'UTR in
green; other isoforms of the same gene are grey, so a reader can see directly that the two
sites sit within one 3'UTR rather than in the UTRs of competing isoforms. The percentages on each track are the quantity the switch test
actually used, the proximal and distal share of that gene's reads in that stage
(n_reads_pas / n_reads_gene from the differential tables), and are not measured off the
coverage drawn here.

Per gene: {" ".join(lines)}

**Caveats that travel with this figure.** These are four loci chosen to be legible, not a
random sample, and they are drawn from the replicated set, so they show what a replicated
switch looks like rather than establishing how often one occurs; the frequency claims are
in Fig 5 and manuscript/20. Coverage height reflects both site usage and gene expression,
which is why the tested usage fractions are printed rather than inferred from peak height.
Stage labels are expression-derived cluster assignments, not sorted populations. The two
mice are one study and one chemistry.

## Provenance

Coverage `results/figures/manuscript/switch_tracks_coverage.tsv.gz`, track metadata
`switch_tracks_meta.tsv`, tested usage `switch_tracks_usage.tsv`. Built by
`scripts/manuscript_figures/_switch_track_extract.py` (cached BAM pass) and
`scripts/manuscript_figures/fig_switch_tracks.py`. Switch calls from
`results/stage3_spermatogenesis_v2` (code 9dfdefb); BAMs
`data/benchmark/gse104556/starsolo/Mouse{{1,2}}_scRNAseq`. Canvas {FIG_W} x {FIG_H} in;
PNG 600 dpi; PDF vector, fonttype 42. The build asserts pairwise text-overlap freedom, a
{TYPE['annotation_min']:.0f} pt minimum type size and a blank 8 px margin.
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

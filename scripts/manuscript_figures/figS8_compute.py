#!/usr/bin/env python3
"""
figS8_compute.py -- manuscript Fig S8 ("figS8_compute"): the compute detail
behind Fig 3d. Fig 3d is the summary (PBMC v1 -> v2); this figure carries the
per-arm numbers, the cohort- and mouse-scale production confirmations, the
competitor runtimes WITH their retry/skip caveats, and where the memory
actually went.

SOURCES OF TRUTH (FIGURES_MANIFEST.md row S8)
  * results/figures/manuscript/fig3_tradeoff_compute.tsv -- the verified
    compute registry (every PBMC arm + competitor, each row carrying its
    runtime_mem.txt / benchmark_headtohead.tsv provenance and its caveat note
    VERBATIM). Read programmatically; nothing retyped.
  * manuscript/19_final_gate_v2.md section 3 (mouse arms, uncontended run,
    "quote peak RSS as the production number and disclose the concurrency
    when quoting wall time") and section 5 (17-library cohort: 9:06:26 ->
    1:06:26 at identical 505,197-site output, 8.2x).
  * manuscript/15_final_gate.md section 5 (+ addendum): the v1 limitation as
    stated (heaviest tool in the panel, ~1.4 cores despite --threads 16, a
    10k-cell PBMC BAM needs a >=300 GB node), the 10 section 5 stop signal
    (RSS > 150 GB -> profile first) that was exceeded, and the profiling
    result -- the peak was NOT peak calling but four dense float64 copies of
    the cells x PAS matrix in _tfidf_signac_method1 (291 GB computed, 287.8 GB
    measured in isolation), fixed by sparse TF-IDF + per-(contig, strand)
    parallel peak calling, byte-identical outputs (18/18 files).

RULES THAT BIND EVERY PANEL (21 section 7 MUST-NOT-CLAIM 5)
  "Never quote wall time without the concurrency disclosure; peak RSS is the
  production number." The v1/v2 benchmark-arm wall times were measured with
  all four arms running concurrently and are labelled so; the uncontended
  27:43 single run is drawn separately. Competitor caveats travel verbatim:
  scAPAtrap's wall is the successful RESUMED run and UNDERSTATES a
  from-scratch run (12:58:50 total machine time); scUTRquant's is two summed
  attempts; SCAPTURE's is two mandatory stages with the annotation prebuild
  untimed.

OUTPUTS
  manuscript/figures/figS8_compute.{png,pdf}   (PNG 600 dpi, PDF fonttype 42)
  manuscript/figures/figS8_compute.caption.md  ('## Legend' = the journal
      legend; the concurrency disclosures and the †/‡ marker explanations
      moved there in the 2026-09-02 submission pass -- the per-bar markers
      themselves stay on the image; '## Provenance' = record-keeping)
  results/figures/manuscript/figS8_compute.tsv          (panels a/b, per arm)
  results/figures/manuscript/figS8_compute_confirm.tsv  (panel c)
  results/figures/manuscript/figS8_compute_memory.tsv   (panel d)

Run:  export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/figS8_compute.py
"""
import os
from pathlib import Path

os.environ.setdefault("LC_ALL", "C")
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import re
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _pubstyle import PAL, TYPE, apply_rc, sentence_case   # shared publication style

apply_rc()   # DESIGN_DIRECTIVES.md item 1: one type scale across every figure
ANN = TYPE["annotation_min"]   # 6 pt floor for on-figure annotation


def sc_title(s):
    """Sentence-case a panel title (directive 6) while keeping the lower-case
    panel letter that prefixes it: 'a   the ...' -> 'a   The ...'."""
    m = re.match(r"^([a-z])(\s+)(.*)$", s, flags=re.S)
    return m.group(1) + m.group(2) + sentence_case(m.group(3)) if m else sentence_case(s)

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
OUTDIR = WD / "results/figures/manuscript"
FIGDIR = WD / "manuscript/figures"
NAME = "figS8_compute"
OUTDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)

INK, MUTED, GRID = "#1B2429", "#5A6B73", "#D8E0E3"
COLOR = {  # Okabe-Ito, same tool assignment as Fig 2/3
    "PeakATail":  "#0072B2",
    "polyApipe":  "#D55E00",
    "SCAPTURE":   "#009E73",
    "Sierra":     "#CC79A7",
    "scAPAtrap":  "#E69F00",
    "scUTRquant": "#56B4E9",
}

# ---------------------------------------------------------------------------
# the verified compute registry (Fig 3d's source TSV)
# ---------------------------------------------------------------------------
REG = OUTDIR / "fig3_tradeoff_compute.tsv"
reg = pd.read_csv(REG, sep="\t")


# Fig 3's compute registry renamed its current-version row 'PeakATail v2' -> 'PeakATail'
# in the 2026-09-03 main-figure design pass (DESIGN_DIRECTIVES.md item 4: the mains carry
# no v1->v2 history).  S8 IS the version-history supplement, so it still names v1 and v2:
# resolve the renamed row here.  Same registry row, same numbers -- the EXPECT asserts
# below still pin 34:37.49 / 12.53 GB, so a wrong row cannot pass silently.
REG_ALIAS = {"PeakATail v2": ("PeakATail v2", "PeakATail")}


def rrow(tool):
    for name in REG_ALIAS.get(tool, (tool,)):
        q = reg[reg.tool == name]
        if len(q) == 1:
            return q.iloc[0]
    raise AssertionError((tool, "no unique row in", str(REG)))


v1 = rrow("PeakATail v1")
v2 = rrow("PeakATail v2")
v2u = rrow("PeakATail v2 uncontended")
shp = rrow("PeakATail shipped")
v1ip = rrow("PeakATail v1 IP arm")
v2ip = rrow("PeakATail v2 IP arm")

# EXPECT vs 19 section 3 / 15 section 5 (the verified headline pairs)
assert v1.wall_str == "3:45:53" and round(v1.peak_rss_gb, 1) == 293.7, (v1.wall_str, v1.peak_rss_gb)
assert v2.wall_str == "34:37.49" and round(v2.peak_rss_gb, 2) == 12.53, (v2.wall_str, v2.peak_rss_gb)
assert v2u.wall_str == "27:43", v2u.wall_str
assert v1ip.wall_str == "3:44:37" and round(v1ip.peak_rss_gb, 1) == 239.5
assert v2ip.wall_str == "32:59.38" and round(v2ip.peak_rss_gb, 2) == 11.12
assert round(shp.peak_rss_gb, 1) == 105.5

# ---------------------------------------------------------------------------
# panel a/b arm list: (display label, registry row, bar style tag, disclosure)
# ---------------------------------------------------------------------------
ARMS = [
    ("SCAPTURE", rrow("SCAPTURE"), "competitor",
     "two mandatory stages summed (PAScall 8:57:58 + PASquant 3:15:44); annotation prebuild not timed"),
    ("scAPAtrap †", rrow("scAPAtrap"), "competitor",
     "wall = the successful RESUMED run (4:04:29) after an OOM-killed first attempt skipped three stages; "
     "UNDERSTATES a from-scratch run (12:58:50 total machine time)"),
    ("PeakATail v1 (no-IP arm)", v1, "v1", "single arm, --threads 16"),
    ("PeakATail v1 IP arm", v1ip, "v1", "single arm, --threads 16 (the v1 default call set)"),
    ("polyApipe", rrow("polyApipe"), "competitor", "single clean run, exit 0"),
    ("PeakATail shipped (coverage-only)", shp, "shipped", "pre-clip-seeding caller, v1-era benchmark run"),
    ("Sierra", rrow("Sierra"), "competitor", "regtools 18:53 + FindPeaks/CountPeaks 2:25:04, one clean run"),
    ("PeakATail v2 (no-IP arm) *", v2, "v2", "CONCURRENT: measured with all four arms running at once"),
    ("PeakATail v2 IP arm *", v2ip, "v2", "CONCURRENT: measured with all four arms running at once"),
    ("scUTRquant ‡", rrow("scUTRquant"), "competitor",
     "25:52.06 snakemake run + 4:53.13 salvage rerun of 4 jobs after a missing-nbformat exit 1; both summed"),
    ("PeakATail v2, uncontended", v2u, "v2u", "single uncontended run (19 section 3; peak RSS not re-measured, v2 value shown)"),
]


def fmt_wall(row):
    if isinstance(row.wall_str, str) and row.wall_str:
        s = row.wall_str.split(".")[0]
        return s if ":" in s else f"{float(row.wall_h):.2f} h"
    h = float(row.wall_h)
    hh = int(h)
    mm = int(round((h - hh) * 60))
    return f"{hh}:{mm:02d} h"


def bar_face(tag, tool_color):
    # v1 arms pale, v2 solid, shipped hatched, uncontended open
    if tag == "v1":
        return dict(color=tool_color, alpha=0.40, edgecolor=tool_color, linewidth=0.8)
    if tag == "shipped":
        return dict(color="white", edgecolor=tool_color, hatch="////", linewidth=0.9)
    if tag == "v2u":
        return dict(color="white", edgecolor=tool_color, linewidth=1.2)
    return dict(color=tool_color, edgecolor=tool_color, linewidth=0.8)


# canvas: the pre-submission 8.7 x 9.6 in canvas carried a 12-line footer
# caption; it now lives in the caption sidecar ('## Legend') and the canvas
# height drops by the freed space.  Width stays 8.7 in (221 mm): panel d's
# multi-line y labels deliberately occupy the inter-column gap, interleaving
# with panel c's bar annotations -- narrowing collides them (legibility
# outranks the 180 mm width target).
fig = plt.figure(figsize=(8.7, 7.6))
gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 0.85],
                      left=0.235, right=0.965, top=0.945, bottom=0.070,
                      hspace=0.42, wspace=0.26)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0])
axD = fig.add_subplot(gs[1, 1])


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED, length=2.5, labelsize=TYPE["tick"])
    ax.set_axisbelow(True)


tsvAB = []
ys = np.arange(len(ARMS))[::-1]
for (lab, row, tag, disclosure), y in zip(ARMS, ys):
    tool = row.tool.split()[0] if row.tool.startswith("PeakATail") else row.tool
    c = COLOR.get(tool, "#999999")
    face = bar_face(tag, c)
    axA.barh(y, float(row.wall_h), height=0.62, zorder=3, **face)
    axA.text(float(row.wall_h) + 0.15, y, fmt_wall(row), va="center", fontsize=ANN, color=INK, zorder=4)
    axB.barh(y, float(row.peak_rss_gb), height=0.62, zorder=3, **face)
    axB.text(float(row.peak_rss_gb) * 1.18, y, f"{row.peak_rss_gb:.1f} GB"
             + (" (v2 value)" if tag == "v2u" else ""), va="center", fontsize=ANN, color=INK, zorder=4,
             bbox=dict(facecolor="white", alpha=0.92, edgecolor="none", pad=0.5))
    tsvAB.append(dict(panel="a,b", arm=lab.replace(" *", "").replace(" †", "").replace(" ‡", ""),
                      tool=tool, wall_h=float(row.wall_h), wall_str=row.wall_str if isinstance(row.wall_str, str) else "",
                      peak_rss_gb=float(row.peak_rss_gb), style_tag=tag,
                      disclosure=disclosure, registry_note_verbatim=row.note, source=row.source))

axA.set_yticks(ys)
axA.set_yticklabels([sentence_case(a[0]) for a in ARMS], fontsize=TYPE["tick"])
axA.set_ylim(-0.7, len(ARMS) + 0.5)   # headroom shrunk with the moved notes
axA.set_xlim(0, 13.8)
axA.set_xlabel(sentence_case("wall time on the PBMC 10k v3 BAM (hours)"), fontsize=TYPE["axis_label"])
# the concurrency disclosure and the †/‡ explanations moved to the legend
# (2026-09-02 submission pass); the per-bar *, †, ‡ markers stay on the image
axA.set_title(sc_title("a   wall time per arm"),
              loc="left", fontweight="bold", fontsize=TYPE["panel_title"])
style(axA)
axA.grid(True, axis="x", color=GRID, lw=0.5, alpha=0.7)

axB.set_yticks(ys)
axB.set_yticklabels([])
axB.set_ylim(-0.7, len(ARMS) + 0.5)   # headroom shrunk with the moved notes
axB.set_xscale("log")
axB.set_xlim(2, 900)
axB.set_xticks([3, 10, 30, 100, 300])
axB.set_xticklabels(["3", "10", "30", "100", "300"])
axB.set_xlabel(sentence_case("peak RSS (GB, log scale) — the production number"),
               fontsize=TYPE["axis_label"])
axB.set_title(sc_title("b   peak memory per arm:\n293.7 GB → 12.53 GB (23.5×)"),
              loc="left", fontweight="bold", fontsize=TYPE["panel_title"])
axB.axvline(150, color=INK, lw=0.8, ls=(0, (4, 2)), zorder=2)
# short line label; "(10 §5) — exceeded by v1, then profiled" moved to the legend
# label beside the reference line, not straddling it (design pass 2026-09-03)
axB.text(150 * 1.07, len(ARMS) - 0.1, "150 GB stop signal",
         fontsize=ANN, color=INK, ha="left", va="bottom", linespacing=1.25)
style(axB)
axB.grid(True, axis="x", color=GRID, lw=0.5, alpha=0.7)
RSS_FOLD = v1.peak_rss_gb / v2.peak_rss_gb
assert round(RSS_FOLD, 1) == 23.5, RSS_FOLD  # arithmetic on the two verified registry values

# ---------------------------------------------------------------------------
# panel c: production confirmations at other scales (19 sections 3/5, 15 section 5)
# ---------------------------------------------------------------------------
def hms(h, m, s):
    return h + m / 60 + s / 3600


confirm = [
    # label, v1 wall h, v2 wall h, v1 str, v2 str, v1 rss, v2 rss, threads note, source
    ("17-library Laughney cohort\n(~224 GB of BAM, one unified PAS space)",
     hms(9, 6, 26), hms(1, 6, 26), "9:06:26", "1:06:26", 13.63, 13.53,
     "identical YAML/flags; v2 --peak-workers 16; unified PAS 505,197 = 505,197 (identical output scale)",
     "manuscript/19_final_gate_v2.md section 5"),
    ("testis mouse 1", hms(1, 1, 9), hms(0, 14, 10), "1:01:09", "14:10", 23.1, 3.07,
     "v1 --threads 12 (15 section 5); v2 concurrent with the other three arms (19 section 3)",
     "manuscript/15_final_gate.md section 5; manuscript/19_final_gate_v2.md section 3"),
    ("testis mouse 2", hms(1, 1, 48), hms(0, 16, 0), "1:01:48", "16:00", 22.4, 3.65,
     "v1 --threads 12 (15 section 5); v2 concurrent with the other three arms (19 section 3)",
     "manuscript/15_final_gate.md section 5; manuscript/19_final_gate_v2.md section 3"),
]
COHORT_FOLD = confirm[0][1] / confirm[0][2]
assert round(COHORT_FOLD, 1) == 8.2, COHORT_FOLD   # 19 section 5's own "8.2x"

yc = np.arange(len(confirm))[::-1] * 1.0
bh = 0.30
for (lab, w1, w2, s1, s2, r1, r2, note, src), y in zip(confirm, yc):
    axC.barh(y + bh / 2 + 0.015, w1, height=bh, color=COLOR["PeakATail"], alpha=0.40,
             edgecolor=COLOR["PeakATail"], linewidth=0.8, zorder=3)
    axC.barh(y - bh / 2 - 0.015, w2, height=bh, color=COLOR["PeakATail"],
             edgecolor=COLOR["PeakATail"], linewidth=0.8, zorder=3)
    axC.text(w1 + 0.12, y + bh / 2 + 0.015, f"{s1}   ({r1:.2f} GB)", va="center", fontsize=ANN, color=INK)
    axC.text(w2 + 0.12, y - bh / 2 - 0.015, f"{s2}   ({r2:.2f} GB)", va="center", fontsize=ANN, color=INK)
axC.text(confirm[0][1] / 2, yc[0] + bh + 0.30, "8.2× — identical output:\n505,197 unified PAS on both codes",
         ha="center", va="bottom", fontsize=ANN, color=INK, linespacing=1.25)
axC.set_yticks(yc)
axC.set_yticklabels([sentence_case(c0[0]) for c0 in confirm], fontsize=TYPE["tick"])
axC.set_xlim(0, 10.6)
axC.set_ylim(yc[-1] - 0.75, yc[0] + 1.30)
# commit hashes moved to the legend sidecar's provenance; the key stays
axC.set_xlabel(sentence_case("wall time (hours; light = v1, solid = v2)"), fontsize=TYPE["axis_label"])
# the '(wall, with peak RSS)' qualifier moved to the Legend sidecar (directive 2, no-loss:
# the '(c) ... fell 9:06:26 -> 1:06:26 (8.2x) ... with peak RSS ...' sentence carries it)
axC.set_title(sc_title("c   the #97 fix confirmed at cohort\nand mouse scale"),
              loc="left", fontweight="bold", fontsize=TYPE["panel_title"])
style(axC)
axC.grid(True, axis="x", color=GRID, lw=0.5, alpha=0.7)

# ---------------------------------------------------------------------------
# panel d: where the memory went (15 section 5 addendum, profiling verified FIXED)
# ---------------------------------------------------------------------------
mem = [
    ("v1 peak RSS, full PBMC run", float(v1.peak_rss_gb), COLOR["PeakATail"], 0.40,
     "results/figures/manuscript/fig3_tradeoff_compute.tsv (runtime_mem.txt)"),
    ("the culprit in isolation:\n4 dense float64 copies of the\ncells × PAS matrix (clustering,\nnot peak calling)", 287.8,
     "#D55E00", 1.0, "manuscript/15_final_gate.md section 5 addendum (profiling, results/perf/)"),
    ("after perf/clip-memory:\nsparse TF-IDF + parallel calling\n(byte-identical outputs, 18/18)", 12.45,
     COLOR["PeakATail"], 1.0, "manuscript/15_final_gate.md section 5 addendum (uncontended PR run)"),
]
yd = np.arange(len(mem))[::-1] * 1.0
for (lab, v, c, a, src), y in zip(mem, yd):
    axD.barh(y, v, height=0.55, color=c, alpha=a, edgecolor=c, linewidth=0.8, zorder=3)
    axD.text(v + 8, y, f"{v:.1f} GB", va="center", fontsize=ANN, color=INK)
axD.set_yticks(yd)
axD.set_yticklabels([sentence_case(m[0]) for m in mem], fontsize=ANN)
axD.set_xlim(0, 390)
axD.set_ylim(yd[-1] - 0.7, yd[0] + 0.8)
# two-line axis label: the single line overran the right canvas edge at the shared
# type scale (design pass 2026-09-03)
axD.set_xlabel(sentence_case("peak RSS (GB)\n291 GB computed = 4 × 23,303 × 390,493 × 8 B"),
               fontsize=TYPE["axis_label"])
# the '(15 s5 addendum)' citation moved to the Legend/Provenance sidecar (directive 2)
axD.set_title(sc_title("d   where the memory went — profiling,\nnot peak calling"),
              loc="left", fontweight="bold", fontsize=TYPE["panel_title"])
style(axD)
axD.grid(True, axis="x", color=GRID, lw=0.5, alpha=0.7)

# ---------------------------------------------------------------------------
# the journal legend (single source: the caption sidecar's '## Legend'; the
# former 12-line on-figure footer plus the moved marker explanations,
# substance verbatim)
# ---------------------------------------------------------------------------
caption = (
    "Figure S8 | Compute detail behind Fig 3d. (a, b) Every PBMC 10k v3 arm from the verified compute registry "
    "(fig3_tradeoff_compute.tsv; each bar's provenance and caveat travel in figS8_compute.tsv verbatim). "
    "Rule in force (21 §7 item 5): never quote wall time without the concurrency disclosure — the v2 arm "
    "wall times (34:37 no-IP, 32:59 IP; marked * on the figure) were measured with all four benchmark arms "
    "running CONCURRENTLY; the "
    "uncontended single run is 27:43 — and peak RSS is the production number: 293.7 GB (v1) → 12.53 GB (v2), "
    "23.5×. v1 exceeded the pre-registered 150 GB stop signal (10 §5; dashed line in b) and was profiled "
    "rather than shipped; "
    "on v1 a 10k-cell PBMC BAM needed a ≥300 GB node and used ~1.4 cores despite --threads 16 (15 §5). "
    "Competitor caveats carried verbatim: scAPAtrap's (†) 4:04:29 is the successful RESUMED run after an "
    "OOM-killed first attempt and UNDERSTATES a from-scratch run (12:58:50 total machine time, disclosed, "
    "not summed); scUTRquant's (‡) 30:45 sums two attempts (salvage rerun after exit 1); SCAPTURE's 12:14 sums "
    "its two mandatory stages with the annotation prebuild untimed; polyApipe and Sierra ran clean. "
    "(c) The #97 performance fix confirmed at two other scales on the same box: the 17-library Laughney "
    "cohort (~224 GB of BAM, identical YAML and flags) fell 9:06:26 → 1:06:26 (8.2×) at IDENTICAL output — "
    "505,197 unified PAS on both codes — and the mouse arms fell 1:01:09/1:01:48 → 14:10/16:00 with peak "
    "RSS 23.1/22.4 → 3.07/3.65 GB (v1 --threads 12; v2 concurrent). (d) Why: profiling located the v1 peak "
    "in clustering (_tfidf_signac_method1), not peak calling — four dense float64 copies of the cells × PAS "
    "matrix (291 GB computed, 287.8 GB measured in isolation); the fix (sparse TF-IDF, per-(contig, strand) "
    "parallel calling, streamed CB filter) is byte-identical on all 18/18 compared outputs."
)

for ext, kw in (("png", dict(dpi=600)), ("pdf", {})):
    p = FIGDIR / f"{NAME}.{ext}"
    fig.savefig(p, **kw)
    print("wrote", p)
p = OUTDIR / f"{NAME}.png"
fig.savefig(p, dpi=600); print("wrote", p)

# ---------------------------------------------------------------------------
# audit TSVs
# ---------------------------------------------------------------------------
p = OUTDIR / f"{NAME}.tsv"
pd.DataFrame(tsvAB).to_csv(p, sep="\t", index=False, float_format="%.6f"); print("wrote", p)

conf_rows = []
for lab, w1, w2, s1, s2, r1, r2, note, src in confirm:
    conf_rows.append(dict(panel="c", scale=lab.replace("\n", " "),
                          v1_wall_h=w1, v1_wall_str=s1, v1_peak_rss_gb=r1,
                          v2_wall_h=w2, v2_wall_str=s2, v2_peak_rss_gb=r2,
                          speedup=w1 / w2, note=note, source=src))
p = OUTDIR / f"{NAME}_confirm.tsv"
pd.DataFrame(conf_rows).to_csv(p, sep="\t", index=False, float_format="%.6f"); print("wrote", p)

mem_rows = [dict(panel="d", quantity=m[0].replace("\n", " "), peak_rss_gb=m[1], source=m[4]) for m in mem]
mem_rows.append(dict(panel="d", quantity="dense-copy arithmetic (computed): 4 x 23,303 x 390,493 x 8 B",
                     peak_rss_gb=291.0, source="manuscript/15_final_gate.md section 5 addendum"))
p = OUTDIR / f"{NAME}_memory.tsv"
pd.DataFrame(mem_rows).to_csv(p, sep="\t", index=False, float_format="%.6f"); print("wrote", p)

# ---------------------------------------------------------------------------
# sidecar caption
# ---------------------------------------------------------------------------
cav = reg[["tool", "kind", "note", "source"]].copy()
cav_lines = "\n".join(f"- **{r.tool}** ({r.kind}): {r.note} — `{r.source}`" for r in cav.itertuples())
cap_md = f"""# Fig S8 — `figS8_compute` caption (generated by `scripts/manuscript_figures/figS8_compute.py`; submission pass 2026-09-02: the on-figure footer, the concurrency-disclosure prose and the †/‡ explanations moved into the Legend below — the per-bar *, †, ‡ markers stay on the image; image at 600 dpi)

## Legend

{caption}
Values and verbatim caveats: `results/figures/manuscript/figS8_compute*.tsv`.

**Binding disclosure rules (legend content).**
- 21 §7 MUST-NOT-CLAIM 5, verbatim: "Never quote wall time without the concurrency disclosure; peak RSS is
  the production number." The starred arms in panel a are the four-arm concurrent measurements of
  19 §3; 27:43 is the uncontended single run from the #97 PR (its peak RSS was not re-measured, so the v2
  concurrent RSS is shown for that bar and labelled).
- The v1 compute row is a *limitation that was disclosed and then fixed*: 15 §5 states it exceeded the
  10 §5 stop signal (RSS > 150 GB → profile first) and was the heaviest tool in the panel; the addendum
  records the profile and the fix, verified byte-identical (18/18 files incl. pasbed.bed, pas_support.tsv,
  matrices, clusters).
- Competitor rows are the tools' own verified runs on the same box (v1-era benchmark); their caveat notes
  travel verbatim from the registry (below) and none is hidden or summed away.

## Provenance

**Registry rows carried verbatim (`fig3_tradeoff_compute.tsv` note column):**
{cav_lines}

**Panel c** numbers: 19 §5 (cohort table: wall 9:06:26 vs 1:06:26, peak RSS 13.63 vs 13.53 GB, unified PAS
505,197 = 505,197, "same output scale at 1/8 the wall time") and 19 §3 / 15 §5 (mouse arms; v1 at
--threads 12, v2 measured with all four arms concurrent); light bars = v1 code 4efeb125, solid = v2 code
9dfdefb. **Panel d** numbers: 15 §5 addendum
(4 × 23,303 × 390,493 × 8 B = 291 GB computed; 287.8 GB measured in isolation; PBMC 12.45 GB / 27m43s
uncontended after the fix; mouse1 3.68 GB / 9m03s).

Sources: every plotted value in `results/figures/manuscript/figS8_compute*.tsv` with a source column;
panels a/b re-read `fig3_tradeoff_compute.tsv` (verified) rather than retyping it.
"""
cap_md += """
**Design pass 2026-09-03** (`manuscript/figures/DESIGN_DIRECTIVES.md`, supplement light pass). Type comes from
the shared style module `scripts/manuscript_figures/_pubstyle.py` (`apply_rc()`), and every on-figure annotation
sits at or above the 6 pt floor. Axis labels, panel titles and prose tick labels are sentence-cased through
`_pubstyle.sentence_case()`, canonical identifiers preserved and the lower-case panel letters kept. Two title
qualifiers left the image for the Legend above: panel c's *(wall, with peak RSS)* — no-loss, the \"(c) ... fell
9:06:26 → 1:06:26 (8.2×) ... with peak RSS 23.1/22.4 → 3.07/3.65 GB\" sentence carries it — and panel d's
*(15 §5 addendum)* citation, which is the panel-d source line in this Provenance. Overlaps fixed: panel b's
150 GB stop-signal label now sits beside its reference line instead of straddling it, panel b's bar values
carry an opaque backing so that line passes behind them, and panel d's axis label is set on two lines (one
line overran the canvas edge). Panel d and the confirm TSV regenerate byte-identical; every plotted value in
`figS8_compute.tsv` is unchanged — only its `registry_note_verbatim` column moved, because it travels verbatim
from `fig3_tradeoff_compute.tsv`, which the main-figure pass re-noted on 2026-09-02. That same pass renamed the
registry's current-version row `PeakATail v2` → `PeakATail` (directive 4: no v1→v2 history in the mains); this
script resolves the renamed row through `REG_ALIAS`, and the unchanged EXPECT asserts (34:37.49 / 12.53 GB)
still pin it.
"""
p = FIGDIR / f"{NAME}.caption.md"; p.write_text(cap_md); print("wrote", p)

# ---------------------------------------------------------------------------
# edge check: outer 8 px must be blank
# ---------------------------------------------------------------------------
from PIL import Image
im = np.asarray(Image.open(FIGDIR / f"{NAME}.png").convert("L"))
edge = 8
border = np.concatenate([im[:edge, :].ravel(), im[-edge:, :].ravel(),
                         im[:, :edge].ravel(), im[:, -edge:].ravel()])
ink = int((border < 250).sum())
print(f"edge check: {im.shape[1]}x{im.shape[0]} px, ink pixels in the outer {edge} px = {ink}")
assert ink == 0, "INK IN THE OUTER 8 PX -- something is clipped"

print(f"figS8_compute: PBMC {v1.wall_str} / {v1.peak_rss_gb:.1f} GB -> {v2.wall_str} / "
      f"{v2.peak_rss_gb:.2f} GB (RSS {RSS_FOLD:.1f}x); cohort 9:06:26 -> 1:06:26 ({COHORT_FOLD:.1f}x) "
      f"at identical 505,197 unified PAS")

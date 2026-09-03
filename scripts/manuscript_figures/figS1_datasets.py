#!/usr/bin/env python3
"""
figS1_datasets.py -- manuscript Fig S1 ("figS1_datasets"): the dataset table
made visual, including the per-BAM poly(A) clip rate. Replaces the retired
`cohort_qc` (21 D5/section 4.3; its caveats transfer to this figure's scope --
none of its lambda_gradient-era cohort numbers are reused here).

SOURCES OF TRUTH (FIGURES_MANIFEST.md row S1)
  * v2 run_config.json files of the four benchmark arms (seq-len, CB length,
    identical parameter blocks) -- read programmatically below.
  * the Stage-3 cohort run manifest results/stage3_laughney_v3/run_manifest.json
    (code 9dfdefb3, 17 GSM datasets = 14 patients) -- read programmatically.
  * manuscript/04_datasets.md (accessions, chemistry, cell counts).
  * manuscript/10_caller_fix_plan.md section R2 (PBMC = CellRanger 3.0.0 /
    10x v3 / 91 bp R2; testis = STARsolo / 10x v2 / 98 bp R2).
  * manuscript/26_second_donor_preregistration.md (pbmc4k library facts,
    verifier verdict FIXED; head-sample clip rates 2.2565% / 3.6965%).
  * manuscript/23_algorithm_roadmap.md section 2/[V] section 5: the CORRECTED
    genome-wide qualifying poly(A) clip rate 0.5730% (3,195,067 / 557,564,408
    accepted CB reads). NEVER 1.152% -- that constant is superseded
    (05_figure_index.md standing correction).
  * manuscript/19_final_gate_v2.md section 5 (~224 GB cohort BAM, 505,197
    unified PAS) and manuscript/20_stage3_replication.md (29,063 curated ->
    18,651 label-confirmed cells; the cohort's own 0.015% clip-rate warning is
    INVALID -- head-of-chr1 sampling -- and is never repeated as fact).

NUMBER POLICY -- a cell with no verified source is printed as an em-dash and
the reason is recorded in the caption + audit TSV; nothing is approximated.
Dropped cells: testis clip rate (10 section R2 required the measurement; no
verified value exists in the manuscript set), cohort clip rate (the only
in-run estimate is the invalid head-sample warning, 20), cohort read length /
upstream pipeline version (not recorded in the verified set).

OUTPUTS
  manuscript/figures/figS1_datasets.{png,pdf}  (PNG 600 dpi, PDF fonttype 42)
  manuscript/figures/figS1_datasets.caption.md ('## Legend' + '## Provenance')
  results/figures/manuscript/figS1_datasets.tsv          (every printed cell)
  results/figures/manuscript/figS1_datasets_cliprate.tsv (panel b values)

Run:  export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/figS1_datasets.py
"""
import json
import os
import sys
import textwrap
from pathlib import Path

os.environ.setdefault("LC_ALL", "C")
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _pubstyle import PAL, TYPE, apply_rc, sentence_case   # shared publication style

apply_rc()   # DESIGN_DIRECTIVES.md item 1: one type scale across every figure
ANN = TYPE["annotation_min"]   # 6 pt floor for on-figure annotation


def sc_title(s):
    """Sentence-case a panel title (directive 6) while keeping the lower-case
    panel letter that prefixes it: 'a   the ...' -> 'a   The ...'."""
    import re
    m = re.match(r"^([a-z])(\s+)(.*)$", s, flags=re.S)
    return m.group(1) + m.group(2) + sentence_case(m.group(3)) if m else sentence_case(s)

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
OUTDIR = WD / "results/figures/manuscript"
FIGDIR = WD / "manuscript/figures"
NAME = "figS1_datasets"
OUTDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)

INK, MUTED, GRID = "#1B2429", "#5A6B73", "#D8E0E3"
BLUE, GREEN, ORANGE, GREY = "#0072B2", "#009E73", "#D55E00", "#999999"

# ---------------------------------------------------------------------------
# machine-readable inputs: the four v2 run_config.json + the cohort manifest
# ---------------------------------------------------------------------------
RUN_CONFIGS = {
    "pbmc_10k_v3": WD / "results/benchmark_tools/pbmc_10k_v3/peakatail_clipseeded_final_v2_ipfilt/run/run_config.json",
    "pbmc4k":      WD / "results/benchmark_tools/pbmc4k_donor2/run_ipfilt/run/run_config.json",
    "testis_m1":   WD / "results/benchmark_tools/gse104556/peakatail_clipseeded_final_v2/mouse1/run/run_config.json",
    "testis_m2":   WD / "results/benchmark_tools/gse104556/peakatail_clipseeded_final_v2/mouse2/run/run_config.json",
}
seqlen, cblen = {}, {}
for k, p in RUN_CONFIGS.items():
    assert p.exists(), p
    v = json.loads(p.read_text())["variables"]
    seqlen[k], cblen[k] = int(v["seqlen"]), int(v["cb_len"])

# EXPECT: run_config seq-lens equal the verified read lengths (26 section 3;
# 10 section R2: PBMC 91 bp R2, testis 98 bp R2; pbmc4k 98 bp).
assert seqlen == {"pbmc_10k_v3": 91, "pbmc4k": 98, "testis_m1": 98, "testis_m2": 98}, seqlen
assert set(cblen.values()) == {16}, cblen

COHORT_MANIFEST = WD / "results/stage3_laughney_v3/run_manifest.json"
cm = json.loads(COHORT_MANIFEST.read_text())
assert cm["tool"]["commit_short"] == "9dfdefb", cm["tool"]["commit_short"]
assert "17 GSMs = 14 patients" in cm["pre_registered_design"]["unit_of_replication"]

# ---------------------------------------------------------------------------
# clip rates (the evidence-channel column) -- all verified, definitions named
# ---------------------------------------------------------------------------
CLIP_GW_NUM, CLIP_GW_DEN = 3_195_067, 557_564_408      # 23 [V] section 5
CLIP_GW_PBMC = 100.0 * CLIP_GW_NUM / CLIP_GW_DEN       # 0.5730 % genome-wide qualifying
assert round(CLIP_GW_PBMC, 4) == 0.5730, CLIP_GW_PBMC  # EXPECT vs 23 section 2
CLIP_SUPERSEDED = 1.152        # the retired docstring constant -- never quote as fact
CLIP_HEAD_PBMC = 2.2565        # head-sample QC estimator, PBMC (26 A1.2, verifier-exact)
CLIP_HEAD_PBMC4K = 3.6965      # head-sample QC estimator, pbmc4k (26 R4.4, verifier-exact)
CLIP_GW_ESTIMATOR_TRUTH = 0.5364  # what the estimator's own criteria measure genome-wide (23 section 3 Step 0)
CLIP_WARN = 0.3                # caller warn threshold (26 R4.4)
assert CLIP_HEAD_PBMC / CLIP_GW_ESTIMATOR_TRUTH > 4.0  # the ~4.2x head bias (05 standing correction)

# ---------------------------------------------------------------------------
# the table -- every cell carries its source; None = verified-source missing
# ---------------------------------------------------------------------------
COLS = ["Library", "Species / accession", "Chemistry", "Upstream pipeline",
        "Cells", "Reads / input", "R2 length\n(--seq-len)", "Poly(A) clip rate\n(% of CB reads)"]
DASH = "—"


def cell(text, source, note=""):
    return dict(text=text, source=source, note=note)


rows = [
    ("pbmc_10k_v3", [
        cell("PBMC 10k v3\n(human donor 1)", "manuscript/04_datasets.md (b)"),
        cell("human\n10x public pbmc_10k_v3", "manuscript/04_datasets.md (b)"),
        cell("10x 3' v3", "10 sR2; 26 s2"),
        cell("CellRanger 3.0.0\n(GRCh38-3.0.0 refdata)", "10 sR2; 26 s2/A1"),
        cell("~10,000 (est.)", "manuscript/04_datasets.md; 26 sR4.2"),
        cell("557,564,408\naccepted CB reads", "23 [V] s5 (clip-rate denominator)"),
        cell(f"{seqlen['pbmc_10k_v3']} bp", "run_config.json; 26 s3 (100,000/100,000 first records)"),
        cell("0.5730 genome-wide qualifying\n(3,195,067 / 557,564,408)\nhead-sample estimator: 2.2565",
             "23 s2/[V] s5; 26 A1.2", "corrects the superseded 1.152%"),
    ]),
    ("pbmc4k", [
        cell("pbmc4k\n(human donor 2)", "manuscript/26 s2"),
        cell("human\n10x public pbmc4k\n(flow cell H53GNBCXY)", "26 s2/V"),
        cell("10x 3' v2", "26 s2 (web summary)"),
        cell("CellRanger 2.1.0\n(GRCh38-1.2.0 refdata)", "26 s2 (web summary + BAM @PG)"),
        cell("4,340 (est.)", "26 s2 (web summary)"),
        cell("87,433 mean reads/cell;\nBAM 33.24 GB", "26 s2"),
        cell(f"{seqlen['pbmc4k']} bp", "run_config.json; 26 s2 (200,000/200,000 first records)"),
        cell("head-sample estimator: 3.6965\n(genome-wide not measured)", "26 sR4.4/V",
             "1.64x donor 1's estimator value; both far above the 0.3% warn threshold"),
    ]),
    ("testis_m1", [
        cell("testis mouse 1\n(GSE104556)", "manuscript/04_datasets.md (b)"),
        cell("mouse\nGEO GSE104556", "manuscript/04_datasets.md (b)"),
        cell("10x 3' v2", "10 sR2"),
        cell("STARsolo\n(version not recorded)", "10 sR2 / 10 s5 gate table"),
        cell("2,500\n(both mice together)", "manuscript/04_datasets.md (b)"),
        cell(DASH, "no verified per-BAM read count in the manuscript set"),
        cell(f"{seqlen['testis_m1']} bp", "run_config.json; 10 sR2; 26 s3"),
        cell(DASH, "10 sR2 required the measurement; no verified value exists"),
    ]),
    ("testis_m2", [
        cell("testis mouse 2\n(GSE104556)", "manuscript/04_datasets.md (b)"),
        cell("mouse\nGEO GSE104556", "manuscript/04_datasets.md (b)"),
        cell("10x 3' v2", "10 sR2"),
        cell("STARsolo\n(version not recorded)", "10 sR2 / 10 s5 gate table"),
        cell("(see mouse 1)", "manuscript/04_datasets.md (b): 2,500 across 2 mice"),
        cell(DASH, "no verified per-BAM read count in the manuscript set"),
        cell(f"{seqlen['testis_m2']} bp", "run_config.json; 10 sR2; 26 s3"),
        cell(DASH, "10 sR2 required the measurement; no verified value exists"),
    ]),
    ("laughney", [
        cell("Laughney LUAD\ncohort (17 libs,\ncohort-level row)", "20; cohort run_manifest.json"),
        cell("human\nGEO GSE123904\n17 GSMs = 14 patients", "13 addendum; 20; cohort run_manifest.json"),
        cell("10x 3' (v2 era;\nHiSeq 2500)", "manuscript/04_datasets.md (c)"),
        cell(DASH, "upstream pipeline/version not recorded in the verified set"),
        cell("29,063 curated;\n18,651 label-\nconfirmed (64.2%)", "20 (verified SOUND)"),
        cell("~224 GB of BAM;\n505,197 unified PAS", "19 s5"),
        cell(DASH, "not recorded in the verified set"),
        cell(DASH, "only in-run estimate is the INVALID head-of-chr1 0.015% warning (20) -- never repeated as fact"),
    ]),
]

# ---------------------------------------------------------------------------
# figure: panel a = the table, panel b = the clip-rate evidence channel
# (two separate gridspecs: the table wants the full width, the bar panel needs
#  a real left margin for its y labels — the single-gridspec version clipped
#  them off the figure edge and failed its own 8-px check)
# ---------------------------------------------------------------------------
# Publication layout: the on-figure caption paragraph lives in the sidecar's
# Legend now, so the canvas keeps only the two panels (height 7.6 -> 6.2 in).
# Width stays 8.9 in: the 8-column table's columns are tuned to it, and a
# 7.7-in trial collided the header cells -- legibility outranks the 180 mm
# width target for this landscape supplement.
fig = plt.figure(figsize=(8.9, 6.2))
gsT = fig.add_gridspec(1, 1, left=0.015, right=0.985, top=0.945, bottom=0.411)
gsB = fig.add_gridspec(1, 1, left=0.130, right=0.985, top=0.330, bottom=0.079)
axT = fig.add_subplot(gsT[0])
axB = fig.add_subplot(gsB[0])

axT.set_xlim(0, 1); axT.set_ylim(0, 1); axT.axis("off")
axT.set_title(sc_title("a   the libraries every headline number runs on"), loc="left",
              fontweight="bold", fontsize=TYPE["panel_title"])

# column x-positions (left edges, axes fraction)
colx = [0.000, 0.125, 0.255, 0.345, 0.480, 0.600, 0.745, 0.815]
colw = [0.125, 0.130, 0.090, 0.135, 0.120, 0.145, 0.070, 0.185]
n_show_rows = len(rows)
header_h = 0.115
row_h = (1.0 - header_h) / n_show_rows

for j, cname in enumerate(COLS):
    axT.text(colx[j] + 0.004, 1.0 - header_h / 2, cname, ha="left", va="center",
             fontsize=ANN + 0.4, color="white", fontweight="bold", linespacing=1.15)
axT.add_patch(matplotlib.patches.Rectangle((0, 1.0 - header_h), 1.0, header_h,
              facecolor="#3D5A6C", edgecolor="none", zorder=0))

tsv_rows = []
for i, (lib, cells) in enumerate(rows):
    y1 = 1.0 - header_h - i * row_h
    y0 = y1 - row_h
    if i % 2 == 0:
        axT.add_patch(matplotlib.patches.Rectangle((0, y0), 1.0, row_h,
                      facecolor="#F0F4F6", edgecolor="none", zorder=0))
    axT.axhline(y0, color=GRID, lw=0.5, zorder=1)
    for j, c in enumerate(cells):
        is_dash = c["text"] == DASH
        axT.text(colx[j] + 0.004, y0 + row_h / 2, c["text"], ha="left", va="center",
                 fontsize=ANN if j else ANN + 0.4, color=MUTED if is_dash else INK,
                 fontweight="bold" if j == 0 else "normal", linespacing=1.2, zorder=2)
        tsv_rows.append(dict(panel="a", library=lib, column=COLS[j].replace("\n", " "),
                             value=c["text"].replace("\n", " "), source=c["source"],
                             note=c["note"]))
axT.axhline(1.0 - header_h, color="#3D5A6C", lw=0.8, zorder=1)
axT.axhline(1.0 - header_h - n_show_rows * row_h, color="#3D5A6C", lw=0.8, zorder=1)

# ---- panel b: clip rate, the channel the method runs on ----
bars = [
    ("PBMC 10k v3\ngenome-wide", CLIP_GW_PBMC, BLUE, "",
     "23 s2/[V] s5: 3,195,067 / 557,564,408 accepted CB reads"),
    ("PBMC 10k v3\nhead-sample est.", CLIP_HEAD_PBMC, BLUE, "////",
     "26 A1.2 (verifier-exact); first 200,000 CB reads = head of chr1"),
    ("pbmc4k\nhead-sample est.", CLIP_HEAD_PBMC4K, GREEN, "////",
     "26 R4.4/V (verifier-exact); genome-wide value not measured"),
]
ys = np.arange(len(bars))[::-1] * 0.9
for (lab, v, c, hatch, src), y in zip(bars, ys):
    axB.barh(y, v, height=0.55, color="white" if hatch else c, edgecolor=c,
             linewidth=1.0, hatch=hatch, zorder=3)
    axB.text(v + 0.05, y, f"{v:.4f}%", va="center", fontsize=ANN + 0.4, color=INK, zorder=4)
    tsv_rows.append(dict(panel="b", library=lab.replace("\n", " "), column="clip_rate_pct",
                         value=f"{v:.4f}", source=src,
                         note="hatched = head-sample estimator (known ~4.2x head bias)" if hatch else
                              "genome-wide qualifying rate (the figure constant)"))
axB.set_yticks(ys)
axB.set_yticklabels([sentence_case(b[0]) for b in bars], fontsize=TYPE["tick"])
axB.axvline(CLIP_WARN, color=INK, lw=0.8, ls=(0, (4, 2)), zorder=2)
axB.text(CLIP_WARN, ys[0] + 0.62, "0.3% caller warn threshold (26)", fontsize=ANN,
         color=INK, ha="left", va="bottom")
axB.axvline(CLIP_SUPERSEDED, color=ORANGE, lw=0.8, ls=(0, (1, 2)), zorder=2)
# label sits in the free space right of the retired-constant line on the TOP bar
# row (the genome-wide bar ends at 0.57, so nothing is drawn there) -- the old
# placement ran the text straight through the pbmc4k bar
axB.text(CLIP_SUPERSEDED + 0.06, ys[0], "1.152% -- SUPERSEDED docstring constant,\ndo not quote (05 standing correction)",
         fontsize=ANN, color=ORANGE, ha="left", va="center")
axB.set_xlim(0, 4.4)
axB.set_ylim(ys[-1] - 0.75, ys[0] + 1.05)
axB.set_xlabel(sentence_case("qualifying poly(A) soft-clip rate (% of CB reads)"),
               fontsize=TYPE["axis_label"])
# panel-level descriptive qualifier '(where a verified value exists)' moved to the
# Legend sidecar (DESIGN_DIRECTIVES.md item 2, no-loss: it is stated in '(b) ...' there)
axB.set_title(sc_title("b   the evidence channel: poly(A) clip rate per BAM"),
              loc="left", fontweight="bold", fontsize=TYPE["panel_title"])
for s in ("top", "right"):
    axB.spines[s].set_visible(False)
for s in ("left", "bottom"):
    axB.spines[s].set_color(GRID)
axB.tick_params(colors=MUTED, length=2.5, labelsize=TYPE["tick"])
axB.grid(True, axis="x", color=GRID, lw=0.5, alpha=0.7)
axB.set_axisbelow(True)

# The legend paragraph -- no longer drawn on the image (surgery pass,
# 2026-09-02): it is written to the caption sidecar's '## Legend' below.
caption = (
    "Figure S1 | The datasets behind every headline number, with the poly(A) clip rate -- the "
    "evidence channel the clip-seeded caller runs on. (a) Per-library table; all four benchmark "
    "arms ran the identical v2 parameter block (108 recorded parameters; the only library "
    "parameter is --seq-len = the R2 read length, read here from each arm's run_config.json). "
    "An em-dash means no verified measurement exists in the manuscript set and the value is "
    "omitted rather than approximated: testis clip rates (10 §R2 required the measurement; "
    "none was recorded), testis/cohort per-BAM read counts, and the cohort clip rate -- the "
    "cohort's own 0.015% warning is invalid (head-of-chr1 sampling; 20) and is never repeated "
    "as fact. (b) Clip rates where a verified value exists. The genome-wide qualifying rate on "
    "PBMC 10k v3 is 0.5730% (3,195,067 / 557,564,408 accepted CB reads; 23 §2) -- the "
    "corrected constant, superseding the 1.152% docstring figure. Hatched bars are the caller's "
    "head-sample QC estimator (first 200,000 CB reads of a coordinate-sorted BAM = head of "
    "chr1), which reads 2.2565% on PBMC where its own criteria measure 0.5364% genome-wide "
    "(~4.2× head bias; 23 §3, 05 standing correction) -- so estimator values are "
    "comparable to each other, not to the genome-wide bar. pbmc4k's estimator value is 1.64× "
    "donor 1's: the evidence channel is richer on donor 2, and both are far above the 0.3% warn "
    "threshold (26 §R4.4). Sources per cell: figS1_datasets.tsv."
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
df = pd.DataFrame(tsv_rows)
p = OUTDIR / f"{NAME}.tsv"; df.to_csv(p, sep="\t", index=False); print("wrote", p)
clip = pd.DataFrame([
    dict(quantity="genome_wide_qualifying_clip_rate_pct", library="pbmc_10k_v3", value=round(CLIP_GW_PBMC, 4),
         numerator=CLIP_GW_NUM, denominator=CLIP_GW_DEN, source="manuscript/23_algorithm_roadmap.md s2/[V] s5"),
    dict(quantity="head_sample_estimator_clip_rate_pct", library="pbmc_10k_v3", value=CLIP_HEAD_PBMC,
         numerator="", denominator="first 200,000 CB reads", source="manuscript/26 A1.2 (verifier-exact)"),
    dict(quantity="head_sample_estimator_clip_rate_pct", library="pbmc4k", value=CLIP_HEAD_PBMC4K,
         numerator="", denominator="first 200,000 CB reads", source="manuscript/26 R4.4/V (verifier-exact)"),
    dict(quantity="estimator_criteria_genome_wide_truth_pct", library="pbmc_10k_v3", value=CLIP_GW_ESTIMATOR_TRUTH,
         numerator="", denominator="", source="manuscript/23 s3 Step 0; 05 standing correction (caption only)"),
    dict(quantity="superseded_docstring_constant_pct", library="pbmc_10k_v3", value=CLIP_SUPERSEDED,
         numerator="", denominator="", source="RETIRED -- 05 standing correction; reference line only, never quote"),
    dict(quantity="warn_threshold_pct", library="all", value=CLIP_WARN,
         numerator="", denominator="", source="manuscript/26 R4.4"),
])
p = OUTDIR / f"{NAME}_cliprate.tsv"; clip.to_csv(p, sep="\t", index=False); print("wrote", p)

# ---------------------------------------------------------------------------
# sidecar caption
# ---------------------------------------------------------------------------
cap_md = f"""# Fig S1 — `figS1_datasets` caption (generated by `scripts/manuscript_figures/figS1_datasets.py`)

## Legend

{caption}

**Dropped cells and why (number policy: never approximate).**
- testis mouse 1/2 poly(A) clip rate: `10` §R2 flagged the STARsolo/98-bp arm as unmeasured and required the
  measurement before Stage 1 Phase 2; no verified value was recorded anywhere in the manuscript set.
- testis / cohort per-BAM read counts and cohort R2 length: not recorded in the verified set.
- cohort upstream pipeline version: not recorded in the verified set.
- cohort clip rate: the only in-run estimate is `clip_rate_pct 0.015` from the head-sampling estimator on a
  coordinate-sorted BAM (head of chr1) — `20` (disclosure 2) rules "**Never repeat the 'no clip evidence' /
  '0.015% clip rate' justification as fact**"; the same warning fired on five libraries.

**Cohort row scope.** 17 GSM libraries = 14 patients (`13` addendum; cohort `run_manifest.json`, code
`9dfdefb3`); the pre-registered primary excludes MetBone (15 GSMs = 12 patients, `20`). Cells are the Stage-3
curated set (29,063) and its label-confirmed subset (18,651, 64.2%; `20`, verified SOUND) — not a CellRanger
cell estimate. The retired `cohort_qc` numbers (55,422 cells / 16,500-site 17/17 intersection) come from the
pre-clip-seeded lambda_gradient run and are not reused (05 §R1).

**Clip-rate definitions (do not mix).** Genome-wide qualifying rate (PBMC 0.5730%; `23` §2) is the figure
constant and supersedes 1.152%. The head-sample QC estimator (first 200,000 CB reads; PBMC 2.2565%, pbmc4k
3.6965%, both verifier-exact in `26`) carries a known ~4.2× head bias (its own criteria measure 0.5364%
genome-wide on PBMC; `23` §3 Step 0) — estimator values are compared only with each other.

## Provenance

Sources: every printed cell's source is a column of `results/figures/manuscript/figS1_datasets.tsv`;
panel b: `figS1_datasets_cliprate.tsv`. Machine-read inputs: the four v2 `run_config.json` files and
`results/stage3_laughney_v3/run_manifest.json`. PNG 600 dpi, PDF vector with subsetted TrueType
(fonttype 42, no Type 3). The working-phase render carried the legend paragraph on the image; at the
2026-09-02 surgery pass it moved here, and the image keeps the table, panel titles, axis labels and
the two reference-line labels only.

Design pass 2026-09-03 (`DESIGN_DIRECTIVES.md`, supplement light pass): type now comes from the shared
`scripts/manuscript_figures/_pubstyle.py` scale (`apply_rc()`), every on-figure annotation sits at or above
the 6 pt floor, axis label and panel titles are sentence-cased through `_pubstyle.sentence_case()`
(canonical identifiers preserved), and panel b's descriptive qualifier *(where a verified value exists)*
left the image for the Legend above (no-loss: it is the "(b) Clip rates where a verified value exists"
sentence). Numbers, panels and audit TSVs are unchanged.
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

print(f"figS1_datasets: PBMC genome-wide clip {CLIP_GW_PBMC:.4f}% (never 1.152%); "
      f"head-sample {CLIP_HEAD_PBMC}% vs pbmc4k {CLIP_HEAD_PBMC4K}%; seq-lens {seqlen}")

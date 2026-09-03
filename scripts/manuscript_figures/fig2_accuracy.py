#!/usr/bin/env python3
"""
fig2_accuracy.py -- manuscript Fig 2 ("fig2_accuracy"): where the shipping tool
sits against the competitor panel, and what it costs at a matched call budget.

2026-09-02 PUBLICATION DESIGN PASS (manuscript/figures/DESIGN_DIRECTIVES.md).
What changed, and why:
  * DIRECTIVE 4 -- development history is OFF the main figures.  The "shipped
    (pre-fix)" point, the two no-IP intermediate points and the arrowed
    development path are NOT PLOTTED any more.  Their rows are still read,
    still asserted and still written to the audit TSV (`plotted = False`); the
    story they told lives in Fig S12 (version history) and Fig S8 (compute),
    with one pointer sentence in the Legend.
  * DIRECTIVE 5 -- exactly ONE precision/recall plane in the paper.  Fig 3's
    trade surface (the molecule-support sweep) FOLDED IN HERE: PeakATail is now
    drawn as its current *operating curve* (support >=1 ... >=10) with exactly
    two named, user-choosable points emphasised -- the pre-registered precision
    default (>=2 molecules) and the >=1-molecule sensitivity arm.  This is the
    same mark the old development path occupied, so the plane did not get
    busier; it got honest.  Fig 3 sharpened to pure robustness.
  * The matched-call-count head-to-head, previously buried in the caption, is
    now panels c/d -- it is the honest comparison and it belongs on the page.
  * The tier-decomposition panel is GONE from Fig 2: Fig S2 panel a already is
    the tier composition of every v2 arm, so keeping it here was a pairwise
    redundancy (directive 5).  Its numbers stay in the Legend and in
    `fig2_accuracy_tiers.tsv`, which is still written for the record.
  * F1 isolines, the superseded two-sided gate and the per-call-set null ticks
    left the plotted plane (crowding; all three are recorded in
    `fig2_accuracy_reference_lines.tsv` and quantified in the Legend).  The
    pre-registered gate line stays, with a short label.
  * All type comes from `_pubstyle` (shared palette / markers / type scale);
    every prose label goes through `sentence_case()` (directive 6).

VERSION SWITCH -- env var FINAL_BENCHMARK_VERSION picks which PeakATail run is
plotted.  The competitor / catalog arms are the same TSVs in both versions.
  v2  (default; THE PAPER NUMBERS)  code 9dfdefb = #93 + #96 (IP-filter
      minus-strand fix) + #97 (performance), arms run 2026-08-21 16:14-16:50.
      Sources of truth: manuscript/19_final_gate_v2.md and
      results/benchmark_tools/final_v2_verify/VERIFIED_v2.md (verifier: FIXED).
  v1  (FINAL_BENCHMARK_VERSION=v1; kept selectable for the record)
      code 4efeb125, manuscript/15_final_gate.md (verified SOUND).
  The operating CURVE (the >=1/2/3/5/10 sweep) exists only for the v2 arms, so
  under v1 the two named operating points are drawn without the connecting
  curve and the figure says so.  Panels c/d are the D3 matched-N record
  (manuscript/25), which is a v2 analysis in both modes -- stated in the Legend.

SOURCE OF TRUTH -- every plotted value is READ from a score TSV written by the
single scoring path scripts/benchmark_tools/score_tool.py (100 bp cutoff,
strand-matched point mode, curated PolyASite 2.0 representative sites,
per-dataset detected-gene recall denominator), from the shared molecule-sweep
module `_molsweep.py` (same scorer), or from the verified matched-N tables of
results/perf_gap/D3_head_to_head/ (manuscript/25, verifier FIXED).  Nothing
numeric is typed by hand except the facts that exist only in the gate write-up
(pre-registration timestamps, code commit, compute) -- those are cited to
19 (v2) / 15 (v1) in the Legend.

PANELS
  a  PBMC 10k v3: atlas-agreement precision@100 (y) vs detected-gene recall@100
     (x).  Competitors at their own operating points (one marker/colour per
     tool, `_pubstyle.TOOL_STYLE`); scUTRquant hollow (catalog-based, not
     ranked); scTail absent (not runnable: R1 = 28 bp).  PeakATail as its
     molecule-support operating curve with the precision default (>=2 mol) and
     the sensitivity arm (>=1 mol) named.  Pre-registered gate P >= 0.50.
  b  Same for GSE104556 testis; competitors mean +- range over the two mice,
     PeakATail curve = mean of the two mice, range bars on the two named points.
  c  PBMC 10k v3 at MATCHED CALL BUDGET: atlas-agreement precision@100 as each
     tool is truncated to a common number of calls (D3 rank sweep); the filled
     marker on each line is that tool's own budget.
  d  Same for testis mouse 1.

NAMING RULES (manuscript/01 number policy)
  * "atlas-agreement precision (100 bp)", never bare "precision".
  * recall shown = detected-gene recall R_det; full-atlas recall in the Legend.
  * .../peakatail_clipseeded_final_v2/pas_tier1_ge2mol_noIP_POSTHOC.bed is NOT
    the pre-registered default and is not plotted.

OUTPUTS
  manuscript/figures/fig2_accuracy.{png,pdf}   (PNG 600 dpi, PDF fonttype 42)
  manuscript/figures/fig2_accuracy.caption.md  (sidecar: '## Legend' = the journal
                                                legend, single source of the
                                                caption; '## Provenance')
  results/figures/manuscript/fig2_accuracy.tsv            (every arm; `plotted` flag)
  results/figures/manuscript/fig2_accuracy_sweep.tsv      (the operating curve, a/b)
  results/figures/manuscript/fig2_accuracy_matchedN.tsv   (panels c/d)
  results/figures/manuscript/fig2_accuracy_tiers.tsv      (NOT PLOTTED since
                                                           2026-09-02; kept for
                                                           the record, see Fig S2)
  results/figures/manuscript/fig2_accuracy_reference_lines.tsv (gates, isolines, nulls)

Run:  export LC_ALL=C; python3 scripts/manuscript_figures/fig2_accuracy.py
      export LC_ALL=C; FINAL_BENCHMARK_VERSION=v1 python3 .../fig2_accuracy.py
      (the env-var name FINAL_BENCHMARK_VERSION is kept through the rename: it
       selects a run, not a stem, and the reproduction commands recorded in
       15/19 use it)
"""
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("LC_ALL", "C")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _pubstyle import PAL, TOOL_STYLE, TYPE, apply_rc, sentence_case   # noqa: E402
from _molsweep import molecule_sweep, PREREG_K, SENSITIVITY_K          # noqa: E402

apply_rc()

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
BT = WD / "results/benchmark_tools"
D3 = WD / "results/perf_gap/D3_head_to_head/tsv"
OUTDIR = WD / "results/figures/manuscript"
FIGDIR = WD / "manuscript/figures"
NAME = "fig2_accuracy"
OUTDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)

CUTOFF = 100.0
INK, MUTED, GRID = PAL["ink"], PAL["muted"], PAL["grid"]
HUMAN, MOUSE = "pbmc_10k_v3", "gse104556"
TOOL_ORDER = ["PeakATail", "polyApipe", "scAPAtrap", "SCAPTURE", "Sierra", "scUTRquant"]

# ---------------------------------------------------------------------------
# Which PeakATail run: v2 (default, the paper numbers) or v1 (kept for the record).
# Only the PeakATail arm directories / score-file prefixes and the hand-cited
# gate facts differ; competitor arms are identical.
# ---------------------------------------------------------------------------
VERSION = os.environ.get("FINAL_BENCHMARK_VERSION", "v2").lower()
assert VERSION in ("v1", "v2"), f"FINAL_BENCHMARK_VERSION must be v1 or v2, got {VERSION!r}"
VERSIONS = {
    "v2": dict(
        h_dir="peakatail_clipseeded_final_v2", h_pfx="pbmc_final_v2",
        h_ip_dir="peakatail_clipseeded_final_v2_ipfilt", h_ip_pfx="pbmc_final_v2_ipfilt",
        m_dir="peakatail_clipseeded_final_v2", m_pfx="testis_m{r}_final_v2",
        commit="9dfdefb", commit_note="#96 + #97 merged",
        gate_doc="19", gate_doc_path="manuscript/19_final_gate_v2.md",
        verdict="verifier verdict FIXED",
        verified_table="results/benchmark_tools/final_v2_verify/VERIFIED_v2.md",
        verified_short="final_v2_verify/VERIFIED_v2.md (FIXED)",
        prereg="gate 01:19, original arms 02:59, v2 re-run 16:14, all 2026-08-21",
        compute="12.5 GB / ~28-35 min on PBMC (was 294 GB / 3:46 h)",
        ip_trade="the IP filter buys +5.3 pp atlas-agreement precision for −1.2 pp R_det (mouse 1)",
        prev=dict(label="v1 run (code 4efeb125, 15)", P=(0.7167, 0.7414, 0.7554)),
        expect_h=(46524, 0.7062, 0.1754), expect_m1=(26255, 0.7450), expect_m2=(26526, 0.7572),
        expect_frac_single=72.2,
        # definition nit, 19 section 4: the same ratio computed on the no-IP arm is 72.1 %
        # (1 - 62,110 / 222,955), which is the verifier recount the manuscript text quotes.
        single_mol_note=("this ratio is 1 − n(default) / n(tier-1) **on the IP arm**; the same ratio on the "
                         "no-IP arm is 72.1% (1 − 62,110 / 222,955) — the verifier's recount and the value the "
                         "manuscript text quotes (19 §4), superseding 15 §4's 72.5% (the v1 IP arm)"),
        sweep=True,
    ),
    "v1": dict(
        h_dir="peakatail_clipseeded_final", h_pfx="pbmc_final",
        h_ip_dir="peakatail_clipseeded_final_ipfilt", h_ip_pfx="pbmc_final_ipfilt",
        m_dir="peakatail_clipseeded_final", m_pfx="testis_m{r}_final",
        commit="4efeb125", commit_note="pre-fix, Stage-1d not merged",
        gate_doc="15", gate_doc_path="manuscript/15_final_gate.md",
        verdict="verified SOUND",
        verified_table="manuscript/15_final_gate.md",
        verified_short="15_final_gate.md (SOUND)",
        prereg="gate committed 01:19, arms started 02:59 (2026-08-21)",
        compute="294 GB / 3:46 h on PBMC",
        ip_trade="the IP filter buys +4.9 pp atlas-agreement precision for −1.4 pp R_det (mouse 1)",
        prev=None,
        expect_h=(44394, 0.7167, 0.1707), expect_m1=(25991, 0.7414), expect_m2=(26164, 0.7554),
        expect_frac_single=72.5,
        single_mol_note=("this ratio is 1 − n(default) / n(tier-1) **on the v1 IP arm** (15 §4); the "
                         "same ratio on the no-IP arm is 72.1% (1 − 62,110 / 222,955), the verifier's recount (19 §4)"),
        # the molecule-support sweep was only ever computed on the v2 arms
        sweep=False,
    ),
}[VERSION]
H_DIR, H_PFX = VERSIONS["h_dir"], VERSIONS["h_pfx"]
H_IP_DIR, H_IP_PFX = VERSIONS["h_ip_dir"], VERSIONS["h_ip_pfx"]
M_DIR, M_PFX = VERSIONS["m_dir"], VERSIONS["m_pfx"]
CODE_COMMIT = VERSIONS["commit"]          # gate write-up header
PREREG = VERSIONS["prereg"]               # gate write-up header
GATE_DOC = VERSIONS["gate_doc"]

# ---------------------------------------------------------------------------
# Arm registry: (dataset, tool, call_set label, short label, replicate, tsv, role)
# role: competitor | catalog | path (PeakATail operating point, in path order)
#
# 2026-09-02: role names path0..path4 are kept so the audit TSV keeps its
# vocabulary, but they are no longer a *path* on the figure -- directive 4
# retired the development arrows.  Only path3 (>=1-molecule sensitivity arm)
# and path4 (precision default) are plotted; path0/1/2 and the "extra" outputs
# are read, asserted and written to the TSV with plotted = False.
# ---------------------------------------------------------------------------
H = BT / HUMAN
M = BT / MOUSE
ARMS = [
    # ---- PBMC competitors --------------------------------------------------
    (HUMAN, "scUTRquant", "scUTRquant (catalog)", "scUTRquant*", "", H / "scutrquant/score_scutrquant.tsv", "catalog"),
    (HUMAN, "SCAPTURE",   "SCAPTURE",  "SCAPTURE",  "", H / "scapture/score_scapture.tsv",   "competitor"),
    (HUMAN, "polyApipe",  "polyApipe", "polyApipe", "", H / "polyapipe/score_polyapipe.tsv", "competitor"),
    (HUMAN, "Sierra",     "Sierra",    "Sierra",    "", H / "sierra/score_sierra.tsv",       "competitor"),
    (HUMAN, "scAPAtrap",  "scAPAtrap", "scAPAtrap", "", H / "scapatrap/score_scapatrap.tsv", "competitor"),
    # ---- PBMC PeakATail outputs --------------------------------------------
    # NOT PLOTTED since 2026-09-02 (directive 4: development history); kept for the record.
    (HUMAN, "PeakATail", "PeakATail shipped (pre-fix caller)", "shipped", "",
     H / "peakatail/score_peakatail.tsv", "path0"),
    (HUMAN, "PeakATail", "PeakATail both tiers, no IP filter", "both tiers\n(no IP)", "",
     H / H_DIR / f"score_{H_PFX}.tsv", "path1"),
    (HUMAN, "PeakATail", "PeakATail tier-1 >=1 molecule, no IP filter", "tier-1 >=1 mol\n(no IP)", "",
     H / H_DIR / f"score_{H_PFX}__tier1.tsv", "path2"),
    # PLOTTED: the two named, user-choosable operating points
    (HUMAN, "PeakATail", "PeakATail tier-1 >=1 molecule, IP filter", "tier-1 >=1 mol (IP)", "",
     H / H_IP_DIR / f"score_{H_IP_PFX}__tier1.tsv", "path3"),
    (HUMAN, "PeakATail", "PeakATail precision default (tier-1, IP, >=2 molecules; pre-registered)",
     "precision default\n(tier-1, IP, >=2 mol)", "",
     H / H_IP_DIR / f"score_{H_IP_PFX}__PRESPEC_precision_default.tsv", "path4"),
    # ---- PBMC extra PeakATail outputs (record only; tier TSV) ---------------
    (HUMAN, "PeakATail", "PeakATail both tiers, IP filter", "both tiers (IP)", "",
     H / H_IP_DIR / f"score_{H_IP_PFX}.tsv", "extra"),
    (HUMAN, "PeakATail", "PeakATail tier-2 (coverage-only), IP filter", "tier-2 (IP)", "",
     H / H_IP_DIR / f"score_{H_IP_PFX}__tier2.tsv", "extra"),
    (HUMAN, "PeakATail", "PeakATail tier-2 (coverage-only), no IP filter", "tier-2 (no IP)", "",
     H / H_DIR / f"score_{H_PFX}__tier2.tsv", "extra"),
]
for rep in ("mouse1", "mouse2"):
    r = rep[-1]
    ARMS += [
        (MOUSE, "scUTRquant", "scUTRquant (catalog)", "scUTRquant*", rep, M / f"scutrquant/score_scutrquant_{rep}.tsv", "catalog"),
        (MOUSE, "polyApipe",  "polyApipe", "polyApipe", rep, M / f"polyapipe/score_polyapipe_{rep}.tsv", "competitor"),
        (MOUSE, "Sierra",     "Sierra",    "Sierra",    rep, M / f"sierra/score_sierra_{rep}.tsv", "competitor"),
        (MOUSE, "scAPAtrap",  "scAPAtrap", "scAPAtrap", rep, M / f"scapatrap/{rep}/score_scapatrap_{rep}.tsv", "competitor"),
        (MOUSE, "PeakATail", "PeakATail shipped (pre-fix caller)", "shipped", rep,
         M / f"peakatail/score_peakatail_{rep}.tsv", "path0"),
        (MOUSE, "PeakATail", "PeakATail both tiers, IP filter", "both tiers (IP)", rep,
         M / M_DIR / rep / f"score_{M_PFX.format(r=r)}.tsv", "path1"),
        (MOUSE, "PeakATail", "PeakATail tier-1 >=1 molecule, IP filter", "tier-1 >=1 mol (IP)", rep,
         M / M_DIR / rep / f"score_{M_PFX.format(r=r)}__tier1.tsv", "path3"),
        (MOUSE, "PeakATail", "PeakATail precision default (tier-1, IP, >=2 molecules; pre-registered)",
         "precision default\n(tier-1, IP, >=2 mol)", rep,
         M / M_DIR / rep / f"score_{M_PFX.format(r=r)}__PRESPEC_precision_default.tsv", "path4"),
        (MOUSE, "PeakATail", "PeakATail tier-2 (coverage-only), IP filter", "tier-2 (IP)", rep,
         M / M_DIR / rep / f"score_{M_PFX.format(r=r)}__tier2.tsv", "extra"),
    ]
# SCAPTURE mouse: mouse1 plotted. The mouse-2 SITE-LEVEL run COMPLETED
# (gse104556/scapture/DONE.mouse2.ok: 24,076 points, P@100 0.672, scored in
# mouse2/score_scapture_mouse2.tsv; 15 §3 reports both mice 0.694 / 0.672);
# only the per-cell PASquant step failed, which site-level benchmarking does
# not use.  The panel keeps the single-mouse point; never describe the mouse-2
# run as invalid (21 §11 stale-caption fix, 2026-09-02).
ARMS.append((MOUSE, "SCAPTURE", "SCAPTURE", "SCAPTURE", "mouse1",
             M / "scapture/mouse1/score_scapture_mouse1.tsv", "competitor"))

PLOTTED_ROLES = ("competitor", "catalog", "path3", "path4")


def read_arm(tsv: Path) -> dict:
    """Pull the cutoff-100 values straight out of one score_tool.py TSV."""
    df = pd.read_csv(tsv, sep="\t")
    at = df[df.cutoff_bp == CUTOFF]

    def one(panel, series, reference, col="value", replicate=None):
        q = at[(at.panel == panel) & (at.series == series) & (at.reference == reference)]
        if replicate is not None:
            q = q[q.replicate == replicate]
        assert len(q) == 1, (tsv, panel, series, reference, replicate, len(q))
        return float(q[col].iloc[0])

    nulls = [one("precision", "null_genic", "atlas_full", replicate=s) for s in (1, 2, 3)]
    out = dict(
        n=int(one("precision", "real", "atlas_full", col="n_query")),
        n_matched=int(one("precision", "real", "atlas_full", col="n_matched")),
        P=one("precision", "real", "atlas_full"),
        R_det=one("recall", "real", "atlas_detected"),
        R_det_denominator=int(one("recall", "real", "atlas_detected", col="n_query")),
        R_full=one("recall", "real", "atlas_full"),
        R_full_denominator=int(one("recall", "real", "atlas_full", col="n_query")),
        F1_det=one("f1", "real", "atlas_detected"),
        null_P_seed1=nulls[0], null_P_seed2=nulls[1], null_P_seed3=nulls[2],
        null_P_mean=float(np.mean(nulls)),
    )
    # self-consistency of the scorer's F1 with its own P and R_det
    f1 = 2 * out["P"] * out["R_det"] / (out["P"] + out["R_det"])
    assert abs(f1 - out["F1_det"]) < 1e-5, (tsv, f1, out["F1_det"])
    return out


rows = []
for ds, tool, call_set, short, rep, tsv, role in ARMS:
    assert tsv.exists(), tsv
    d = read_arm(tsv)
    d.update(dataset=ds, tool=tool, call_set=call_set, short=short,
             replicate=rep or "single", role=role,
             source=str(tsv.relative_to(WD)), cutoff_bp=int(CUTOFF))
    rows.append(d)
pts = pd.DataFrame(rows)

# ---- consistency with the verified gate table for this VERSION ------------
def pick(ds, role, rep="single"):
    q = pts[(pts.dataset == ds) & (pts.role == role) & (pts.replicate == rep)]
    assert len(q) == 1
    return q.iloc[0]


_d = pick(HUMAN, "path4");  assert (_d.n, round(_d.P, 4), round(_d.R_det, 4)) == VERSIONS["expect_h"], _d
_m1 = pick(MOUSE, "path4", "mouse1"); assert (_m1.n, round(_m1.P, 4)) == VERSIONS["expect_m1"], _m1
_m2 = pick(MOUSE, "path4", "mouse2"); assert (_m2.n, round(_m2.P, 4)) == VERSIONS["expect_m2"], _m2
assert _d.R_det_denominator == 285136 and _m1.R_det_denominator == 126686
GATE_P = 0.50          # 13 section 1, pre-registered precision-first gate -- PLOTTED
ORIG_P, ORIG_F1 = 0.38, 0.261   # 10 section 5 original two-sided gate -- NOT PLOTTED since 2026-09-02
F1_ISO = [0.1, 0.2, 0.3, 0.4]   # NOT PLOTTED since 2026-09-02 (kept in the reference-line TSV)
for r in pts[pts.role == "path4"].itertuples():
    assert r.P >= GATE_P, r
# SIDECAR-GENERATION CHECK (the value left the plot on 2026-09-02; the assert did not):
# the Legend sentence "no >=1-molecule output clears the original P >= 0.38 floor" is
# guarded here, not by a drawn line.
for r in pts[(pts.dataset == HUMAN) & (pts.role.isin(["path1", "path2", "path3"]))].itertuples():
    assert r.P < ORIG_P, r
# SIDECAR-GENERATION CHECK: the per-call-set genic-shuffle nulls used to be drawn as
# tick marks under each point; they are now a Legend sentence bounded by this assert.
NULL_MAX = float(pts.null_P_mean.max())
assert NULL_MAX < 0.030, NULL_MAX

# ---------------------------------------------------------------------------
# the operating curve (DIRECTIVE 5: Fig 3's trade surface folds in here)
# ---------------------------------------------------------------------------
if VERSIONS["sweep"]:
    sweep = molecule_sweep()
    # the curve's two named points must BE the two plotted arms, to 4 dp
    _sw_h2 = sweep[(sweep.dataset == "pbmc") & (sweep.min_molecules == PREREG_K)].iloc[0]
    _sw_h1 = sweep[(sweep.dataset == "pbmc") & (sweep.min_molecules == SENSITIVITY_K)].iloc[0]
    _h1 = pick(HUMAN, "path3")
    assert (int(_sw_h2.n), round(_sw_h2.atlas_agreement_precision, 4)) == (int(_d.n), round(_d.P, 4))
    assert (int(_sw_h1.n), round(_sw_h1.atlas_agreement_precision, 4)) == (int(_h1.n), round(_h1.P, 4))
else:
    sweep = None
    print("VERSION=v1: the molecule-support sweep exists only for the v2 arms -- "
          "the two named operating points are drawn without the connecting curve")

# ---------------------------------------------------------------------------
# panels c/d: matched call budget (D3 rank sweep, manuscript/25, verifier FIXED)
# ---------------------------------------------------------------------------
MN_COL = {"PeakATail": "PeakATail", "polyApipe": "polyApipe", "scAPAtrap": "scAPAtrap",
          "Sierra": "Sierra", "SCAPTURE": "SCAPTURE*", "scUTRquant": "scUTRquant(catalog)"}
MN_SRC = {HUMAN: D3 / "01d_matched_n_headline_pbmc.tsv",
          MOUSE: D3 / "01e_matched_n_headline_mouse1.tsv"}
mn_rows = []
for ds, src in MN_SRC.items():
    assert src.exists(), src
    t = pd.read_csv(src, sep="\t")
    for tool in TOOL_ORDER:
        c = MN_COL[tool]
        for i in t.index:
            v = t.at[i, c + "_P100"]
            if pd.isna(v):
                continue
            mn_rows.append(dict(panel="c" if ds == HUMAN else "d", dataset=ds, tool=tool,
                                N=int(t.at[i, "N"]),
                                atlas_agreement_precision_100bp=float(v),
                                recall_detected_genes_100bp=float(t.at[i, c + "_Rdet"]),
                                F1_detected_genes_100bp=float(t.at[i, c + "_F1"]),
                                plotted=True, source=str(src.relative_to(WD))))
mn = pd.DataFrame(mn_rows)

# each tool's OWN call budget (the filled marker on its line).  The D3 rank sweep is a
# v2 analysis (manuscript/25) in BOTH render modes, so PeakATail's budget here is the v2
# arm's, not the selected VERSION's -- stated in the Legend.
D3_PA_NATIVE = {HUMAN: 46524, MOUSE: 26255}
NATIVE_N = {}
for ds in (HUMAN, MOUSE):
    for tool in TOOL_ORDER:
        if tool == "PeakATail":
            NATIVE_N[(ds, tool)] = D3_PA_NATIVE[ds]
        else:
            q = pts[(pts.dataset == ds) & (pts.tool == tool)]
            if ds == MOUSE:
                q = q[q.replicate == "mouse1"]
            NATIVE_N[(ds, tool)] = int(q.n.iloc[0]) if len(q) else None
mn["is_native_budget"] = [NATIVE_N[(r.dataset, r.tool)] == r.N for r in mn.itertuples()]

# EXPECT asserts on the plotted matched-N values (manuscript/25 §2, verified)
def _mn(ds, tool, N):
    q = mn[(mn.dataset == ds) & (mn.tool == tool) & (mn.N == N)]
    assert len(q) == 1, (ds, tool, N)
    return round(float(q.atlas_agreement_precision_100bp.iloc[0]), 4)


assert _mn(HUMAN, "PeakATail", 120916) == 0.4036 and _mn(HUMAN, "polyApipe", 120916) == 0.3800
assert _mn(HUMAN, "PeakATail", 46524) == 0.7062 and _mn(HUMAN, "polyApipe", 46524) == 0.6049
assert _mn(HUMAN, "PeakATail", 9456) == 0.9527 and _mn(HUMAN, "polyApipe", 9456) == 0.9671
assert _mn(MOUSE, "PeakATail", 71983) == 0.4488 and _mn(MOUSE, "polyApipe", 71983) == 0.4353
if VERSION == "v2":
    # the D3 native row IS the score-TSV default (one number, two independent pipelines)
    assert D3_PA_NATIVE[HUMAN] == int(_d.n) and D3_PA_NATIVE[MOUSE] == int(_m1.n)
    assert _mn(HUMAN, "PeakATail", int(_d.n)) == round(_d.P, 4)
    assert _mn(MOUSE, "PeakATail", int(_m1.n)) == round(_m1.P, 4)
# SIDECAR-GENERATION CHECK: the recall half of the matched-N claim is a Legend
# sentence, not a plotted series, so its bound is asserted here.
_mnr = mn[(mn.dataset == HUMAN) & (mn.tool.isin(["PeakATail", "polyApipe"]))]
for N in sorted(set(_mnr[_mnr.tool == "polyApipe"].N) & set(_mnr[_mnr.tool == "PeakATail"].N)):
    a = float(_mnr[(_mnr.tool == "PeakATail") & (_mnr.N == N)].recall_detected_genes_100bp.iloc[0])
    b = float(_mnr[(_mnr.tool == "polyApipe") & (_mnr.N == N)].recall_detected_genes_100bp.iloc[0])
    assert a > b, (N, a, b)   # "recall at every matched N tested", 25 §2

# ---------------------------------------------------------------------------
# tier decomposition -- NOT PLOTTED since 2026-09-02.
# DIRECTIVE 5 (pairwise redundancy): Fig S2 panel a already shows the tier
# composition of every final v2 arm, so the Fig 2 panel duplicated a supplement.
# The table is still computed, still asserted and still written to
# fig2_accuracy_tiers.tsv; its numbers travel in the Legend.
# ---------------------------------------------------------------------------
t2 = pts[(pts.dataset == HUMAN) & (pts.short == "tier-2 (IP)")].iloc[0]
t1 = pick(HUMAN, "path3")           # tier-1 >=1 mol, IP
dflt = pick(HUMAN, "path4")         # tier-1 >=2 mol, IP
both = pts[(pts.dataset == HUMAN) & (pts.short == "both tiers (IP)")].iloc[0]
assert t2.n + t1.n == both.n, (t2.n, t1.n, both.n)
n_single = int(t1.n - dflt.n)
frac_single = n_single / t1.n        # v2 0.7223 -> "72.2 %"; v1 0.7253 -> "72.5 %"
# SIDECAR-GENERATION CHECK (the panel went to Fig S2; the assert stayed here)
assert round(frac_single * 100, 1) == VERSIONS["expect_frac_single"], frac_single
m_single = {}
for rep in ("mouse1", "mouse2"):
    a, b = pick(MOUSE, "path3", rep), pick(MOUSE, "path4", rep)
    m_single[rep] = (int(a.n - b.n), (a.n - b.n) / a.n)

tiers = pd.DataFrame([
    dict(slice="tier-2 (coverage-only)", n=int(t2.n), P_of_output=t2.P, output_scored="tier-2 (IP)",
         source=t2.source),
    dict(slice="tier-1 single-molecule (discarded by default)", n=n_single, P_of_output=t1.P,
         output_scored="tier-1 >=1 mol (IP), i.e. singletons + >=2 mol together",
         source=f"{t1.source} minus {dflt.source} (n difference)"),
    dict(slice="tier-1 >=2 molecules (= precision default)", n=int(dflt.n), P_of_output=dflt.P,
         output_scored="precision default", source=dflt.source),
])
tiers["dataset"] = HUMAN
tiers["arm"] = H_IP_DIR
tiers["plotted"] = False        # not plotted since 2026-09-02 -- see Fig S2 panel a
tiers["frac_of_tier1_single_molecule"] = frac_single
tiers["mouse1_frac_tier1_single_molecule"] = m_single["mouse1"][1]
tiers["mouse2_frac_tier1_single_molecule"] = m_single["mouse2"][1]

# ===========================================================================
# FIGURE
# ===========================================================================
# Canvas at final print width: 180 mm = 7.09 in double column.  Two square-ish
# rows of two panels with a single shared tool key between them -- the key that
# used to be repeated per panel now serves all four.
FIG_W, FIG_H = 7.09, 6.45
fig = plt.figure(figsize=(FIG_W, FIG_H))
gs_top = fig.add_gridspec(1, 2, left=0.086, right=0.988, top=0.9535, bottom=0.6124, wspace=0.215)
gs_bot = fig.add_gridspec(1, 2, left=0.086, right=0.988, top=0.4200, bottom=0.0790, wspace=0.215)
axA = fig.add_subplot(gs_top[0, 0])
axB = fig.add_subplot(gs_top[0, 1])
axC = fig.add_subplot(gs_bot[0, 0])
axD = fig.add_subplot(gs_bot[0, 1])
ref_lines = []
callouts = {"a": 0, "b": 0, "c": 0, "d": 0}   # design law: <= 3 short callouts per panel


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED, length=2.4, width=0.7, labelsize=TYPE["tick"])
    ax.grid(True, color=GRID, lw=0.5, alpha=0.85)
    ax.set_axisbelow(True)


def panel_tag(ax, letter, title):
    ax.text(-0.005, 1.035, letter, transform=ax.transAxes, fontsize=TYPE["panel_letter"],
            fontweight="bold", va="bottom", ha="left", color=INK)
    ax.text(0.058, 1.038, sentence_case(title), transform=ax.transAxes,
            fontsize=TYPE["panel_title"], va="bottom", ha="left", color=INK)


OBSTACLES = {}          # panel key -> [(x, y, pad_pt), ...] everything already drawn


def occupy(key, x, y, pad=4.0):
    OBSTACLES.setdefault(key, []).append((float(x), float(y), pad))


# candidate label slots, in preference order: (dx pt, dy pt, ha, va)
SLOTS = [(-7, -7, "right", "top"), (8, 8, "left", "bottom"), (8, -7, "left", "top"),
         (-7, 8, "right", "bottom"), (0, -10, "center", "top"), (0, 10, "center", "bottom"),
         (11, 0, "left", "center"), (-11, 0, "right", "center")]


def callout(ax, key, s, xy, color=INK, weight="normal", slots=SLOTS, literal=False):
    """A short on-image note, auto-placed into the first slot that touches nothing.

    Bounding-box discipline (DESIGN_DIRECTIVES.md item 1): the text box is tested
    against every marker and curve sample already drawn in this panel and against
    every text already placed, so no label can land on a mark or on another label.
    """
    callouts[key] += 1
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    boxes = []
    for (ox, oy, pad) in OBSTACLES.get(key, []):
        px, py = ax.transData.transform((ox, oy))
        r = pad * fig.dpi / 72.0
        boxes.append(matplotlib.transforms.Bbox.from_extents(px - r, py - r, px + r, py + r))
    for t in ax.texts:
        boxes.append(t.get_window_extent(renderer=rend).expanded(1.06, 1.14))
    # `literal` labels start with a relational operator, not prose: directive 6's
    # initial-capital rule does not apply to them ("≥2 mol", never "≥2 Mol").
    txt = s if literal else sentence_case(s)
    for dx, dy, ha, va in slots:
        a = ax.annotate(txt, xy=xy, xytext=(dx, dy), textcoords="offset points",
                        fontsize=TYPE["annotation"], color=color, ha=ha, va=va,
                        fontweight=weight, zorder=9)
        bb = a.get_window_extent(renderer=rend).expanded(1.08, 1.16)
        ax_bb = ax.get_window_extent(renderer=rend)
        if any(bb.overlaps(b) for b in boxes) or not ax_bb.fully_contains(bb.x0, bb.y0) \
                or not ax_bb.fully_contains(bb.x1, bb.y1):
            a.remove()
            continue
        return a
    raise AssertionError(f"panel {key}: no free slot for label {s!r}")


def gate_line(ax, key):
    ax.axhline(GATE_P, color=INK, lw=0.9, ls=(0, (4, 2)), zorder=2)
    ax.text(0.012, GATE_P - 0.018, sentence_case("pre-registered gate P ≥ 0.50"),
            ha="left", va="top", fontsize=TYPE["annotation"], color=INK, zorder=9)
    callouts[key] += 1
    ref_lines.append(dict(panel=key, kind="gate", label="pre-registered gate (13 section 1)",
                          P=GATE_P, F1="", plotted=True))


def draw_competitors(ax, sub, paired):
    """One marker/colour identity per tool (_pubstyle.TOOL_STYLE); mouse = mean +- range."""
    for tool in ["scUTRquant", "SCAPTURE", "polyApipe", "Sierra", "scAPAtrap"]:
        q = sub[sub.tool == tool]
        if q.empty:
            continue
        st = TOOL_STYLE[tool]
        mk = dict(marker=st["marker"], ms=6.0, mec=st["color"], mew=1.2,
                  mfc=st.get("mfc", st["color"]), ls="none", zorder=6)
        if paired and len(q) == 2:
            x, y = q.R_det.mean(), q.P.mean()
            ax.errorbar(x, y, xerr=[[x - q.R_det.min()], [q.R_det.max() - x]],
                        yerr=[[y - q.P.min()], [q.P.max() - y]],
                        color=st["color"], lw=0.7, capsize=1.4, zorder=5, ls="none")
        else:
            x, y = q.R_det.iloc[0], q.P.iloc[0]
        # a surface ring keeps overlapping markers legible
        ax.plot(x, y, mfc="none", mec="white", mew=2.4, marker=st["marker"], ms=6.0,
                ls="none", zorder=5)
        ax.plot(x, y, **mk)
        occupy(ax.get_label(), x, y, pad=8.5)   # a label must not sit beside a competitor marker (judge pass)
        ref_lines.append(dict(panel=ax.get_label(), kind="null_genic mean (3 seeds)", label=tool,
                              P=q.null_P_mean.mean(), F1="", plotted=False))


def draw_operating_curve(ax, ds_keys, sub, paired):
    """PeakATail's CURRENT operating curve + the two named, user-choosable points."""
    c = TOOL_STYLE["PeakATail"]["color"]
    key = ax.get_label()
    if sweep is not None:
        s = sweep[sweep.dataset.isin(ds_keys)].groupby("min_molecules").agg(
            P=("atlas_agreement_precision", "mean"), R=("recall_detected_genes", "mean")
        ).sort_index()
        ax.plot(s.R.to_numpy(), s.P.to_numpy(), "-", color=c, lw=1.1, alpha=0.55, zorder=4)
        mid = s.drop(index=[k for k in (SENSITIVITY_K, PREREG_K) if k in s.index])
        ax.plot(mid.R.to_numpy(), mid.P.to_numpy(), ls="none", marker="o", ms=3.6,
                mfc=c, mec="white", mew=0.8, alpha=0.85, zorder=5)
        # the whole curve is an obstacle, not just its vertices
        for i in range(len(s) - 1):
            for f in np.linspace(0, 1, 9):
                occupy(key, s.R.iloc[i] + f * (s.R.iloc[i + 1] - s.R.iloc[i]),
                       s.P.iloc[i] + f * (s.P.iloc[i + 1] - s.P.iloc[i]), pad=2.6)
    named = []
    for role, k, filled in (("path4", PREREG_K, True), ("path3", SENSITIVITY_K, False)):
        q = sub[sub.role == role]
        if paired:
            x, y = q.R_det.mean(), q.P.mean()
            ax.errorbar(x, y, xerr=[[x - q.R_det.min()], [q.R_det.max() - x]],
                        yerr=[[y - q.P.min()], [q.P.max() - y]],
                        color=c, lw=0.7, capsize=1.4, zorder=6, ls="none")
        else:
            x, y = q.R_det.iloc[0], q.P.iloc[0]
        ax.plot(x, y, mfc="none", mec="white", mew=2.6, marker="D",
                ms=8.0 if filled else 7.0, ls="none", zorder=6)
        ax.plot(x, y, marker="D", ms=8.0 if filled else 7.0, mfc=c if filled else "white",
                mec=c, mew=1.4, ls="none", zorder=7)
        occupy(key, x, y, pad=6.5)
        named.append((x, y))
        ref_lines.append(dict(panel=ax.get_label(), kind="null_genic mean (3 seeds)",
                              label=f"PeakATail >={k} molecules", P=q.null_P_mean.mean(),
                              F1="", plotted=False))
    return named


# ---- panels a / b: the one precision/recall plane -------------------------
XLIM, YLIM = (0.0, 0.345), (0.0, 1.0)
subH = pts[pts.dataset == HUMAN]
subM = pts[pts.dataset == MOUSE]

for ax, key, sub, ds_keys, paired, title in (
        (axA, "a", subH, ["pbmc"], False, "PBMC 10k v3"),
        (axB, "b", subM, ["mouse1", "mouse2"], True, "GSE104556 testis")):
    ax.set_label(key)
    ax.set_xlim(*XLIM); ax.set_ylim(*YLIM)
    ax.set_xticks([0.0, 0.1, 0.2, 0.3])
    draw_competitors(ax, sub[sub.role.isin(["competitor", "catalog"])], paired=paired)
    (p_default, p_sens) = draw_operating_curve(ax, ds_keys, sub, paired=paired)
    gate_line(ax, key)
    # the two named, user-choosable operating points; the key below spells them out
    callout(ax, key, "≥2 mol", p_default, color=TOOL_STYLE["PeakATail"]["color"],
            weight="bold", literal=True)
    callout(ax, key, "≥1 mol", p_sens, color=TOOL_STYLE["PeakATail"]["color"], literal=True)
    # plain "R_det", not mathtext R$_{det}$: a mathtext subscript renders at 0.7x the
    # label size (7.5 -> 5.25 pt), under the 6 pt floor, and the type guard below
    # measures Text artists, not the glyph runs inside mathtext (judge finding 3)
    ax.set_xlabel(sentence_case("detected-gene recall R_det @100 bp"))
    ax.set_ylabel(sentence_case("atlas-agreement precision @100 bp"))
    panel_tag(ax, key, title)
    style(ax)

# ---- panels c / d: matched call budget ------------------------------------
for ax, ds, key, title in ((axC, HUMAN, "c", "PBMC 10k v3"),
                           (axD, MOUSE, "d", "GSE104556 testis, mouse 1")):
    ax.set_label(key)
    sub = mn[mn.dataset == ds]
    for tool in TOOL_ORDER:
        q = sub[sub.tool == tool].sort_values("N")
        if q.empty:
            continue
        st = TOOL_STYLE[tool]
        lead = tool == "PeakATail"
        if len(q) > 1:
            ax.plot(q.N.to_numpy(), q.atlas_agreement_precision_100bp.to_numpy(),
                    ls=(0, (3, 1.6)) if tool == "scUTRquant" else "-",
                    color=st["color"], lw=1.6 if lead else 1.0,
                    alpha=1.0 if lead else 0.9, zorder=6 if lead else 4)
        nat = q[q.is_native_budget]
        if len(nat):
            x, y = float(nat.N.iloc[0]), float(nat.atlas_agreement_precision_100bp.iloc[0])
            ax.plot(x, y, marker=st["marker"], ms=7.0, mfc="none", mec="white", mew=2.4,
                    ls="none", zorder=6)
            ax.plot(x, y, marker=st["marker"], ms=7.0 if lead else 5.8,
                    mfc=st.get("mfc", st["color"]), mec=st["color"], mew=1.3,
                    ls="none", zorder=7)
        elif len(q) == 1:   # SCAPTURE has no non-circular ranking column: points only
            ax.plot(q.N, q.atlas_agreement_precision_100bp, marker=st["marker"], ms=5.0,
                    mfc=st["color"], mec="white", mew=0.8, ls="none", alpha=0.9, zorder=5)
    ax.set_xscale("log")
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel(sentence_case("calls kept per tool (matched budget, log)"))
    ax.set_ylabel(sentence_case("atlas-agreement precision @100 bp"))
    panel_tag(ax, key, title)
    style(ax)
axC.set_xlim(7.5e3, 1.05e6)
axC.set_xticks([1e4, 1e5, 1e6])
axC.set_xticklabels(["10k", "100k", "1M"])
axD.set_xlim(1.8e4, 1.8e5)
axD.set_xticks([2e4, 5e4, 1e5])
axD.set_xticklabels(["20k", "50k", "100k"])
for ax in (axC, axD):
    ax.minorticks_off()

# ---- one shared key for all four panels: our operating points, then the field
PA = TOOL_STYLE["PeakATail"]["color"]
row1 = [   # labels open with a relational operator, so directive 6's capital does not apply
    Line2D([], [], ls="none", marker="D", ms=6.4, mfc=PA, mec=PA, mew=1.4,
           label="≥2 molecules — precision default"),
    Line2D([], [], ls="none", marker="D", ms=5.8, mfc="white", mec=PA, mew=1.4,
           label="≥1 molecule — sensitivity arm"),
    Line2D([], [], color=PA, lw=1.1, alpha=0.6, marker="o", ms=3.6, mfc=PA, mec="white",
           label=sentence_case("PeakATail operating curve")),
]
row2 = []
for tool in ["polyApipe", "scAPAtrap", "SCAPTURE", "Sierra", "scUTRquant"]:
    st = TOOL_STYLE[tool]
    row2.append(Line2D([], [], ls="none", color=st["color"], marker=st["marker"], ms=5.6,
                       mfc=st.get("mfc", st["color"]), mec=st["color"], mew=1.2,
                       label="scUTRquant (catalog)" if tool == "scUTRquant" else tool))
lg1 = fig.legend(handles=row1, loc="center", bbox_to_anchor=(0.537, 0.531), ncol=3, frameon=False,
                 handletextpad=0.4, columnspacing=1.7, borderpad=0.0, fontsize=TYPE["annotation"])
# NB: fig.legend() already registers the legend on the figure; the fig.add_artist(lg1)
# that used to stand here registered the SAME artist a second time, so row 1 was drawn
# twice -- invisible in the raster (identical pixels) but two stacked text/marker objects
# in the vector PDF (design-judge finding 2, 2026-09-03).
fig.legend(handles=row2, loc="center", bbox_to_anchor=(0.537, 0.4985), ncol=5, frameon=False,
           handletextpad=0.4, columnspacing=1.7, borderpad=0.0, fontsize=TYPE["annotation"])

for k, v in callouts.items():
    assert v <= 3, f"panel {k}: {v} callouts (design law: at most 3)"

# ---------------------------------------------------------------------------
# bounding-box discipline: no two pieces of text may overlap
# ---------------------------------------------------------------------------
def text_overlaps(fig):
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    items = []
    for ax in fig.axes:
        for t in ax.texts:
            items.append((ax.get_label(), t))
        for t in (ax.xaxis.label, ax.yaxis.label, ax.title):
            if t.get_text():
                items.append((ax.get_label(), t))
        for t in (ax.get_xaxis().get_ticklabels() + ax.get_yaxis().get_ticklabels()):
            if t.get_text():
                items.append((ax.get_label(), t))
        if ax.get_legend() is not None:
            items += [(ax.get_label(), t) for t in ax.get_legend().get_texts()]
    items += [("fig", t) for t in fig.texts]
    for lg in fig.legends:
        items += [("fig", t) for t in lg.get_texts()]
    bad = []
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            a, b = items[i][1], items[j][1]
            if not (a.get_visible() and b.get_visible()):
                continue
            ba = a.get_window_extent(renderer=rend).expanded(1.02, 1.06)
            bb = b.get_window_extent(renderer=rend).expanded(1.02, 1.06)
            if ba.overlaps(bb):
                bad.append((items[i][0], repr(a.get_text())[:34], items[j][0],
                            repr(b.get_text())[:34]))
    return bad


_bad = text_overlaps(fig)
for row in _bad:
    print("TEXT OVERLAP:", row)
assert not _bad, f"{len(_bad)} text/text overlaps -- fix before shipping"

# design law (DESIGN_DIRECTIVES.md item 1, added at the 2026-09-03 no-loss audit):
# no annotation may be set below the minimum print size.  The overlap assert above
# and the 8-px edge assert below already cover collisions and clipping; this closes
# the third law so all six mains carry the same three checks.
_small = sorted({(round(float(t.get_fontsize()), 2), t.get_text()[:34])
                 for ax in fig.axes
                 for t in (list(ax.texts) + [ax.xaxis.label, ax.yaxis.label, ax.title]
                           + list(ax.get_xaxis().get_ticklabels())
                           + list(ax.get_yaxis().get_ticklabels())
                           + (list(ax.get_legend().get_texts()) if ax.get_legend() else []))
                 if t.get_text().strip() and t.get_visible()
                 and float(t.get_fontsize()) < TYPE["annotation_min"]}
                | {(round(float(t.get_fontsize()), 2), t.get_text()[:34]) for t in fig.texts
                   if t.get_text().strip() and float(t.get_fontsize()) < TYPE["annotation_min"]}
                | {(round(float(t.get_fontsize()), 2), t.get_text()[:34])
                   for lg in fig.legends for t in lg.get_texts()
                   if t.get_text().strip() and float(t.get_fontsize()) < TYPE["annotation_min"]})
assert not _small, ("annotation below the minimum print size", _small)
print(f"design check: no text overlaps, no text below {TYPE['annotation_min']} pt")


for ext, kw in (("png", dict(dpi=600)), ("pdf", {})):
    p = FIGDIR / f"{NAME}.{ext}"
    fig.savefig(p, **kw)
    print("wrote", p)
p = OUTDIR / f"{NAME}.png"
fig.savefig(p, dpi=600); print("wrote", p)

# ---------------------------------------------------------------------------
# TSVs: every value, plotted or kept for the record
# ---------------------------------------------------------------------------
PANEL_OF = {"competitor": {HUMAN: "a", MOUSE: "b"}, "catalog": {HUMAN: "a", MOUSE: "b"},
            "path3": {HUMAN: "a", MOUSE: "b"}, "path4": {HUMAN: "a", MOUSE: "b"}}
pts_out = pts.copy()
pts_out["short"] = pts_out["short"].str.replace("\n", " ")
pts_out["plotted"] = pts_out.role.isin(PLOTTED_ROLES)
pts_out["panel"] = [PANEL_OF.get(r.role, {}).get(r.dataset, "") for r in pts_out.itertuples()]
cols = ["panel", "plotted", "dataset", "tool", "call_set", "short", "replicate", "role", "cutoff_bp",
        "n", "n_matched", "P", "R_det", "R_det_denominator", "R_full", "R_full_denominator", "F1_det",
        "null_P_seed1", "null_P_seed2", "null_P_seed3", "null_P_mean", "source"]
pts_out = pts_out[cols].rename(columns={"P": "atlas_agreement_precision_100bp",
                                        "R_det": "recall_detected_genes_100bp",
                                        "R_full": "recall_full_atlas_100bp",
                                        "F1_det": "F1_detected_genes_100bp"})
p = OUTDIR / f"{NAME}.tsv"; pts_out.to_csv(p, sep="\t", index=False, float_format="%.6f"); print("wrote", p)

if sweep is not None:
    sw = sweep.copy()
    sw.insert(0, "panel", ["a" if d == "pbmc" else "b" for d in sw.dataset])
    sw.insert(1, "plotted", True)
    sw.insert(2, "tool", "PeakATail v2 (IP arm)")
    p = OUTDIR / f"{NAME}_sweep.tsv"
    sw.to_csv(p, sep="\t", index=False, float_format="%.6f"); print("wrote", p)

p = OUTDIR / f"{NAME}_matchedN.tsv"
mn.to_csv(p, sep="\t", index=False, float_format="%.6f"); print("wrote", p)

# NOT PLOTTED since 2026-09-02 -- the tier panel moved to Fig S2 (directive 5);
# the table is still emitted so the numbers stay pinned to a reproducible file.
p = OUTDIR / f"{NAME}_tiers.tsv"; tiers.to_csv(p, sep="\t", index=False, float_format="%.6f"); print("wrote", p)

# NOT PLOTTED since 2026-09-02 (except the pre-registered gate): the F1 isolines,
# the superseded two-sided gate and the per-call-set nulls keep their rows here.
#
# NO-LOSS (2026-09-03 audit): the loop above records a null row only for the call
# sets the redesigned plane still draws.  The arms the design pass retired --
# "shipped (pre-fix)", the no-IP intermediates, both-tiers and tier-2 -- had a
# per-call-set null row in the pre-redesign file, and the Legend's bound
# (NULL_MAX) and its PBMC/mouse ranges are read from THOSE rows.  They are
# re-emitted here with plotted=False so this file stays the complete per-call-set
# null record the Legend points at, and nothing that was in it was dropped.
# (matched on the null value, not the label: a drawn row is labelled by tool or by
#  operating point, a retired row by its call-set short name, but both are the mean
#  over the panel's replicates, so the value is the identity of the call set here.)
_have = {(r["panel"], round(float(r["P"]), 9)) for r in ref_lines
         if r["kind"].startswith("null_genic")}
for _ds, _panel in ((HUMAN, "a"), (MOUSE, "b")):
    for _short, _q in pts[pts.dataset == _ds].groupby("short", sort=False):
        _v = float(_q.null_P_mean.mean())
        if (_panel, round(_v, 9)) in _have:
            continue
        ref_lines.append(dict(panel=_panel, kind="null_genic mean (3 seeds)",
                              label=" ".join(str(_short).split()),   # no newline in a TSV cell
                              P=_v, F1="", plotted=False))
        _have.add((_panel, round(_v, 9)))
ref_lines.append(dict(panel="a,b", kind="gate", label="original two-sided gate (10 section 5)",
                      P=ORIG_P, F1=ORIG_F1, plotted=False))
for f in F1_ISO:
    ref_lines.append(dict(panel="a,b", kind="F1 isoline", label=f"F1 = {f}", P="", F1=f, plotted=False))
p = OUTDIR / f"{NAME}_reference_lines.tsv"
_rl = pd.DataFrame(ref_lines)
_rl.to_csv(p, sep="\t", index=False, float_format="%.6f"); print("wrote", p)

# SIDECAR-GENERATION CHECK: the Legend sentence "every call set's 3-seed
# gene-body-shuffled null P@100 is at most NULL_MAX ... per-call-set values are in
# fig2_accuracy_reference_lines.tsv" is only true while that file really holds one
# row per scored call set and its extremes are the ones quoted.
_nl = _rl[_rl.kind.str.startswith("null_genic")]
_want = pts[pts.dataset == HUMAN].short.nunique() + pts[pts.dataset == MOUSE].short.nunique()
assert len(_nl) == _want, \
    ("per-call-set null rows missing from the reference-lines record", len(_nl), _want)
assert round(float(_nl.P.astype(float).max()), 4) == round(NULL_MAX, 4), \
    (float(_nl.P.astype(float).max()), NULL_MAX)
assert round(float(_nl[_nl.panel == "a"].P.astype(float).max()), 4) == \
    round(float(subH.null_P_mean.max()), 4)                       # the quoted PBMC top

# ---------------------------------------------------------------------------
# sidecar caption -- the single source of the journal legend.  Everything the
# 2026-09-02 design pass took off the image lands here (no-loss rule).
# ---------------------------------------------------------------------------
poly_h = subH[subH.tool == "polyApipe"].iloc[0]
poly_m = subM[subM.tool == "polyApipe"]
rel_h = (dflt.R_det / poly_h.R_det - 1) * 100
rel_m = (subM[subM.role == "path4"].R_det.mean() / poly_m.R_det.mean() - 1) * 100
dF1_h = dflt.F1_det - poly_h.F1_det
dF1_m = subM[subM.role == "path4"].F1_det.mean() - poly_m.F1_det.mean()
m_full = subM[subM.role == "path4"].sort_values("replicate").R_full.tolist()
h1 = pick(HUMAN, "path3")
m1_1 = pick(MOUSE, "path3", "mouse1")
m2_1 = pick(MOUSE, "path3", "mouse2")

_prev = VERSIONS["prev"]
prev_note = "" if not _prev else (
    f"The {_prev['label']} gave {_prev['P'][0]:.3f} / {_prev['P'][1]:.3f} / {_prev['P'][2]:.3f}; the "
    f"−1.1 pp PBMC move is the corrected minus-strand IP window (#96), reconstructed key-for-key. ")
_ip = VERSIONS["ip_trade"]
ip_note = _ip[0].upper() + _ip[1:]

if sweep is not None:
    _sw = sweep.copy()
    _swh = _sw[_sw.dataset == "pbmc"].set_index("min_molecules")
    _swm = _sw[_sw.dataset != "pbmc"].groupby("min_molecules").agg(
        P=("atlas_agreement_precision", "mean"), R=("recall_detected_genes", "mean"),
        n=("n", "mean"))
    curve_md = (
        f"The light line is the tool's molecule-support operating curve on the IP-filtered arms: "
        f"≥1 / ≥2 / ≥3 / ≥5 / ≥10 distinct clip molecules give P@100 "
        f"{_swh.loc[1].atlas_agreement_precision:.3f} / {_swh.loc[2].atlas_agreement_precision:.3f} / "
        f"{_swh.loc[3].atlas_agreement_precision:.3f} / {_swh.loc[5].atlas_agreement_precision:.3f} / "
        f"{_swh.loc[10].atlas_agreement_precision:.3f} at R_det "
        f"{_swh.loc[1].recall_detected_genes:.3f} / {_swh.loc[2].recall_detected_genes:.3f} / "
        f"{_swh.loc[3].recall_detected_genes:.3f} / {_swh.loc[5].recall_detected_genes:.3f} / "
        f"{_swh.loc[10].recall_detected_genes:.3f} on PBMC (n {int(_swh.loc[1].n):,} → "
        f"{int(_swh.loc[10].n):,}), and P@100 {_swm.loc[1].P:.3f} → {_swm.loc[10].P:.3f} at R_det "
        f"{_swm.loc[1].R:.3f} → {_swm.loc[10].R:.3f} as the two-mouse mean. **Only the ≥2-molecule "
        f"point was pre-registered** (13 §1); ≥1 is the shipped sensitivity arm and ≥3/≥5/≥10 "
        f"are descriptive points, never gated, drawn small and unlabelled for that reason. (The post-hoc sweep "
        f"that informed the ≥2 threshold is disclosed in 12 CORRECTION and is never cited as a result.) "
        f"Every curve point was recomputed by `score_tool.py` from `pas.bed` column 5 with the reference "
        f"arguments of `scripts/benchmark_tools/stage2_final_launch.sh`; the ≥1 and ≥2 points reproduce "
        f"{GATE_DOC} §1 to 4 dp on all three arms (asserted). ")
else:
    curve_md = ("The molecule-support operating curve was computed only on the v2 arms, so this v1 render "
                "draws the two named operating points without the connecting curve. ")

caption = (
    "Figure 2 | Where PeakATail sits against the field, and what it costs at a matched call budget "
    f"(PeakATail code {CODE_COMMIT}, {VERSIONS['commit_note']}). "
    "Definitions: atlas-agreement precision = fraction of calls within 100 bp (point, strand-matched) of a "
    "PolyASite 2.0 representative site; it is agreement with a curated atlas, not ground truth (atlas-novel "
    "true sites count as false positives). R_det denominator = atlas sites in genes detected in the dataset "
    f"({int(dflt.R_det_denominator):,} human / {int(_m1.R_det_denominator):,} mouse); full-atlas recall of the "
    f"default {dflt.R_full:.3f} (PBMC; {int(dflt.R_full_denominator):,}) and {m_full[0]:.3f} / {m_full[1]:.3f} "
    f"(mice; {int(_m1.R_full_denominator):,}). One scorer (score_tool.py), same denominators for every tool; "
    "3-seed gene-body-shuffled nulls. scUTRquant is catalog-based (hollow marker, dashed line): shown, not "
    "ranked with the de novo tools. "
    f"Pre-registration: the precision-first default and its P@100 ≥ 0.50 gate were committed before the "
    f"ORIGINAL run and are unchanged here ({PREREG}; 13 §1, {GATE_DOC}); the ≥2-molecule threshold was "
    "pre-registered before that run (12 CORRECTION). "
    f"Gate: PASS on all three ({dflt.P:.3f} / {_m1.P:.3f} / {_m2.P:.3f}); no ≥1-molecule output clears the "
    "superseded original P ≥ 0.38 floor (10 §5), which is why that line is recorded in "
    "`fig2_accuracy_reference_lines.tsv` rather than drawn. "
    + prev_note +
    f"{ip_note}; compute {VERSIONS['compute']}; {GATE_DOC} §1/§3. "
    # CORRECTED 2026-08-21 (verifier pass): the pre-correction clause ("recall at or below polyApipe's ...
    # the F1_det lead is not a meaningful margin") stated a real limitation but compared 46,524 of our calls
    # with 120,916 of polyApipe's. The builder's first replacement dropped the limitation entirely and kept
    # only the favourable half; both halves must be here. So: the single-point deltas STAY (rel_h / rel_m /
    # dF1_h / dF1_m, data-driven), and the matched-N context is added and bounded to a specific N.
    # Verified source: manuscript/25_competitive_position.md §2, §8.1.
    f"Caveats: single PBMC donor / BAM; the two mice are one study and chemistry; at its own default the tool's "
    f"recall is below polyApipe's ({rel_h:+.0f}% / {rel_m:+.0f}%), F1_det lead {dF1_h:+.3f} / {dF1_m:+.3f} — "
    f"at 2.6x fewer calls, which is exactly the confound panels c/d remove. "
    "The non-IP PBMC ≥2-molecule file (pas_tier1_ge2mol_noIP_POSTHOC) is not the default and is not shown."
)

removed_md = f"""**What this figure deliberately does not show** (2026-09-02 publication design pass,
`manuscript/figures/DESIGN_DIRECTIVES.md` items 4 and 5; nothing was deleted, everything below is still
computed, asserted and written to the audit TSVs).
*Development history* — the "shipped (pre-fix)" caller and the two no-IP intermediate outputs are no longer
drawn: end users choose between shipping operating points, not between our development stages. Their values are
in `fig2_accuracy.tsv` with `plotted = False` (PBMC shipped P@100 {pick(HUMAN, 'path0').P:.3f} at
R_det {pick(HUMAN, 'path0').R_det:.3f}, n {int(pick(HUMAN, 'path0').n):,}; both tiers no IP
{pick(HUMAN, 'path1').P:.3f} / {pick(HUMAN, 'path1').R_det:.3f}; tier-1 ≥1 mol no IP
{pick(HUMAN, 'path2').P:.3f} / {pick(HUMAN, 'path2').R_det:.3f}), and the version-to-version story is **Fig S12**
(caller versions) and **Fig S8** (compute).
*Tier decomposition* — the former panel d duplicated **Fig S2 panel a**, which is the tier composition of every
final v2 arm. On the PBMC IP arm all calls scored are n = {int(both.n):,} (P@100 {both.P:.3f}), splitting
tier-2 {int(t2.n):,} (P@100 {t2.P:.3f}) /
tier-1 single-molecule {n_single:,} / tier-1 ≥2 molecules {int(dflt.n):,} (P@100 {dflt.P:.3f}), i.e.
{frac_single*100:.1f}% of tier-1 sites are single-molecule and are dropped by the default
(mouse {m_single['mouse1'][1]*100:.1f}% / {m_single['mouse2'][1]*100:.1f}%); the ≥1-molecule arm scores
{t1.P:.3f}. Definition matters here: {VERSIONS['single_mol_note']}. Table: `fig2_accuracy_tiers.tsv`.
*F1 isolines and the superseded gate* — F1_det isolines 0.1–0.4 and the original two-sided gate
(P ≥ 0.38 & F1_det > 0.261, 10 §5) were chrome on a plane that now carries an operating curve; both are
rows in `fig2_accuracy_reference_lines.tsv`.
*Genic-shuffle nulls* — the per-call-set null ticks are gone from the plane; every call set's 3-seed
gene-body-shuffled null P@100 is at most {NULL_MAX:.4f} (PBMC call sets {float(subH.null_P_mean.min()):.4f}–{float(subH.null_P_mean.max()):.4f},
mouse {float(subM.null_P_mean.min()):.4f}–{float(subM.null_P_mean.max()):.4f}); per-call-set values are in
`fig2_accuracy_reference_lines.tsv` and the null design is **Fig S3**."""

framing_md = f"""**Competitive framing (`manuscript/25_competitive_position.md`, corrected 2026-08-21).** Panels c/d
are the drawn form of this: at **matched call count** PeakATail leads polyApipe on recall at every N tested and on
precision at every N from 20,320 upward, on both datasets — at polyApipe's own N = 120,916, 0.4036 / 0.2336 vs
0.3800 / 0.1988; at our own N = 46,524, 0.7062 / 0.1754 vs 0.6049 / 0.1375; mouse 1 at N = 71,983,
0.4488 / 0.3034 vs 0.4353 / 0.2231 — and on **both** mice the ≥1-molecule arm beats polyApipe on precision,
recall **and** F1 simultaneously ({m1_1.P:.4f} / {m1_1.R_det:.4f} / {m1_1.F1_det:.4f} and {m2_1.P:.4f} /
{m2_1.R_det:.4f} / {m2_1.F1_det:.4f} vs 0.4005 / 0.2499 / 0.3078 and 0.4119 / 0.2508 / 0.3118). The
≥2-molecule default is also **not** the F1 optimum ({h1.F1_det:.4f} / {m1_1.F1_det:.4f} / {m2_1.F1_det:.4f}
for the ≥1-molecule arm vs the default's {dflt.F1_det:.4f} / {_m1.F1_det:.4f} / {_m2.F1_det:.4f}) — a
deliberate reliability choice, not a ceiling. Qualifiers that travel with panels c/d: the claim holds among
**de novo** tools only (catalog-based scUTRquant is above PeakATail on both axes at several matched N — visible
as the dashed line above the blue one at small budgets), the precision half is not established for
20,320 ≤ N ≤ 35,759 because of competitor-side ties and is **false below N ≈ 15,000, where polyApipe is
higher (0.9671 vs 0.9527 at N = 9,456)** while the recall half holds at every N, "PeakATail at scAPAtrap's N" is
not reachable (787,138 > our maximum 333,920, so the blue line stops), and SCAPTURE has no non-circular ranking
column so it cannot be truncated and appears as isolated markers. Full grid, tie-break brackets and provenance:
`manuscript/25_competitive_position.md` §2, §6, §8; every plotted point:
`results/figures/manuscript/fig2_accuracy_matchedN.tsv`."""

cap_md = f"""# Fig 2 — `fig2_accuracy` caption (generated by `scripts/manuscript_figures/fig2_accuracy.py`)

## Legend

{caption}

**Panel a** — PBMC 10k v3: atlas-agreement precision @100 bp (y) against detected-gene recall R_det @100 bp (x).
Competitors sit at their own operating points, one marker and colour per tool throughout the paper
(`scripts/manuscript_figures/_pubstyle.py`). {curve_md}Two points are named because they are the two the user
chooses between: the filled diamond is the **pre-registered precision default** (tier-1 ∩ IP-pass ∩ ≥2
distinct clip molecules) and the open diamond the **≥1-molecule sensitivity arm**. Dashed line: the
pre-registered gate P@100 ≥ 0.50 (13 §1). scTail is absent: not runnable on this BAM (R1 = 28 bp).
**Panel b** — GSE104556 testis; competitor markers and the two named points are the mean of the two mice with
range bars, and the curve is the two-mouse mean (per-mouse values in the audit TSVs). All mouse arms ran with the
IP filter. SCAPTURE is plotted for mouse 1 only, but its mouse-2 site-level run COMPLETED (24,076 points, P@100
0.672 — `gse104556/scapture/mouse2/score_scapture_mouse2.tsv`; 15 §3 reports both mice, 0.694 / 0.672); only
the per-cell PASquant step failed, which site-level benchmarking does not use — never describe the mouse-2 run as
invalid (21 §11 stale-caption fix, 2026-09-02).
**Panels c and d** — the same tools compared at a **matched call budget**: each ranked call set is truncated to a
common number of calls and re-scored (D3 rank sweep, `results/perf_gap/D3_head_to_head/`), so precision is read at
equal N rather than at each tool's own N. The filled marker on each line is that tool's own shipped budget
(PeakATail {D3_PA_NATIVE[HUMAN]:,} human / {D3_PA_NATIVE[MOUSE]:,} mouse 1; polyApipe {NATIVE_N[(HUMAN, 'polyApipe')]:,} /
{NATIVE_N[(MOUSE, 'polyApipe')]:,}; scAPAtrap {NATIVE_N[(HUMAN, 'scAPAtrap')]:,} /
{NATIVE_N[(MOUSE, 'scAPAtrap')]:,}). This is the honest head-to-head: the recall gap at each tool's own default
(see Caveats above) is a call-budget effect, and it reverses when the budget is equalised. Panels c/d are the
manuscript/25 record, which is a v2 analysis in both render modes of this script.

{framing_md}

**Atlas-independent corroboration (C2; 16 §v2, verified FIXED).** The same precision default is corroborated by
truth that owes nothing to the atlas: 0.7647 [Wilson 95% CI 0.7608–0.7685] of its 46,524 sites (35,575) lie within
25 bp (strand-matched) of a poly(A)-verified Kinnex x3p 3′ end at ≥5 UMI, vs a gene-body-shuffled null of 0.0072
(10 seeds; 106.5×); the atlas-known hexamer-pass complement scores 0.8941 on the same truth, calibrating the
ceiling. Caveats travel (21 §5.1-C2): the Kinnex truth comes from different donors, and its "≥5 UMI" support is an
alignment-record count that was never de-duplicated, so the support thresholds are slightly optimistic. Expanded in
Fig S7 (the pre-registered trusted-novel NEGATIVE result lives there); every value:
`results/figures/manuscript/figS7_novelfunnel_concordance.tsv` (set CAL_default_output_all).

{removed_md}

## Provenance

Verified: {VERSIONS['verified_short']}. Every plotted value: `results/figures/manuscript/fig2_accuracy*.tsv`
(`fig2_accuracy.tsv` = the per-arm points with a `plotted` flag;
{'`_sweep.tsv` = the operating curve; ' if sweep is not None else '`_sweep.tsv` is NOT written by this v1 render — the operating curve exists only for the v2 arms; '}
`_matchedN.tsv` = panels c/d; `_tiers.tsv` = the tier table, not plotted since 2026-09-02;
`_reference_lines.tsv` = gates, F1 isolines and per-call-set nulls, with a `plotted` flag).
Sources: `{VERSIONS["gate_doc_path"]}` + `{VERSIONS["verified_table"]}` ({VERSIONS["verdict"]}), the score TSVs
listed in the `source` column of `fig2_accuracy.tsv`, the shared sweep module
`scripts/manuscript_figures/_molsweep.py`, and `results/perf_gap/D3_head_to_head/tsv/01d,01e`
(`manuscript/25_competitive_position.md`). Script: `scripts/manuscript_figures/fig2_accuracy.py`
(env switch `FINAL_BENCHMARK_VERSION={VERSION}` selected this run; the data-rich pre-2026-09-02 layout is
reachable from git history — this pass changed the design, not a number).
"""
p = FIGDIR / f"{NAME}.caption.md"; p.write_text(cap_md); print("wrote", p)

# ---------------------------------------------------------------------------
# edge check: nothing may be clipped; the outer 8 px of the PNG must be blank
# ---------------------------------------------------------------------------
try:
    from PIL import Image
    im = np.asarray(Image.open(FIGDIR / f"{NAME}.png").convert("L"))
    edge = 8
    border = np.concatenate([im[:edge, :].ravel(), im[-edge:, :].ravel(),
                             im[:, :edge].ravel(), im[:, -edge:].ravel()])
    ink = int((border < 250).sum())
    print(f"edge check: {im.shape[1]}x{im.shape[0]} px, ink pixels in the outer {edge} px = {ink}")
    assert ink == 0, "INK IN THE OUTER 8 PX -- something is clipped"
except ImportError:
    print("edge check skipped (no PIL)")

print(f"VERSION={VERSION} code={CODE_COMMIT}: precision default P@100 "
      f"{_d.P:.4f} / {_m1.P:.4f} / {_m2.P:.4f} (n {int(_d.n):,} / {int(_m1.n):,} / {int(_m2.n):,}); "
      f"sensitivity arm {h1.P:.4f} / {m1_1.P:.4f} / {m2_1.P:.4f}; "
      f"callouts per panel {callouts}")

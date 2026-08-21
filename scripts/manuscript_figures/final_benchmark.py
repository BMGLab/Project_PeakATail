#!/usr/bin/env python3
"""
final_benchmark.py -- manuscript Fig 2 ("final_benchmark"): the final Stage-2
run head-to-head against the competitor panel.

VERSION SWITCH -- env var FINAL_BENCHMARK_VERSION picks which PeakATail run is
plotted.  The competitor / catalog arms are the same TSVs in both versions.
  v2  (default; THE PAPER NUMBERS)  code 9dfdefb = #93 + #96 (IP-filter
      minus-strand fix) + #97 (performance), arms run 2026-08-21 16:14-16:50.
      Sources of truth: manuscript/19_final_gate_v2.md and
      results/benchmark_tools/final_v2_verify/VERIFIED_v2.md (verifier: FIXED).
  v1  (FINAL_BENCHMARK_VERSION=v1; kept selectable for the record)
      code 4efeb125, manuscript/15_final_gate.md (verified SOUND).

SOURCE OF TRUTH -- every plotted value is READ from a score TSV written by the
single scoring path scripts/benchmark_tools/score_tool.py (100 bp cutoff,
strand-matched point mode, curated PolyASite 2.0 representative sites,
per-dataset detected-gene recall denominator).  Nothing numeric is typed by
hand except the facts that exist only in the gate write-up (pre-registration
timestamps, code commit, compute, the IP-filter precision/recall trade) --
those are cited to 19 (v2) / 15 (v1) on the figure.  A consistency assert
below checks the headline ones against the verified gate table.

PANELS
  a  PBMC 10k v3: atlas-agreement precision@100 (y) vs detected-gene recall@100
     (x), every tool; scUTRquant hollow (catalog-based, not ranked); scTail
     absent (not runnable: R1 = 28 bp); PeakATail as a connected path of
     operating points shipped -> both tiers (no IP) -> tier-1 >=1 mol (no IP)
     -> tier-1 >=1 mol (IP) -> precision default (tier-1, IP, >=2 molecules).
     F1 isolines 0.1..0.4; pre-registered gate P >= 0.50 (13 section 1);
     original two-sided gate P >= 0.38 & F1_det > 0.261 (10 section 5);
     3-seed genic-shuffle null (mean) as small markers at the bottom.
  b  Same for GSE104556 testis, two mice as mean +- range (both mice in TSV).
     The mouse arms all ran with the IP filter, so the path has no "no IP" step.
  c  n sites per call set, log scale, same tools / outputs, both datasets.
  d  PeakATail tier decomposition on the PBMC IP arm: tier-2 / tier-1
     single-molecule / tier-1 >=2 molecules (= default), with each output's
     P@100.  72.2 % of tier-1 sites are single-molecule on v2 (derived from the
     two score-TSV n's: 1 - 46,524 / 167,565; 72.5 % on v1).

NAMING RULES (manuscript/01 number policy)
  * "atlas-agreement precision (100 bp)", never bare "precision".
  * recall shown = detected-gene recall R_det; full-atlas recall in caption.
  * .../peakatail_clipseeded_final_v2/pas_tier1_ge2mol_noIP_POSTHOC.bed is NOT
    the pre-registered default and is not plotted.

OUTPUTS
  manuscript/figures/final_benchmark.{png,pdf}   (PNG 300 dpi, PDF fonttype 42)
  manuscript/figures/final_benchmark.caption.md  (sidecar caption)
  results/figures/manuscript/final_benchmark.tsv            (every point, a/b/c)
  results/figures/manuscript/final_benchmark_tiers.tsv      (panel d)
  results/figures/manuscript/final_benchmark_reference_lines.tsv (gates, isolines, nulls)

Run:  export LC_ALL=C; python3 scripts/manuscript_figures/final_benchmark.py
      export LC_ALL=C; FINAL_BENCHMARK_VERSION=v1 python3 .../final_benchmark.py
"""
import os
import textwrap
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "DejaVu Sans"
matplotlib.rcParams["font.size"] = 7.5
matplotlib.rcParams["axes.titlesize"] = 8
matplotlib.rcParams["axes.labelsize"] = 8
matplotlib.rcParams["legend.fontsize"] = 6.3

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
BT = WD / "results/benchmark_tools"
OUTDIR = WD / "results/figures/manuscript"
FIGDIR = WD / "manuscript/figures"
NAME = "final_benchmark"
OUTDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)

CUTOFF = 100.0
INK, MUTED, GRID = "#1B2429", "#5A6B73", "#D8E0E3"

# Okabe-Ito colourblind-safe palette (figures/README.md convention).
COLOR = {
    "PeakATail":  "#0072B2",
    "polyApipe":  "#D55E00",
    "SCAPTURE":   "#009E73",
    "Sierra":     "#CC79A7",
    "scAPAtrap":  "#E69F00",
    "scUTRquant": "#56B4E9",
}
HUMAN, MOUSE = "pbmc_10k_v3", "gse104556"

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
        ip_trade="the IP filter buys +5.3 pp precision for \u22121.2 pp recall (mouse 1)",
        prev=dict(label="v1 run (code 4efeb125, 15)", P=(0.7167, 0.7414, 0.7554)),
        expect_h=(46524, 0.7062, 0.1754), expect_m1=(26255, 0.7450), expect_m2=(26526, 0.7572),
        expect_frac_single=72.2,
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
        ip_trade="the IP filter buys +4.9 pp precision for \u22121.4 pp recall (mouse 1)",
        prev=None,
        expect_h=(44394, 0.7167, 0.1707), expect_m1=(25991, 0.7414), expect_m2=(26164, 0.7554),
        expect_frac_single=72.5,
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
    # ---- PBMC PeakATail path (order matters) --------------------------------
    (HUMAN, "PeakATail", "PeakATail shipped (pre-fix caller)", "shipped", "",
     H / "peakatail/score_peakatail.tsv", "path0"),
    (HUMAN, "PeakATail", "PeakATail both tiers, no IP filter", "both tiers\n(no IP)", "",
     H / H_DIR / f"score_{H_PFX}.tsv", "path1"),
    (HUMAN, "PeakATail", "PeakATail tier-1 >=1 molecule, no IP filter", "tier-1 >=1 mol\n(no IP)", "",
     H / H_DIR / f"score_{H_PFX}__tier1.tsv", "path2"),
    (HUMAN, "PeakATail", "PeakATail tier-1 >=1 molecule, IP filter", "tier-1 >=1 mol (IP)", "",
     H / H_IP_DIR / f"score_{H_IP_PFX}__tier1.tsv", "path3"),
    (HUMAN, "PeakATail", "PeakATail precision default (tier-1, IP, >=2 molecules; pre-registered)",
     "precision default\n(tier-1, IP, >=2 mol)", "",
     H / H_IP_DIR / f"score_{H_IP_PFX}__PRESPEC_precision_default.tsv", "path4"),
    # ---- PBMC extra PeakATail outputs (panel c/d only, not on the path) -----
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
# SCAPTURE mouse: mouse1 only (mouse2 run truncated by disk exhaustion; 09).
ARMS.append((MOUSE, "SCAPTURE", "SCAPTURE", "SCAPTURE", "mouse1",
             M / "scapture/mouse1/score_scapture_mouse1.tsv", "competitor"))


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
GATE_P = 0.50          # 13 section 1, pre-registered precision-first gate
ORIG_P, ORIG_F1 = 0.38, 0.261   # 10 section 5 original two-sided gate (PBMC)
for r in pts[pts.role == "path4"].itertuples():
    assert r.P >= GATE_P, r
for r in pts[(pts.dataset == HUMAN) & (pts.role.isin(["path1", "path2", "path3"]))].itertuples():
    assert r.P < ORIG_P, r   # "no >=1-molecule output clears the 0.38 floor"

# ---------------------------------------------------------------------------
# panel d: tier decomposition on the PBMC IP arm (all n from the score TSVs)
# ---------------------------------------------------------------------------
t2 = pts[(pts.dataset == HUMAN) & (pts.short == "tier-2 (IP)")].iloc[0]
t1 = pick(HUMAN, "path3")           # tier-1 >=1 mol, IP
dflt = pick(HUMAN, "path4")         # tier-1 >=2 mol, IP
both = pts[(pts.dataset == HUMAN) & (pts.short == "both tiers (IP)")].iloc[0]
assert t2.n + t1.n == both.n, (t2.n, t1.n, both.n)
n_single = int(t1.n - dflt.n)
frac_single = n_single / t1.n        # v2 0.7223 -> "72.2 %"; v1 0.7253 -> "72.5 %"
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
tiers["frac_of_tier1_single_molecule"] = frac_single
tiers["mouse1_frac_tier1_single_molecule"] = m_single["mouse1"][1]
tiers["mouse2_frac_tier1_single_molecule"] = m_single["mouse2"][1]

# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(8.3, 10.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.82], width_ratios=[1, 1],
                      left=0.15, right=0.975, top=0.965, bottom=0.18,
                      hspace=0.38, wspace=0.30)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0])
axD = fig.add_subplot(gs[1, 1])

F1_ISO = [0.1, 0.2, 0.3, 0.4]
ref_lines = []


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED, length=2.5, labelsize=7)
    ax.grid(True, color=GRID, lw=0.5, alpha=0.7)
    ax.set_axisbelow(True)


def f1_iso(ax, xmax):
    R = np.linspace(0.005, xmax, 400)
    for f in F1_ISO:
        with np.errstate(divide="ignore", invalid="ignore"):
            P = f * R / (2 * R - f)
        ok = (2 * R - f > 0) & (P <= 1.0) & (P >= 0)
        ax.plot(R[ok], P[ok], color="#C9D1D5", lw=0.7, zorder=1)
        # label where the isoline meets the right edge or the top
        xr = R[ok][-1]; yr = P[ok][-1]
        ax.text(xr, yr, f"F1 {f:.1f}", color="#9AA7AD", fontsize=5.8,
                ha="right", va="bottom", zorder=1)


def gates(ax, xmax, original=True, original_note=""):
    ax.axhline(GATE_P, color=INK, lw=0.9, ls=(0, (4, 2)), zorder=2)
    ax.text(xmax * 0.995, GATE_P + 0.008, "pre-registered gate P ≥ 0.50 (13 §1)",
            ha="right", va="bottom", fontsize=5.9, color=INK)
    ref_lines.append(dict(panel=ax.get_label(), kind="gate", label="pre-registered gate (13 section 1)",
                          P=GATE_P, F1=""))
    if original:
        ax.axhline(ORIG_P, color=MUTED, lw=0.8, ls=(0, (1, 2)), zorder=2)
        R = np.linspace(ORIG_F1 / 2 + 0.002, xmax, 400)
        P = ORIG_F1 * R / (2 * R - ORIG_F1)
        ok = (P <= 1.0)
        ax.plot(R[ok], P[ok], color=MUTED, lw=0.8, ls=(0, (1, 2)), zorder=2)
        ax.text(0.004, ORIG_P - 0.012,
                "original gate P ≥ 0.38 & F1 > 0.261 (10 §5)" + original_note,
                ha="left", va="top", fontsize=5.8, color=MUTED)
        ref_lines.append(dict(panel=ax.get_label(), kind="gate", label="original two-sided gate (10 section 5)",
                              P=ORIG_P, F1=ORIG_F1))
    for f in F1_ISO:
        ref_lines.append(dict(panel=ax.get_label(), kind="F1 isoline", label=f"F1 = {f}", P="", F1=f))


PATH_ROLES_H = ["path0", "path1", "path2", "path3", "path4"]
PATH_ROLES_M = ["path0", "path1", "path3", "path4"]


def draw_competitors(ax, sub, paired):
    """sub: rows of one dataset, competitor+catalog roles.  paired: mouse."""
    handles = []
    for tool in ["scUTRquant", "SCAPTURE", "polyApipe", "Sierra", "scAPAtrap"]:
        q = sub[sub.tool == tool]
        if q.empty:
            continue
        c = COLOR[tool]
        hollow = tool == "scUTRquant"
        mk = dict(marker="o", ms=6.5, mec=c, mew=1.3, mfc="white" if hollow else c,
                  ls="none", zorder=5)
        if paired and len(q) == 2:
            x, y = q.R_det.mean(), q.P.mean()
            ax.errorbar(x, y, xerr=[[x - q.R_det.min()], [q.R_det.max() - x]],
                        yerr=[[y - q.P.min()], [q.P.max() - y]],
                        color=c, lw=0.8, capsize=1.5, zorder=4, ls="none")
            h, = ax.plot(x, y, **mk)
        else:
            x, y = q.R_det.iloc[0], q.P.iloc[0]
            h, = ax.plot(x, y, **mk)
        nm = "scUTRquant* (catalog)" if hollow else tool
        if paired and len(q) == 1:
            nm += " (m1 only)"
        handles.append((h, nm))
        # null markers at the bottom, same x
        ax.plot(x, q.null_P_mean.mean(), marker="|", ms=5, mec=c, mew=1.0, ls="none", zorder=3)
        ref_lines.append(dict(panel=ax.get_label(), kind="null_genic mean (3 seeds)", label=tool,
                              P=q.null_P_mean.mean(), F1=""))
    return handles


def draw_path(ax, sub, roles, paired, labels):
    c = COLOR["PeakATail"]
    xs, ys = [], []
    for i, role in enumerate(roles):
        q = sub[sub.role == role]
        if paired:
            x, y = q.R_det.mean(), q.P.mean()
            ax.errorbar(x, y, xerr=[[x - q.R_det.min()], [q.R_det.max() - x]],
                        yerr=[[y - q.P.min()], [q.P.max() - y]],
                        color=c, lw=0.8, capsize=1.5, zorder=6, ls="none")
        else:
            x, y = q.R_det.iloc[0], q.P.iloc[0]
        xs.append(x); ys.append(y)
        last = role == "path4"
        first = role == "path0"
        ax.plot(x, y, marker="D" if last else "o", ms=8 if last else 5.5,
                mfc=c if last else ("white" if first else c),
                mec=c, mew=1.3, ls="none", zorder=7)
        ax.plot(x, q.null_P_mean.mean(), marker="|", ms=5, mec=c, mew=1.0, ls="none", zorder=3)
        ref_lines.append(dict(panel=ax.get_label(), kind="null_genic mean (3 seeds)",
                              label=f"PeakATail {q.short.iloc[0]}".replace("\n", " "),
                              P=q.null_P_mean.mean(), F1=""))
        lab, dx, dy, ha, va = labels[role]
        ax.annotate(lab, (x, y), xytext=(dx, dy), textcoords="offset points",
                    fontsize=6.2, color=c, ha=ha, va=va, zorder=8,
                    fontweight="bold" if last else "normal")
    ax.plot(xs, ys, color=c, lw=1.1, alpha=0.75, zorder=6)
    # arrow heads along the path
    for (x0, y0), (x1, y1) in zip(zip(xs[:-1], ys[:-1]), zip(xs[1:], ys[1:])):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", color=c, lw=0.9, alpha=0.75,
                                    shrinkA=4, shrinkB=5), zorder=6)


# ---- panel a -------------------------------------------------------------
axA.set_label("a")
subH = pts[pts.dataset == HUMAN]
XMAX_H = 0.42
f1_iso(axA, XMAX_H)
gates(axA, XMAX_H)
hH = draw_competitors(axA, subH[subH.role.isin(["competitor", "catalog"])], paired=False)
labels_H = {
    "path0": ("shipped\n(pre-fix)", -4, 7, "right", "bottom"),
    "path1": ("both tiers (no IP)", 0, -8, "center", "top"),
    "path2": ("tier-1 ≥1 mol\n(no IP)", 6, 0, "left", "center"),
    "path3": ("tier-1 ≥1 mol (IP)", 7, 13, "left", "center"),
    "path4": ("precision default\n(tier-1, IP, ≥2 mol)", 9, 0, "left", "center"),
}
draw_path(axA, subH, PATH_ROLES_H, paired=False, labels=labels_H)
axA.set_xlim(0, XMAX_H); axA.set_ylim(0, 0.9)
axA.set_xlabel("detected-gene recall R_det @100 bp (n = 285,136 atlas sites)")
axA.set_ylabel("atlas-agreement precision @100 bp")
axA.set_title("a   PBMC 10k v3 (human, single donor)", loc="left", fontweight="bold")
axA.text(0.985, 0.02, "scTail: not runnable on this BAM (R1 = 28 bp)\n"
         "| = 3-seed genic-shuffle null (mean) per call set",
         transform=axA.transAxes, ha="right", va="bottom", fontsize=5.8, color=MUTED)
hH.append((Line2D([], [], marker="o", ms=5.5, color=COLOR["PeakATail"], mec=COLOR["PeakATail"], lw=1.1),
           "PeakATail operating points (final run, arrows = filter steps)"))
hH.append((Line2D([], [], marker="D", ms=7, color="none", mfc=COLOR["PeakATail"], mec=COLOR["PeakATail"]),
           "PeakATail pre-registered precision default"))
style(axA)

# ---- panel b -------------------------------------------------------------
axB.set_label("b")
subM = pts[pts.dataset == MOUSE]
XMAX_M = 0.42
f1_iso(axB, XMAX_M)
gates(axB, XMAX_M, original=True, original_note=" — PBMC reference")
hM = draw_competitors(axB, subM[subM.role.isin(["competitor", "catalog"])], paired=True)
labels_M = {
    "path0": ("shipped (pre-fix)", 8, -13, "left", "top"),
    "path1": ("both tiers (IP)", 7, -2, "left", "top"),
    "path3": ("tier-1 ≥1 mol (IP)", 7, 0, "left", "center"),
    "path4": ("precision default\n(tier-1, IP, ≥2 mol)", 9, 0, "left", "center"),
}
draw_path(axB, subM, PATH_ROLES_M, paired=True, labels=labels_M)
axB.set_xlim(0, XMAX_M); axB.set_ylim(0, 0.9)
axB.set_xlabel("detected-gene recall R_det @100 bp (n = 126,686)")
axB.set_ylabel("atlas-agreement precision @100 bp")
axB.set_title("b   GSE104556 testis (2 mice, mean ± range)", loc="left", fontweight="bold")
axB.text(0.985, 0.02, "all mouse arms ran with the IP filter (no ‘no IP’ step)\n"
         "SCAPTURE: mouse 1 only (mouse-2 run invalid)\n"
         "| = 3-seed genic-shuffle null (mean) per call set",
         transform=axB.transAxes, ha="right", va="bottom", fontsize=5.8, color=MUTED)
style(axB)
fig.legend([h for h, _ in hH], [n for _, n in hH], loc="center", bbox_to_anchor=(0.56, 0.528),
           ncol=4, frameon=False, handletextpad=0.4, columnspacing=1.4, labelspacing=0.4)

# ---- panel c: n sites per call set (log) ---------------------------------
axC.set_label("c")
order = [  # (tool, short label in pts) top -> bottom
    ("scUTRquant", "scUTRquant*"), ("SCAPTURE", "SCAPTURE"), ("polyApipe", "polyApipe"),
    ("Sierra", "Sierra"), ("scAPAtrap", "scAPAtrap"),
    ("PeakATail", "shipped"), ("PeakATail", "both tiers\n(no IP)"), ("PeakATail", "tier-1 >=1 mol\n(no IP)"),
    ("PeakATail", "both tiers (IP)"), ("PeakATail", "tier-1 >=1 mol (IP)"), ("PeakATail", "tier-2 (IP)"),
    ("PeakATail", "precision default\n(tier-1, IP, >=2 mol)"),
]
ylabels = {  # PT = PeakATail (final run)
    "scUTRquant*": "scUTRquant* (catalog)", "shipped": "PT shipped (pre-fix)",
    "both tiers\n(no IP)": "PT both tiers (no IP)", "tier-1 >=1 mol\n(no IP)": "PT tier-1 ≥1 mol (no IP)",
    "both tiers (IP)": "PT both tiers (IP)", "tier-1 >=1 mol (IP)": "PT tier-1 ≥1 mol (IP)",
    "tier-2 (IP)": "PT tier-2 only (IP)",
    "precision default\n(tier-1, IP, >=2 mol)": "PT precision default",
}
bh = 0.26
yt, ytl = [], []
for i, (tool, short) in enumerate(order):
    y = len(order) - 1 - i
    c = COLOR[tool]
    qh = pts[(pts.dataset == HUMAN) & (pts.tool == tool) & (pts.short == short)]
    q1 = pts[(pts.dataset == MOUSE) & (pts.tool == tool) & (pts.short == short) & (pts.replicate == "mouse1")]
    q2 = pts[(pts.dataset == MOUSE) & (pts.tool == tool) & (pts.short == short) & (pts.replicate == "mouse2")]
    hollow = tool == "scUTRquant"
    kw = dict(height=bh, edgecolor=c, linewidth=0.8, zorder=3)
    if not qh.empty:
        axC.barh(y + bh, qh.n.iloc[0], color="white" if hollow else c, **kw)
        axC.text(qh.n.iloc[0] * 1.12, y + bh, f"{int(qh.n.iloc[0]):,}", va="center", fontsize=5.3, color=INK)
    if not q1.empty:
        axC.barh(y, q1.n.iloc[0], color="white" if hollow else c, alpha=1 if hollow else 0.55, hatch="////", **kw)
        axC.text(q1.n.iloc[0] * 1.12, y, f"{int(q1.n.iloc[0]):,}", va="center", fontsize=5.3, color=INK)
    if not q2.empty:
        axC.barh(y - bh, q2.n.iloc[0], color="white" if hollow else c, alpha=1 if hollow else 0.3, hatch="....", **kw)
        axC.text(q2.n.iloc[0] * 1.12, y - bh, f"{int(q2.n.iloc[0]):,}", va="center", fontsize=5.3, color=INK)
    elif not q1.empty:
        axC.text(1.2e4, y - bh, "m2: no valid run", va="center", fontsize=5.3, color=MUTED)
    elif qh.empty is False and q1.empty:
        axC.text(1.2e4, y, "not run on mouse", va="center", fontsize=5.3, color=MUTED)
    yt.append(y); ytl.append(ylabels.get(short, short))
axC.set_yticks(yt); axC.set_yticklabels(ytl, fontsize=5.9)
axC.set_xscale("log"); axC.set_xlim(1e4, 6e6)
axC.set_xlabel("n sites scored per call set (log scale)")
axC.set_title("c   Call-set size", loc="left", fontweight="bold")
axC.legend([matplotlib.patches.Patch(facecolor="#888888", edgecolor="#888888"),
            matplotlib.patches.Patch(facecolor="#888888", alpha=0.55, hatch="////", edgecolor="#888888"),
            matplotlib.patches.Patch(facecolor="#888888", alpha=0.3, hatch="....", edgecolor="#888888")],
           ["PBMC 10k v3", "testis mouse 1", "testis mouse 2"], loc="lower right", frameon=False,
           handlelength=1.4, labelspacing=0.3, fontsize=6.0)
style(axC)
axC.grid(True, axis="x", color=GRID, lw=0.5); axC.grid(False, axis="y")
axC.set_ylim(-0.6, len(order) - 0.4)

# ---- panel d: tier decomposition -----------------------------------------
axD.set_label("d")
c = COLOR["PeakATail"]
slices = [
    ("tier-2\n(coverage-only)", int(t2.n), "#C9D1D5", t2.P, "tier-2 only"),
    ("tier-1\nsingle-molecule", n_single, "#8FBBDD", t1.P, "tier-1 ≥1 mol (all)"),
    ("tier-1 ≥2 molecules\n= precision default", int(dflt.n), c, dflt.P, "default"),
]
left = 0
for k, (lab, n, col, P, _) in enumerate(slices):
    axD.barh(1.0, n, left=left, height=0.55, color=col, edgecolor="white", linewidth=1.0, zorder=3)
    if k < 2:
        axD.text(left + n / 2, 1.0, f"{lab}\nn = {n:,}", ha="center", va="center", fontsize=5.7,
                 color=INK, zorder=4)
    else:  # too narrow for an in-bar label: put it above the slice, right-aligned
        axD.annotate(f"{lab}\nn = {n:,}", xy=(left + n / 2, 1.28), xytext=(left + n, 1.42),
                     ha="right", va="bottom", fontsize=5.7, color=col, zorder=4,
                     arrowprops=dict(arrowstyle="-", color=col, lw=0.6))
    left += n
axD.text(0, 1.36, f"PBMC IP arm, all calls scored:\nn = {int(both.n):,}  (P@100 {both.P:.3f})",
         ha="left", va="bottom", fontsize=6.0, color=MUTED)
# precision of each scored output, as a dot strip below
outs = [
    ("tier-2 only", t2.P, t2.n, "#C9D1D5"),
    ("tier-1 ≥1 mol\n(singletons + ≥2 mol)", t1.P, t1.n, "#8FBBDD"),
    ("precision default\n(tier-1 ≥2 mol)", dflt.P, dflt.n, c),
]
x_scale = both.n  # map P in [0,1] onto the bar width so the two rows share an x axis
for k, (lab, P, n, col) in enumerate(outs):
    yk = 0.30 - 0.22 * k
    axD.plot([0, P * x_scale], [yk, yk], color=col, lw=3.2, solid_capstyle="butt", zorder=3)
    axD.text(P * x_scale + x_scale * 0.012, yk, f"{P:.3f}  {lab.split(chr(10))[0]}",
             va="center", ha="left", fontsize=5.8, color=INK, zorder=5,
             bbox=dict(facecolor="white", edgecolor="none", pad=0.6))
axD.plot([GATE_P * x_scale] * 2, [-0.28, 0.42], color=INK, lw=0.8, ls=(0, (4, 2)), zorder=2)
axD.text(GATE_P * x_scale, 0.44, "gate 0.50", ha="center", va="bottom", fontsize=5.6, color=INK)
axD.text(0, -0.46, "atlas-agreement precision @100 bp of each scored output\n"
         "(x axis: 0 → 1 spans the width of the top bar)",
         fontsize=5.7, color=MUTED, va="top")
axD.text(0, 1.95, textwrap.fill(f"{frac_single*100:.1f}% of tier-1 sites are single-molecule ({n_single:,} / {int(t1.n):,}) "
         f"and are dropped by the default; mouse {m_single['mouse1'][1]*100:.0f}% / {m_single['mouse2'][1]*100:.0f}%.", 62),
         fontsize=6.0, color=INK, va="bottom")
axD.set_xlim(0, both.n * 1.30); axD.set_ylim(-0.80, 2.55)
axD.set_yticks([]); axD.set_xticks([0, 1e5, 2e5, 3e5])
axD.set_xticklabels(["0", "100k", "200k", "300k"])
axD.set_xlabel("n sites (top bar)")
axD.set_title("d   PeakATail tier decomposition, PBMC IP arm", loc="left", fontweight="bold")
for s in ("top", "right", "left"):
    axD.spines[s].set_visible(False)
axD.spines["bottom"].set_color(GRID)
axD.tick_params(colors=MUTED, length=2.5, labelsize=7)

# ---- footer caption -------------------------------------------------------
pA_pp = pick(MOUSE, "path4", "mouse1")
poly_h = subH[subH.tool == "polyApipe"].iloc[0]
poly_m = subM[subM.tool == "polyApipe"]
rel_h = (dflt.R_det / poly_h.R_det - 1) * 100
rel_m = (subM[subM.role == "path4"].R_det.mean() / poly_m.R_det.mean() - 1) * 100
dF1_h = dflt.F1_det - poly_h.F1_det
dF1_m = subM[subM.role == "path4"].F1_det.mean() - poly_m.F1_det.mean()
m_full = subM[subM.role == "path4"].sort_values("replicate").R_full.tolist()

_prev = VERSIONS["prev"]
prev_note = "" if not _prev else (
    f"The {_prev['label']} gave {_prev['P'][0]:.3f} / {_prev['P'][1]:.3f} / {_prev['P'][2]:.3f}; the "
    f"\u22121.1 pp PBMC move is the corrected minus-strand IP window (#96), reconstructed key-for-key. ")
_ip = VERSIONS["ip_trade"]
ip_note = _ip[0].upper() + _ip[1:]

caption = (
    "Fig 2 | Trustworthy PAS detection: the final Stage-2 run (PeakATail code "
    f"{CODE_COMMIT}, {VERSIONS['commit_note']}) against the competitor panel. "
    "Definitions: atlas-agreement precision = fraction of calls within 100 bp (point, strand-matched) of a "
    "PolyASite 2.0 representative site; it is agreement with a curated atlas, not ground truth (atlas-novel "
    f"true sites count as false positives). R_det denominator = atlas sites in genes detected in the dataset "
    f"({int(dflt.R_det_denominator):,} human / {int(_m1.R_det_denominator):,} mouse); full-atlas recall of the "
    f"default {dflt.R_full:.3f} (PBMC; {int(dflt.R_full_denominator):,}) and {m_full[0]:.3f} / {m_full[1]:.3f} "
    f"(mice; {int(_m1.R_full_denominator):,}). One scorer (score_tool.py), same denominators for every tool; "
    "3-seed gene-body-shuffled nulls. scUTRquant is catalog-based (hollow, not ranked). "
    f"Pre-registration: the precision-first default and its P@100 ≥ 0.50 gate were committed before the "
    f"ORIGINAL run and are unchanged here ({PREREG}; 13 §1, {GATE_DOC}); the ≥2-molecule threshold was "
    "pre-registered before that run (12 CORRECTION). "
    f"Gate: PASS on all three ({dflt.P:.3f} / {_m1.P:.3f} / {_m2.P:.3f}); no ≥1-molecule output clears the "
    "original P ≥ 0.38 floor. "
    + prev_note +
    f"{ip_note}; compute {VERSIONS['compute']}; {GATE_DOC} §1/§3. "
    f"Caveats: single PBMC donor / BAM; the two mice are one study and chemistry; the default's recall is at "
    f"or below polyApipe's ({rel_h:+.0f}% PBMC, {rel_m:+.0f}% mouse) and the F1_det lead "
    f"({dF1_h:+.3f} / {dF1_m:+.3f}) is not a meaningful margin — the claim is precision, not overall accuracy. "
    "The non-IP PBMC ≥2-molecule file (pas_tier1_ge2mol_noIP_POSTHOC) is not the default and is not shown. "
    f"Verified: {VERSIONS['verified_short']}. Every plotted value: "
    "results/figures/manuscript/final_benchmark*.tsv."
)
_cap_wrapped = textwrap.fill(caption, 165)
_n_cap_lines = _cap_wrapped.count("\n") + 1
# approved layout: the footer block must stay clear of the panel c/d x-axis labels
assert _n_cap_lines <= 12, f"caption is {_n_cap_lines} wrapped lines; >12 collides with the panel x labels"
fig.text(0.03, 0.012, _cap_wrapped, fontsize=5.8, color=INK, va="bottom", ha="left",
         linespacing=1.35)

for ext, kw in (("png", dict(dpi=300)), ("pdf", {})):
    p = FIGDIR / f"{NAME}.{ext}"
    fig.savefig(p, **kw)
    print("wrote", p)
p = OUTDIR / f"{NAME}.png"
fig.savefig(p, dpi=300); print("wrote", p)

# ---------------------------------------------------------------------------
# TSVs: every plotted value
# ---------------------------------------------------------------------------
pts_out = pts.copy()
pts_out["short"] = pts_out["short"].str.replace("\n", " ")
pts_out["panel"] = pts_out.apply(
    lambda r: ("a" if r.dataset == HUMAN else "b") + ",c" if r.role != "extra" else "c" + (",d" if r.dataset == HUMAN else ""),
    axis=1)
cols = ["panel", "dataset", "tool", "call_set", "short", "replicate", "role", "cutoff_bp", "n", "n_matched",
        "P", "R_det", "R_det_denominator", "R_full", "R_full_denominator", "F1_det",
        "null_P_seed1", "null_P_seed2", "null_P_seed3", "null_P_mean", "source"]
pts_out = pts_out[cols].rename(columns={"P": "atlas_agreement_precision_100bp",
                                        "R_det": "recall_detected_genes_100bp",
                                        "R_full": "recall_full_atlas_100bp",
                                        "F1_det": "F1_detected_genes_100bp"})
p = OUTDIR / f"{NAME}.tsv"; pts_out.to_csv(p, sep="\t", index=False, float_format="%.6f"); print("wrote", p)
p = OUTDIR / f"{NAME}_tiers.tsv"; tiers.to_csv(p, sep="\t", index=False, float_format="%.6f"); print("wrote", p)
p = OUTDIR / f"{NAME}_reference_lines.tsv"
pd.DataFrame(ref_lines).to_csv(p, sep="\t", index=False, float_format="%.6f"); print("wrote", p)

# ---------------------------------------------------------------------------
# sidecar caption
# ---------------------------------------------------------------------------
cap_md = f"""# Fig 2 — `final_benchmark` caption (generated by `scripts/manuscript_figures/final_benchmark.py`)

{caption}

**Panel a** — PBMC 10k v3: atlas-agreement precision @100 bp (y) vs detected-gene recall R_det @100 bp (x) for
every tool; PeakATail as a connected path of operating points (shipped → both tiers, no IP → tier-1 ≥1 mol, no IP
→ tier-1 ≥1 mol, IP → precision default = tier-1 ∩ IP-pass ∩ ≥2 distinct clip molecules). Light-grey F1_det
isolines 0.1–0.4; dashed line = pre-registered gate P ≥ 0.50 (13 §1); dotted line + dotted isoline = original
two-sided gate P ≥ 0.38 & F1_det > 0.261 (10 §5). Small ticks at the bottom = mean of the 3-seed gene-body-shuffled
null for each call set. scTail is absent: not runnable on this BAM (R1 = 28 bp).
**Panel b** — GSE104556 testis, two mice as mean ± range (both mice in the TSV); all mouse arms ran with the IP
filter so the path has no "no IP" step; SCAPTURE is mouse 1 only. **Panel c** — n sites scored per call set, log
scale. **Panel d** — PeakATail tier decomposition on the PBMC IP arm: tier-2 / tier-1 single-molecule / tier-1 ≥2
molecules (= default), with the atlas-agreement precision @100 bp of each scored output (tier-2 {t2.P:.3f},
tier-1 ≥1 mol {t1.P:.3f}, default {dflt.P:.3f}); {frac_single*100:.1f}% of tier-1 sites are single-molecule
({n_single:,} / {int(t1.n):,}; mouse {m_single['mouse1'][1]*100:.1f}% / {m_single['mouse2'][1]*100:.1f}%).

Sources: `{VERSIONS["gate_doc_path"]}` + `{VERSIONS["verified_table"]}` ({VERSIONS["verdict"]}) and the score TSVs listed in the `source` column of
`results/figures/manuscript/final_benchmark.tsv`. Reference lines: `final_benchmark_reference_lines.tsv`;
panel d: `final_benchmark_tiers.tsv`.
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
      f"tier-1 single-molecule {frac_single*100:.1f}%")

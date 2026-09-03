#!/usr/bin/env python3
"""
figS11_gatehistory.py -- manuscript Fig S11 ("figS11_gatehistory"): the
pre-registration timeline as a figure, plus every gate outcome on every arm
(the T3 gate table drawn). This is the figure that makes the pre-registration
AUDITABLE rather than asserted (21 section 4.3 row S9 of the historical plan;
FIGURES_MANIFEST row S11).

SOURCES OF TRUTH (every date, sha and clock time from the named documents)
  * 10 section 5 / "Definition of done" (2026-08-19): the ORIGINAL two-sided
    gate -- P@100 >= 0.38 AND F1_det > 0.261 on PBMC.
  * 12 (run 2026-08-20) + its CORRECTION (2026-08-21): the stale Stage-2 FAIL
    ("under the full pre-registered gate, NO arm passes"), and the post-hoc
    >=2-clip sweep -- directory created 12 minutes after the gate FAIL was
    recorded; "It is post hoc, full stop."
  * 13 section 1 + 15 (verified from git): the precision-first default and its
    P@100 >= 0.50 gate committed in `0e27b1a` at 2026-08-21 01:18:59; the
    final arms started 02:59:36 (code 4efeb125; run 02:59-06:46; SOUND).
  * 19 (run 2026-08-21 16:14-16:50, code 9dfdefb): the v2 re-run against
    UNCHANGED gates (verifier FIXED).
  * 24 (written 2026-08-21, last write 23:34 per 26's cross-reference) + its
    Amendment 3 (2026-08-23): the prime adoption criterion, and its FAIL at
    delta 0.000000 -- an identity, not a regression.
  * 26 (written 2026-08-21 23:47-23:53; driver frozen 23:48:30; run
    2026-08-22 00:07:45-00:23:26; verifier FIXED 2026-08-22): the second-donor
    test, gate PASS 0.8279.
  * 27 (written 2026-08-21; clock time not recorded in the document): the
    3'UTR singleton promotion pre-registration -- confirmatory evaluation NOT
    RUN; PENDING; ships flag-off until its section 3 passes.

NUMBER POLICY -- gate VALUES are read programmatically where score TSVs exist
(v1 arms from results/benchmark_tools/.../score_*.tsv, v2 from the verified
fig2_accuracy.tsv, pbmc4k from its run's score TSV) and EXPECT-asserted
against the documents; Stage-2 (2026-08-20) values exist only in 12 (its
count matrices were voided, its BEDs superseded) and are cited constants.
Every timestamp is quoted from a named document; where a document records a
date but no clock time, the marker's position is a DISPLAY position and the
audit TSV says so.

OUTPUTS
  manuscript/figures/figS11_gatehistory.{png,pdf}  (PNG 600 dpi, PDF fonttype 42)
  manuscript/figures/figS11_gatehistory.caption.md (single source of the legend;
      the image itself carries only the timeline/table data, no caption prose)
  results/figures/manuscript/figS11_gatehistory.tsv        (timeline events)
  results/figures/manuscript/figS11_gatehistory_gates.tsv  (panel b: T3 rows)

Run:  export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/figS11_gatehistory.py
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
import textwrap
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
BT = WD / "results/benchmark_tools"
OUTDIR = WD / "results/figures/manuscript"
FIGDIR = WD / "manuscript/figures"
NAME = "figS11_gatehistory"
OUTDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)

INK, MUTED, GRID = "#1B2429", "#5A6B73", "#D8E0E3"
C_PREREG, C_PASS, C_FAIL, C_IDENT, C_PEND = "#0072B2", "#009E73", "#D55E00", "#CC79A7", "#999999"
BOXFACE = "#F7FAFB"

# ---------------------------------------------------------------------------
# programmatic gate values + EXPECT asserts against the documents
# ---------------------------------------------------------------------------
def score100(tsv: Path):
    df = pd.read_csv(tsv, sep="\t")
    at = df[df.cutoff_bp == 100.0]
    P = float(at[(at.panel == "precision") & (at.series == "real") & (at.reference == "atlas_full")].value.iloc[0])
    F1 = float(at[(at.panel == "f1") & (at.series == "real") & (at.reference == "atlas_detected")].value.iloc[0])
    n = int(at[(at.panel == "precision") & (at.series == "real") & (at.reference == "atlas_full")].n_query.iloc[0])
    return n, P, F1


# v1 (code 4efeb125; 15 sections 1-2)
v1h = score100(BT / "pbmc_10k_v3/peakatail_clipseeded_final_ipfilt/score_pbmc_final_ipfilt__PRESPEC_precision_default.tsv")
v1m1 = score100(BT / "gse104556/peakatail_clipseeded_final/mouse1/score_testis_m1_final__PRESPEC_precision_default.tsv")
v1m2 = score100(BT / "gse104556/peakatail_clipseeded_final/mouse2/score_testis_m2_final__PRESPEC_precision_default.tsv")
v1t1_noip = score100(BT / "pbmc_10k_v3/peakatail_clipseeded_final/score_pbmc_final__tier1.tsv")
v1t1_ip = score100(BT / "pbmc_10k_v3/peakatail_clipseeded_final_ipfilt/score_pbmc_final_ipfilt__tier1.tsv")
assert (v1h[0], round(v1h[1], 4)) == (44394, 0.7167), v1h          # 15 section 1
assert (v1m1[0], round(v1m1[1], 4)) == (25991, 0.7414), v1m1
assert (v1m2[0], round(v1m2[1], 4)) == (26164, 0.7554), v1m2
assert (v1t1_noip[0], round(v1t1_noip[1], 4), round(v1t1_noip[2], 4)) == (222955, 0.3032, 0.3001), v1t1_noip
assert (v1t1_ip[0], round(v1t1_ip[1], 4), round(v1t1_ip[2], 4)) == (161595, 0.3546, 0.3013), v1t1_ip

# v2 (code 9dfdefb; 19) from the verified Fig 2 audit TSV
f2 = pd.read_csv(OUTDIR / "fig2_accuracy.tsv", sep="\t")


def f2v(dataset, role, rep="single"):
    q = f2[(f2.dataset == dataset) & (f2.role == role) & (f2.replicate == rep)]
    assert len(q) == 1
    r = q.iloc[0]
    return int(r.n), float(r.atlas_agreement_precision_100bp), float(r.F1_detected_genes_100bp)


v2h, v2m1, v2m2 = f2v("pbmc_10k_v3", "path4"), f2v("gse104556", "path4", "mouse1"), f2v("gse104556", "path4", "mouse2")
v2t1_noip, v2t1_ip = f2v("pbmc_10k_v3", "path2"), f2v("pbmc_10k_v3", "path3")
assert (v2h[0], round(v2h[1], 4)) == (46524, 0.7062), v2h          # 19 section 1
assert (v2m1[0], round(v2m1[1], 4)) == (26255, 0.7450), v2m1
assert (v2m2[0], round(v2m2[1], 4)) == (26526, 0.7572), v2m2
assert (v2t1_noip[0], round(v2t1_noip[1], 4), round(v2t1_noip[2], 4)) == (222955, 0.3032, 0.3001), v2t1_noip
assert (v2t1_ip[0], round(v2t1_ip[1], 4), round(v2t1_ip[2], 4)) == (167565, 0.3520, 0.3046), v2t1_ip

# pbmc4k (26) from the run's own score TSV
d2 = score100(BT / "pbmc4k_donor2/run_ipfilt/score_pbmc4k_final_v2_ipfilt__PRESPEC_precision_default.tsv")
assert (d2[0], round(d2[1], 4)) == (20672, 0.8279), d2             # 26 R2/R4

# Stage-2 (2026-08-20): cited constants -- the run's matrices were voided and
# its call sets superseded, so 12 is the only source (verified to 6 decimals
# there by its own verifier).
S2 = dict(both=(377236, 0.188, 0.235), t1=(197410, 0.308, 0.290), t1ip=(142400, 0.362, 0.290))

# gate arithmetic, asserted rather than narrated
ORIG_P, ORIG_F1, GATE_P = 0.38, 0.261, 0.50
assert S2["t1ip"][1] < ORIG_P and S2["t1"][1] < ORIG_P and S2["both"][2] < ORIG_F1
for n, P, F1 in (v1t1_noip, v1t1_ip, v2t1_noip, v2t1_ip):
    assert P < ORIG_P and F1 > ORIG_F1        # F1 pass / P FAIL, every time
for n, P, F1 in (v1h, v1m1, v1m2, v2h, v2m1, v2m2, d2):
    assert P >= GATE_P                         # precision gate PASS, every time

# the pre-registration margin: default committed 1:40:37 before the arms started
_commit_s = 1 * 3600 + 18 * 60 + 59
_arms_s = 2 * 3600 + 59 * 60 + 36
_delta = _arms_s - _commit_s
assert _delta == 6037                          # 1 h 40 m 37 s
DELTA_STR = "1:40:37"

# ---------------------------------------------------------------------------
# timeline events. t = days since 2026-08-19 00:00 (all clock times from docs;
# display positions flagged where the doc records a date but no clock time)
# ---------------------------------------------------------------------------
def T(day, h=0, m=0, s=0):
    return day + (h + m / 60 + s / 3600) / 24


EVENTS = [
    dict(id="orig_gate", t=T(0, 12), t_end=None, kind="prereg",
         ts="2026-08-19 (date; no clock time in 10)", display_note="display position 12:00",
         box=(0.52, 2.05, 0.98),
         text="ORIGINAL two-sided gate\npre-registered (10 §5):\nP@100 ≥ 0.38 AND\nF1_det > 0.261 (PBMC)",
         verdict="", source="manuscript/10_caller_fix_plan.md section 5 / 'Definition of done' (dated 2026-08-19)"),
    dict(id="stage2_fail", t=T(1, 12), t_end=None, kind="fail",
         ts="2026-08-20 (date; no clock time in 12)", display_note="display position 12:00",
         box=(1.28, -2.0, 0.98),
         text="stale Stage-2 PBMC run (12):\nunder the FULL gate\nNO arm passes\n(best P 0.362 < 0.38) — FAIL",
         verdict="FAIL", source="manuscript/12_stage2_gate.md (run 2026-08-20; CORRECTED 2026-08-21)"),
    dict(id="sweep", t=T(2, 0, 30), t_end=None, kind="disclosure",
         ts="sweep dir created 12 min after the FAIL was recorded (12 CORRECTION §1); section added 2026-08-21",
         display_note="display position 00:30 — clock not recorded; ordered before the 01:18:59 commit per 15",
         box=(1.82, -3.55, 1.14),
         text="post-hoc ≥2-clip sweep — created 12 min\nafter the FAIL was recorded; 'post hoc,\nfull stop' (12 CORRECTION). Informed the\nthreshold; DISCLOSED, never a result",
         verdict="", source="manuscript/12_stage2_gate.md 'Post-hoc sensitivity' + CORRECTION section 1"),
    dict(id="default_commit", t=T(2, 1, 18, 59), t_end=None, kind="prereg",
         ts="2026-08-21 01:18:59 (verified from git, 15)", display_note="",
         box=(1.55, 2.35, 1.06),
         text="precision-first default\n(tier-1 ∩ IP ∩ ≥2 mol) + gate\nP@100 ≥ 0.50 ×3 committed\n`0e27b1a` 01:18:59 (13 §1) —\narms started 1:40:37 LATER",
         verdict="", source="manuscript/13_reliability_positioning.md section 1; manuscript/15_final_gate.md (git-verified)"),
    dict(id="v1_run", t=T(2, 2, 59, 36), t_end=T(2, 6, 46), kind="run",
         ts="2026-08-21 02:59:36 - 06:46", display_note="",
         box=(2.42, -1.95, 1.02),
         text="v1 final arms, code 4efeb125\n(15, SOUND): 0.50 gate PASS ×3\n(0.7167/0.7414/0.7554);\noriginal 0.38 floor FAIL (all\n≥1-mol arms; F1 half passes)",
         verdict="PASS+FAIL", source="manuscript/15_final_gate.md sections 1-2"),
    dict(id="v2_run", t=T(2, 16, 14), t_end=T(2, 16, 50), kind="run",
         ts="2026-08-21 16:14 - 16:50", display_note="",
         box=(3.10, -3.5, 1.08),
         text="v2 re-run, code 9dfdefb (#96+#97),\ngates UNCHANGED (19, FIXED):\n0.50 gate PASS ×3 (0.7062/0.7450/\n0.7572); 0.38 floor still FAIL on P",
         verdict="PASS+FAIL", source="manuscript/19_final_gate_v2.md sections 1-2"),
    dict(id="utr_prereg", t=T(2, 20), t_end=None, kind="prereg",
         ts="2026-08-21 (clock time not recorded in 27)", display_note="display position 20:00",
         box=(2.52, 3.55, 1.06),
         text="3'UTR singleton promotion\npre-registered (27): rule, arms,\nfalsifiers fixed; discovery data\ndisqualified — still UNTESTED",
         verdict="PENDING", source="manuscript/27_utr_singleton_prereg.md (status header)"),
    dict(id="prime_prereg", t=T(2, 23, 34), t_end=None, kind="prereg",
         ts="2026-08-21, last write 23:34 (26's cross-reference)", display_note="",
         box=(3.42, 2.3, 1.06),
         text="prime adoption criterion\npre-registered (24 §3.1):\nΔP ≥ −0.005 & ΔR ≥ +0.010\n& ΔF1 > 0, on all three",
         verdict="", source="manuscript/24_prime_preregistration.md; write time per manuscript/26 preamble"),
    dict(id="donor2_prereg", t=T(2, 23, 47), t_end=None, kind="prereg",
         ts="2026-08-21 23:47-23:53; driver frozen 23:48:30", display_note="",
         box=(4.02, 3.6, 1.10),
         text="second-donor test pre-registered\n(26): same gate/scorer/nulls,\nnothing retunable — BEFORE any\npbmc4k number existed",
         verdict="", source="manuscript/26_second_donor_preregistration.md (preamble + A1 + verifier record)"),
    dict(id="donor2_run", t=T(3, 0, 7, 45), t_end=T(3, 0, 23, 26), kind="run",
         ts="2026-08-22 00:07:45 - 00:23:26; verifier FIXED same day", display_note="",
         box=(4.00, -1.95, 1.06),
         text="pbmc4k run, frozen 9dfdefb3\n(26): gate P@100 ≥ 0.50\nPASS at 0.8279 (38× null) —\n4th library; verifier FIXED",
         verdict="PASS", source="manuscript/26_second_donor_preregistration.md R1-R4 + V"),
    dict(id="amendment3", t=T(4, 12), t_end=None, kind="identity",
         ts="2026-08-23 (date; no clock time in 24)", display_note="display position 12:00",
         box=(4.42, 1.85, 1.10),
         text="Amendment 3 (24): prime default\nbyte-identical to v2's arm →\ncriterion FAILS at Δ 0.000000 —\nan identity, not a regression;\ncriterion targeted a v2 default\nthat never existed",
         verdict="FAIL (identity)", source="manuscript/24_prime_preregistration.md Amendment 3 (2026-08-23)"),
]

KIND_COLOR = dict(prereg=C_PREREG, run=INK, fail=C_FAIL, disclosure=MUTED, identity=C_IDENT)
VERD_COLOR = {"PASS": C_PASS, "FAIL": C_FAIL, "PASS+FAIL": C_PASS, "PENDING": C_PEND,
              "FAIL (identity)": C_IDENT, "": C_PREREG}

# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
# canvas: the footer caption band (bottom 0.135 of an 11.4 in canvas) moved to
# the legend sidecar; height shrinks by the freed space, both panels keep size
fig = plt.figure(figsize=(8.7, 10.15))
gsA = fig.add_gridspec(1, 1, left=0.025, right=0.985, top=0.961, bottom=0.585)
gsB = fig.add_gridspec(1, 1, left=0.025, right=0.985, top=0.523, bottom=0.029)
axA = fig.add_subplot(gsA[0])
axB = fig.add_subplot(gsB[0])

axA.set_xlim(-0.05, 5.05)
axA.set_ylim(-4.6, 4.6)
axA.axis("off")
axA.set_title(sc_title("a   the pre-registration timeline — every gate fixed "
                       "before the number it judges"),
              loc="left", fontweight="bold", fontsize=TYPE["panel_title"])

# day bands + labels
DAYS = ["Aug 19", "Aug 20", "Aug 21", "Aug 22", "Aug 23"]
for d in range(5):
    if d % 2 == 0:
        axA.axvspan(d, d + 1, color="#F0F4F6", zorder=0)
    axA.text(d + 0.5, -4.45, f"{DAYS[d]}, 2026", ha="center", va="bottom", fontsize=TYPE["tick"],
             color=MUTED, fontweight="bold")
axA.axhline(0, color=INK, lw=1.2, zorder=2)

# hour ticks on Aug 21 (the crowded day)
for h in (6, 12, 18):
    axA.plot([2 + h / 24] * 2, [-0.10, 0.10], color=MUTED, lw=0.6, zorder=2)
    axA.text(2 + h / 24, 0.14, f"{h:02d}h", ha="center", va="bottom", fontsize=ANN, color=MUTED)

tl_rows = []
for ev in EVENTS:
    c = KIND_COLOR[ev["kind"]]
    vc = VERD_COLOR[ev["verdict"]]
    bx, by, bw = ev["box"]
    # marker on the axis
    if ev["t_end"]:
        axA.add_patch(matplotlib.patches.Rectangle((ev["t"], -0.09), max(ev["t_end"] - ev["t"], 0.012),
                      0.18, facecolor=c, edgecolor="none", zorder=4))
        anchor = (ev["t"] + ev["t_end"]) / 2
    else:
        axA.plot(ev["t"], 0, marker="D" if ev["kind"] == "prereg" else "o", ms=5.5,
                 mfc=c if ev["kind"] != "disclosure" else "white", mec=c, mew=1.2, ls="none", zorder=5)
        anchor = ev["t"]
    # verdict glyph next to the axis marker (PASS+FAIL drawn as two glyphs,
    # each in its own colour -- the two gates genuinely disagree by design)
    gy = -0.42 if by < 0 else 0.30
    gva = "top" if by < 0 else "bottom"
    # glyphs sit on the side AWAY from the label box: centred on the marker they were
    # crossed by their own leader line (design pass 2026-09-03)
    gside = -0.09 if bx >= anchor else 0.09
    if ev["verdict"] == "PASS+FAIL":
        axA.text(anchor + gside - 0.035, gy, "✓", ha="right", va=gva, fontsize=7.5,
                 color=C_PASS, fontweight="bold", zorder=6)
        axA.text(anchor + gside + 0.035, gy, "✗", ha="left", va=gva, fontsize=7.5,
                 color=C_FAIL, fontweight="bold", zorder=6)
    else:
        glyph = {"PASS": "✓", "FAIL": "✗", "PENDING": "?", "FAIL (identity)": "=0", "": ""}[ev["verdict"]]
        if glyph:
            gx = anchor + (-0.055 if bx >= anchor else 0.055)
            axA.text(gx, gy, glyph, ha="right" if bx >= anchor else "left", va=gva,
                     fontsize=7.5, color=vc, fontweight="bold", zorder=6)
    # leader line + label box
    ytop = by - 0.55 * np.sign(by) * 0  # box center
    axA.plot([anchor, bx], [0.22 * np.sign(by), by - np.sign(by) * 0.72], color=c, lw=0.6,
             alpha=0.75, zorder=3)
    axA.text(bx, by, ev["text"], ha="center", va="center", fontsize=ANN, color=INK,
             linespacing=1.32, zorder=7,
             bbox=dict(boxstyle="round,pad=0.42", facecolor=BOXFACE, edgecolor=vc if ev["verdict"] else c,
                       linewidth=1.0))
    tl_rows.append(dict(event=ev["id"], timestamp=ev["ts"], kind=ev["kind"],
                        verdict=ev["verdict"], label=ev["text"].replace("\n", " "),
                        display_note=ev["display_note"], source=ev["source"]))

# tiny legend
axA.plot(0.06, -2.9, marker="D", ms=5.5, mfc=C_PREREG, mec=C_PREREG, ls="none")
axA.text(0.12, -2.9, sentence_case("pre-registration"), fontsize=ANN, color=INK, va="center")
axA.add_patch(matplotlib.patches.Rectangle((0.035, -3.35), 0.05, 0.16, facecolor=INK))
axA.text(0.12, -3.27, sentence_case("run (span)"), fontsize=ANN, color=INK, va="center")
axA.text(0.045, -3.68, "✓", fontsize=7, color=C_PASS, fontweight="bold", va="center")
axA.text(0.12, -3.68, sentence_case("gate PASS"), fontsize=ANN, color=INK, va="center")
axA.text(0.045, -4.06, "✗", fontsize=7, color=C_FAIL, fontweight="bold", va="center")
axA.text(0.12, -4.06, sentence_case("gate FAIL"), fontsize=ANN, color=INK, va="center")

# ---------------------------------------------------------------------------
# panel b: the gate table (T3 drawn) -- every gate, every arm, every verdict
# ---------------------------------------------------------------------------
axB.set_xlim(0, 1)
axB.set_ylim(0, 1)
axB.axis("off")
# the '(table T3 drawn)' provenance aside moved to the Legend sidecar (directive 2,
# no-loss: it is the '(b) Table T3 drawn: every pre-registered gate ...' sentence)
axB.set_title(sc_title("b   every pre-registered gate and every outcome"),
              loc="left", fontweight="bold", fontsize=TYPE["panel_title"])

GATES = [
    ("ORIGINAL two-sided gate — P@100 ≥ 0.38 AND F1_det > 0.261, PBMC ≥1-molecule output "
     "(10 §5, 2026-08-19)", [
        ("Stage-2 run 2026-08-20 (stale; 12)",
         "both tiers n 377,236: P 0.188 / F1 0.235; tier-1 n 197,410: P 0.308 / F1 0.290;\n"
         "tier-1+IP n 142,400: P 0.362 / F1 0.290 — 'under the full gate, NO arm passes'",
         "FAIL", "manuscript/12_stage2_gate.md (verified to 6 dp there; CORRECTED framing)"),
        ("v1 final run 2026-08-21 (15 §2)",
         f"tier-1 ≥1 mol no-IP n {v1t1_noip[0]:,}: P {v1t1_noip[1]:.4f} / F1 {v1t1_noip[2]:.4f}; "
         f"IP arm n {v1t1_ip[0]:,}: P {v1t1_ip[1]:.4f} / F1 {v1t1_ip[2]:.4f}",
         "F1 pass / P FAIL", "score_pbmc_final__tier1.tsv; score_pbmc_final_ipfilt__tier1.tsv (15 SOUND)"),
        ("v2 re-run 2026-08-21 (19 §2)",
         f"tier-1 ≥1 mol no-IP n {v2t1_noip[0]:,}: P {v2t1_noip[1]:.4f} (byte-identical to v1); "
         f"IP arm n {v2t1_ip[0]:,}: P {v2t1_ip[1]:.4f} / F1 {v2t1_ip[2]:.4f}",
         "F1 pass / P FAIL", "fig2_accuracy.tsv (19 FIXED) — the ≥1-mol outputs remain the sensitivity arm"),
    ]),
    ("PRECISION-FIRST default gate — P@100 ≥ 0.50 on PBMC AND both mice (13 §1; commit `0e27b1a` "
     f"01:18:59, {DELTA_STR} before the arms)", [
        ("v1 final run 2026-08-21 (15 §1)",
         f"{v1h[1]:.4f} / {v1m1[1]:.4f} / {v1m2[1]:.4f}  (n {v1h[0]:,} / {v1m1[0]:,} / {v1m2[0]:,})",
         "PASS ×3", "score_*__PRESPEC_precision_default.tsv, v1 arms (15 SOUND)"),
        ("v2 re-run 2026-08-21, gates unchanged (19 §1)",
         f"{v2h[1]:.4f} / {v2m1[1]:.4f} / {v2m2[1]:.4f}  (n {v2h[0]:,} / {v2m1[0]:,} / {v2m2[0]:,})",
         "PASS ×3", "fig2_accuracy.tsv (19 FIXED)"),
        ("pbmc4k second donor 2026-08-22 (26 §6, same gate)",
         f"{d2[1]:.4f}  (n {d2[0]:,}; 38× its genic-shuffle null) — the 4th library",
         "PASS", "pbmc4k_donor2/run_ipfilt/score_*__PRESPEC_precision_default.tsv (26 FIXED)"),
    ]),
    ("PRIME adoption criterion — ΔP ≥ −0.005 & ΔR ≥ +0.010 & ΔF1_det > 0 on all three datasets, "
     "default vs default (24 §3.1, 2026-08-21; all prime numbers EXPLORATORY per §3.4)", [
        ("prime default vs v2 default (Amendment 3, 2026-08-23)",
         "prime's default output byte-identical to v2's manuscript arm on all four libraries →\n"
         "every Δ exactly 0.000000 — an identity, not a regression; the criterion was written\n"
         "against a v2 default that never existed (--ip-filter-mode default is `annotate`)",
         "FAIL (identity)", "manuscript/24_prime_preregistration.md Amendment 3; results/prime_bench/"),
    ]),
    ("3'UTR SINGLETON promotion — R2 arm; no precision loss with a recall gain, deciding evidence = "
     "testis mouse 2 + non-discovery data (27 §1/§3, 2026-08-21; Rule T)", [
        ("confirmatory evaluation",
         "NOT RUN — discovery-set numbers are labelled 'DISCOVERY SET — not evidence';\n"
         "the promotion ships behind a flag, default OFF, until 27 §3 passes",
         "PENDING", "manuscript/27_utr_singleton_prereg.md (status header + §1.3)"),
    ]),
]

# layout: group header rows + data rows
row_specs = []          # (kind, payload)
for g, (gtitle, rows) in enumerate(GATES):
    row_specs.append(("g", gtitle))
    for r in rows:
        row_specs.append(("r", r))

# group headers are wrapped for DRAWING only (the audit TSV keeps the full string):
# the longest ran past the canvas edge at the shared 6 pt floor (design pass 2026-09-03)
GHEAD_WRAP = 148
H_G, H_R = 0.062, 0.082  # heights (axes fraction) for a 1-line group header / data row
H_G2 = 0.088             # ... and for a header that wraps to two lines


def _gh(payload):
    return textwrap.fill(payload, GHEAD_WRAP)


def _hg(payload):
    return H_G if _gh(payload).count("\n") == 0 else H_G2


total = sum(_hg(p) if k == "g" else H_R for k, p in row_specs)
scale = 0.995 / total
y = 1.0
gate_rows = []
for kind, payload in row_specs:
    h = (_hg(payload) if kind == "g" else H_R) * scale
    y0 = y - h
    if kind == "g":
        axB.add_patch(matplotlib.patches.Rectangle((0, y0), 1.0, h, facecolor="#3D5A6C",
                      edgecolor="none", zorder=1))
        axB.text(0.006, y0 + h / 2, _gh(payload), ha="left", va="center", fontsize=ANN,
                 color="white", fontweight="bold", zorder=2, linespacing=1.25)
        current_gate = payload
    else:
        arm, values, verdict, src = payload
        axB.add_patch(matplotlib.patches.Rectangle((0, y0), 1.0, h, facecolor="white",
                      edgecolor="none", zorder=1))
        axB.axhline(y0, color=GRID, lw=0.5, zorder=2)
        axB.text(0.006, y0 + h / 2, sentence_case(arm), ha="left", va="center", fontsize=ANN,
                 color=INK, zorder=3, linespacing=1.25)
        axB.text(0.290, y0 + h / 2, values, ha="left", va="center", fontsize=ANN, color=INK,
                 zorder=3, linespacing=1.3)
        vc = VERD_COLOR.get(verdict.split()[0], C_FAIL) if verdict not in VERD_COLOR else VERD_COLOR[verdict]
        if verdict.startswith("F1 pass"):
            vc = C_FAIL
        axB.add_patch(matplotlib.patches.FancyBboxPatch((0.872, y0 + h / 2 - 0.016), 0.120, 0.032,
                      boxstyle="round,pad=0.004", facecolor=vc, edgecolor="none", zorder=3))
        axB.text(0.932, y0 + h / 2, verdict, ha="center", va="center", fontsize=ANN,
                 color="white", fontweight="bold", zorder=4)
        gate_rows.append(dict(gate=current_gate, evaluated_on=arm,
                              values=values.replace("\n", " "), verdict=verdict, source=src))
    y = y0
axB.axhline(y, color="#3D5A6C", lw=0.8, zorder=2)

# ---------------------------------------------------------------------------
# legend (moved OFF the image; the .caption.md sidecar is its single source)
# ---------------------------------------------------------------------------
legend = (
    "Figure S11 | The pre-registration record, drawn so it can be audited rather than taken on trust. "
    "(a) Timeline. The original two-sided gate (P@100 ≥ 0.38 AND F1_det > 0.261) was fixed in 10 §5 on "
    "2026-08-19, one day before the run it judged; the stale Stage-2 run FAILED it in full (12: 'NO arm "
    "passes'). The ≥2-molecule threshold was then examined in a sweep created 12 minutes after that FAIL "
    "was recorded — post hoc, full stop (12 CORRECTION) — and is disclosed as such: it informed the choice "
    "and is never cited as a result. The precision-first default and its P@100 ≥ 0.50 gate were committed "
    f"(`0e27b1a`, 01:18:59, git-verified) {DELTA_STR} before the final arms started (02:59:36); the v1 run "
    "passed on all three call sets, and the v2 re-run (code 9dfdefb, 16:14, gates UNCHANGED) passed again "
    "while the original 0.38 precision floor kept failing on every ≥1-molecule output — both facts shown, "
    "neither hidden. Three further pre-registrations followed on 2026-08-21, each written before its "
    "number existed: the prime adoption criterion (24, last write 23:34), the second-donor test (26, "
    "23:47, driver frozen 19 min before the caller started), and the 3'UTR singleton promotion (27, still "
    "untested). Outcomes: pbmc4k PASS 0.8279 (2026-08-22, verifier FIXED); prime FAIL at Δ 0.000000 — a "
    "byte-identity, not a regression, because the criterion was drafted against a v2 default that never "
    "existed (Amendment 3, 2026-08-23); singleton PENDING, flag default OFF. Markers whose documents "
    "record a date but no clock time sit at display positions, flagged in the audit TSV. (b) Table T3 "
    "drawn: every pre-registered gate, every arm it was evaluated on, every value and every verdict — "
    "including both FAILs and the untested pre-registration. Gate values are re-read from the score TSVs "
    "where they exist (v1 arms, v2 via fig2_accuracy.tsv, pbmc4k) and EXPECT-asserted against 15/19/26; "
    "Stage-2 values are quoted from 12 (its matrices were voided; its own verifier reproduced them to 6 "
    "decimals)."
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
pd.DataFrame(tl_rows).to_csv(p, sep="\t", index=False); print("wrote", p)
p = OUTDIR / f"{NAME}_gates.tsv"
pd.DataFrame(gate_rows).to_csv(p, sep="\t", index=False); print("wrote", p)

# ---------------------------------------------------------------------------
# sidecar caption
# ---------------------------------------------------------------------------
cap_md = f"""# Fig S11 — `figS11_gatehistory` sidecar (generated by `scripts/manuscript_figures/figS11_gatehistory.py`; single source of the legend)

## Legend

{legend}

**What this figure is for.** A reviewer asking "were the gates really fixed in advance?" should be able to
check every claim here against a named document and, for the timing that matters most, against git: the
default and its 0.50 gate are in commit `0e27b1a` (2026-08-21 01:18:59) and the arms started 02:59:36
({DELTA_STR} later; both timestamps verified from git by the 15 verifier). The disclosed deviations are ON
the figure, not in a footnote: the ≥2-molecule threshold was informed by a post-hoc sweep (created 12
minutes after the Stage-2 FAIL was recorded — 12 CORRECTION §1 withdraws the earlier 'pre-identified in
the plan' defence), and the second pre-registration ADDED a gate rather than replacing the failed one —
the original 0.38 precision floor is still reported, still failing, on every ≥1-molecule output in 15 §2
and 19 §2.

**Timestamp honesty.** Documents 10, 12, 24 (Amendment 3) and 27 record dates but no clock times; their
markers sit at display positions (noon, or an ordered slot within the day) and the audit TSV's
`display_note` column says exactly which. The 27 marker is ordered within 2026-08-21 only by its
citations (it cites 24 §3.1 and 25); no clock claim is made. The prime pre-registration's 23:34 is the
last-write time recorded in 26's preamble, not a git timestamp (24 was not committed at write time — the
26 verifier records the same limitation for 26 itself, closed there by the frozen driver's mtime chain).

**Verdict colour code (also encoded as glyphs, never colour alone).** PASS = ✓ green (#009E73);
FAIL = ✗ vermillion (#D55E00); the v1/v2 runs carry ✓✗ because the two gates disagree by design (the
0.50 default gate passes while the original 0.38 floor keeps failing on the sensitivity arms — the paper
reports both); prime's "=0" (#CC79A7) is a FAIL by identity (every Δ exactly 0.000000, byte-identical
output); PENDING = ? grey (untested pre-registration, flag default OFF).

## Provenance

Sources: every timeline event and every table row carries its document/TSV source in
`results/figures/manuscript/figS11_gatehistory.tsv` and `figS11_gatehistory_gates.tsv`. Gate values are
programmatic reads of the v1 score TSVs, the verified `fig2_accuracy.tsv` (v2) and the pbmc4k run's score
TSV, EXPECT-asserted against 15 §1-2, 19 §1-2 and 26 R2; Stage-2 values are cited from 12.
Every event and row: `results/figures/manuscript/figS11_gatehistory.tsv` / `figS11_gatehistory_gates.tsv`.
Script: `scripts/manuscript_figures/figS11_gatehistory.py`; rendered at PNG 600 dpi / PDF fonttype 42.
"""
cap_md += """
**Design pass 2026-09-03** (`manuscript/figures/DESIGN_DIRECTIVES.md`, supplement light pass). Type comes from
the shared style module `scripts/manuscript_figures/_pubstyle.py` (`apply_rc()`), and every on-figure label —
timeline event boxes, hour ticks, key, and all three table columns — now sits at or above the 6 pt floor
(the hour ticks were 4.8 pt). Panel titles and the prose parts of the table's arm column are sentence-cased
through `_pubstyle.sentence_case()`, canonical identifiers preserved and the lower-case panel letters kept.
Panel b's *(table T3 drawn)* aside left the image for the Legend above (no-loss: it is the \"(b) Table T3 drawn:
every pre-registered gate ...\" sentence). Two overlaps were fixed: each verdict glyph now sits on the side of
its marker away from its label box, where its own leader line used to run through it; and the gate group
headers are wrapped for drawing only — the PRIME header reached the canvas edge at the larger type — while the
audit TSV keeps each header as one unwrapped string. No panel, number or audit TSV changed; both TSVs
regenerate byte-identical.

*Flagged for the PI, not changed here (it would re-encode the verdict colours):* the PENDING neutral
`#999999` gives low contrast twice — the `?` glyph over the light day band, and white pill text on the grey.
The other verdict pills (PASS, FAIL, identity) use the same white-on-colour treatment, so a change would have
to be made to the set, not to one member.
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

print(f"figS11_gatehistory: 0.50 gate PASS x7 (v1 {v1h[1]:.4f}/{v1m1[1]:.4f}/{v1m2[1]:.4f}, "
      f"v2 {v2h[1]:.4f}/{v2m1[1]:.4f}/{v2m2[1]:.4f}, pbmc4k {d2[1]:.4f}); original 0.38 floor FAIL "
      f"throughout; prime FAIL at delta 0.000000 (identity); 27 PENDING; commit->arms {DELTA_STR}")

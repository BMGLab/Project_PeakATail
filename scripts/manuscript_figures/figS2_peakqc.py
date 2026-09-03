#!/usr/bin/env python3
"""
figS2_peakqc.py -- manuscript Fig S2 ("figS2_peakqc"): v2 peak-call QC over the
existing v2 outputs (READ-ONLY) plus the tier-2 cleavage-offset census and its
negative disposition.

Per FIGURES_MANIFEST.md row S2 and manuscript/21 (plan): widths, per-gene PAS
counts, tier composition, molecule-count distributions (arm always named,
MUST-NOT-CLAIM 8), and the ~90-105 nt cleavage-offset analysis (07 section 3).
Cited from R1.

SOURCES OF TRUTH (verified; nothing is re-run, nothing dump-derived)
  * Tier composition / headline n and P: the score_tool.py TSVs of the final v2
    arms (the same files Fig 2 reads; verified in 19 + final_v2_verify).
  * Width / per-gene / molecule-support distributions: the v2 run outputs
    themselves (run/pasbed.bed, run/pas_gene.tsv, pas_*.bed) -- read-only QC of
    files already behind the verified 19-gate numbers.
  * Cleavage-offset census: results/paramsweep/leadB/ (dist_* signed-distance
    files, offsetsweep_* post-hoc sweeps) and results/paramsweep/arms_all.tsv,
    with per-lead verdicts in results/paramsweep/VERDICTS.md section 3
    (executing manuscript/28 section 3 step 1; acceptance criteria 28 section 2,
    pre-registered).  The original ~90-105 nt observation is 07 section 3.
    The census ran on the development slices (PBMC chr19+21, mouse1 chr18+19)
    -- labelled as such on the figure; slice caveat: paramsweep README section 1.

NUMBER POLICY (manuscript/01, 21 section 7)
  * every count names its arm; "atlas-agreement precision", never bare
    "precision"; Kinnex t5 = ">=5 long-read records at the terminus", never a
    de-duplicated UMI count (MUST-NOT-CLAIM 8).
  * the 72.2% single-molecule share is the IP arm's (scored); the no-IP arm's
    72.1% is what the manuscript text quotes (19 section 4) -- both stated.
  * the offset correction did NOT rescue tier-2: stated on the figure
    (P@100 0.085 vs genic null 0.028 after correction; the internal-priming
    decoy rate more than triples).  This supplement documents the diagnosis and
    the negative disposition, not a fix.

PANELS
  a  Tier composition of every final v2 arm (scored n from the score TSVs):
     tier-2 / tier-1 single-molecule / tier-1 >=2 molecules (= the
     pre-registered default on the IP arms).
  b  Tier-2 call-width distributions (run/pasbed.bed intervals; tier-1 calls
     are 1-bp cleavage points by construction).
  c  Calls per gene of the precision default, per library.
  d  Tier-1 distinct-clip-molecule distribution per library, >=2 threshold
     marked.
  e  Post-hoc offset sweep (slice): tier-2 P@10 vs shift on both truths, the
     four agreeing offset estimates, tier-1 control peaking at +0.
  f  Signed-distance histograms to the nearest atlas site, tier-2 before /
     after correction, both species (slice).
  g  Tier-1 control histogram (+0 offset; shifting it by +95 destroys it) and
     the boxed negative disposition.

OUTPUTS
  manuscript/figures/figS2_peakqc.{png,pdf}    600 dpi PNG / fonttype-42 PDF
  manuscript/figures/figS2_peakqc.caption.md   script-written sidecar
                                               ('## Legend' + '## Provenance')
  results/figures/manuscript/figS2_peakqc_composition.tsv    (a)
  results/figures/manuscript/figS2_peakqc_widths.tsv         (b)
  results/figures/manuscript/figS2_peakqc_pergene.tsv        (c)
  results/figures/manuscript/figS2_peakqc_support.tsv        (d)
  results/figures/manuscript/figS2_peakqc_offset_sweep.tsv   (e)
  results/figures/manuscript/figS2_peakqc_offset_hist.tsv    (f, g)
  results/figures/manuscript/figS2_peakqc_disposition.tsv    (e-g verdict rows)

Run:  export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/figS2_peakqc.py
"""
import os
import re
import textwrap
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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
BT = WD / "results/benchmark_tools"
PS = WD / "results/paramsweep"
OUTDIR = WD / "results/figures/manuscript"
FIGDIR = WD / "manuscript/figures"
NAME = "figS2_peakqc"
OUTDIR.mkdir(parents=True, exist_ok=True)

INK, MUTED, GRID = "#1B2429", "#5A6B73", "#D8E0E3"
# Okabe-Ito, colour follows the entity (library / tier slice), fixed order.
C_PBMC, C_M1, C_M2 = "#0072B2", "#009E73", "#CC79A7"
C_TIER2, C_SINGLE, C_DEFAULT = "#C9D1D5", "#8FBBDD", "#0072B2"   # fig2 panel-d colours
C_KINNEX = "#D55E00"

# ---------------------------------------------------------------------------
# arm registry -- every count on this figure names its arm
# ---------------------------------------------------------------------------
H = BT / "pbmc_10k_v3"
M = BT / "gse104556"
ARMS = {
    "pbmc_noip": dict(
        lib="PBMC 10k v3", arm="final v2, no IP filter (peakatail_clipseeded_final_v2)",
        short="PBMC 10k v3\nno-IP arm", ip=False, color=C_PBMC,
        d=H / "peakatail_clipseeded_final_v2", pfx="pbmc_final_v2", default_scored=False),
    "pbmc_ip": dict(
        lib="PBMC 10k v3", arm="final v2, IP filter (peakatail_clipseeded_final_v2_ipfilt)",
        short="PBMC 10k v3\nIP arm", ip=True, color=C_PBMC,
        d=H / "peakatail_clipseeded_final_v2_ipfilt", pfx="pbmc_final_v2_ipfilt", default_scored=True),
    "mouse1": dict(
        lib="GSE104556 testis mouse 1", arm="final v2, IP filter (peakatail_clipseeded_final_v2/mouse1)",
        short="testis mouse 1\nIP arm", ip=True, color=C_M1,
        d=M / "peakatail_clipseeded_final_v2/mouse1", pfx="testis_m1_final_v2", default_scored=True),
    "mouse2": dict(
        lib="GSE104556 testis mouse 2", arm="final v2, IP filter (peakatail_clipseeded_final_v2/mouse2)",
        short="testis mouse 2\nIP arm", ip=True, color=C_M2,
        d=M / "peakatail_clipseeded_final_v2/mouse2", pfx="testis_m2_final_v2", default_scored=True),
}
IP_ARMS = ["pbmc_ip", "mouse1", "mouse2"]   # the arms behind the pre-registered default


def score_np(tsv: Path):
    """(n, atlas-agreement P@100) straight out of one score_tool.py TSV."""
    df = pd.read_csv(tsv, sep="\t")
    q = df[(df.cutoff_bp == 100.0) & (df.panel == "precision") & (df.series == "real")
           & (df.reference == "atlas_full")]
    assert len(q) == 1, tsv
    return int(q.n_query.iloc[0]), float(q.value.iloc[0])


# ---------------------------------------------------------------------------
# panel a data: tier composition, scored n (verified path: same TSVs as Fig 2)
# ---------------------------------------------------------------------------
comp_rows = []
comp = {}
for key, a in ARMS.items():
    d, pfx = a["d"], a["pfx"]
    n_both, P_both = score_np(d / f"score_{pfx}.tsv")
    n_t1, P_t1 = score_np(d / f"score_{pfx}__tier1.tsv")
    n_t2, P_t2 = score_np(d / f"score_{pfx}__tier2.tsv")
    assert n_t1 + n_t2 == n_both, (key, n_t1, n_t2, n_both)
    if a["default_scored"]:
        n_def, P_def = score_np(d / f"score_{pfx}__PRESPEC_precision_default.tsv")
        def_label = "pre-registered precision default (tier-1, IP, >=2 molecules)"
    else:
        # the no-IP >=2-molecule file is POSTHOC and is NOT the default (21 sec.7);
        # scored n only, to show the funnel -- never presented as an operating point.
        n_def, P_def = score_np(d / f"score_{pfx}__tier1_ge2mol_noIP_POSTHOC.tsv")
        def_label = "tier-1 >=2 molecules, no IP -- POSTHOC, NOT the default"
    comp[key] = dict(n_both=n_both, n_t1=n_t1, n_t2=n_t2, n_def=n_def,
                     n_single=n_t1 - n_def, P_both=P_both, P_t1=P_t1, P_t2=P_t2, P_def=P_def)
    for slice_name, n, P, src in [
            ("tier-2 (coverage-only)", n_t2, P_t2, f"score_{pfx}__tier2.tsv"),
            ("tier-1 single-molecule", n_t1 - n_def, np.nan, f"score_{pfx}__tier1.tsv minus default (n difference)"),
            (def_label, n_def, P_def,
             f"score_{pfx}__PRESPEC_precision_default.tsv" if a["default_scored"]
             else f"score_{pfx}__tier1_ge2mol_noIP_POSTHOC.tsv")]:
        comp_rows.append(dict(panel="a", library=a["lib"], arm=a["arm"], slice=slice_name,
                              n_scored=n, atlas_agreement_precision_100bp=P,
                              n_all_calls_scored=n_both, tier1_scored=n_t1,
                              source=str((d / src).relative_to(WD))))

# EXPECT: the verified 19 / final_v2_verify headline numbers (same as Fig 2 asserts)
assert comp["pbmc_ip"]["n_def"] == 46524 and round(comp["pbmc_ip"]["P_def"], 4) == 0.7062
assert comp["pbmc_ip"]["n_t1"] == 167565 and round(comp["pbmc_ip"]["P_t1"], 4) == 0.3520
assert comp["pbmc_ip"]["n_t2"] == 166355 and comp["pbmc_ip"]["n_both"] == 333920
assert comp["pbmc_noip"]["n_t1"] == 222955 and comp["pbmc_noip"]["n_def"] == 62110
assert (comp["mouse1"]["n_def"], round(comp["mouse1"]["P_def"], 4)) == (26255, 0.7450)
assert (comp["mouse2"]["n_def"], round(comp["mouse2"]["P_def"], 4)) == (26526, 0.7572)
FRAC_SINGLE_IP = comp["pbmc_ip"]["n_single"] / comp["pbmc_ip"]["n_t1"]
assert round(FRAC_SINGLE_IP * 100, 1) == 72.2   # IP arm; the no-IP arm's is 72.1% (19 sec.4)
# cross-check against the Fig 2 audit TSV (same quantities, independent read)
_f2 = pd.read_csv(OUTDIR / "fig2_accuracy_tiers.tsv", sep="\t")
assert int(_f2[_f2.slice.str.startswith("tier-2")].n.iloc[0]) == comp["pbmc_ip"]["n_t2"]
# the audit TSV stores %.6f-rounded values, so compare at that precision
assert abs(float(_f2.frac_of_tier1_single_molecule.iloc[0]) - FRAC_SINGLE_IP) < 5e-7

# ---------------------------------------------------------------------------
# panels b-d data: read the v2 run outputs (read-only)
# ---------------------------------------------------------------------------
BED_COLS = ["chrom", "start", "end", "pas_id", "score", "strand"]
width_rows, pergene_rows, support_rows = [], [], []
WBINS = [1, 50, 100, 150, 200, 300, 400, 600, 800, 1200, 2000, 10**9]
WLAB = ["1-49", "50-99", "100-149", "150-199", "200-299", "300-399", "400-599",
        "600-799", "800-1199", "1200-1999", ">=2000"]
SBINS = [1, 2, 3, 4, 5, 10, 20, 10**9]
SLAB = ["1", "2", "3", "4", "5-9", "10-19", ">=20"]
GBINS = [1, 2, 3, 4, 5, 6, 7, 8, 10**9]
GLAB = ["1", "2", "3", "4", "5", "6", "7", ">=8"]
qc = {}
for key in IP_ARMS:
    a = ARMS[key]
    pasbed = pd.read_csv(a["d"] / "run/pasbed.bed", sep="\t", header=None, names=BED_COLS)
    t2w = (pasbed.loc[pasbed.score == 0, "end"] - pasbed.loc[pasbed.score == 0, "start"]).to_numpy()
    t1 = pasbed[pasbed.score > 0]
    t1w = (t1.end - t1.start).to_numpy()
    assert (t1w == 1).all(), f"{key}: tier-1 calls are expected to be 1-bp points"
    # tier-1 distinct clip molecules = the BED score column: selecting score>=2
    # reproduces the pre-registered default BED exactly (checked here).
    dflt = pd.read_csv(a["d"] / "pas_PRESPEC_precision_default.bed", sep="\t", header=None,
                       names=BED_COLS)
    assert len(dflt) == (t1.score >= 2).sum(), key
    assert set(dflt.pas_id) == set(t1.loc[t1.score >= 2, "pas_id"]), key
    # per-gene calls of the default
    pg = pd.read_csv(a["d"] / "run/pas_gene.tsv", sep="\t", header=None,
                     names=["pas_id", "gene_id"]).set_index("pas_id")["gene_id"]
    genes = dflt.pas_id.map(pg)
    assert genes.notna().all(), f"{key}: default calls missing from pas_gene.tsv"
    per_gene = genes.value_counts()
    qc[key] = dict(t2w=t2w, t1_support=t1.score.to_numpy(), per_gene=per_gene,
                   n_t1_bed=len(t1), n_t2_bed=len(t2w), n_def_bed=len(dflt))
    # BED counts differ from scored n only by non-primary-contig rows (<0.05%)
    for bed_n, sc_n in [(len(t1), comp[key]["n_t1"]), (len(t2w), comp[key]["n_t2"]),
                        (len(dflt), comp[key]["n_def"])]:
        assert abs(bed_n - sc_n) / sc_n < 5e-4, (key, bed_n, sc_n)
    src_bed = str((a["d"] / "run/pasbed.bed").relative_to(WD))
    wh = np.histogram(t2w, bins=WBINS)[0]
    for lab, c in zip(WLAB, wh):
        width_rows.append(dict(panel="b", library=a["lib"], arm=a["arm"], tier="tier-2",
                               width_bin_bp=lab, n_calls=int(c), frac=c / len(t2w),
                               median_width_bp=float(np.median(t2w)),
                               p90_width_bp=float(np.quantile(t2w, 0.90)), source=src_bed))
    sh = np.histogram(t1.score, bins=SBINS)[0]
    for lab, c in zip(SLAB, sh):
        support_rows.append(dict(panel="d", library=a["lib"], arm=a["arm"], tier="tier-1",
                                 clip_molecules_bin=lab, n_calls=int(c), frac=c / len(t1),
                                 frac_single_molecule=float((t1.score == 1).mean()),
                                 source=src_bed))
    gh = np.histogram(per_gene.to_numpy(), bins=GBINS)[0]
    for lab, c in zip(GLAB, gh):
        pergene_rows.append(dict(panel="c", library=a["lib"], arm=a["arm"],
                                 output="pre-registered precision default",
                                 calls_per_gene_bin=lab, n_genes=int(c), frac_genes=c / len(per_gene),
                                 n_genes_total=int(len(per_gene)),
                                 mean_calls_per_gene=float(per_gene.mean()),
                                 source=f"{src_bed} + run/pas_gene.tsv + pas_PRESPEC_precision_default.bed"))

# ---------------------------------------------------------------------------
# panels e-g data: the lead-B offset census (slice; results/paramsweep/leadB)
# ---------------------------------------------------------------------------
def read_sweep(name):
    df = pd.read_csv(PS / "leadB" / name, sep="\t")
    df["source"] = f"results/paramsweep/leadB/{name}"
    return df

sw_pbmc_t2 = read_sweep("offsetsweep_pbmc_base_tier2.tsv")
sw_m1_t2 = read_sweep("offsetsweep_m1_base_tier2.tsv")
sw_pbmc_def = read_sweep("offsetsweep_pbmc_base_default.tsv")
# the four agreeing estimates (VERDICTS.md sec.3): two tool --auto values +
# three post-hoc argmaxes recomputed here from the sweep TSVs
AUTO_PBMC, AUTO_M1 = 103, 87    # tool's own --auto-cleavage-offset (VERDICTS sec.3; A-crest 94)
am_atlas_pbmc = int(sw_pbmc_t2.loc[sw_pbmc_t2.atlas_P10.idxmax(), "shift"])
am_kin_pbmc = int(sw_pbmc_t2.loc[sw_pbmc_t2.kin5_P10.idxmax(), "shift"])
am_atlas_m1 = int(sw_m1_t2.loc[sw_m1_t2.atlas_P10.idxmax(), "shift"])
assert am_atlas_pbmc == am_kin_pbmc == am_atlas_m1 == 95   # VERDICTS sec.3 table
assert int(sw_pbmc_def.loc[sw_pbmc_def.atlas_P10.idxmax(), "shift"]) == 0  # tier-1 control


def read_dist_txt(fname):
    """Parse a leadB dist_*.txt into rows: (set, truth, shift, bin_lo, pct) + summary."""
    txt = (PS / "leadB" / fname).read_text()
    out = []
    header = re.compile(r"\[(?P<name>[^|]+)\|(?P<truth>[^|]+)\| post-hoc shift \+(?P<shift>\d+)\]"
                        r"\s+n=(?P<n>\d+)\s+median=(?P<med>[+-]?\d+)")
    blocks = txt.split("\n\n")
    for b in blocks:
        m = header.search(b)
        if not m:
            continue
        lines = [l for l in b.splitlines()]
        bins = vals = None
        for l in lines:
            if l.strip().startswith("bin :"):
                bins = [int(x) for x in l.split(":", 1)[1].split()]
            if l.strip().startswith("%%"):
                vals = [float(x) for x in l.split(":", 1)[1].split()]
        assert bins is not None and vals is not None and len(bins) == len(vals), fname
        mm = re.search(r"\|d\|<=10 ([0-9.]+)\s+\|d\|<=25 ([0-9.]+)", b)
        mode = re.search(r"modal 10-bp bin=\[([+-]?\d+),", b)
        for lo, v in zip(bins, vals):
            out.append(dict(set=m["name"].strip(), truth=m["truth"].strip(),
                            shift=int(m["shift"]), n=int(m["n"]), median_d=int(m["med"]),
                            modal_10bp_bin_lo=int(mode.group(1)),
                            frac_within_10=float(mm.group(1)), frac_within_25=float(mm.group(2)),
                            bin_lo=lo, bin_width_bp=20, pct_of_calls=v,
                            source=f"results/paramsweep/leadB/{fname}"))
    return pd.DataFrame(out)

d_pbmc = read_dist_txt("dist_pbmc_tier2.txt")
d_m1 = read_dist_txt("dist_m1_tier2.txt")
d_ctrl = read_dist_txt("dist_pbmc_default.txt")
d_tool95 = read_dist_txt("dist_pbmc_tier2_tool95.txt")
d_tool87 = read_dist_txt("dist_m1_tier2_tool87.txt")


def block(df, truth, shift):
    # NB "shift" must be bracket-accessed: df.shift is the pandas method
    q = df[(df.truth == truth) & (df["shift"] == shift)].sort_values("bin_lo")
    assert len(q) == 25, (truth, shift, len(q))
    return q

pb_before = block(d_pbmc, "atlas", 0)
pb_after = block(d_pbmc, "atlas", 95)
m1_before = block(d_m1, "atlas", 0)
m1_after = block(d_m1, "atlas", 87)
ct_before = block(d_ctrl, "atlas", 0)
ct_after = block(d_ctrl, "atlas", 95)
# EXPECT: VERDICTS.md sec.3 signed-distance table
assert (pb_before.frac_within_10.iloc[0], pb_before.modal_10bp_bin_lo.iloc[0],
        pb_before.median_d.iloc[0]) == (0.0072, 90, 64)
assert (pb_after.frac_within_10.iloc[0], pb_after.modal_10bp_bin_lo.iloc[0],
        pb_after.median_d.iloc[0]) == (0.0200, 0, 13)
assert (m1_before.frac_within_10.iloc[0], m1_before.modal_10bp_bin_lo.iloc[0],
        m1_before.median_d.iloc[0]) == (0.0106, 100, 89)
assert (m1_after.frac_within_10.iloc[0], m1_after.modal_10bp_bin_lo.iloc[0]) == (0.0370, 10)
assert (ct_before.frac_within_10.iloc[0], ct_before.modal_10bp_bin_lo.iloc[0],
        ct_before.median_d.iloc[0]) == (0.6316, 0, 0)
assert (ct_after.modal_10bp_bin_lo.iloc[0], ct_after.median_d.iloc[0]) == (-100, -91)
# real corrected tool runs (annotation values)
assert (d_tool95[d_tool95.truth == "atlas"].median_d.iloc[0]) == 3
assert (d_tool87[d_tool87.truth == "atlas"].median_d.iloc[0]) == 16

# "does not rescue the tier" -- real tool runs, arms_all.tsv (VERDICTS sec.3)
arms_all = pd.read_csv(PS / "arms_all.tsv", sep="\t")


def arm_row(arm, view="tier2"):
    q = arms_all[(arms_all.arm == arm) & (arms_all.view == view)]
    assert len(q) == 1, (arm, view)
    return q.iloc[0]

r_base = arm_row("pbmc_base")
r_off95 = arm_row("pbmc_B_off95")
r_m1b = arm_row("m1_base")
r_m1auto = arm_row("m1_B_auto")
assert round(r_off95["P@100"], 4) == 0.0851 and round(r_off95["null_P@100_mean"], 4) == 0.0283
assert round(r_base["kin_decoy@25"], 4) == 0.0200 and round(r_off95["kin_decoy@25"], 4) == 0.0646
assert round(r_m1auto["P@100"], 4) == 0.1789 and round(r_m1auto["null_P@100_mean"], 4) == 0.0157
assert round(r_base["R_det@100"], 4) == 0.0365 and round(r_off95["R_det@100"], 4) == 0.0369
assert round(r_m1b["R_det@100"], 4) == 0.0218 and round(r_m1auto["R_det@100"], 4) == 0.0193
disposition = pd.DataFrame([
    dict(panel="g", arm="pbmc_base (slice)", view="tier-2", n=int(r_base.n_scored),
         P100=r_base["P@100"], null_P100=r_base["null_P@100_mean"],
         R_det100=r_base["R_det@100"], kin_decoy25=r_base["kin_decoy@25"],
         source="results/paramsweep/arms_all.tsv"),
    dict(panel="g", arm="pbmc_B_off95 = --cleavage-offset 95 (slice)", view="tier-2",
         n=int(r_off95.n_scored), P100=r_off95["P@100"], null_P100=r_off95["null_P@100_mean"],
         R_det100=r_off95["R_det@100"], kin_decoy25=r_off95["kin_decoy@25"],
         source="results/paramsweep/arms_all.tsv"),
    dict(panel="g", arm="m1_base (slice)", view="tier-2", n=int(r_m1b.n_scored),
         P100=r_m1b["P@100"], null_P100=r_m1b["null_P@100_mean"],
         R_det100=r_m1b["R_det@100"], kin_decoy25=np.nan, source="results/paramsweep/arms_all.tsv"),
    dict(panel="g", arm="m1_B_auto = --auto-cleavage-offset (87 bp; slice)", view="tier-2",
         n=int(r_m1auto.n_scored), P100=r_m1auto["P@100"], null_P100=r_m1auto["null_P@100_mean"],
         R_det100=r_m1auto["R_det@100"], kin_decoy25=np.nan, source="results/paramsweep/arms_all.tsv"),
])

# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
# Publication layout: the footer caption and the boxed disposition text live in
# the caption sidecar's Legend now (surgery pass, 2026-09-02), so the canvas
# keeps only the panel band (12.4 -> 10.0 in tall; 8.3 -> 7.9 in wide -- the
# three dense bottom panels tolerate that much and no more at this type size).
fig = plt.figure(figsize=(7.9, 10.0))
gs = fig.add_gridspec(3, 6, left=0.118, right=0.972, top=0.963, bottom=0.058,
                      hspace=0.48, wspace=1.25)
axA = fig.add_subplot(gs[0, 0:3])
axB = fig.add_subplot(gs[0, 3:6])
axC = fig.add_subplot(gs[1, 0:3])
axD = fig.add_subplot(gs[1, 3:6])
axE = fig.add_subplot(gs[2, 0:2])
axF = fig.add_subplot(gs[2, 2:4])
axG = fig.add_subplot(gs[2, 4:6])


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED, length=2.5, labelsize=TYPE["tick"])
    ax.grid(True, color=GRID, lw=0.5, alpha=0.7)
    ax.set_axisbelow(True)


# ---- panel a: tier composition -------------------------------------------
order = ["pbmc_noip", "pbmc_ip", "mouse1", "mouse2"]
for i, key in enumerate(order):
    c, a = comp[key], ARMS[key]
    y = len(order) - 1 - i
    left = 0
    slices = [(c["n_t2"], C_TIER2), (c["n_single"], C_SINGLE), (c["n_def"], C_DEFAULT)]
    for n, col in slices:
        axA.barh(y, n, left=left, height=0.62, color=col, edgecolor="white", lw=0.8, zorder=3)
        left += n
    axA.text(left + 6000, y, f"{c['n_both']:,}", va="center", fontsize=ANN, color=INK)
    # second line goes under the bar; the widest bar's line starts further left and
    # sits lower so it clears both the bar above and the axes edge on the right
    x2, y2 = (1.85e5, y - 0.46) if key == "pbmc_noip" else (left + 6000, y - 0.27)
    if a["default_scored"]:
        axA.text(x2, y2, f"default {c['n_def']:,}  (P@100 {c['P_def']:.3f})",
                 va="center", fontsize=ANN, color=MUTED)
    else:
        axA.text(x2, y2, f">=2 mol {c['n_def']:,} (POSTHOC, not the default)",
                 va="center", fontsize=ANN, color=MUTED)
axA.set_yticks(range(len(order)))
axA.set_yticklabels([sentence_case(ARMS[k]["short"]) for k in reversed(order)],
                    fontsize=TYPE["tick"])
# floor lowered so the frameless key clears the bottom arm's per-arm data line
# (they touched at the shared 6 pt annotation floor; design pass 2026-09-03)
axA.set_ylim(-1.15, len(order) - 0.5)
axA.set_xlim(0, 6.15e5)
axA.set_xticks([0, 1e5, 2e5, 3e5, 4e5])
axA.set_xticklabels(["0", "100k", "200k", "300k", "400k"])
axA.set_xlabel(sentence_case("calls scored (score_tool.py n, cutoff 100 bp)"))
axA.set_title(sc_title("a   tier composition, all final v2 arms"), loc="left", fontweight="bold")
axA.legend([matplotlib.patches.Patch(facecolor=c, edgecolor="white") for c in (C_TIER2, C_SINGLE, C_DEFAULT)],
           [sentence_case(s) for s in ("tier-2 (coverage-only)", "tier-1 single-molecule",
                                       "tier-1 >=2 mol (default, IP arms)")],
           loc="lower right", bbox_to_anchor=(1.0, 0.0), frameon=False,
           handlelength=1.2, labelspacing=0.3, fontsize=ANN)
# (the single-molecule-share / arm-naming note moved to the caption Legend)
style(axA)
axA.grid(False, axis="y")

# ---- panel b: tier-2 widths ----------------------------------------------
x = np.arange(len(WLAB))
wdf = pd.DataFrame(width_rows)
for k, key in enumerate(IP_ARMS):
    a = ARMS[key]
    q = wdf[wdf.library == a["lib"]].sort_values("width_bin_bp", key=lambda s: [WLAB.index(v) for v in s])
    axB.bar(x + (k - 1) * 0.27, q.frac * 100, width=0.25, color=a["color"], zorder=3,
            label=f"{a['short'].splitlines()[0]} (median {q.median_width_bp.iloc[0]:.0f} bp)")
axB.set_xticks(x)
axB.set_xticklabels(WLAB, rotation=45, ha="right", fontsize=ANN)
axB.set_xlabel(sentence_case("tier-2 call width (bp; run/pasbed.bed interval)"))
axB.set_ylabel(sentence_case("% of tier-2 calls"))
axB.set_title(sc_title("b   tier-2 call widths (IP arms)"), loc="left", fontweight="bold")
axB.legend(loc="upper right", frameon=False, handlelength=1.2, labelspacing=0.3)
# (the tier-1-width/tier-2-interval explanation moved to the caption Legend)
style(axB)

# ---- panel c: calls per gene of the default ------------------------------
gdf = pd.DataFrame(pergene_rows)
x = np.arange(len(GLAB))
for k, key in enumerate(IP_ARMS):
    a = ARMS[key]
    q = gdf[gdf.library == a["lib"]]
    q = q.set_index("calls_per_gene_bin").loc[GLAB]
    axC.bar(x + (k - 1) * 0.27, q.frac_genes * 100, width=0.25, color=a["color"], zorder=3,
            label=(f"{a['short'].splitlines()[0]}: {q.n_genes_total.iloc[0]:,} genes, "
                   f"mean {q.mean_calls_per_gene.iloc[0]:.2f}"))
axC.set_xticks(x)
axC.set_xticklabels(GLAB)
axC.set_xlabel(sentence_case("calls per gene, pre-registered precision default"))
axC.set_ylabel(sentence_case("% of genes with >=1 default call"))
axC.set_title(sc_title("c   calls per gene (precision default)"), loc="left", fontweight="bold")
axC.legend(loc="upper right", frameon=False, handlelength=1.2, labelspacing=0.3)
style(axC)

# ---- panel d: tier-1 molecule support ------------------------------------
sdf = pd.DataFrame(support_rows)
x = np.arange(len(SLAB))
for k, key in enumerate(IP_ARMS):
    a = ARMS[key]
    q = sdf[sdf.library == a["lib"]].set_index("clip_molecules_bin").loc[SLAB]
    axD.bar(x + (k - 1) * 0.27, q.frac * 100, width=0.25, color=a["color"], zorder=3,
            label=f"{a['short'].splitlines()[0]}: {q.frac_single_molecule.iloc[0]*100:.1f}% single-molecule")
axD.axvline(0.5, color=INK, lw=0.9, ls=(0, (4, 2)), zorder=4)
axD.text(0.56, axD.get_ylim()[1] * 0.02 + 66, ">=2 molecules kept\nby the default", fontsize=ANN,
         color=INK, ha="left", va="top")
axD.set_xticks(x)
axD.set_xticklabels(SLAB)
axD.set_xlabel(sentence_case("distinct clip molecules per tier-1 call (BED score)"))
axD.set_ylabel(sentence_case("% of tier-1 calls"))
axD.set_ylim(0, 78)
axD.set_title(sc_title("d   tier-1 molecule support (IP arms)"), loc="left", fontweight="bold")
axD.legend(loc="upper right", frameon=False, handlelength=1.2, labelspacing=0.3,
           bbox_to_anchor=(1.0, 0.78))
style(axD)

# ---- panel e: offset sweep (slice) ---------------------------------------
def norm(v):
    return v / v.max()

axE.plot(sw_pbmc_t2["shift"], norm(sw_pbmc_t2.atlas_P10), color=C_PBMC, lw=1.3, zorder=5,
         label=f"PBMC tier-2 vs atlas (max P@10 {sw_pbmc_t2.atlas_P10.max():.3f})")
axE.plot(sw_pbmc_t2["shift"], norm(sw_pbmc_t2.kin5_P10), color=C_KINNEX, lw=1.1, ls=(0, (4, 2)), zorder=5,
         label=f"PBMC tier-2 vs Kinnex t5 (max {sw_pbmc_t2.kin5_P10.max():.3f})")
axE.plot(sw_m1_t2["shift"], norm(sw_m1_t2.atlas_P10), color=C_M1, lw=1.3, zorder=5,
         label=f"mouse tier-2 vs atlas (max {sw_m1_t2.atlas_P10.max():.3f})")
axE.plot(sw_pbmc_def["shift"], norm(sw_pbmc_def.atlas_P10), color=INK, lw=1.0, ls=(0, (1, 1.5)), zorder=4,
         label=f"control: tier-1 default (max {sw_pbmc_def.atlas_P10.max():.3f} at +0)")
# reference line stops below the annotation band: drawn full-height it ran through
# the legend entries and the --auto-cleavage-offset note (design pass 2026-09-03)
axE.axvline(95, ymax=0.63, color=MUTED, lw=0.8, zorder=2)
axE.text(98, 0.16, "+95", ha="left", fontsize=ANN, color=MUTED)
for xo, lab in ((AUTO_PBMC, "auto 103 (PBMC)"), (AUTO_M1, "auto 87 (mouse)")):
    axE.plot([xo], [1.045], marker="v", ms=4, color=MUTED, zorder=6, clip_on=False)
axE.text(0.0, 0.78, "--auto-cleavage-offset:\n87 (mouse), 103 (PBMC);\nA-crest 94",
         transform=axE.transAxes, fontsize=ANN, color=MUTED, ha="left", va="top", linespacing=1.3)
axE.set_xlim(0, 200)
axE.set_ylim(0, 1.62)
axE.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
axE.set_xlabel(sentence_case("post-hoc 3' shift applied to calls (nt)"))
axE.set_ylabel(sentence_case("P@10, relative to each curve's max"))
axE.set_title(sc_title("e   offset census (slices)"), loc="left", fontweight="bold")
axE.legend(loc="upper left", bbox_to_anchor=(-0.02, 1.02), frameon=False,
           handlelength=1.3, labelspacing=0.3, fontsize=ANN)
style(axE)

# ---- panel f: tier-2 signed distance before/after ------------------------
def steps(ax, blockdf, color, ls, lw, label):
    e = np.append(blockdf.bin_lo.to_numpy(), blockdf.bin_lo.iloc[-1] + 20)
    ax.stairs(blockdf.pct_of_calls.to_numpy(), e, color=color, ls=ls, lw=lw, label=label, zorder=4)

steps(axF, pb_before, C_PBMC, "-", 1.3, "PBMC before (mode +90..+100)")
steps(axF, pb_after, C_PBMC, (0, (4, 2)), 1.1, "PBMC shifted +95 (mode 0..+10)")
steps(axF, m1_before, C_M1, "-", 1.3, "mouse before (mode +100..+110)")
steps(axF, m1_after, C_M1, (0, (4, 2)), 1.1, "mouse shifted +87 (mode +10..+20)")
axF.axvspan(90, 105, color=GRID, alpha=0.55, zorder=1)
axF.axvline(0, ymax=0.73, color=MUTED, lw=0.6, zorder=2)   # stops below the legend band
axF.set_xlim(-200, 300)
axF.set_ylim(0, 7.0)
axF.set_yticks([0, 1, 2, 3, 4, 5, 6])
axF.set_xlabel(sentence_case("signed distance to nearest atlas site (bp)\n"
                             "d = truth − call; + = truth downstream"))
axF.set_ylabel(sentence_case("% of tier-2 calls per 20-bp bin"))
axF.set_title(sc_title("f   tier-2 offset, before / after"), loc="left", fontweight="bold")
axF.legend(loc="upper left", frameon=False, handlelength=1.5, labelspacing=0.3, fontsize=ANN)
# data annotations stay; the heavy-tails caveat sentence moved to the caption Legend
axF.text(0.98, 0.73, "shaded: ~90–105 nt (07 §3)\nmedians +64 → +13 (PBMC),\n"
         "+89 → +13 (mouse); real\ncorrected runs: +3 / +16",
         transform=axF.transAxes, fontsize=ANN, color=MUTED, ha="right", va="top", linespacing=1.3,
         bbox=dict(facecolor="white", alpha=0.92, edgecolor="none", pad=1.0), zorder=6)
style(axF)

# ---- panel g: tier-1 control + disposition -------------------------------
steps(axG, ct_before, C_DEFAULT, "-", 1.3, "tier-1 default arm, before: mode [0,+10), median +0")
steps(axG, ct_after, MUTED, (0, (4, 2)), 1.1, "same calls shifted +95: mode [−100,−90)")
axG.axvline(0, ymax=0.86, color=MUTED, lw=0.6, zorder=2)   # stops below the legend band
axG.set_xlim(-200, 300)
axG.set_ylim(0, 50)
# two-line axis label: the single line overran the right canvas edge at the shared
# type scale (design pass 2026-09-03)
axG.set_xlabel(sentence_case("signed distance to nearest\natlas site (bp)"))
axG.set_ylabel(sentence_case("% of calls per 20-bp bin"))
axG.set_title(sc_title("g   tier-1 control: no offset"), loc="left", fontweight="bold")
axG.legend(loc="upper right", frameon=False, handlelength=1.5, labelspacing=0.3, fontsize=ANN)
# short data label stays; "the offset is tier-2-specific" moved to the caption Legend
axG.text(0.97, 0.55,
         "clip-anchored tier-1\nsits at +0\n(|d|≤10: 63.2%)",
         transform=axG.transAxes, fontsize=ANN, color=INK, ha="right", va="top", linespacing=1.35)
# The boxed negative-disposition verdict is no longer drawn on the image
# (surgery pass, 2026-09-02): its every sentence lives in the caption Legend's
# "Negative disposition" paragraph, and its numbers in figS2_peakqc_disposition.tsv.
style(axG)

# ---- legend paragraph (written to the caption sidecar, not drawn) ---------
caption = (
    "Figure S2 | v2 peak-call QC and the tier-2 cleavage-offset census. "
    "(a) Tier composition of every final v2 arm (code 9dfdefb), scored n from the single scoring path "
    "(score_tool.py, cutoff 100 bp): the pre-registered precision default (tier-1, IP filter, >=2 distinct "
    f"clip molecules) keeps {comp['pbmc_ip']['n_def']:,} of {comp['pbmc_ip']['n_t1']:,} tier-1 PBMC calls "
    f"(P@100 {comp['pbmc_ip']['P_def']:.4f}); the no-IP >=2-molecule file is POSTHOC and is not the default. "
    "(b) Tier-2 calls are wide intervals (whole peak regions); tier-1 calls are 1-bp clip-anchored cleavage "
    "points. (c) Default calls per gene. (d) Tier-1 distinct-clip-molecule distributions; the >=2 threshold "
    "is the pre-registered default cut (12 CORRECTION). "
    "(e-g) The ~90-105 nt tier-2 coverage offset (predicted in 07 §3: 10x R2 coverage runs out before the "
    "poly(A) junction; measured in the Lead-B census on the development slices, PBMC chr19+21 / mouse1 "
    "chr18+19). Four independent estimates agree: tool --auto-cleavage-offset 103 (PBMC) / 87 (mouse), and "
    "post-hoc argmax +95 on atlas P@10 (both species) and on Kinnex t5 P@10 (PBMC; t5 = >=5 long-read records "
    "at the terminus, not de-duplicated UMIs). The signed-distance mode moves to ~0 when corrected (f) and the "
    "clip-anchored tier-1 control sits at +0 (g) -- the offset is real and tier-2-specific. It is also a "
    "negative result: correction fails every pre-registered acceptance criterion of 28 §2 (negative "
    "disposition below). "
    "Sources: 19 + final_v2_verify (FIXED) score TSVs; v2 run outputs (read-only); results/paramsweep/ "
    "arms_all.tsv, leadB/, VERDICTS.md §2-§3 (slice sweep executing 28 §3 step 1; provisional "
    "pending verifier). Slice caveat: the PBMC slice is clip-richer than the genome (paramsweep README §1); "
    "whether the offset holds for all 166,355 full-BAM tier-2 sites is an inference from the slices, not a "
    "measurement. Every plotted value: results/figures/manuscript/figS2_peakqc_*.tsv."
)

for ext, kw in (("png", dict(dpi=600)), ("pdf", {})):
    p = FIGDIR / f"{NAME}.{ext}"
    fig.savefig(p, **kw)
    print("wrote", p)
p = OUTDIR / f"{NAME}.png"
fig.savefig(p, dpi=600)
print("wrote", p)

# ---------------------------------------------------------------------------
# audit TSVs (every plotted value, with a source column)
# ---------------------------------------------------------------------------
sweep_out = pd.concat([sw_pbmc_t2.assign(panel="e", curve="pbmc_tier2"),
                       sw_m1_t2.assign(panel="e", curve="m1_tier2"),
                       sw_pbmc_def.assign(panel="e", curve="pbmc_default_control")],
                      ignore_index=True)
est = pd.DataFrame([
    dict(panel="e", curve="estimate", label="tool --auto-cleavage-offset, PBMC (aataaa_spacing, A-crest 94)",
         shift=AUTO_PBMC, source="results/paramsweep/VERDICTS.md sec.3"),
    dict(panel="e", curve="estimate", label="tool --auto-cleavage-offset, mouse 1 (A-crest 94)",
         shift=AUTO_M1, source="results/paramsweep/VERDICTS.md sec.3"),
    dict(panel="e", curve="estimate", label="post-hoc argmax tier-2 atlas P@10, PBMC",
         shift=am_atlas_pbmc, source="results/paramsweep/leadB/offsetsweep_pbmc_base_tier2.tsv"),
    dict(panel="e", curve="estimate", label="post-hoc argmax tier-2 atlas P@10, mouse",
         shift=am_atlas_m1, source="results/paramsweep/leadB/offsetsweep_m1_base_tier2.tsv"),
    dict(panel="e", curve="estimate", label="post-hoc argmax tier-2 Kinnex t5 P@10, PBMC",
         shift=am_kin_pbmc, source="results/paramsweep/leadB/offsetsweep_pbmc_base_tier2.tsv"),
])
hist_out = pd.concat([pb_before.assign(panel="f"), pb_after.assign(panel="f"),
                      m1_before.assign(panel="f"), m1_after.assign(panel="f"),
                      ct_before.assign(panel="g"), ct_after.assign(panel="g"),
                      d_tool95[d_tool95.truth == "atlas"].assign(panel="f_annotation"),
                      d_tool87[d_tool87.truth == "atlas"].assign(panel="f_annotation")],
                     ignore_index=True)
outs = [("composition", pd.DataFrame(comp_rows)), ("widths", pd.DataFrame(width_rows)),
        ("pergene", pd.DataFrame(pergene_rows)), ("support", pd.DataFrame(support_rows)),
        ("offset_sweep", pd.concat([sweep_out, est], ignore_index=True)),
        ("offset_hist", hist_out), ("disposition", disposition)]
for suffix, df in outs:
    p = OUTDIR / f"{NAME}_{suffix}.tsv"
    df.to_csv(p, sep="\t", index=False, float_format="%.6f")
    print("wrote", p)

# ---------------------------------------------------------------------------
# sidecar caption
# ---------------------------------------------------------------------------
cap_md = f"""# Fig S2 — `figS2_peakqc` caption (generated by `scripts/manuscript_figures/figS2_peakqc.py`)

## Legend

{caption}

**Panel a** — every count names its arm. PBMC no-IP arm: {comp['pbmc_noip']['n_both']:,} calls
({comp['pbmc_noip']['n_t1']:,} tier-1 / {comp['pbmc_noip']['n_t2']:,} tier-2); PBMC IP arm:
{comp['pbmc_ip']['n_both']:,} ({comp['pbmc_ip']['n_t1']:,} / {comp['pbmc_ip']['n_t2']:,}), default
{comp['pbmc_ip']['n_def']:,}; mice (IP): {comp['mouse1']['n_both']:,} / {comp['mouse2']['n_both']:,} calls,
defaults {comp['mouse1']['n_def']:,} (P@100 {comp['mouse1']['P_def']:.4f}) and {comp['mouse2']['n_def']:,}
(P@100 {comp['mouse2']['P_def']:.4f}). Single-molecule share of tier-1: {FRAC_SINGLE_IP*100:.1f}% on the PBMC
IP arm ({comp['pbmc_ip']['n_single']:,} / {comp['pbmc_ip']['n_t1']:,}); the no-IP arm's 72.1%
(1 − 62,110 / 222,955) is the verifier's recount and the value the manuscript text quotes (19 §4) —
per 21 §7 MUST-NOT-CLAIM 8, the arm is always named. Atlas-agreement precision @100 bp of each scored
output is in `figS2_peakqc_composition.tsv`; tier-2 scores ~{comp['pbmc_ip']['P_t2']:.3f} (PBMC IP arm),
which is what motivates the offset census below.
**Panel b** — tier-2 interval widths from `run/pasbed.bed` (medians: PBMC IP
{float(np.median(qc['pbmc_ip']['t2w'])):.0f} bp, mouse1 {float(np.median(qc['mouse1']['t2w'])):.0f} bp,
mouse2 {float(np.median(qc['mouse2']['t2w'])):.0f} bp). Tier-1 calls are 1-bp points (asserted). The Lead-A
sweep measured why tier-2 is coarse: `lambda_gradient`'s gradient candidates never survive its own filters, so
the fallback emits one whole-peak interval per peak — slice medians 334 bp (PBMC) / 255 bp (mouse), 0 of
16,414 / 0 of 2,722 candidate intervals partitioned (VERDICTS §2; a detector finding, not a `--max-pas` ceiling).
**Panel c** — calls per gene of the pre-registered default. **Panel d** — distinct clip molecules per tier-1
call (the BED score; selecting score >= 2 reproduces the default BED exactly, asserted per arm).
**Panels e–g** — the Lead-B cleavage-offset census on the development slices (PBMC chr19+21, mouse1 chr18+19;
baselines reproduce the reference runs to the call, paramsweep README §5). e: post-hoc shift sweep, each curve
normalised to its own max (absolute maxima printed in the legend); the four agreeing estimates are marked.
f: signed-distance histograms (d = truth − call, transcript orientation), tier-2 before vs after the shift;
the [+90,+105] band is 07 §3's predicted offset. The tails are heavy: only ~7% of tier-2 calls lie within
100 bp of any atlas site, so the mode describes the near-atlas minority. g: the tier-1 control (default arm) peaks at [0,+10) with
median +0 and is destroyed by the same +95 shift — the offset is tier-2-specific, which is why v2 already
exempts clip-supported rows (`cleavage_offset.py::rewrite_bed_3prime_offset(skip_supported=True)`).
**Negative disposition (VERDICTS §3; offset correction does NOT rescue tier-2)** — correcting the offset
triples tier-2 P@10 but leaves the tier at {r_off95['P@100']:.3f} atlas-agreement P@100 vs a
{r_off95['null_P100'] if 'null_P100' in r_off95 else r_off95['null_P@100_mean']:.3f} genic-shuffle null on PBMC
(3.0× the null; was 2.2× before correction; {r_m1auto['P@100']:.3f} vs {r_m1auto['null_P@100_mean']:.3f} on
mouse), buys no detected-gene recall (R_det@100 {r_base['R_det@100']:.4f} → {r_off95['R_det@100']:.4f}), and more
than triples the internal-priming decoy rate ({r_base['kin_decoy@25']:.3f} → {r_off95['kin_decoy@25']:.3f} at
matched IP-veto; ×12.8 with no veto) — the failure mode 24 §3.2.2 named in advance. FAIL on every 28 §2
criterion; nothing is adopted; the default arm is byte-identical across every offset arm (tier-1 exempt by
design), and the prime branch's TASK-E disposition (report an `inferred_cleavage` column,
do not move `pasbed.bed`) stands. Kinnex "t5" support is an alignment-record count, not de-duplicated UMIs
(21 §7 MUST-NOT-CLAIM 8).

Caveats: the census ran on development slices — no headline claim may rest on them (24 §2), and whether the
offset holds for all {comp['pbmc_ip']['n_t2']:,} full-BAM tier-2 sites (IP arm) is an inference; the paramsweep
sweep is PROVISIONAL (unseen by the verifier). Width/per-gene/support distributions are computed from the v2
run outputs read-only; BED-level counts differ from scored n by <0.05% (non-primary contigs).

## Provenance

Sources: `manuscript/19_final_gate_v2.md` + `results/benchmark_tools/final_v2_verify/VERIFIED_v2.md` (FIXED)
score TSVs; v2 run outputs; `results/paramsweep/{{arms_all.tsv,leadB/,VERDICTS.md,README.md}}` (executing
`manuscript/28` §3 step 1 under the pre-registered criteria of 28 §2); `manuscript/07` §3. Every plotted value:
`results/figures/manuscript/figS2_peakqc_*.tsv` (source column per row). PNG 600 dpi, PDF vector with
subsetted TrueType (fonttype 42, no Type 3). The working-phase render carried the legend paragraph, the
boxed disposition verdict and several in-panel explanations on the image; at the 2026-09-02 surgery pass
they moved into the Legend above, and the image keeps panel letters, short titles, axis labels, legends
and data annotations only.
"""
cap_md += """
**Design pass 2026-09-03** (`manuscript/figures/DESIGN_DIRECTIVES.md`, supplement light pass). Type comes from
the shared style module `scripts/manuscript_figures/_pubstyle.py` (`apply_rc()`), and every on-figure annotation,
key entry and tick label now sits at or above the 6 pt floor. Axis labels, panel titles, key entries and prose
tick labels are sentence-cased through `_pubstyle.sentence_case()`, canonical identifiers preserved and the
lower-case panel letters kept. Four overlaps were fixed at that larger type: the +95 marker in **e** and the
zero references in **f** and **g** now stop below their key bands instead of running through the key entries and
the `--auto-cleavage-offset` note; **f**'s shaded-band note carries an opaque backing so the zero reference
passes behind it; **a**'s floor was lowered so its key clears the bottom arm's data line; **g**'s axis label is
set on two lines (one line overran the canvas edge); and the left margin was widened for the **a** arm labels.
No panel, number or audit TSV changed; all seven TSVs regenerate byte-identical.
"""
p = FIGDIR / f"{NAME}.caption.md"
p.write_text(cap_md)
print("wrote", p)

# ---------------------------------------------------------------------------
# edge check: outer 8 px of the PNG must be blank
# ---------------------------------------------------------------------------
from PIL import Image
im = np.asarray(Image.open(FIGDIR / f"{NAME}.png").convert("L"))
edge = 8
border = np.concatenate([im[:edge, :].ravel(), im[-edge:, :].ravel(),
                         im[:, :edge].ravel(), im[:, -edge:].ravel()])
ink = int((border < 250).sum())
print(f"edge check: {im.shape[1]}x{im.shape[0]} px, ink pixels in the outer {edge} px = {ink}")
assert ink == 0, "INK IN THE OUTER 8 PX -- something is clipped"

print(f"figS2_peakqc: default {comp['pbmc_ip']['n_def']:,}/{comp['mouse1']['n_def']:,}/"
      f"{comp['mouse2']['n_def']:,}; tier-1 single-molecule {FRAC_SINGLE_IP*100:.1f}% (IP arm); "
      f"offset estimates {AUTO_PBMC}/{AUTO_M1}/+95; post-correction tier-2 P@100 "
      f"{r_off95['P@100']:.4f} vs null {r_off95['null_P@100_mean']:.4f} -- negative disposition")

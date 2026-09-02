#!/usr/bin/env python3
"""
fig3_tradeoff.py -- manuscript Fig 3 ("fig3_tradeoff"): the operating-point
trade surface, the matching-window sensitivity check, biological-replicate reproducibility and
the compute footprint of the final v2 caller (code 9dfdefb, after PRs #93/#96/#97).

PANELS
  a  TRADE SURFACE.  Atlas-agreement precision @100 bp vs detected-gene recall @100 bp as the
     molecule-support threshold sweeps >=1, 2, 3, 5, 10 distinct clip molecules on the IP-filtered
     v2 arms (PBMC 10k v3, testis mouse 1, testis mouse 2).  ONLY the >=2 point is the
     pre-registered default (13 section 1); >=3/>=5/>=10 and >=1 are descriptive.  Every point is
     recomputed here by score_tool.py with the reference arguments of stage2_final_launch.sh --
     nothing is interpolated.  The >=1 and >=2 points reproduce manuscript/19 to 4 dp (asserted).
  b  MATCHING-WINDOW SENSITIVITY.  Precision and detected-gene recall at 10 / 25 / 50 / 100 bp for
     the v2 default, the v2 >=1-molecule sensitivity arm and the five competitors, PBMC.  Answers
     "is the precision lead an artefact of a loose window?".  Source: the verified score TSVs.
  c  REPLICATE REPRODUCIBILITY.  The two GSE104556 testis mice as biological replicates:
     strand-matched site agreement (bedtools closest -s -d -t first, both directions) at 25 and
     100 bp plus a strand-aware Jaccard, computed here for the v2 default and the v2 >=1-molecule
     arm; competitor and shipped-caller values are the v1-era head-to-head run.
  d  COMPUTE.  Wall time vs peak RSS on the PBMC BAM: PeakATail v1 -> v2 (the resolved Stage-1d
     item) against the five competitors.

INPUTS (all primary or verified)
  results/benchmark_tools/pbmc_10k_v3/peakatail_clipseeded_final_v2_ipfilt/{pas.bed,runtime_mem.txt,score_*.tsv}
  results/benchmark_tools/gse104556/peakatail_clipseeded_final_v2/mouse{1,2}/{pas.bed,runtime_mem.txt,score_*.tsv,.work_*/}
  results/benchmark_tools/{pbmc_10k_v3,gse104556}/<competitor>/score_*.tsv
  results/benchmark_tools/gse104556/depth_concordance.tsv        (verified, manuscript/05 + 09)
  results/figures/manuscript/benchmark_headtohead.tsv            (verified, manuscript/05)
  manuscript/19_final_gate_v2.md, manuscript/18_spermatogenesis_final.md   (quoted facts)

OUTPUTS
  manuscript/figures/fig3_tradeoff.{png,pdf}   (PNG 600 dpi, PDF fonttype 42)
  manuscript/figures/fig3_tradeoff.caption.md  ('## Legend' = the journal legend,
                                                single source of the caption;
                                                '## Provenance')
  results/figures/manuscript/fig3_tradeoff.tsv               (panels a + b)
  results/figures/manuscript/fig3_tradeoff_reproducibility.tsv (panel c)
  results/figures/manuscript/fig3_tradeoff_compute.tsv         (panel d)
  results/figures/manuscript/fig3_tradeoff.png                 (copy)

RUN
  cd /mnt/ssd1/Projects/PeakATail_wd
  export LC_ALL=C            # REQUIRED (machine locale tr_TR corrupts sort/bedtools)
  python3 scripts/manuscript_figures/fig3_tradeoff.py

The heavy steps (15 score_tool runs, 16 bedtools concordance runs) cache into
$FIG3_WORKDIR (default /mnt/ssd0/emaout/fig3_trade_repro_work); delete it to force a recompute.
"""
import os
import subprocess
import textwrap
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ["LC_ALL"] = "C"

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "DejaVu Sans"
matplotlib.rcParams["font.size"] = 7.0
matplotlib.rcParams["axes.titlesize"] = 8.0
matplotlib.rcParams["axes.labelsize"] = 7.4
matplotlib.rcParams["legend.fontsize"] = 6.0
matplotlib.rcParams["xtick.labelsize"] = 6.6
matplotlib.rcParams["ytick.labelsize"] = 6.6

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
BT = WD / "results/benchmark_tools"
OUTDIR = WD / "results/figures/manuscript"
FIGDIR = WD / "manuscript/figures"
NAME = "fig3_tradeoff"
WORK = Path(os.environ.get("FIG3_WORKDIR", "/mnt/ssd0/emaout/fig3_trade_repro_work"))
OUTDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)

CODE = "9dfdefb"
INK, MUTED, GRID = "#1B2429", "#5A6B73", "#D8E0E3"
# Okabe-Ito colourblind-safe palette (manuscript/figures/README.md convention)
BLUE, GREEN, VERM, PURPLE, ORANGE, SKY, GREY = (
    "#0072B2", "#009E73", "#D55E00", "#CC79A7", "#E69F00", "#56B4E9", "#999999")
DSCOL = {"pbmc": BLUE, "mouse1": GREEN, "mouse2": PURPLE}
DSLAB = {"pbmc": "PBMC 10k v3", "mouse1": "testis mouse 1", "mouse2": "testis mouse 2"}
TOOLCOL = {"PeakATail v2 default": BLUE, "PeakATail v2 tier-1 ≥1 mol": SKY,
           "polyApipe": VERM, "SCAPTURE": GREEN, "Sierra": PURPLE,
           "scAPAtrap": ORANGE, "scUTRquant*": GREY}

PY = str(WD / "tools/PeakATail/.venv/bin/python")
SCORE = str(WD / "scripts/benchmark_tools/score_tool.py")
R = WD / "data/references"
G = BT / "gse104556"
MOUSE_REFS = ["--atlas", str(R / "atlases/polyasite2.GRCm38.96.rep_sites.bed6"),
              "--tes", str(R / "atlases/tes.protein_coding.GRCm38.102.bed6"),
              "--genome", str(R / "mouse/chrom.sizes.filt"),
              "--genebodies", str(G / "shared_refs/genebodies.merged.bed"),
              "--detected-atlas", str(G / "shared_refs/pas2.in_detected_genes.bed")]
HUMAN_REFS = ["--detected-atlas", str(BT / "shared_refs_pbmc/pas2.in_detected_genes.bed")]

ARM = {  # dataset -> (species, IP-filtered v2 arm directory, verified score-file prefix)
    "pbmc":   ("human", BT / "pbmc_10k_v3/peakatail_clipseeded_final_v2_ipfilt", "pbmc_final_v2_ipfilt"),
    "mouse1": ("mouse", G / "peakatail_clipseeded_final_v2/mouse1", "testis_m1_final_v2"),
    "mouse2": ("mouse", G / "peakatail_clipseeded_final_v2/mouse2", "testis_m2_final_v2"),
}
THRESHOLDS = [1, 2, 3, 5, 10]
PREREG_K = 2
CUTOFFS = [10, 25, 50, 100]
GATE_P, ORIG_P = 0.50, 0.38


def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, executable="/bin/bash",
                          env=dict(os.environ, LC_ALL="C"), capture_output=True, text=True).stdout


# ---------------------------------------------------------------------------
# panel a -- molecule-support sweep, recomputed with score_tool.py
# ---------------------------------------------------------------------------
def sweep_scores():
    rows = []
    for ds, (sp, arm, _pfx) in ARM.items():
        out = WORK / ds
        out.mkdir(parents=True, exist_ok=True)
        for k in THRESHOLDS:
            lab = f"{ds}_ge{k}mol"
            tsv = out / f"score_{lab}.tsv"
            if not tsv.exists():
                bed = out / f"pas_ge{k}mol.bed"
                sh(f"awk -F'\\t' -v K={k} '$5>=K' {arm/'pas.bed'} > {bed}")
                refs = HUMAN_REFS if sp == "human" else MOUSE_REFS
                subprocess.run([PY, SCORE, str(bed), lab, "--outdir", str(out)] + refs,
                               check=True, env=dict(os.environ, LC_ALL="C"),
                               stdout=subprocess.DEVNULL)
                print(f"  computed {tsv}")
            d = pd.read_csv(tsv, sep="\t")
            def pick(panel, ref, series="real"):
                q = d[(d.panel == panel) & (d.series == series) & (d.reference == ref)
                      & (d.cutoff_bp == 100.0)]
                return q
            pr = pick("precision", "atlas_full").iloc[0]
            rd = pick("recall", "atlas_detected").iloc[0]
            rf = pick("recall", "atlas_full").iloc[0]
            f1 = pick("f1", "atlas_detected").iloc[0]
            nulls = pick("precision", "atlas_full", "null_genic").sort_values("replicate")["value"].tolist()
            rows.append(dict(panel="a", dataset=ds, arm=f"v2 IP arm, >={k} molecules",
                             min_molecules=k, pre_registered=(k == PREREG_K),
                             cutoff_bp=100, n=int(pr.n_query), n_matched=int(pr.n_matched),
                             atlas_agreement_precision=float(pr.value),
                             recall_detected_genes=float(rd.value),
                             recall_detected_denominator=int(rd.n_query),
                             recall_full_atlas=float(rf.value),
                             recall_full_denominator=int(rf.n_query),
                             F1_detected_genes=float(f1.value),
                             null_P_seed1=nulls[0], null_P_seed2=nulls[1], null_P_seed3=nulls[2],
                             null_P_mean=float(np.mean(nulls)), plotted=True,
                             source=str(tsv.relative_to(WORK.parent)) + " (recomputed by this script)"))
    return pd.DataFrame(rows)


sweep = sweep_scores()

# Cross-check against manuscript/19 section 1 (verified FIXED). If these fail the pipeline drifted.
VERIFIED_19 = {  # (dataset, k) -> (n, P@100, R_det@100)
    ("pbmc", 2):   (46524, 0.7062, 0.1754), ("pbmc", 1):   (167565, 0.3520, 0.2685),
    ("mouse1", 2): (26255, 0.7450, 0.2048), ("mouse1", 1): (52792, 0.5686, 0.2806),
    ("mouse2", 2): (26526, 0.7572, 0.2080), ("mouse2", 1): (51842, 0.5880, 0.2826),
}
for (ds, k), (n, p, r) in VERIFIED_19.items():
    q = sweep[(sweep.dataset == ds) & (sweep.min_molecules == k)].iloc[0]
    assert int(q.n) == n, (ds, k, q.n, n)
    assert abs(q.atlas_agreement_precision - p) < 5e-5, (ds, k, q.atlas_agreement_precision, p)
    assert abs(q.recall_detected_genes - r) < 5e-5, (ds, k, q.recall_detected_genes, r)
print("panel a: >=1 and >=2 points reproduce manuscript/19 section 1 to 4 dp on all three arms")

# ---------------------------------------------------------------------------
# panel b -- matching-window sensitivity from the verified score TSVs
# ---------------------------------------------------------------------------
CUTOFF_ARMS = [
    # label, dataset, score TSV
    ("PeakATail v2 default", "pbmc", BT / "pbmc_10k_v3/peakatail_clipseeded_final_v2_ipfilt"
     / "score_pbmc_final_v2_ipfilt__PRESPEC_precision_default.tsv"),
    ("PeakATail v2 tier-1 ≥1 mol", "pbmc", BT / "pbmc_10k_v3/peakatail_clipseeded_final_v2_ipfilt"
     / "score_pbmc_final_v2_ipfilt__tier1.tsv"),
    ("scUTRquant*", "pbmc", BT / "pbmc_10k_v3/scutrquant/score_scutrquant.tsv"),
    ("SCAPTURE", "pbmc", BT / "pbmc_10k_v3/scapture/score_scapture.tsv"),
    ("polyApipe", "pbmc", BT / "pbmc_10k_v3/polyapipe/score_polyapipe.tsv"),
    ("Sierra", "pbmc", BT / "pbmc_10k_v3/sierra/score_sierra.tsv"),
    ("scAPAtrap", "pbmc", BT / "pbmc_10k_v3/scapatrap/score_scapatrap.tsv"),
    # mouse rows: audit TSV only (not plotted), so the window check is auditable on both datasets
    ("PeakATail v2 default", "mouse1", G / "peakatail_clipseeded_final_v2/mouse1"
     / "score_testis_m1_final_v2__PRESPEC_precision_default.tsv"),
    ("PeakATail v2 default", "mouse2", G / "peakatail_clipseeded_final_v2/mouse2"
     / "score_testis_m2_final_v2__PRESPEC_precision_default.tsv"),
    ("PeakATail v2 tier-1 ≥1 mol", "mouse1", G / "peakatail_clipseeded_final_v2/mouse1"
     / "score_testis_m1_final_v2__tier1.tsv"),
    ("PeakATail v2 tier-1 ≥1 mol", "mouse2", G / "peakatail_clipseeded_final_v2/mouse2"
     / "score_testis_m2_final_v2__tier1.tsv"),
    ("scUTRquant*", "mouse1", G / "scutrquant/score_scutrquant_mouse1.tsv"),
    ("scUTRquant*", "mouse2", G / "scutrquant/score_scutrquant_mouse2.tsv"),
    ("SCAPTURE", "mouse1", G / "scapture/mouse1/score_scapture_mouse1.tsv"),
    ("polyApipe", "mouse1", G / "polyapipe/score_polyapipe_mouse1.tsv"),
    ("polyApipe", "mouse2", G / "polyapipe/score_polyapipe_mouse2.tsv"),
    ("Sierra", "mouse1", G / "sierra/score_sierra_mouse1.tsv"),
    ("Sierra", "mouse2", G / "sierra/score_sierra_mouse2.tsv"),
    ("scAPAtrap", "mouse1", G / "scapatrap/mouse1/score_scapatrap_mouse1.tsv"),
    ("scAPAtrap", "mouse2", G / "scapatrap/mouse2/score_scapatrap_mouse2.tsv"),
]
crows = []
for tool, ds, tsv in CUTOFF_ARMS:
    d = pd.read_csv(tsv, sep="\t")
    for c in CUTOFFS + [200]:
        pr = d[(d.panel == "precision") & (d.series == "real") & (d.reference == "atlas_full")
               & (d.cutoff_bp == float(c))].iloc[0]
        rd = d[(d.panel == "recall") & (d.series == "real") & (d.reference == "atlas_detected")
               & (d.cutoff_bp == float(c))].iloc[0]
        crows.append(dict(panel="b", dataset=ds, tool=tool, cutoff_bp=c,
                          n=int(pr.n_query), n_matched=int(pr.n_matched),
                          atlas_agreement_precision=float(pr.value),
                          recall_detected_genes=float(rd.value),
                          recall_detected_denominator=int(rd.n_query),
                          plotted=(ds == "pbmc" and c in CUTOFFS),
                          source=str(tsv.relative_to(WD))))
cut = pd.DataFrame(crows)

# ---------------------------------------------------------------------------
# panel c -- cross-mouse replicate agreement, computed here with bedtools
# ---------------------------------------------------------------------------
def pt_bed(mouse, suffix):
    n = mouse[-1]
    stem = f"testis_m{n}_final_v2__{suffix}"
    return G / f"peakatail_clipseeded_final_v2/{mouse}/.work_{stem}/{stem}"


def repro_table():
    cache = WORK / "repro_v2.tsv"
    if cache.exists():
        return pd.read_csv(cache, sep="\t")
    rows = []
    tmp = WORK / "repro"
    tmp.mkdir(parents=True, exist_ok=True)
    for arm, lab in (("PRESPEC_precision_default", "PeakATail v2 default"),
                     ("tier1", "PeakATail v2 tier-1 ≥1 mol")):
        A1, A2 = f"{pt_bed('mouse1', arm)}.pt.bed", f"{pt_bed('mouse2', arm)}.pt.bed"
        for w in (25, 100):
            inter = union = 0
            for s, tag in (("+", "p"), ("-", "m")):
                paths = []
                for i, src in ((1, A1), (2, A2)):
                    p = tmp / f"j{i}_{tag}.bed"
                    sh(f"""awk -F'\\t' -v S='{s}' -v W={w} 'BEGIN{{OFS="\\t"}} $6==S{{st=$2-W; if(st<0)st=0; print $1,st,$3+W}}' {src} """
                       f"| sort -k1,1 -k2,2n | bedtools merge -i - > {p}")
                    paths.append(p)
                j = sh(f"bedtools jaccard -a {paths[0]} -b {paths[1]} | tail -1").split()
                inter += int(j[0]); union += int(j[1])
            jac = inter / union
            for a, b in (("mouse1", "mouse2"), ("mouse2", "mouse1")):
                qa = A1 if a == "mouse1" else A2
                qb = A2 if a == "mouse1" else A1
                nb = f"{pt_bed(b, arm)}.null"
                n, m = (int(x) for x in sh(
                    f"bedtools closest -s -d -t first -a {qa} -b {qb} 2>/dev/null | "
                    f"awk -F'\\t' -v W={w} '{{n++; if($NF>=0&&$NF<=W)m++}} END{{print n, m+0}}'").split())
                ch = np.mean([float(sh(
                    f"bedtools closest -s -d -t first -a {qa} -b {nb}.s{sd}.pt.bed 2>/dev/null | "
                    f"awk -F'\\t' -v W={w} '{{n++; if($NF>=0&&$NF<=W)m++}} END{{print m/n}}'"))
                    for sd in (1, 2, 3)])
                rows.append(dict(tool=lab, direction=f"{a}->{b}", window_bp=w, n_query=n,
                                 n_agree=m, concordance=m / n, chance_mean=float(ch),
                                 jaccard_strand_aware=jac, jaccard_intersect_bp=inter,
                                 jaccard_union_bp=union, provenance="recomputed by this script"))
    df = pd.DataFrame(rows)
    df.to_csv(cache, sep="\t", index=False, float_format="%.6f")
    print(f"  computed {cache}")
    return df


rep_v2 = repro_table()

# v1-era per-tool reproducibility (verified: manuscript/05 adversarial pass, manuscript/09 table)
dc = pd.read_csv(G / "depth_concordance.tsv", sep="\t")
DC_LAB = {"PeakATail": "PeakATail shipped (v1)", "polyApipe": "polyApipe (v1)",
          "Sierra": "Sierra (v1)", "scAPAtrap": "scAPAtrap (v1)"}
rep_v1 = pd.DataFrame([dict(tool=DC_LAB[r.tool], direction=r.direction, window_bp=100,
                            n_query=int(r.n_query), n_agree=int(round(r.conc_all * r.n_query)),
                            concordance=float(r.conc_all), chance_mean=float(r.conc_null),
                            jaccard_strand_aware=np.nan, jaccard_intersect_bp=np.nan,
                            jaccard_union_bp=np.nan,
                            provenance="results/benchmark_tools/gse104556/depth_concordance.tsv "
                                       "(v1-era head-to-head run, verified in manuscript/05)")
                       for r in dc.itertuples() if r.tool in DC_LAB])
repro = pd.concat([rep_v2, rep_v1], ignore_index=True)
repro.insert(0, "panel", "c")

# ---------------------------------------------------------------------------
# panel d -- compute (PBMC BAM)
# ---------------------------------------------------------------------------
def rt(path):
    """(wall hours, peak RSS 'GB' = ru_maxrss kbytes / 1e6, the convention of 19 sec.3 / 15 sec.5)"""
    txt = Path(path).read_text()
    el = [l for l in txt.splitlines() if "Elapsed (wall clock)" in l][0].split()[-1]
    parts = [float(x) for x in el.split(":")]
    hours = (parts[0] * 3600 + parts[1] * 60 + parts[2]) / 3600 if len(parts) == 3 \
        else (parts[0] * 60 + parts[1]) / 3600
    kb = float([l for l in txt.splitlines() if "Maximum resident set size" in l][0].split()[-1])
    return hours, kb / 1e6, el


hh = pd.read_csv(OUTDIR / "benchmark_headtohead.tsv", sep="\t")
hh = hh[(hh.panel == "e_resources") & (hh.dataset == "pbmc_10k_v3")]
res = {r.tool: {} for r in hh.itertuples()}
for r in hh.itertuples():
    res[r.tool][r.series] = r.y

crows = []
w_v2, g_v2, el_v2 = rt(BT / "pbmc_10k_v3/peakatail_clipseeded_final_v2/runtime_mem.txt")
w_v1, g_v1, el_v1 = rt(BT / "pbmc_10k_v3/peakatail_clipseeded_final/runtime_mem.txt")
crows.append(dict(tool="PeakATail v1", label="PeakATail v1", wall_h=w_v1, wall_str=el_v1,
                  peak_rss_gb=g_v1, kind="peakatail",
                  note="final Stage-2 clip-seeded arm WITHOUT --ip-filter, code 4efeb125, --threads 16; "
                       "paired with the v2 no-IP arm below so the Stage-1d arrow is like for like",
                  source="results/benchmark_tools/pbmc_10k_v3/peakatail_clipseeded_final/runtime_mem.txt"))
crows.append(dict(tool="PeakATail v2", label="PeakATail v2", wall_h=w_v2, wall_str=el_v2,
                  peak_rss_gb=g_v2, kind="peakatail",
                  note="code 9dfdefb, clip-seeded arm WITHOUT --ip-filter, --threads 16, measured with all "
                       "four arms running concurrently; this is the 34:37 / 12.53 GB pair quoted in 19 sec.3. "
                       "The IP-filtered default arm used in panels a-c was cheaper on both codes (audit rows below)",
                  source="results/benchmark_tools/pbmc_10k_v3/peakatail_clipseeded_final_v2/runtime_mem.txt"))
crows.append(dict(tool="PeakATail v2 uncontended", label="", wall_h=27 / 60 + 43 / 3600, wall_str="27:43",
                  peak_rss_gb=g_v2, kind="peakatail_uncontended",
                  note="single uncontended run quoted in 19 sec.3 (peak RSS not re-measured; v2 value reused)",
                  source="manuscript/19_final_gate_v2.md section 3"))
crows.append(dict(tool="PeakATail shipped", label="PeakATail shipped (coverage-only, v1)",
                  wall_h=res["PeakATail"]["runtime_s"], wall_str="3:27:05",
                  peak_rss_gb=res["PeakATail"]["rss_gib"] * 1048576 / 1e6, kind="shipped",
                  note="pre-clip-seeding caller, v1-era benchmark run",
                  source="results/figures/manuscript/benchmark_headtohead.tsv (verified, manuscript/05)"))
# the IP-filtered arms -- the call sets panels a-c actually use.  NOT plotted (panel d shows the
# no-IP pair, which is the like-for-like v1 -> v2 comparison and the pair 19 sec.3 / 15 sec.5 quote),
# but recorded here so the arm identity of every compute number is auditable.
w_v2i, g_v2i, el_v2i = rt(BT / "pbmc_10k_v3/peakatail_clipseeded_final_v2_ipfilt/runtime_mem.txt")
w_v1i, g_v1i, el_v1i = rt(BT / "pbmc_10k_v3/peakatail_clipseeded_final_ipfilt/runtime_mem.txt")
crows.append(dict(tool="PeakATail v1 IP arm", label="", wall_h=w_v1i, wall_str=el_v1i,
                  peak_rss_gb=g_v1i, kind="audit_only",
                  note="AUDIT ROW, not plotted: the IP-filtered v1 arm (the v1 default call set); "
                       "matches 15 sec.5 (3:44:37 / 239.5 GB)",
                  source="results/benchmark_tools/pbmc_10k_v3/peakatail_clipseeded_final_ipfilt/runtime_mem.txt"))
crows.append(dict(tool="PeakATail v2 IP arm", label="", wall_h=w_v2i, wall_str=el_v2i,
                  peak_rss_gb=g_v2i, kind="audit_only",
                  note="AUDIT ROW, not plotted: the IP-filtered v2 arm whose calls are in panels a-c; "
                       "the 32:59 / 11.12 GB pair of 19 sec.3, also concurrent",
                  source="results/benchmark_tools/pbmc_10k_v3/peakatail_clipseeded_final_v2_ipfilt/runtime_mem.txt"))
# caveats carried over verbatim from the RUNTIME registry of benchmark_headtohead.py, so the
# compute audit trail does not lose the retry / skipped-stage disclosures of the source runs
COMP_NOTE = {
    "polyApipe": "single clean run, exit 0",
    "scAPAtrap": "wall_h = the successful RESUMED run (4:04:29), which skipped three stages already "
                 "done by an OOM-killed first attempt, so it UNDERSTATES a from-scratch run "
                 "(12:58:50 total machine time across both); disclosed, not summed",
    "SCAPTURE": "two mandatory stages of one run (PAScall 8:57:58 + PASquant 3:15:44); annotation "
                "prebuild not timed",
    "scUTRquant": "25:52.06 snakemake run + 4:53.13 salvage rerun of 4 jobs after a missing-nbformat "
                  "exit 1; both attempts summed, none hidden",
    "Sierra": "regtools junctions extract 18:53.30 + FindPeaks/CountPeaks 2:25:04, one clean run",
}
for t, lab in (("polyApipe", "polyApipe"), ("scAPAtrap", "scAPAtrap"), ("SCAPTURE", "SCAPTURE"),
               ("scUTRquant", "scUTRquant*"), ("Sierra", "Sierra")):
    crows.append(dict(tool=t, label=lab, wall_h=res[t]["runtime_s"], wall_str="",
                      peak_rss_gb=res[t]["rss_gib"] * 1048576 / 1e6, kind="competitor",
                      note="competitor's own verified run on the same box (v1-era benchmark). "
                           + COMP_NOTE[t],
                      source="results/figures/manuscript/benchmark_headtohead.tsv (verified, manuscript/05) "
                             "-> its RUNTIME registry entry, sourced from the tool's own run log"))
comp = pd.DataFrame(crows)
comp.insert(0, "panel", "d")

# ---------------------------------------------------------------------------
# FIGURE
# ---------------------------------------------------------------------------
# 2026-09-02 submission pass: the on-figure title and footer moved to the
# sidecar Legend; the canvas height shrank by the freed bands (10.0 -> 8.1 in,
# same axes sizes in inches) and the width narrowed to the double-column norm
# (7.4 -> 7.09 in = 180 mm).
FIG_W, FIG_H = 7.09, 8.1
fig = plt.figure(figsize=(FIG_W, FIG_H))
gs_top = fig.add_gridspec(2, 2, left=0.084, right=0.988, top=0.969, bottom=0.506,
                          width_ratios=[1.02, 0.98], hspace=0.46, wspace=0.235)
GS_BOT_BOTTOM = 0.068
gs_bot = fig.add_gridspec(1, 2, left=0.157, right=0.988, top=0.416, bottom=GS_BOT_BOTTOM,
                          width_ratios=[1.0, 0.90], wspace=0.34)
axA = fig.add_subplot(gs_top[0:2, 0])
axB1 = fig.add_subplot(gs_top[0, 1])
axB2 = fig.add_subplot(gs_top[1, 1])
axC = fig.add_subplot(gs_bot[0, 0])
axD = fig.add_subplot(gs_bot[0, 1])
for ax in (axA, axB1, axB2, axC, axD):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.spines["left"].set_color(MUTED)
    ax.spines["bottom"].set_color(MUTED)
    ax.tick_params(colors=MUTED, length=2.4, width=0.7)


def panel_tag(ax, letter, title, pad=1.015):
    x = 0.0
    if letter:
        ax.text(-0.004, pad, letter, transform=ax.transAxes, fontsize=10.5, fontweight="bold",
                va="bottom", ha="left", color=INK)
        x = 0.052
    ax.text(x, pad, title, transform=ax.transAxes, fontsize=7.7,
            va="bottom", ha="left", color=INK)


# ---- panel a: trade surface ------------------------------------------------
axA.grid(True, color=GRID, lw=0.5, zorder=0)
axA.set_axisbelow(True)
nmin, nmax = sweep.n.min(), sweep.n.max()


def msize(n):
    return 20 + 175 * (np.log10(n) - np.log10(nmin)) / (np.log10(nmax) - np.log10(nmin))


for ds in ("pbmc", "mouse1", "mouse2"):
    s = sweep[sweep.dataset == ds].sort_values("min_molecules")
    c = DSCOL[ds]
    axA.plot(s.recall_detected_genes.to_numpy(), s.atlas_agreement_precision.to_numpy(),
             "-", color=c, lw=1.15, zorder=3, alpha=0.85)
    axA.scatter(s.recall_detected_genes, s.atlas_agreement_precision, s=[msize(n) for n in s.n],
                facecolor=c, edgecolor="white", linewidth=0.7, zorder=4,
                label=f"{DSLAB[ds]} (IP arm)")
    d = s[s.min_molecules == PREREG_K].iloc[0]
    axA.scatter([d.recall_detected_genes], [d.atlas_agreement_precision], s=msize(d.n) + 200,
                facecolor="none", edgecolor=INK, linewidth=1.2, zorder=5)

# threshold labels on the PBMC path only (the mouse paths carry the same five, same order).
# Two labels move out of the default below-right slot: >=5 would sit on the mouse-2 curve, and
# >=3's slot is crossed diagonally by the PBMC >=3 -> >=2 segment.
PBMC_LBL = {5: (+0.0055, +0.008, "left", "bottom"), 3: (-0.005, -0.008, "right", "top")}
for r in sweep[sweep.dataset == "pbmc"].itertuples():
    dx, dy, ha, va = PBMC_LBL.get(int(r.min_molecules), (+0.0055, -0.024, "left", "top"))
    axA.text(r.recall_detected_genes + dx, r.atlas_agreement_precision + dy,
             f"≥{r.min_molecules}", fontsize=6.4, color=BLUE, ha=ha, va=va, zorder=6)
for _k in (1, 10):   # the mouse paths carry the same five thresholds; label their endpoints
    _r = sweep[(sweep.dataset == "mouse1") & (sweep.min_molecules == _k)].iloc[0]
    axA.text(_r.recall_detected_genes - 0.005, _r.atlas_agreement_precision, f"≥{_k}",
             fontsize=6.4, color=GREEN, ha="right", va="center", zorder=6)

axA.axhline(GATE_P, color=INK, lw=0.95, ls=(0, (4, 2)), zorder=2)
axA.text(0.303, GATE_P + 0.010, "pre-registered gate  P ≥ 0.50  (13 §1)",
         ha="right", va="bottom", fontsize=6.2, color=INK)
axA.axhline(ORIG_P, color=MUTED, lw=0.85, ls=(0, (1, 2)), zorder=2)
axA.text(0.303, ORIG_P + 0.010, "original two-sided gate  P ≥ 0.38",
         ha="right", va="bottom", fontsize=6.2, color=MUTED)

d = sweep[(sweep.dataset == "pbmc") & (sweep.min_molecules == PREREG_K)].iloc[0]
# (the "only point on any path that was pre-specified; the rest are
#  descriptive" qualification moved to the sidecar Legend, panel-a paragraph;
#  the ring keeps its short label)
axA.annotate("pre-registered default:\n≥2 molecules (13 §1)",
             xy=(d.recall_detected_genes, d.atlas_agreement_precision),
             xytext=(0.106, 0.588), fontsize=6.1, color=INK, ha="center", va="center",
             bbox=dict(boxstyle="round,pad=0.32", facecolor="white", edgecolor=MUTED, lw=0.55),
             arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8,
                             connectionstyle="arc3,rad=-0.20"), zorder=7)
_nl = {ds: float(sweep[(sweep.dataset == ds) & (sweep.min_molecules == PREREG_K)].null_P_mean.iloc[0])
       for ds in ("pbmc", "mouse1", "mouse2")}
axA.text(0.0575, 0.307,
         "3-seed gene-body-shuffled null P@100 at the default:\n"
         f"{_nl['pbmc']:.4f} (PBMC), {_nl['mouse1']:.4f} / {_nl['mouse2']:.4f} (mice)",
         fontsize=6.0, color=MUTED, ha="left", va="bottom")
axA.set_xlim(0.055, 0.305)
axA.set_ylim(0.30, 1.0)
axA.set_xlabel("detected-gene recall R$_{det}$ @100 bp")
axA.set_ylabel("atlas-agreement precision @100 bp")
leg = axA.legend(loc="upper right", bbox_to_anchor=(1.005, 1.005), frameon=False,
                 handletextpad=0.35, borderpad=0.2, labelspacing=0.30, scatterpoints=1)
for h in leg.legend_handles:
    h.set_sizes([26])
sz = [axA.scatter([], [], s=msize(n), facecolor=GREY, edgecolor="white", linewidth=0.6,
                  label=f"{n // 1000}k") for n in (10000, 50000, 168000)]
leg2 = axA.legend(handles=sz, loc="upper right", bbox_to_anchor=(1.005, 0.845), frameon=False,
                  handletextpad=0.5, borderpad=0.2, labelspacing=0.62, scatterpoints=1,
                  title="n sites called", title_fontsize=6.0, alignment="left")
axA.add_artist(leg)
panel_tag(axA, "a", "trade surface: molecule-support sweep (v2 IP arms)")

# ---- panels b1 / b2: matching-window sensitivity ---------------------------
plotted = cut[(cut.dataset == "pbmc") & (cut.cutoff_bp.isin(CUTOFFS))]
ORDER = ["PeakATail v2 default", "scUTRquant*", "SCAPTURE", "PeakATail v2 tier-1 ≥1 mol",
         "polyApipe", "Sierra", "scAPAtrap"]
SHORT = {"PeakATail v2 default": "PeakATail ≥2 mol (default)",
         "PeakATail v2 tier-1 ≥1 mol": "PeakATail ≥1 mol",
         "scUTRquant*": "scUTRquant* (catalog)", "SCAPTURE": "SCAPTURE",
         "polyApipe": "polyApipe", "Sierra": "Sierra", "scAPAtrap": "scAPAtrap"}
for ax, col, ylab, title, letter in (
        (axB1, "atlas_agreement_precision", "atlas-agreement precision",
         "precision vs matching window (PBMC)", "b"),
        (axB2, "recall_detected_genes", "detected-gene recall R$_{det}$",
         "detected-gene recall vs matching window (PBMC)", "")):
    ax.grid(True, color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for tool in ORDER:
        s = plotted[plotted.tool == tool].sort_values("cutoff_bp")
        ls = (0, (3, 1.6)) if tool == "scUTRquant*" else "-"
        wide = tool == "PeakATail v2 default"
        ax.plot(s.cutoff_bp.to_numpy(), s[col].to_numpy(), linestyle=ls, color=TOOLCOL[tool],
                lw=1.7 if wide else 1.0, marker="o", markersize=3.2 if wide else 2.4,
                markeredgecolor="white", markeredgewidth=0.5, zorder=5 if wide else 3,
                label=SHORT[tool])
    ax.set_xscale("log")
    ax.set_xticks(CUTOFFS)
    ax.set_xticklabels([str(c) for c in CUTOFFS])
    ax.minorticks_off()
    ax.set_xlabel("matching window (bp, strand-matched)")
    ax.set_ylabel(ylab)
    panel_tag(ax, letter, title)
axB1.set_ylim(0.0, 1.20)
axB1.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
axB2.set_ylim(0.0, 0.335)   # headroom shrank with the P@10/P@100 block's move to the legend
axB2.set_yticks([0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30])
axB1.legend(loc="upper left", bbox_to_anchor=(-0.012, 1.012), frameon=False, ncol=2,
            handlelength=1.5, handletextpad=0.4, labelspacing=0.20, columnspacing=0.9,
            borderpad=0.1, fontsize=5.9)
_p10 = {t: float(plotted[(plotted.tool == t) & (plotted.cutoff_bp == 10)].atlas_agreement_precision.iloc[0])
        / float(plotted[(plotted.tool == t) & (plotted.cutoff_bp == 100)].atlas_agreement_precision.iloc[0])
        for t in ORDER}
# (the P@10/P@100 window-dependence block and its "not a loose-window artefact"
#  conclusion moved to the sidecar Legend, panel-b paragraph)

# ---- panel c: replicate reproducibility ------------------------------------
CROWS = [("scAPAtrap (v1)", ORANGE, "scAPAtrap"),
         ("Sierra (v1)", PURPLE, "Sierra"),
         ("PeakATail v2 default", BLUE, "PeakATail v2\ndefault (≥2 mol)"),
         ("PeakATail shipped (v1)", GREY, "PeakATail\nshipped"),
         ("PeakATail v2 tier-1 ≥1 mol", SKY, "PeakATail v2\n≥1 mol"),
         ("polyApipe (v1)", VERM, "polyApipe")]
axC.grid(True, axis="x", color=GRID, lw=0.5, zorder=0)
axC.set_axisbelow(True)
bh = 0.30
for i, (tool, c, _lab) in enumerate(CROWS):
    y = len(CROWS) - 1 - i
    s = repro[(repro.tool == tool) & (repro.window_bp == 100)]
    for j, dr in enumerate(("mouse1->mouse2", "mouse2->mouse1")):
        v = float(s[s.direction == dr].concordance.iloc[0])
        axC.barh(y + (0.5 - j) * bh, v, height=bh * 0.94, color=c,
                 alpha=1.0 if j == 0 else 0.55, edgecolor=c, linewidth=0.5, zorder=3)
        axC.text(v + 0.010, y + (0.5 - j) * bh, f"{v:.3f}", va="center", ha="left",
                 fontsize=5.9, color=INK, zorder=4)
    q = repro[(repro.tool == tool) & (repro.window_bp == 25)]
    for j, dr in enumerate(("mouse1->mouse2", "mouse2->mouse1")):
        if len(q):
            axC.plot([float(q[q.direction == dr].concordance.iloc[0])], [y + (0.5 - j) * bh],
                     marker="D", markersize=3.2, color=INK, markeredgecolor="white",
                     markeredgewidth=0.5, zorder=6)
axC.set_yticks(range(len(CROWS)))
axC.set_yticklabels([lab for _t, _c, lab in CROWS][::-1], fontsize=6.3)
axC.set_xlim(0, 1.10)
axC.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
# the band below the bars that held the Jaccard / single-donor text block is
# gone (the block moved to the sidecar Legend), so the y range holds bars only
axC.set_ylim(-0.70, len(CROWS) + 0.85)
axC.set_xlabel("replicate agreement, mouse 1 vs mouse 2\n(fraction of query sites with a strand-matched call ≤100 bp in the other mouse)",
               fontsize=6.2)
axC.plot([0.0095, 0.0095], [-0.36, len(CROWS) - 0.58], color=MUTED, lw=0.8, ls=(0, (1, 2)),
         zorder=2)   # same data span the old axes-fraction axvline covered
axC.legend(handles=[Patch(facecolor=GREY, edgecolor=GREY, label="mouse 1 → mouse 2"),
                    Patch(facecolor=GREY, alpha=0.55, edgecolor=GREY, label="mouse 2 → mouse 1"),
                    Line2D([], [], marker="D", ls="none", color=INK, markersize=3.2,
                           label="same pair, ≤25 bp window"),
                    Line2D([], [], ls=(0, (1, 2)), color=MUTED, lw=0.8,
                           label="chance ≤0.010 (genic shuffle)")],
           loc="upper right", bbox_to_anchor=(1.015, 1.008), frameon=False, ncol=1,
           handlelength=1.3, handletextpad=0.45, labelspacing=0.28, borderpad=0.1, fontsize=5.9)
jd = repro[repro.tool == "PeakATail v2 default"]
j25 = float(jd[jd.window_bp == 25].jaccard_strand_aware.iloc[0])
j100 = float(jd[jd.window_bp == 100].jaccard_strand_aware.iloc[0])
# (the strand-aware-Jaccard values and the "PBMC is a single donor — no human
#  biological replicate exists here" caveat moved to the sidecar Legend,
#  panel-c paragraph)
panel_tag(axC, "c", "biological-replicate reproducibility (GSE104556 testis)")

# ---- panel d: compute ------------------------------------------------------
axD.grid(True, color=GRID, lw=0.5, zorder=0)
axD.set_axisbelow(True)
v1 = comp[comp.tool == "PeakATail v1"].iloc[0]
v2 = comp[comp.tool == "PeakATail v2"].iloc[0]
un = comp[comp.tool == "PeakATail v2 uncontended"].iloc[0]
axD.annotate("", xy=(v2.wall_h, v2.peak_rss_gb), xytext=(v1.wall_h, v1.peak_rss_gb),
             arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=1.3, shrinkA=5, shrinkB=6,
                             connectionstyle="arc3,rad=0.22"), zorder=3)
axD.plot([un.wall_h, v2.wall_h], [un.peak_rss_gb, v2.peak_rss_gb], "-", color=BLUE, lw=0.9,
         alpha=0.6, zorder=3)
axD.scatter([un.wall_h], [un.peak_rss_gb], s=32, facecolor="white", edgecolor=BLUE,
            linewidth=1.1, zorder=4)
# (x-factor, y-factor, ha, va) for each label
LBL = {"PeakATail v1": (0.88, 1.0, "right", "center"),
       "PeakATail v2": (1.0, 0.62, "center", "top"),
       "PeakATail shipped (coverage-only, v1)": (1.0, 0.60, "center", "top"),
       "polyApipe": (1.0, 0.68, "center", "top"),
       "scAPAtrap": (1.12, 1.25, "left", "bottom"),
       "SCAPTURE": (1.0, 1.50, "center", "bottom"),
       "scUTRquant*": (1.14, 0.80, "left", "center"),
       "Sierra": (0.80, 1.0, "right", "center")}
for r in comp.itertuples():
    if r.kind in ("peakatail_uncontended", "audit_only"):
        continue
    key = r.label.split("\n")[0]
    col = BLUE if r.kind == "peakatail" else (INK if r.kind == "shipped" else TOOLCOL[r.label])
    axD.scatter([r.wall_h], [r.peak_rss_gb], s=44, facecolor=col, edgecolor="white",
                linewidth=0.7, zorder=5)
    fx, fy, ha, va = LBL[key]
    axD.text(r.wall_h * fx, r.peak_rss_gb * fy, r.label.replace("\n", " "), fontsize=6.1,
             color=col, ha=ha, va=va, zorder=6)
axD.set_xscale("log")
axD.set_yscale("log")
axD.set_xlim(0.25, 26)
axD.set_ylim(2.6, 1500)
axD.set_xticks([0.5, 1, 2, 4, 8, 16])
axD.set_xticklabels(["0.5", "1", "2", "4", "8", "16"])
axD.set_yticks([4, 10, 30, 100, 300])
axD.set_yticklabels(["4", "10", "30", "100", "300"])
axD.minorticks_off()
axD.set_xlabel("wall time on the PBMC 10k v3 BAM (h, log)")
axD.set_ylabel("peak RSS (GB, log)")
_v1i = comp[comp.tool == "PeakATail v1 IP arm"].iloc[0]
_v2i = comp[comp.tool == "PeakATail v2 IP arm"].iloc[0]
# (the Stage-1d prose — arm identity, "no ≥300 GB node needed", the IP-arm
#  audit pair — moved to the sidecar Legend, panel-d paragraph; the arrow keeps
#  a short value label and the open circle keeps its key)
axD.text(0.015, 0.985,
         f"v1 → v2: {v1.peak_rss_gb:.0f} → {v2.peak_rss_gb:.1f} GB "
         f"({v1.peak_rss_gb / v2.peak_rss_gb:.1f}×), {v1.wall_str} → {v2.wall_str}\n"
         "open circle = uncontended 27:43 (19 §3)",
         transform=axD.transAxes, fontsize=6.0, color=INK, ha="left", va="top", linespacing=1.4)
panel_tag(axD, "d", "compute footprint, PBMC BAM")

# ---- legend text (sidecar only) --------------------------------------------
# The figure-level title line and the footer below both moved to the sidecar
# '## Legend' at the 2026-09-02 submission pass (journal style: the legend
# opens 'Figure 3 | <title>.'); nothing is drawn on the canvas here.  The
# footer string is kept verbatim as `footer` and written into the Legend.
footer = (
    f"Code {CODE} (PRs #93/#96/#97 merged); every PeakATail v2 value is from the v2 run recorded in "
    "manuscript/19_final_gate_v2.md and results/benchmark_tools/final_v2_verify/VERIFIED_v2.md (verifier verdict FIXED). "
    "PeakATail points labelled v1 or 'shipped' (c and d) are the 4efeb125 records (15 §5, 05). "
    "DEFINITIONS — atlas-agreement precision @W = fraction of called sites whose nearest same-strand curated PolyASite 2.0 "
    "representative site lies ≤W bp away (bedtools closest -s -d -t first, point mode); detected-gene recall R_det@W = fraction of "
    "the detected-gene-restricted atlas (285,136 human / 126,686 mouse sites) with a called site ≤W bp on the same strand; "
    "full-atlas recall (569,005 / 301,006 sites) is in the audit TSV. Molecule support = distinct UMI-deduplicated poly(A)-clip "
    "molecules at a site (pas.bed column 5). "
    "PRE-REGISTERED vs DESCRIPTIVE — in panel a only the ≥2-molecule point and its P ≥ 0.50 gate were pre-registered (13 §1, "
    "committed 01:18:59 on 2026-08-21, before the arms started at 02:59:36); ≥1, ≥3, ≥5 and ≥10 molecules are descriptive points "
    "computed for this figure and were never gated. Every panel-a point and both v2 rows of panel c were computed here with "
    "score_tool.py / bedtools using the reference arguments of scripts/benchmark_tools/stage2_final_launch.sh; the ≥1 and ≥2 points "
    "reproduce 19 §1 to 4 dp on all three arms (asserted in the script). "
    "REPLICATION (panel c) — cross-mouse per-gene switch effects ρ = 0.655 (n = 917; null 0.001 ± 0.032) and 6,049 / 9,469 / 11,263 "
    "same-direction replicated PAS per stage pair with 0 in all 15 null pairings come from 18 (testis record on the v1 code). "
    "COMPETITORS — competitor precision, recall, reproducibility and compute are from their own verified runs (15 §3, 09, 05 and "
    "results/benchmark_tools/gse104556/depth_concordance.tsv); they were not re-run for v2, so panels c and d compare v2 PeakATail "
    "against v1-era competitor measurements. scUTRquant* is catalog-based: shown, not ranked with the de novo tools. scTail is absent "
    "(unrunnable on this BAM, R1 = 28 bp); SCAPTURE: mouse 1 plotted (its mouse-2 site-level run completed and scored 0.672 — "
    "15 §3; only the per-cell PASquant step failed, so it has no cross-mouse reproducibility row in c). "
    "SCOPE — PBMC 10k v3 is one donor and one CellRanger BAM; the two mice are one study and one chemistry. Peak RSS is GNU time -v "
    "maximum resident set size as kbytes/1e6, the convention of 19 §3 and 15 §5; the v2 wall time was measured with all four arms "
    "running concurrently. PANEL-d ARM — both ends of the v1→v2 arrow are the arm run WITHOUT --ip-filter (like for like, the pair "
    f"15 §5 / 19 §3 quote); the IP-filtered default arm of panels a–c was cheaper on both codes ({_v1i.wall_str} / "
    f"{_v1i.peak_rss_gb:.1f} GB → {_v2i.wall_str} / {_v2i.peak_rss_gb:.2f} GB, own logs, audit rows in the compute TSV), which also "
    "carries each competitor run's retry/skipped-stage caveats. "
    "Every plotted value: results/figures/manuscript/fig3_tradeoff{,_reproducibility,_compute}.tsv."
)

for ext, kw in (("png", dict(dpi=600)), ("pdf", {})):
    p = FIGDIR / f"{NAME}.{ext}"
    fig.savefig(p, **kw)
    print("wrote", p)
fig.savefig(OUTDIR / f"{NAME}.png", dpi=600)
print("wrote", OUTDIR / f"{NAME}.png")

# ---------------------------------------------------------------------------
# audit TSVs -- every plotted value
# ---------------------------------------------------------------------------
sweep_out = sweep.copy()
sweep_out.insert(3, "tool", "PeakATail v2 (IP arm)")
p = OUTDIR / f"{NAME}.tsv"
pd.concat([sweep_out, cut], ignore_index=True).to_csv(p, sep="\t", index=False,
                                                      float_format="%.6f")
print("wrote", p)
p = OUTDIR / f"{NAME}_reproducibility.tsv"
repro.to_csv(p, sep="\t", index=False, float_format="%.6f")
print("wrote", p)
p = OUTDIR / f"{NAME}_compute.tsv"
comp.to_csv(p, sep="\t", index=False, float_format="%.6f")
print("wrote", p)

# ---------------------------------------------------------------------------
# sidecar caption
# ---------------------------------------------------------------------------
rr = repro[(repro.tool == "PeakATail v2 default") & (repro.window_bp == 100)].concordance
r1 = repro[(repro.tool == "PeakATail v2 tier-1 ≥1 mol") & (repro.window_bp == 100)].concordance
_d12 = (float(repro[(repro.tool == "PeakATail v2 default") & (repro.window_bp == 100) & (repro.direction == "mouse1->mouse2")].concordance.iloc[0])
        - float(repro[(repro.tool == "PeakATail v2 tier-1 ≥1 mol") & (repro.window_bp == 100) & (repro.direction == "mouse1->mouse2")].concordance.iloc[0]))
_d21 = (float(repro[(repro.tool == "PeakATail v2 default") & (repro.window_bp == 100) & (repro.direction == "mouse2->mouse1")].concordance.iloc[0])
        - float(repro[(repro.tool == "PeakATail v2 tier-1 ≥1 mol") & (repro.window_bp == 100) & (repro.direction == "mouse2->mouse1")].concordance.iloc[0]))
_ch = repro[(repro.window_bp == 100) & (repro.tool.isin([t for t, _c, _l in CROWS]))].chance_mean
_ch_lo, _ch_hi = float(_ch.min()), float(_ch.max())
_ch_pt = float(repro[(repro.tool == "PeakATail v2 default") & (repro.window_bp == 100)].chance_mean.max())
_ch_sa = float(repro[(repro.tool == "scAPAtrap (v1)") & (repro.window_bp == 100)].chance_mean.max())
dfl = {ds: sweep[(sweep.dataset == ds) & (sweep.min_molecules == 2)].iloc[0]
       for ds in ("pbmc", "mouse1", "mouse2")}
k10 = {ds: sweep[(sweep.dataset == ds) & (sweep.min_molecules == 10)].iloc[0]
       for ds in ("pbmc", "mouse1", "mouse2")}
cap_md = f"""# Fig 3 — `fig3_tradeoff` caption (generated by `scripts/manuscript_figures/fig3_tradeoff.py`)

## Legend

Figure 3 | Operating-point trade surface, window sensitivity, replicate reproducibility and compute.

**Panel a — trade surface.** Atlas-agreement precision @100 bp (y) against detected-gene recall
R_det @100 bp (x) for the IP-filtered v2 arms as the molecule-support threshold sweeps ≥1, ≥2, ≥3, ≥5,
≥10 distinct clip molecules, on PBMC 10k v3 and both GSE104556 testis mice. Marker area scales with
the number of sites called (n {int(sweep.n.min()):,}–{int(sweep.n.max()):,}). Thresholds are labelled on the PBMC path and at both ends of
the mouse 1 path; all three paths carry the same five thresholds in the same order. Black rings mark the **pre-registered
default (≥2 molecules, `13_reliability_positioning.md` §1)** — **only that point was pre-registered;
≥1, ≥3, ≥5 and ≥10 are descriptive points computed for this figure and were never gated.** (The
post-hoc sweep that informed the ≥2-molecule threshold is disclosed in 12 CORRECTION and is never
cited as a result.) Dashed
line: the pre-registered gate P ≥ 0.50. Dotted line: the original two-sided gate's precision side,
P ≥ 0.38. Every point was recomputed with `scripts/benchmark_tools/score_tool.py` from `pas.bed`
column 5 using the reference arguments of `scripts/benchmark_tools/stage2_final_launch.sh`; the ≥1 and
≥2 points reproduce `19_final_gate_v2.md` §1 to 4 decimal places on all three arms (asserted in the
script). Defaults: PBMC {dfl['pbmc'].atlas_agreement_precision:.4f} / {dfl['pbmc'].recall_detected_genes:.4f}
(n {int(dfl['pbmc'].n):,}); mouse 1 {dfl['mouse1'].atlas_agreement_precision:.4f} / {dfl['mouse1'].recall_detected_genes:.4f};
mouse 2 {dfl['mouse2'].atlas_agreement_precision:.4f} / {dfl['mouse2'].recall_detected_genes:.4f}. Raising the threshold to
≥10 molecules would reach P {k10['pbmc'].atlas_agreement_precision:.3f} / {k10['mouse1'].atlas_agreement_precision:.3f} /
{k10['mouse2'].atlas_agreement_precision:.3f} at R_det {k10['pbmc'].recall_detected_genes:.3f} / {k10['mouse1'].recall_detected_genes:.3f} /
{k10['mouse2'].recall_detected_genes:.3f} — reported as the shape of the surface, not as a proposal.

**Panel b — matching-window sensitivity.** Atlas-agreement precision (top) and detected-gene recall
(bottom) at 10 / 25 / 50 / 100 bp for the v2 default, the v2 ≥1-molecule sensitivity arm and the five
competitors on PBMC, from the verified per-tool score TSVs (values at 200 bp and the mouse arms are in
the audit TSV). The check the panel exists for: the de novo precision ordering does not change when the
window is tightened to 10 bp, and the ratio P@10 / P@100 is {_p10['PeakATail v2 default']:.2f} for the PeakATail default and
{_p10['polyApipe']:.2f} for polyApipe against {_p10['SCAPTURE']:.2f} (SCAPTURE), {_p10['Sierra']:.2f} (Sierra) and {_p10['scAPAtrap']:.2f} (scAPAtrap) —
the precision lead is not an artefact of a loose window. Catalog-based scUTRquant* is above the
PeakATail default at every window ({_p10['scUTRquant*']:.2f}) and is shown but not ranked with the de novo tools.

**Panel c — biological-replicate reproducibility.** The two testis mice as biological replicates:
the fraction of one mouse's call set with a strand-matched call in the other within 100 bp
(`bedtools closest -s -d -t first`), both directions, with the ≤25 bp value as a diamond. Computed
here for the v2 default ({rr.min():.3f}–{rr.max():.3f} at 100 bp, {float(repro[(repro.tool == 'PeakATail v2 default') & (repro.window_bp == 25)].concordance.min()):.3f}–{float(repro[(repro.tool == 'PeakATail v2 default') & (repro.window_bp == 25)].concordance.max()):.3f} at 25 bp;
strand-aware Jaccard of ±window site sets {j25:.3f} at ±25 bp and {j100:.3f} at ±100 bp) and for the v2
≥1-molecule arm ({r1.min():.3f}–{r1.max():.3f}) — the ≥2-molecule threshold buys {_d12 * 100:.1f} / {_d21 * 100:.1f} points of replicate
agreement over the ≥1-molecule arm.
Competitor and shipped-caller rows are the v1-era head-to-head run
(`results/benchmark_tools/gse104556/depth_concordance.tsv`, verified in `05_figure_index.md`), so this
panel compares a v2 PeakATail arm against v1-era competitor measurements. Chance ≤0.010 (each query set
against the other replicate's three gene-body-shuffled nulls). The dotted line sits at 0.0095, just above
the **highest** per-tool chance level; the six values span {_ch_lo:.4f}–{_ch_hi:.4f} and are in the audit TSV, because a
denser call set has a higher chance of a nearby match — the PeakATail default's own chance level is
{_ch_pt:.4f} against scAPAtrap's {_ch_sa:.4f}, so the raw ranking is the conservative reading of this panel.
**Honest scope: the v2 default is above
polyApipe and above PeakATail's own shipped caller, but below Sierra (0.790 / 0.834) and scAPAtrap
(0.884 / 0.793) in both directions — and scAPAtrap's lead is measured on a set its own
`reducePeaks(min.cells=10, min.count=10)` step has already depth-cleaned.** Supporting replication evidence from
`18_spermatogenesis_final.md` (**still the `4efeb125` record — its v2 re-run after #96/#97 has not
happened**): cross-mouse per-gene switch effects ρ = 0.655 (n = 917; null
0.001 ± 0.032) and 6,049 / 9,469 / 11,263 same-direction replicated PAS per stage pair, 0 in all 15 null
pairings. **PBMC is a single donor and a single CellRanger BAM, so no human biological replicate exists
in this benchmark.**

**Panel d — compute.** Wall time against peak RSS on the PBMC 10k v3 BAM, log–log. The arrow is the
resolved Stage-1d item: PeakATail v1 ({v1.wall_str}, {v1.peak_rss_gb:.1f} GB, code 4efeb125) → v2
({v2.wall_str}, {v2.peak_rss_gb:.2f} GB, code 9dfdefb) — {v1.peak_rss_gb / v2.peak_rss_gb:.1f}× less peak RSS, so the ≥300 GB node
requirement is gone. **Both ends of that arrow are the clip-seeded arm run *without* `--ip-filter`** —
the like-for-like pair, and the pair `15` §5 and `19` §3 quote. The IP-filtered default arm whose call
sets panels a–c show was cheaper on both codes ({_v1i.wall_str} / {_v1i.peak_rss_gb:.1f} GB →
{_v2i.wall_str} / {_v2i.peak_rss_gb:.2f} GB, read here from the arms' own `runtime_mem.txt`); those two
measurements are audit rows in `fig3_tradeoff_compute.tsv` and are not plotted. The Stage-1d fixes cut peak RSS
{v1.peak_rss_gb:.0f} → {v2.peak_rss_gb:.1f} GB on the PBMC BAM — no ≥300 GB node is needed.
The open circle is the uncontended single-run wall time (27:43, `19` §3); the
filled v2 point was measured with all four arms running concurrently, so peak RSS is the production
number and the wall time must be quoted with the concurrency disclosed (~28–35 min). Competitor and
shipped-caller resources are their own verified runs on the same box
(`results/figures/manuscript/benchmark_headtohead.tsv`, verified in `05_figure_index.md`), and each
competitor row carries that run's retry / skipped-stage caveat in the compute TSV's `note` column —
scAPAtrap's 4.07 h is a resumed run that skipped three completed stages and so understates a
from-scratch run (12:58:50 of total machine time across both attempts). Peak RSS is
GNU `time -v` maximum resident set size reported as kbytes/1e6, the convention of `19` §3 and `15` §5.

**Definitions, pre-registration and scope** (the former on-figure footer, verbatim — it travels with the
legend): {footer}

## Provenance

**Index paragraph.** `fig3_tradeoff` — Fig 3: the reliability trade surface, its window
sensitivity, replicate reproducibility and compute. Sweeping the molecule-support threshold on the
IP-filtered v2 arms moves atlas-agreement precision @100 bp from
{dfl['pbmc'].atlas_agreement_precision:.3f} to {k10['pbmc'].atlas_agreement_precision:.3f} (PBMC) while detected-gene
recall falls {dfl['pbmc'].recall_detected_genes:.3f} → {k10['pbmc'].recall_detected_genes:.3f}; only the pre-registered
≥2-molecule point is a claim. The de novo precision ordering is unchanged at a 10 bp matching window,
the two testis mice agree on {rr.min():.3f}–{rr.max():.3f} of default sites at 100 bp against a ≤0.010 chance level, and
the Stage-1d fixes cut peak RSS {v1.peak_rss_gb:.0f} → {v2.peak_rss_gb:.1f} GB on the PBMC BAM.

**Not drawn — what Fig 3 still lacks.** The long-read re-ranking panel (Spearman ρ between atlas and
Kinnex precision orderings, 0.90 at 7 of 8 arms) is quarantined and was not regenerated for v2, so it is
not in this figure; the fairness table (≥95% gene-proximal for every tool) is a table, not a panel; and no
competitor was re-run on the v2 code, so no v2-vs-v2 reproducibility or compute comparison exists.

Sources: `manuscript/19_final_gate_v2.md` + `results/benchmark_tools/final_v2_verify/VERIFIED_v2.md`
(FIXED), `manuscript/18_spermatogenesis_final.md`, `manuscript/15_final_gate.md` §3/§5,
`manuscript/09_headtohead_results.md`, `manuscript/05_figure_index.md`, and the per-file `source`
columns of `results/figures/manuscript/fig3_tradeoff.tsv`,
`fig3_tradeoff_reproducibility.tsv` and `fig3_tradeoff_compute.tsv`.
"""
p = FIGDIR / f"{NAME}.caption.md"
p.write_text(cap_md)
print("wrote", p)

# ---------------------------------------------------------------------------
# edge check: nothing may be clipped; the outer 8 px of the PNG must be blank
# ---------------------------------------------------------------------------
from PIL import Image
im = np.asarray(Image.open(FIGDIR / f"{NAME}.png").convert("L"))
edge = 8
border = np.concatenate([im[:edge, :].ravel(), im[-edge:, :].ravel(),
                         im[:, :edge].ravel(), im[:, -edge:].ravel()])
ink = int((border < 250).sum())
print(f"edge check: {im.shape[1]}x{im.shape[0]} px, ink pixels in the outer {edge} px = {ink}")
assert ink == 0, "INK IN THE OUTER 8 PX -- something is clipped"

print(f"code {CODE}: PBMC default P@100 {dfl['pbmc'].atlas_agreement_precision:.4f} / "
      f"R_det {dfl['pbmc'].recall_detected_genes:.4f}; cross-mouse agreement@100 "
      f"{rr.min():.3f}-{rr.max():.3f}")

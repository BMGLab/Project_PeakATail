#!/usr/bin/env python3
"""
benchmark_headtohead.py -- the MANUSCRIPT head-to-head figure (supersedes
scripts/manuscript_figures/benchmark_tools_running.py, which only ever covered
polyApipe + Sierra on pbmc_10k_v3 plus provisional Laughney PeakATail rows).

Six tools x two datasets, one scoring path
(scripts/benchmark_tools/score_tool.py), per-dataset detected-gene
denominators, 3-seed gene-body-constrained nulls.  Written-up conclusions:
manuscript/09_headtohead_results.md -- this figure must agree with it.

STEP 1 -- source of truth
    results/figures/manuscript/benchmark_consolidated.tsv
    (dataset tool arm metric cutoff value reference n_called runtime_s rss_gib notes)
    harvested from every results/benchmark_tools/<dataset>/<tool>/score_*.tsv
    plus wall-clock/peak-RSS parsed out of each run log by hand-verified
    RUNTIME registry entries below (each carries its log path in `notes`).
    Retries are NEVER silently summed: `runtime_s` is the CLEAN run and the
    failed/OOM-killed attempts are spelled out in `notes`.

STEP 2 -- panels
    (a) precision vs matching cutoff, faceted human | mouse, genic-shuffle null
        band drawn
    (b) precision-recall @100 bp, marker per dataset, annotated with n_called
        -- the trade surface
    (c) F1@100 bars per tool per dataset; annotation-based scUTRquant is
        visually separated and footnoted as near-tautological, NOT ranked among
        the de novo tools
    (d) replicate reproducibility (mouse only): concordance@100 overall and
        stratified depth<=1 vs depth>=2
    (e) resources: clean-run wall time and peak RSS (separate panels, no dual axis)

    PeakATail gets no visual privilege (colours are assigned alphabetically);
    provisional/caveated arms are hatched and footnoted; the title states the
    finding.

DEPTH SOURCES for panel (d), all mouse:
    polyApipe   GFF peakdepth (carried in the pas.bed score column)
    PeakATail   per-PAS UMI totals = row sums of run/annotated_matrix.mtx
                (rows are 1:1 with run/pasbed.bed lines -- the same mapping the
                top-N ranking arms used, see manuscript/09)
    Sierra      row sums of counts/matrix.mtx.gz joined on counts/sitenames.tsv.gz
    scAPAtrap   per-peak UMI totals from counts.tsv.gz
    SCAPTURE    NONE.  Its BED score column is 0 for every evaluated peak
                (162,346 human / 81,138 mouse1; verified) and this run produced
                no per-PAS count matrix -- stated on the panel, never invented.

Every plotted value is written to results/figures/manuscript/benchmark_headtohead.tsv.
Re-runnable end to end; the depth/concordance step caches into
results/benchmark_tools/gse104556/depth_concordance.tsv.
"""
import os
import subprocess
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
BT = WD / "results/benchmark_tools"
OUTDIR = WD / "results/figures/manuscript"
FIGDIR = WD / "manuscript/figures"
NAME = "benchmark_headtohead"
for d in (OUTDIR, FIGDIR):
    d.mkdir(parents=True, exist_ok=True)

ENV = dict(os.environ, LC_ALL="C")  # tr_TR locale silently corrupts BED sorting
CUTOFFS = [10, 25, 50, 100, 200]

INK, MUTED, GRID, SURFACE = "#1B2429", "#5A6B73", "#D8E0E3", "#FFFFFF"
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9"]


def save_manuscript(fig, name, **kw):
    """PNG 300 dpi to results/figures/manuscript/ + PNG&PDF to manuscript/figures/."""
    p = OUTDIR / f"{name}.png"
    fig.savefig(p, dpi=300, **kw)
    print("wrote", p)
    for ext in ("png", "pdf"):
        p = FIGDIR / f"{name}.{ext}"
        fig.savefig(p, **({"dpi": 300} if ext == "png" else {}), **kw)
        print("wrote", p)


def sh_out(cmd):
    return subprocess.run(cmd, shell=True, check=True, env=ENV,
                          executable="/bin/bash", capture_output=True,
                          text=True).stdout


# ---------------------------------------------------------------------------
# tools, arms, and the score TSVs that back them
# ---------------------------------------------------------------------------
# Colour is assigned by ALPHABETICAL tool order -- a rule chosen so that
# PeakATail receives no visual privilege of any kind.
TOOLS = ["PeakATail", "polyApipe", "scAPAtrap", "SCAPTURE", "scUTRquant", "Sierra"]
COLOR = dict(zip(TOOLS, PALETTE))
DE_NOVO = [t for t in TOOLS if t != "scUTRquant"]

HUMAN, MOUSE = "pbmc_10k_v3", "gse104556"
DS_LABEL = {HUMAN: "human PBMC (10x pbmc_10k_v3)",
            MOUSE: "mouse testis (GSE104556)"}

# arm -> (dataset, tool, score tsv, replicate label, primary?, caveat)
ARMS = {
    # ---------------- human, one CellRanger BAM ----------------
    "peakatail":       (HUMAN, "PeakATail",  "peakatail/score_peakatail.tsv", "", True, ""),
    "polyapipe":       (HUMAN, "polyApipe",  "polyapipe/score_polyapipe.tsv", "", True, ""),
    "scapatrap":       (HUMAN, "scAPAtrap",  "scapatrap/score_scapatrap.tsv", "", True, ""),
    "scapture":        (HUMAN, "SCAPTURE",   "scapture/score_scapture.tsv", "", True, ""),
    "scutrquant":      (HUMAN, "scUTRquant", "scutrquant/score_scutrquant.tsv", "", True,
                        "annotation-based: fixed UTRome catalog filtered by detection"),
    "sierra":          (HUMAN, "Sierra",     "sierra/score_sierra.tsv", "", True, ""),
    # human sensitivity / diagnostic variants (kept in the TSV, not plotted)
    "polyapipe_all":   (HUMAN, "polyApipe",  "polyapipe/score_polyapipe_all.tsv", "", False,
                        "sensitivity variant: includes misprime=True peaks"),
    "scapture_all":    (HUMAN, "SCAPTURE",   "scapture/score_scapture_all.tsv", "", False,
                        "sensitivity variant: all evaluated peaks, no DeepPASS filter"),
    "sierra_summit":   (HUMAN, "Sierra",     "sierra/score_sierra_summit.tsv", "", False,
                        "sensitivity variant: MaxPosition coverage summit instead of 3' of fit; "
                        "DENOMINATOR CAVEAT: this arm alone was scored against the older "
                        "Laughney-derived detected-gene atlas (285,220 sites / 14,851 genes), "
                        "not shared_refs_pbmc (285,136 / 14,949) used by every other human arm, "
                        "so its atlas_detected recall and F1 are NOT comparable to them "
                        "(recall@100 0.064610 as scored; 0.064818 on the PBMC-native set)"),
    "peakatail_top22000":  (HUMAN, "PeakATail", "peakatail/topN/score_peakatail_top22000.tsv", "", False,
                            "diagnostic: top 22k PAS by total UMI (matched-N vs Sierra/SCAPTURE)"),
    "peakatail_top36000":  (HUMAN, "PeakATail", "peakatail/topN/score_peakatail_top36000.tsv", "", False,
                            "diagnostic: top 36k PAS by total UMI"),
    "peakatail_top106000": (HUMAN, "PeakATail", "peakatail/topN/score_peakatail_top106000.tsv", "", False,
                            "diagnostic: top 106k PAS by total UMI"),
    # ---------------- mouse, two biological replicates ----------------
    "peakatail_mouse1":  (MOUSE, "PeakATail",  "peakatail/score_peakatail_mouse1.tsv", "mouse1", True, ""),
    "peakatail_mouse2":  (MOUSE, "PeakATail",  "peakatail/score_peakatail_mouse2.tsv", "mouse2", True, ""),
    "polyapipe_mouse1":  (MOUSE, "polyApipe",  "polyapipe/score_polyapipe_mouse1.tsv", "mouse1", True, ""),
    "polyapipe_mouse2":  (MOUSE, "polyApipe",  "polyapipe/score_polyapipe_mouse2.tsv", "mouse2", True, ""),
    "scapatrap_mouse1":  (MOUSE, "scAPAtrap",  "scapatrap/mouse1/score_scapatrap_mouse1.tsv", "mouse1", True, ""),
    "scapatrap_mouse2":  (MOUSE, "scAPAtrap",  "scapatrap/mouse2/score_scapatrap_mouse2.tsv", "mouse2", True, ""),
    "scapture_mouse1":   (MOUSE, "SCAPTURE",   "scapture/mouse1/score_scapture_mouse1.tsv", "mouse1", True,
                          "CAVEAT: mouse1 only -- the mouse2 run was truncated by disk exhaustion "
                          "(SCAPTURE exited 0 anyway); no replicate pair, so no reproducibility arm"),
    "scutrquant_mouse1": (MOUSE, "scUTRquant", "scutrquant/score_scutrquant_mouse1.tsv", "mouse1", True,
                          "annotation-based: fixed mm10 UTRome catalog filtered by detection"),
    "scutrquant_mouse2": (MOUSE, "scUTRquant", "scutrquant/score_scutrquant_mouse2.tsv", "mouse2", True,
                          "annotation-based: fixed mm10 UTRome catalog filtered by detection"),
    "sierra_mouse1":     (MOUSE, "Sierra",     "sierra/score_sierra_mouse1.tsv", "mouse1", True, ""),
    "sierra_mouse2":     (MOUSE, "Sierra",     "sierra/score_sierra_mouse2.tsv", "mouse2", True, ""),
    "polyapipe_mouse1_all": (MOUSE, "polyApipe", "polyapipe/score_polyapipe_mouse1_all.tsv", "mouse1", False,
                             "sensitivity variant: includes misprime=True peaks"),
    "polyapipe_mouse2_all": (MOUSE, "polyApipe", "polyapipe/score_polyapipe_mouse2_all.tsv", "mouse2", False,
                             "sensitivity variant: includes misprime=True peaks"),
    "sierra_pooled":     (MOUSE, "Sierra",     "sierra/score_sierra_pooled.tsv", "pooled", False,
                          "sensitivity variant: mouse1+mouse2 peak sets pooled"),
}

# ---------------------------------------------------------------------------
# RUNTIME / PEAK-RSS registry
# ---------------------------------------------------------------------------
# Every entry was read out of the run log named in `src`.  `wall_s` is the
# CLEAN (successful) end-to-end wall clock for that arm, INCLUDING mandatory
# input-prep stages this benchmark had to run for the tool; failed attempts,
# OOM kills and resumed stages are disclosed in `note` and are never summed
# into `wall_s`.  `rss_gib` is the maximum peak RSS over the arm's stages.
KIB = 1024 ** 2  # KiB -> GiB
RUNTIME = {
    "peakatail": dict(
        wall_s=3 * 3600 + 27 * 60 + 6, rss_gib=105546020 / KIB,
        src="pbmc_10k_v3/peakatail/runtime_mem.txt + run.log",
        note="4 attempts; attempts 1-3 failed (1:33:21, 1:52:42, 4:52:58) on three separate "
             "CellRanger-input bugs (CB '-1' suffix length check, unmapped-read reference_end=None, "
             "RG-id underscore split; see manuscript/github/issue8_cellranger_cb_suffix.md). "
             "wall_s = attempt 4, a full run from scratch, exit 0. Failed attempts NOT summed in."),
    "polyapipe": dict(
        wall_s=3 * 3600 + 26 * 60 + 37, rss_gib=13723976 / KIB,
        src="pbmc_10k_v3/polyapipe/run.log", note="single clean run, exit 0"),
    "sierra": dict(
        wall_s=2 * 3600 + 25 * 60 + 4 + 1133.30, rss_gib=10050780 / KIB,
        src="pbmc_10k_v3/sierra/run.log + junctions.log",
        note="end-to-end = regtools junctions extract 18:53.30 (2.25 GiB) + "
             "FindPeaks/CountPeaks 2:25:04 (9.59 GiB); one clean run"),
    "scapatrap": dict(
        wall_s=4 * 3600 + 4 * 60 + 29, rss_gib=14999828 / KIB,
        src="pbmc_10k_v3/scapatrap/run_resume.log (clean) + run.log (killed attempt)",
        note="RETRY: attempt 1 ran 8:54:21 and was OOM-killed (exit 137, 13.39 GiB peak). "
             "wall_s = the successful resumed run (4:04:29), which SKIPPED three stages "
             "(findUniqueMap / index / dedupByPos) already completed by the killed attempt, "
             "so it understates a from-scratch run; total machine time across both was 12:58:50. "
             "Not summed -- disclosed."),
    "scapture": dict(
        wall_s=8 * 3600 + 57 * 60 + 58 + 3 * 3600 + 15 * 60 + 44, rss_gib=24160272 / KIB,
        src="pbmc_10k_v3/scapture/run.log",
        note="two mandatory stages of one run: PAScall 8:57:58 (23.04 GiB) + PASquant 3:15:44 "
             "(3.29 GiB); annotation prebuild (anno.log) not timed. PASquant's cell x PAS matrix "
             "came out empty (36,204 PAS rows, zero cells) -- the PAS call set is unaffected but "
             "this tool has no usable per-PAS depth on this dataset"),
    "scutrquant": dict(
        wall_s=1552.06 + 293.13, rss_gib=4197760 / KIB,
        src="pbmc_10k_v3/scutrquant/run.log + salvage_rerun.log",
        note="RETRY: the 25:52.06 snakemake run did all the compute but exited 1 on a missing "
             "nbformat in the report rule; wall_s adds the 4:53.13 salvage rerun of the remaining "
             "4 jobs after installing it. Both attempts disclosed, no failed compute hidden."),
    "peakatail_mouse1": dict(wall_s=2 * 3600 + 6 * 60 + 48, rss_gib=14721828 / KIB,
                             src="gse104556/peakatail/mouse1/runtime_mem.txt", note="clean run"),
    "peakatail_mouse2": dict(wall_s=1 * 3600 + 26 * 60 + 46, rss_gib=14284300 / KIB,
                             src="gse104556/peakatail/mouse2/runtime_mem.txt", note="clean run"),
    "polyapipe_mouse1": dict(wall_s=1 * 3600 + 18 * 60 + 36, rss_gib=7095552 / KIB,
                             src="gse104556/polyapipe/run.log", note="clean run"),
    "polyapipe_mouse2": dict(wall_s=1 * 3600 + 2 * 60 + 59, rss_gib=7090768 / KIB,
                             src="gse104556/polyapipe/run.log", note="clean run"),
    "sierra_mouse1": dict(wall_s=1 * 3600 + 12 * 60 + 48 + 400.81, rss_gib=4601992 / KIB,
                          src="gse104556/sierra/mouse1/{run,junctions}.log",
                          note="regtools junctions extract 6:40.81 + FindPeaks/CountPeaks 1:12:48"),
    "sierra_mouse2": dict(wall_s=57 * 60 + 49.83 + 321.94, rss_gib=4569940 / KIB,
                          src="gse104556/sierra/mouse2/{run,junctions}.log",
                          note="regtools junctions extract 5:21.94 + FindPeaks/CountPeaks 57:49.83"),
    "scapatrap_mouse1": dict(wall_s=3 * 3600 + 41 * 60 + 54, rss_gib=8257072 / KIB,
                             src="gse104556/scapatrap/run_mouse1.log", note="clean run"),
    "scapatrap_mouse2": dict(wall_s=2 * 3600 + 46 * 60 + 17, rss_gib=8143052 / KIB,
                             src="gse104556/scapatrap/run_mouse2.log", note="clean run"),
    "scapture_mouse1": dict(wall_s=8 * 3600 + 9 * 60 + 53 + 19 * 60 + 29.25, rss_gib=4305280 / KIB,
                            src="gse104556/scapture/run.log",
                            note="PAScall 8:09:53 + PASquant 19:29.25; shared annotation prebuild "
                                 "1:02.84 excluded (one-off, both mice)"),
    "scutrquant_mouse1": dict(wall_s=686.0 + 3872.0, rss_gib=4196864 / KIB,
                              src="gse104556/scutrquant/run.log",
                              note="pipeline 11:26 + 1:04:32 STARsolo->CellRanger-tag BAM prep that "
                                   "THIS BENCHMARK imposed (patched kallisto segfaults on raw "
                                   "STARsolo BAMs); runtime is prep-dominated, not pipeline cost"),
    "scutrquant_mouse2": dict(wall_s=615.97 + 3214.29, rss_gib=4197760 / KIB,
                              src="gse104556/scutrquant/run.log",
                              note="pipeline 10:15.97 + 53:34.29 STARsolo->CellRanger-tag BAM prep "
                                   "imposed by this benchmark; prep-dominated"),
}
# arms with no run of their own (derived from another arm's output)
RUNTIME_DERIVED = {
    "polyapipe_all": "polyapipe", "polyapipe_mouse1_all": "polyapipe_mouse1",
    "polyapipe_mouse2_all": "polyapipe_mouse2", "scapture_all": "scapture",
    "sierra_summit": "sierra", "peakatail_top22000": "peakatail",
    "peakatail_top36000": "peakatail", "peakatail_top106000": "peakatail",
}

# Two DIFFERENT kinds of caveat, deliberately not conflated:
#   HATCHED_DATA   the arm's ACCURACY numbers carry a hard caveat
#   HATCHED_RUN    only the arm's RUNTIME carries a caveat (retries / resumes)
HATCHED_DATA = {"scapture_mouse1"}
HATCHED_RUN = {"peakatail", "scapatrap", "scutrquant"}

# ---------------------------------------------------------------------------
# STEP 1 -- harvest every score TSV into one tidy source-of-truth table
# ---------------------------------------------------------------------------
def harvest():
    rows = []
    for arm, (ds, tool, rel, rep, primary, caveat) in ARMS.items():
        p = BT / ds / rel
        if not p.exists():
            raise SystemExit(f"missing score TSV: {p}")
        s = pd.read_csv(p, sep="\t")
        n_called = int(s.loc[(s.panel == "precision") & (s.series == "real") &
                             (s.reference == "atlas_full"), "n_query"].iloc[0])
        rt = RUNTIME.get(arm) or RUNTIME.get(RUNTIME_DERIVED.get(arm, ""), None)
        derived = arm in RUNTIME_DERIVED
        wall = rt["wall_s"] if rt else np.nan
        rss = rt["rss_gib"] if rt else np.nan
        rt_note = ("re-derived from the " + RUNTIME_DERIVED[arm] + " run output, no separate run; "
                   + rt["note"]) if derived and rt else (rt["note"] if rt else "no timed run")
        base_note = "; ".join(x for x in (caveat, f"runtime source: {rt['src']}" if rt else "",
                                          rt_note) if x)

        def emit(metric, ref, cut, val, extra=""):
            rows.append(dict(dataset=ds, tool=tool, arm=arm, metric=metric,
                             cutoff=cut, value=val, reference=ref,
                             n_called=n_called, runtime_s=wall, rss_gib=rss,
                             notes="; ".join(x for x in (base_note, extra) if x)))

        real = s[(s.series == "real")]
        for _, r in real[real.panel == "precision"].iterrows():
            emit("precision", r.reference, int(r.cutoff_bp), float(r.value))
        for _, r in real[real.panel == "recall"].iterrows():
            emit("recall", r.reference, int(r.cutoff_bp), float(r.value),
                 f"denominator={int(r.n_query)}")
        for _, r in real[real.panel == "f1"].iterrows():
            emit("f1", r.reference, int(r.cutoff_bp), float(r.value))
        nul = s[(s.panel == "precision") & (s.series == "null_genic")]
        for cut, g in nul.groupby("cutoff_bp"):
            emit("precision_null", "atlas_full", int(cut), float(g.value.mean()),
                 f"mean of {len(g)} genic-shuffle seeds (min {g.value.min():.4f}, "
                 f"max {g.value.max():.4f})")
        rows.append(dict(dataset=ds, tool=tool, arm=arm, metric="replicate",
                         cutoff=np.nan, value=np.nan, reference=rep or "single-sample",
                         n_called=n_called, runtime_s=wall, rss_gib=rss,
                         notes=("PRIMARY arm" if primary else "SENSITIVITY/DIAGNOSTIC arm, "
                                "not plotted") + ("; " + base_note if base_note else "")))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# STEP 1b -- mouse replicate reproducibility, depth-stratified (panel d)
# ---------------------------------------------------------------------------
# point BEDs are the genome-filtered, sorted, 3'-base-reduced sets score_tool.py
# already wrote; reusing them guarantees this panel matches the scored n exactly
# (verified: PeakATail m1->m2 @100 reproduces replicate_concordance.tsv to 6 dp).
REPRO = {
    "PeakATail": dict(
        work="peakatail/.work_peakatail_{m}/peakatail_{m}.pt.bed",
        depth="peakatail_mtx", depth_unit="total UMIs (annotated_matrix.mtx row sums)"),
    "polyApipe": dict(
        work="polyapipe/.work_polyapipe_{m}/polyapipe_{m}.pt.bed",
        depth="score_col", depth_unit="poly(A) reads (GFF peakdepth)"),
    "Sierra": dict(
        work="sierra/.work_sierra_{m}/sierra_{m}.pt.bed",
        depth="sierra_counts", depth_unit="total UMIs (CountPeaks matrix row sums)"),
    "scAPAtrap": dict(
        work="scapatrap/{m}/.work_scapatrap_{m}/scapatrap_{m}.pt.bed",
        depth="scapatrap_counts", depth_unit="total UMIs (counts.tsv.gz per-peak sums)"),
}
CACHE = BT / MOUSE / "depth_concordance.tsv"


def _mtx_rowsums(path, gz=False):
    """Row sums of a MatrixMarket coordinate file (comment lines start with %)."""
    cat = "zcat" if gz else "cat"
    out = sh_out(f"{cat} {path} | awk '/^%/{{next}} !d{{d=$1; next}} "
                 f"{{s[$1]+=$3}} END{{for(i=1;i<=d;i++) print s[i]+0}}'")
    return [int(x) for x in out.split()]


def depth_map(tool, m):
    """name -> depth for one tool/replicate."""
    kind = REPRO[tool]["depth"]
    if kind == "score_col":
        return None  # depth already in the BED score column
    if kind == "peakatail_mtx":
        names = [l.split("\t")[3] for l in
                 (BT / MOUSE / f"peakatail/{m}/run/pasbed.bed").read_text().splitlines()]
        sums = _mtx_rowsums(BT / MOUSE / f"peakatail/{m}/run/annotated_matrix.mtx")
        assert len(names) == len(sums), (len(names), len(sums))
        return dict(zip(names, sums))
    if kind == "sierra_counts":
        c = BT / MOUSE / f"sierra/{m}/counts"
        names = sh_out(f"zcat {c/'sitenames.tsv.gz'}").split()
        sums = _mtx_rowsums(c / "matrix.mtx.gz", gz=True)
        assert len(names) == len(sums), (len(names), len(sums))
        return dict(zip(names, sums))
    if kind == "scapatrap_counts":
        out = sh_out(f"zcat {BT/MOUSE/f'scapatrap/{m}/counts.tsv.gz'} | "
                     f"awk -F'\\t' 'NR>1{{s[$1]+=$3}} END{{for(k in s) print k\"\\t\"s[k]}}'")
        return {l.split("\t")[0]: int(l.split("\t")[1]) for l in out.splitlines()}
    raise ValueError(kind)


def repro_table(force=False):
    if CACHE.exists() and not force:
        return pd.read_csv(CACHE, sep="\t")
    tmp = Path(os.environ.get("TMPDIR", "/mnt/ssd0/emaout/peakatail_benchmark/scratch"))
    tmp.mkdir(parents=True, exist_ok=True)
    beds = {}
    for tool, cfg in REPRO.items():
        for m in ("mouse1", "mouse2"):
            src = BT / MOUSE / cfg["work"].format(m=m)
            dst = tmp / f"hh_{tool}_{m}.bed"
            dm = depth_map(tool, m)
            if dm is None:
                dst.write_text(src.read_text())
            else:
                with open(src) as fh, open(dst, "w") as out:
                    for line in fh:
                        f = line.rstrip("\n").split("\t")
                        out.write("\t".join([f[0], f[1], f[2], f[3],
                                             str(dm.get(f[3], -1)), f[5]]) + "\n")
            beds[(tool, m)] = dst

    rows = []
    for tool, cfg in REPRO.items():
        for a, b in (("mouse1", "mouse2"), ("mouse2", "mouse1")):
            # observed concordance, overall + depth-stratified
            out = sh_out(
                f"bedtools closest -s -d -t first -a {beds[(tool,a)]} -b {beds[(tool,b)]} "
                f"2>/dev/null | awk -F'\\t' '{{dep=$5; dd=$NF; n++; hit=(dd>=0&&dd<=100); "
                f"m+=hit; if(dep<0){{u++; next}} if(dep<=1){{n1++; m1+=hit}} "
                f"else {{n2++; m2+=hit}}}} END{{print n,m,n1+0,m1+0,n2+0,m2+0,u+0}}'")
            n, mm, n1, m1, n2, m2, u = (int(x) for x in out.split())
            # chance level: same query vs the 3 genic-shuffle nulls of the OTHER replicate
            wdir = (BT / MOUSE / cfg["work"].format(m=b)).parent
            stem = wdir.name.replace(".work_", "")
            null = []
            for seed in (1, 2, 3):
                o = sh_out(f"bedtools closest -s -d -t first -a {beds[(tool,a)]} "
                           f"-b {wdir}/{stem}.null.s{seed}.pt.bed 2>/dev/null | "
                           f"awk -F'\\t' '{{n++; if($NF>=0&&$NF<=100)m++}} END{{print m/n}}'")
                null.append(float(o))
            rows.append(dict(tool=tool, direction=f"{a}->{b}", n_query=n,
                             conc_all=mm / n, n_depth1=n1, conc_depth1=(m1 / n1 if n1 else np.nan),
                             n_depth2plus=n2, conc_depth2plus=(m2 / n2 if n2 else np.nan),
                             n_no_depth=u, frac_depth1=(n1 / (n1 + n2) if n1 + n2 else np.nan),
                             conc_null=float(np.mean(null)), depth_unit=cfg["depth_unit"]))
    df = pd.DataFrame(rows)
    df.to_csv(CACHE, sep="\t", index=False, float_format="%.6f")
    print("wrote", CACHE)
    return df


# ---------------------------------------------------------------------------
def hhmm(sec):
    h, r = divmod(int(round(sec)), 3600)
    return f"{h}:{r//60:02d}"


def main():
    cons = harvest()
    cons_path = OUTDIR / "benchmark_consolidated.tsv"
    cons.to_csv(cons_path, sep="\t", index=False, float_format="%.6f")
    print("wrote", cons_path, len(cons), "rows")

    rep = repro_table()

    prim = {a for a, v in ARMS.items() if v[4]}
    plotted = []            # every value drawn, dumped to benchmark_headtohead.tsv

    def get(arm, metric, ref, cut):
        q = cons[(cons.arm == arm) & (cons.metric == metric) &
                 (cons.reference == ref) & (cons.cutoff == cut)]
        return float(q.value.iloc[0]) if len(q) else np.nan

    def ncall(arm):
        return int(cons[cons.arm == arm].n_called.iloc[0])

    def arms_of(ds, tool=None):
        return [a for a in ARMS if a in prim and ARMS[a][0] == ds
                and (tool is None or ARMS[a][1] == tool)]

    fig = plt.figure(figsize=(16.6, 19.2), facecolor=SURFACE)
    gs = fig.add_gridspec(3, 4, left=0.055, right=0.985, top=0.855, bottom=0.215,
                          hspace=0.42, wspace=0.40, height_ratios=[1.0, 1.06, 1.0])
    axA = [fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4])]
    axB = fig.add_subplot(gs[1, 0:2])
    axC = fig.add_subplot(gs[1, 2:4])
    axD = fig.add_subplot(gs[2, 0:2])
    axE1 = fig.add_subplot(gs[2, 2])
    axE2 = fig.add_subplot(gs[2, 3])
    for ax in axA + [axB, axC, axD, axE1, axE2]:
        ax.set_facecolor(SURFACE)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(MUTED)
        ax.tick_params(colors=MUTED, labelsize=9)
        ax.grid(True, color=GRID, lw=0.7, zorder=0)
        ax.set_axisbelow(True)

    # ------------------------------------------------------------------ (a)
    for i, (ax, ds) in enumerate(zip(axA, (HUMAN, MOUSE))):
        nul = []
        for arm in arms_of(ds):
            tool = ARMS[arm][1]
            ys = [get(arm, "precision", "atlas_full", c) for c in CUTOFFS]
            nul.append([get(arm, "precision_null", "atlas_full", c) for c in CUTOFFS])
            rep_lbl = ARMS[arm][3]
            ls = "--" if rep_lbl == "mouse2" else "-"
            ax.plot(CUTOFFS, ys, ls, color=COLOR[tool], lw=2.1,
                    marker=("s" if rep_lbl == "mouse2" else "o"), ms=5.5,
                    mfc="none" if tool == "scUTRquant" else COLOR[tool],
                    zorder=3, alpha=0.95)
            for c, y in zip(CUTOFFS, ys):
                plotted.append(dict(panel="a_precision_vs_cutoff", dataset=ds, tool=tool,
                                    arm=arm, x=c, y=y, series="precision_vs_atlas_full"))
        nul = np.array(nul, float)
        lo, hi = np.nanmin(nul, 0), np.nanmax(nul, 0)
        ax.fill_between(CUTOFFS, lo, hi, color=MUTED, alpha=0.35, zorder=1, lw=0)
        ax.plot(CUTOFFS, np.nanmean(nul, 0), color=MUTED, lw=1.2, ls=":", zorder=2)
        for c, l, h in zip(CUTOFFS, lo, hi):
            plotted.append(dict(panel="a_precision_vs_cutoff", dataset=ds, tool="(null)",
                                arm="genic_shuffle_null", x=c, y=(l + h) / 2,
                                series=f"null band {l:.4f}-{h:.4f} (envelope of per-arm 3-seed means)"))
        # the null band is labelled by a legend, not an arrow, so it never
        # crosses a tool's curve; each panel uses whichever corner is empty
        ax.legend(handles=[Patch(facecolor=MUTED, alpha=0.35, edgecolor="none",
                                 label="genic-shuffle null (envelope of per-arm 3-seed means)")],
                  loc="upper left" if ds == HUMAN else "lower right",
                  frameon=False, fontsize=8.4, labelcolor=MUTED)
        ax.set_xscale("log")
        ax.set_xticks(CUTOFFS)
        ax.set_xticklabels([str(c) for c in CUTOFFS])
        ax.set_xlabel("matching cutoff to nearest PolyASite 2.0 rep site (bp, log)", fontsize=9.5)
        ax.set_ylabel("precision", fontsize=9.5)
        ax.set_ylim(0, 0.90)
        ax.set_title(f"(a{i+1}) precision vs cutoff — {DS_LABEL[ds]}",
                     fontsize=11.5, color=INK, loc="left", pad=8)
    axA[1].annotate("solid = mouse1, dashed = mouse2\nSCAPTURE: mouse1 only (mouse2 run invalid)",
                    xy=(0.02, 0.99), xycoords="axes fraction", va="top", ha="left",
                    fontsize=8.4, color=MUTED)

    # ------------------------------------------------------------------ (b)
    MK = {HUMAN: "o", MOUSE: "^"}
    # label offsets in points, tuned by hand so nothing overlaps
    OFF_H = {"PeakATail": (10, -2, "left"), "polyApipe": (-9, 11, "right"),
             "scAPAtrap": (-11, -11, "right"), "SCAPTURE": (-11, 1, "right"),
             "scUTRquant": (-11, 2, "right"), "Sierra": (10, -3, "left")}
    OFF_M = {"PeakATail": (-11, -12, "right"), "polyApipe": (10, 8, "left"),
             "scAPAtrap": (12, 1, "left"), "SCAPTURE": (10, 3, "left"),
             "scUTRquant": (11, 0, "left"), "Sierra": (-10, 8, "right")}
    for ds in (HUMAN, MOUSE):
        for tool in TOOLS:
            arms = arms_of(ds, tool)
            xs = [get(a, "recall", "atlas_detected", 100) for a in arms]
            ys = [get(a, "precision", "atlas_full", 100) for a in arms]
            ns = [ncall(a) for a in arms]
            for a, x, y, n in zip(arms, xs, ys, ns):
                axB.scatter([x], [y], s=170, marker=MK[ds], zorder=4,
                            facecolor="none" if tool == "scUTRquant" else COLOR[tool],
                            edgecolor=COLOR[tool] if tool == "scUTRquant" else INK,
                            linewidth=2.0 if tool == "scUTRquant" else 0.8,
                            hatch="////" if a in HATCHED_DATA else None)
                plotted.append(dict(panel="b_pr_at100", dataset=ds, tool=tool, arm=a,
                                    x=x, y=y, series=f"n_called={n}"))
            if len(arms) == 2:           # join the mouse replicate pair, label once
                axB.plot(xs, ys, color=COLOR[tool], lw=1.1, alpha=0.6, zorder=3)
            lo_n, hi_n = min(ns) / 1000, max(ns) / 1000
            ntxt = f"{lo_n:.0f}k" if abs(hi_n - lo_n) < 1 else f"{lo_n:.0f}–{hi_n:.0f}k"
            dx, dy, ha = (OFF_H if ds == HUMAN else OFF_M)[tool]
            axB.annotate(f"{tool}\n{ntxt}", (float(np.mean(xs)), float(np.mean(ys))),
                         textcoords="offset points", xytext=(dx, dy), ha=ha,
                         fontsize=8.6, color=INK, va="center")
    axB.set_xlabel("recall @100 bp vs the detected-gene-restricted atlas\n"
                   "(human 285,136 sites | mouse 126,686 sites)", fontsize=9.5)
    axB.set_ylabel("precision @100 bp vs PolyASite 2.0", fontsize=9.5)
    axB.set_xlim(0.055, 0.345)
    axB.set_ylim(0, 0.92)
    axB.set_title("(b) the trade surface: precision–recall @100 bp, labelled with n called",
                  fontsize=11.5, color=INK, loc="left", pad=8)
    axB.legend(handles=[Line2D([], [], ls="", marker="o", mfc=MUTED, mec=INK, ms=9,
                               label="human PBMC"),
                        Line2D([], [], ls="", marker="^", mfc=MUTED, mec=INK, ms=9,
                               label="mouse testis (line joins m1–m2)"),
                        Line2D([], [], ls="", marker="o", mfc="none", mec=MUTED, mew=2.0,
                               ms=9, label="annotation-based (scUTRquant)"),
                        Patch(facecolor="white", edgecolor=INK, hatch="////",
                              label="caveated arm (see footnotes)")],
               loc="lower left", frameon=False, fontsize=8.8)

    # ------------------------------------------------------------------ (c)
    order = ["polyApipe", "SCAPTURE", "Sierra", "scAPAtrap", "PeakATail"]  # PBMC F1 rank
    slots = [(t, float(i)) for i, t in enumerate(order)]
    slots.append(("scUTRquant", float(len(order)) + 0.85))
    labels = order + ["scUTRquant*"]
    W = 0.28
    for tool, xc in slots:
        entries = [("human PBMC", -W, arms_of(HUMAN, tool)[0])]
        marms = arms_of(MOUSE, tool)
        for i, a in enumerate(marms):
            entries.append((f"mouse{ARMS[a][3][-1]}", (0.0 if len(marms) == 1 else i * W), a))
        for lbl, off, a in entries:
            v = get(a, "f1", "atlas_detected", 100)
            axC.bar(xc + off, v, W * 0.92, color=COLOR[tool], zorder=3,
                    edgecolor=INK, linewidth=0.6,
                    alpha=1.0 if "human" in lbl else 0.62,
                    hatch="////" if a in HATCHED_DATA else None)
            axC.text(xc + off, v + 0.007, f"{v:.3f}", ha="center", va="bottom",
                     fontsize=7.8, color=INK, rotation=90)
            plotted.append(dict(panel="c_f1_at100", dataset=ARMS[a][0], tool=tool, arm=a,
                                x=xc + off, y=v, series=lbl))
    axC.axvline(slots[-1][1] - 0.72, color=MUTED, lw=1.0, ls=":")
    axC.set_xticks([s[1] for s in slots])
    axC.set_xticklabels(labels, fontsize=9.5, color=INK)
    axC.set_xlim(-0.72, slots[-1][1] + 0.72)
    axC.set_ylabel("F1 @100 bp (detected-gene denominator)", fontsize=9.5)
    axC.set_ylim(0, 0.60)
    axC.set_title("(c) F1 @100 bp — de novo tools ranked by PBMC F1; scUTRquant* separated, not ranked",
                  fontsize=11.5, color=INK, loc="left", pad=8)
    axC.annotate("* annotation-based: a fixed UTRome\ncatalog filtered by detection. Precision\n"
                 "against an annotation-derived atlas\nis near-tautological, so scUTRquant\n"
                 "is shown but never ranked.",
                 xy=(slots[-1][1] + 0.62, 0.585), ha="right", va="top",
                 fontsize=8.0, color=MUTED)
    axC.legend(handles=[Patch(facecolor=MUTED, edgecolor=INK, label="human PBMC"),
                        Patch(facecolor=MUTED, edgecolor=INK, alpha=0.62,
                              label="mouse testis (m1, m2)"),
                        Patch(facecolor="white", edgecolor=INK, hatch="////",
                              label="caveated arm")],
               loc="upper left", frameon=False, fontsize=8.6)

    # ------------------------------------------------------------------ (d)
    dorder = ["scAPAtrap", "Sierra", "PeakATail", "polyApipe"]
    strata = [("conc_all", "all calls", 1.00), ("conc_depth2plus", "depth $\\geq$ 2", 0.72),
              ("conc_depth1", "depth $\\leq$ 1", 0.42)]
    bw = 0.25
    ticklabels = []
    for i, tool in enumerate(dorder):
        g = rep[rep.tool == tool]
        for j, (col, lbl, alpha) in enumerate(strata):
            vs = g[col].to_numpy(float)
            xpos = i + (j - 1) * bw
            ncol = {"conc_all": "n_query", "conc_depth2plus": "n_depth2plus",
                    "conc_depth1": "n_depth1"}[col]
            if np.all(~np.isfinite(vs)):
                axD.text(xpos, 0.02, "no calls\nat this\ndepth", ha="center", va="bottom",
                         fontsize=7.4, color=MUTED)
            else:
                v = float(np.nanmean(vs))
                axD.bar(xpos, v, bw * 0.9, color=COLOR[tool], alpha=alpha, zorder=3,
                        edgecolor=INK, linewidth=0.6)
                axD.plot([xpos, xpos], [np.nanmin(vs), np.nanmax(vs)], color=INK, lw=1.2, zorder=5)
                axD.text(xpos, np.nanmax(vs) + 0.017, f"{v:.2f}", ha="center",
                         fontsize=8.0, color=INK)
            for _, r in g.iterrows():
                plotted.append(dict(panel="d_reproducibility", dataset=MOUSE, tool=tool,
                                    arm=f"{tool}_{r.direction}", x=xpos, y=r[col],
                                    series=f"concordance@100, {lbl}, n={int(r[ncol])}"))
        f1frac = float(g.frac_depth1.mean())
        ticklabels.append(f"{tool}\n{f1frac*100:.0f}% depth $\\leq$ 1")
        plotted.append(dict(panel="d_reproducibility", dataset=MOUSE, tool=tool, arm=tool,
                            x=float(i), y=f1frac, series="fraction of calls with depth<=1"))
    nullmax = float(rep.conc_null.max())
    axD.axhspan(0, nullmax, color=MUTED, alpha=0.4, lw=0, zorder=1)
    axD.text(4.30, 0.30,
             f"Chance level $\\leq$ {nullmax:.3f}\n(shaded strip at the axis floor):\n"
             f"the same query set scored\nagainst a genic-shuffled\ncopy of the other replicate.",
             fontsize=8.0, color=MUTED, ha="left", va="top")
    plotted.append(dict(panel="d_reproducibility", dataset=MOUSE, tool="(null)",
                        arm="genic_shuffle_null", x=np.nan, y=nullmax,
                        series="max chance concordance@100 over tools/directions"))
    axD.set_xticks(range(len(dorder)))
    axD.set_xticklabels(ticklabels, fontsize=8.8, color=INK)
    axD.set_xlim(-0.55, 6.05)
    axD.set_ylim(0, 1.02)
    axD.set_ylabel("fraction of calls with a strand-matched\ncounterpart $\\leq$ 100 bp in the other mouse",
                   fontsize=9.5)
    axD.set_title("(d) replicate reproducibility, mouse testis — and the call depth it rests on",
                  fontsize=11.5, color=INK, loc="left", pad=8)
    axD.legend(handles=[Patch(facecolor=MUTED, edgecolor=INK, alpha=a, label=l)
                        for _, l, a in strata] +
                       [Line2D([], [], color=INK, lw=1.2, label="m1$\\to$m2 / m2$\\to$m1 range")],
               loc="upper left", bbox_to_anchor=(0.735, 0.995), frameon=False,
               fontsize=8.4, ncol=1)
    axD.text(4.30, 0.72,
             "Not shown, and why\n"
             "SCAPTURE — BED score column is 0 for\n"
             "every evaluated peak (162,346 human /\n"
             "81,138 mouse1, verified), so no depth;\n"
             "and its mouse2 run was truncated by\n"
             "disk exhaustion, so no replicate pair.\n"
             "scUTRquant — cross-replicate agreement\n"
             "of a fixed catalog is tautological.",
             fontsize=7.7, color=MUTED, ha="left", va="top")

    # ------------------------------------------------------------------ (e)
    eorder = ["PeakATail", "polyApipe", "scAPAtrap", "SCAPTURE", "scUTRquant", "Sierra"]
    for ax, col, lab, title in ((axE1, "runtime_s", "clean-run wall time (h)", "(e1) wall time"),
                                (axE2, "rss_gib", "peak RSS (GiB)", "(e2) peak memory")):
        for i, tool in enumerate(eorder):
            for ds, off, alpha in ((HUMAN, -0.19, 1.0), (MOUSE, 0.19, 0.62)):
                arms = arms_of(ds, tool)
                vals = [cons[cons.arm == a][col].iloc[0] for a in arms]
                if col == "runtime_s":
                    vals = [v / 3600 for v in vals]
                if not vals or not np.isfinite(np.nanmean(vals)):
                    continue
                v = float(np.nanmean(vals))
                # panel (e) hatching is RUNTIME-only by design; accuracy caveats
                # (HATCHED_DATA) are hatched in (b)/(c), never here.
                hatched = any(a in HATCHED_RUN for a in arms)
                ax.bar(i + off, v, 0.34, color=COLOR[tool], alpha=alpha, zorder=3,
                       edgecolor=INK, linewidth=0.6, hatch="////" if hatched else None)
                top = max(vals)
                if len(vals) > 1:
                    ax.plot([i + off, i + off], [min(vals), top], color=INK, lw=1.2, zorder=5)
                ax.text(i + off, top, "  " + (hhmm(v * 3600) if col == "runtime_s" else f"{v:.0f}"),
                        ha="center", va="bottom", fontsize=7.6, color=INK, rotation=90)
                for a, raw in zip(arms, vals):
                    plotted.append(dict(panel="e_resources", dataset=ds, tool=tool, arm=a,
                                        x=i + off, y=raw, series=col))
        ax.set_xticks(range(len(eorder)))
        ax.set_xticklabels(eorder, rotation=40, ha="right", fontsize=8.8, color=INK)
        ax.set_ylabel(lab, fontsize=9.5)
        ax.set_title(title, fontsize=11.5, color=INK, loc="left", pad=8)
    axE1.set_ylim(0, 16.5)
    axE2.set_ylim(0, 128)
    axE1.legend(handles=[Patch(facecolor=MUTED, edgecolor=INK, label="human"),
                         Patch(facecolor=MUTED, edgecolor=INK, alpha=0.62,
                               label="mouse (bar = mean of m1,m2)"),
                         Patch(facecolor="white", edgecolor=INK, hatch="////",
                               label="retried / resumed run (runtime only)")],
                loc="upper right", frameon=False, fontsize=7.6)
    axE2.annotate("PeakATail's 101 GiB on the human\nBAM is the largest single resource\ngap anywhere in this figure",
                  xy=(0.97, 0.90), xycoords="axes fraction", fontsize=7.8, color=MUTED,
                  ha="right", va="top")

    # ------------------------------------------------------------------ titles
    fig.text(0.055, 0.982,
             "Six PAS callers, two datasets: a trade surface, not a ranking",
             fontsize=19, color=INK, ha="left", va="top", fontweight="bold")
    fig.text(0.055, 0.960,
             "PeakATail ranks LAST of five de novo tools on human PBMC (F1 0.134) and second on mouse testis (0.258–0.260)",
             fontsize=13.5, color=INK, ha="left", va="top")
    fig.text(0.055, 0.938,
             "SCAPTURE is the precision arm among the de novo tools on BOTH datasets (P@100 0.65 human / 0.69 mouse) on the fewest-to-middling calls; polyApipe takes the best\n"
             "de novo F1 on both (0.261 human, 0.308–0.312 mouse) but 62% of its mouse calls are depth-1 singletons that reproduce across biological replicates only 0.31 of the\n"
             "time, against 0.81 for its deeper calls; the fewest-calls arm is Sierra on mouse and SCAPTURE on human, and both buy precision that way; scAPAtrap calls the most\n"
             "and is least precise. Every rank order changes between the two datasets — no single accuracy number describes any tool in this field.",
             fontsize=10.6, color=MUTED, ha="left", va="top", linespacing=1.45)

    paras = [
        "SCORING — One path for every arm (scripts/benchmark_tools/score_tool.py): each tool's output reduced to 1-bp inferred cleavage points, strand-matched with "
        "`bedtools closest -s -d -t first` against curated PolyASite 2.0 representative sites. Recall and F1 use the per-dataset detected-gene-restricted atlas (human "
        "285,136 sites / 14,949 genes; mouse 126,686 sites / 14,014 genes) so denominators are identical across tools. Null = 3 seeds of width-preserving `bedtools shuffle` "
        "inside merged annotated gene bodies. Full harvest: benchmark_consolidated.tsv; every plotted value: benchmark_headtohead.tsv; written-up conclusions: "
        "manuscript/09_headtohead_results.md.",
        "CAVEATED ARMS — SCAPTURE mouse (hatched in b, c) is mouse1 only: the mouse2 run was truncated by disk exhaustion and SCAPTURE exited 0 regardless, so no mouse "
        "replicate pair exists and its mouse numbers are single-replicate. Hatching in panel (e) marks RUNTIME caveats only, never accuracy: PeakATail human is attempt 4 of 4 "
        "(attempts 1-3 ran 1:33:21, 1:52:42 and 4:52:58 before dying on three separate CellRanger-input bugs; those are disclosed, not summed into the bar); scAPAtrap human is "
        "a resumed run (4:04:29) after attempt 1 was OOM-killed at 8:54:21, and it skipped three already-finished stages, so its bar UNDERSTATES a from-scratch run (12:58:50 "
        "of machine time in total); scUTRquant human needed a 4:53 salvage rerun after a missing-nbformat failure, and scUTRquant mouse wall time is dominated by a "
        "STARsolo->CellRanger-tag BAM conversion this benchmark imposed (64 and 54 min of the 76 and 64 min totals), not by the pipeline itself.",
        "PANEL (d) — Depth is per-tool and NOT interconvertible: polyApipe = poly(A)-read peakdepth from its GFF, PeakATail = total UMIs (annotated_matrix.mtx row sums), "
        "Sierra = total UMIs (CountPeaks matrix row sums), scAPAtrap = total UMIs (counts.tsv.gz). scAPAtrap has no depth<=1 bar because its own "
        "reducePeaks(min.cells=10, min.count=10) step removes such peaks before it reports them, so its high concordance is measured on an already depth-filtered set and is "
        "not comparable to the others at face value. PeakATail's depth<=1 stratum includes 715/45,921 (mouse1) and 555/46,672 (mouse2) sites carrying zero UMIs in the "
        "filtered-cell matrix. The two directions differ only through their denominators; the vertical line is their range.",
        "READ WITH CARE — scUTRquant is annotation-based and its numbers are near-tautological; it is shown for context and never ranked with the de novo tools. n called "
        "spans 21k-787k (37-fold), so precision and recall must always be read together with panel (b). SCAPTURE is annotation-GUIDED (peaks called inside gene models, then "
        "filtered by the DeepPASS sequence classifier), the likeliest source of its precision lead. PeakATail's PBMC deficit is not a threshold artefact (its top 22k / 36k / "
        "106k calls by UMI give precision 0.119 / 0.120 / 0.118, identical to the full set) and not a coordinate offset (relaxing to 200 bp lifts it only 0.118 -> 0.166 while "
        "polyApipe goes 0.380 -> 0.410).",
    ]
    wrapped = "\n".join(textwrap.fill(p, 218) for p in paras)
    fig.text(0.055, 0.185, wrapped, fontsize=7.6, color=MUTED, ha="left", va="top",
             linespacing=1.55)

    fig.legend(handles=[Line2D([], [], color=COLOR[t], lw=3.4, label=t) for t in TOOLS],
               loc="upper left", bbox_to_anchor=(0.055, 0.898), frameon=False,
               fontsize=10, ncol=6, columnspacing=1.6, handlelength=1.6,
               title="colour = tool, assigned in alphabetical order so that no tool is visually privileged",
               title_fontsize=8.6, alignment="left")

    save_manuscript(fig, NAME, facecolor=SURFACE)
    plt.close(fig)

    pl = pd.DataFrame(plotted)
    pl.to_csv(OUTDIR / f"{NAME}.tsv", sep="\t", index=False, float_format="%.6f")
    print("wrote", OUTDIR / f"{NAME}.tsv", len(pl), "rows")


if __name__ == "__main__":
    main()

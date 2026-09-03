#!/usr/bin/env python3
"""
fig3_tradeoff.py -- manuscript Fig 3 (stem `fig3_tradeoff`, kept from the
2026-09-02 rename pass; FIGURE_MAP.tsv and FIGURES_MANIFEST.md key off it).

TITLE AFTER THE 2026-09-02 PUBLICATION DESIGN PASS:
    "Robustness: matching-window resolution, replicate reproducibility and compute"
The stem still says "tradeoff" for file-identity reasons only -- the trade curve
itself left this figure (see below).  The Legend carries the current title.

What the design pass changed (manuscript/figures/DESIGN_DIRECTIVES.md):
  * DIRECTIVE 5 -- the trade surface (molecule-support sweep) MOVED to Fig 2's
    precision/recall plane, because Fig 1e / Fig 3a / Fig 2a were three views of
    one plane and the paper may hold exactly one.  The sweep is still COMPUTED
    here (shared module `_molsweep.py`, one implementation for both scripts) and
    still written to `fig3_tradeoff.tsv` with `plotted = False`, so this figure's
    audit TSV keeps every number it had.  Fig 3 is now pure robustness.
  * DIRECTIVE 4 -- development history is off the main figures.  Panel d (compute)
    shows the CURRENT caller beside the competitors ONLY: no v1 point, no v1->v2
    arrow, no "shipped" caller.  Those rows stay in `fig3_tradeoff_compute.tsv`
    and the improvement story is Fig S8 / Fig S12, with one Legend pointer.
    Panel c likewise drops the shipped-caller row (kept in the reproducibility TSV).
  * CROSS-DONOR reproducibility joins cross-mouse in panel c: the human half of
    the replication claim used to be missing from the mains entirely ("PBMC is a
    single donor" was a caveat in the old legend), and the pre-registered second
    donor (26; Fig S4) supplies it.  Fig S4 keeps the full protocol -- gate,
    per-direction nulls, arithmetic ceiling, matched-N control; the main figure
    shows the headline bars only.
  * `_pubstyle` supplies the palette, the per-tool marker identity and the type
    scale; every prose label goes through `sentence_case()` (directive 6).

PANELS
  a  MATCHING-WINDOW SENSITIVITY, precision.  Atlas-agreement precision at
     10 / 25 / 50 / 100 bp for the v2 default, the v2 >=1-molecule sensitivity arm
     and the five competitors, PBMC.  Answers "is the precision lead an artefact
     of a loose window?".  Source: the verified score TSVs.
  b  The same for detected-gene recall.
  c  REPLICATE REPRODUCIBILITY.  Strand-matched site agreement (bedtools closest
     -s -d -t first, both directions) at 100 bp with the 25 bp value as a diamond:
     the two GSE104556 testis mice as biological replicates (all tools), and the
     two human PBMC donors (PeakATail arms only -- no competitor was run on donor 2).
  d  COMPUTE.  Wall time vs peak RSS of the current caller on the PBMC BAM,
     against the five competitors.

INPUTS (all primary or verified)
  results/benchmark_tools/pbmc_10k_v3/peakatail_clipseeded_final_v2*/{runtime_mem.txt,score_*.tsv}
  results/benchmark_tools/gse104556/peakatail_clipseeded_final_v2/mouse{1,2}/{score_*.tsv,.work_*/}
  results/benchmark_tools/{pbmc_10k_v3,gse104556}/<competitor>/score_*.tsv
  results/benchmark_tools/gse104556/depth_concordance.tsv        (verified, manuscript/05 + 09)
  results/benchmark_tools/pbmc4k_donor2/concordance.txt          (verified, manuscript/26, FIXED)
  results/figures/manuscript/benchmark_headtohead.tsv            (verified, manuscript/05)
  manuscript/19_final_gate_v2.md, manuscript/18_spermatogenesis_final.md   (quoted facts)

OUTPUTS
  manuscript/figures/fig3_tradeoff.{png,pdf}   (PNG 600 dpi, PDF fonttype 42)
  manuscript/figures/fig3_tradeoff.caption.md  ('## Legend' = the journal legend,
                                                single source of the caption;
                                                '## Provenance')
  results/figures/manuscript/fig3_tradeoff.tsv               (window sweep + the
                                                              molecule sweep, the
                                                              latter not plotted)
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
import sys
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _pubstyle import PAL, TOOL_STYLE, TYPE, apply_rc, sentence_case   # noqa: E402
from _molsweep import molecule_sweep                                   # noqa: E402

apply_rc()

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
INK, MUTED, GRID = PAL["ink"], PAL["muted"], PAL["grid"]
PA = TOOL_STYLE["PeakATail"]["color"]
G = BT / "gse104556"
D2 = BT / "pbmc4k_donor2"
CUTOFFS = [10, 25, 50, 100]
GATE_P, ORIG_P = 0.50, 0.38     # NOT PLOTTED here since 2026-09-02 (Fig 2 owns the gate line)

# canonical arm labels used by every table in this script
DEFAULT_ARM = "PeakATail v2 default"
SENS_ARM = "PeakATail v2 tier-1 ≥1 mol"


def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, executable="/bin/bash",
                          env=dict(os.environ, LC_ALL="C"), capture_output=True, text=True).stdout


# ---------------------------------------------------------------------------
# molecule-support sweep -- NOT PLOTTED since 2026-09-02.
# DIRECTIVE 5 moved this curve into Fig 2's precision/recall plane (the paper
# holds exactly one such plane).  It is still computed and still written to
# fig3_tradeoff.tsv so this figure's audit trail keeps every value it had, and
# so the >=1 / >=2 reproduction assert against manuscript/19 §1 still runs here.
# ---------------------------------------------------------------------------
sweep = molecule_sweep()
sweep.insert(0, "panel", "not_plotted_moved_to_fig2")
sweep.insert(1, "plotted", False)
sweep.insert(2, "tool", "PeakATail v2 (IP arm)")

# ---------------------------------------------------------------------------
# panels a/b -- matching-window sensitivity from the verified score TSVs
# ---------------------------------------------------------------------------
CUTOFF_ARMS = [
    # label, dataset, score TSV
    (DEFAULT_ARM, "pbmc", BT / "pbmc_10k_v3/peakatail_clipseeded_final_v2_ipfilt"
     / "score_pbmc_final_v2_ipfilt__PRESPEC_precision_default.tsv"),
    (SENS_ARM, "pbmc", BT / "pbmc_10k_v3/peakatail_clipseeded_final_v2_ipfilt"
     / "score_pbmc_final_v2_ipfilt__tier1.tsv"),
    ("scUTRquant*", "pbmc", BT / "pbmc_10k_v3/scutrquant/score_scutrquant.tsv"),
    ("SCAPTURE", "pbmc", BT / "pbmc_10k_v3/scapture/score_scapture.tsv"),
    ("polyApipe", "pbmc", BT / "pbmc_10k_v3/polyapipe/score_polyapipe.tsv"),
    ("Sierra", "pbmc", BT / "pbmc_10k_v3/sierra/score_sierra.tsv"),
    ("scAPAtrap", "pbmc", BT / "pbmc_10k_v3/scapatrap/score_scapatrap.tsv"),
    # mouse rows: audit TSV only (not plotted), so the window check is auditable on both datasets
    (DEFAULT_ARM, "mouse1", G / "peakatail_clipseeded_final_v2/mouse1"
     / "score_testis_m1_final_v2__PRESPEC_precision_default.tsv"),
    (DEFAULT_ARM, "mouse2", G / "peakatail_clipseeded_final_v2/mouse2"
     / "score_testis_m2_final_v2__PRESPEC_precision_default.tsv"),
    (SENS_ARM, "mouse1", G / "peakatail_clipseeded_final_v2/mouse1"
     / "score_testis_m1_final_v2__tier1.tsv"),
    (SENS_ARM, "mouse2", G / "peakatail_clipseeded_final_v2/mouse2"
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
        crows.append(dict(panel="a,b", dataset=ds, tool=tool, cutoff_bp=c,
                          n=int(pr.n_query), n_matched=int(pr.n_matched),
                          atlas_agreement_precision=float(pr.value),
                          recall_detected_genes=float(rd.value),
                          recall_detected_denominator=int(rd.n_query),
                          plotted=(ds == "pbmc" and c in CUTOFFS),
                          source=str(tsv.relative_to(WD))))
cut = pd.DataFrame(crows)

# ---------------------------------------------------------------------------
# panel c -- replicate agreement: cross-mouse (bedtools, computed here) and
# cross-donor (the pre-registered second-donor protocol output, manuscript/26)
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
    for arm, lab in (("PRESPEC_precision_default", DEFAULT_ARM), ("tier1", SENS_ARM)):
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
rep_v2["comparison"] = "testis mouse 1 vs mouse 2"

# v1-era per-tool reproducibility (verified: manuscript/05 adversarial pass, manuscript/09 table).
# The "PeakATail shipped (v1)" row of this table is NOT PLOTTED since 2026-09-02 (directive 4:
# development history off the main figures); it stays in fig3_tradeoff_reproducibility.tsv and
# its range is quoted in the Legend.  The competitor rows ARE plotted.
dc = pd.read_csv(G / "depth_concordance.tsv", sep="\t")
DC_LAB = {"PeakATail": "PeakATail shipped (v1)", "polyApipe": "polyApipe (v1)",
          "Sierra": "Sierra (v1)", "scAPAtrap": "scAPAtrap (v1)"}
rep_v1 = pd.DataFrame([dict(tool=DC_LAB[r.tool], direction=r.direction, window_bp=100,
                            n_query=int(r.n_query), n_agree=int(round(r.conc_all * r.n_query)),
                            concordance=float(r.conc_all), chance_mean=float(r.conc_null),
                            jaccard_strand_aware=np.nan, jaccard_intersect_bp=np.nan,
                            jaccard_union_bp=np.nan,
                            comparison="testis mouse 1 vs mouse 2",
                            provenance="results/benchmark_tools/gse104556/depth_concordance.tsv "
                                       "(v1-era head-to-head run, verified in manuscript/05)")
                       for r in dc.itertuples() if r.tool in DC_LAB])


# ---- cross-donor: the pre-registered second-donor protocol (26 §5, FIXED) ---
def donor_table():
    """Parse results/benchmark_tools/pbmc4k_donor2/concordance.txt into repro rows."""
    txt = (D2 / "concordance.txt").read_text().splitlines()
    arm, nq = None, {}
    rows = []
    ARMLAB = {"PRESPEC_precision_default": DEFAULT_ARM, "tier1": SENS_ARM}
    DIRLAB = {"pbmc4k->pbmc10k": "donor2->donor1", "pbmc10k->pbmc4k": "donor1->donor2"}
    real, nulls = {}, {}
    for line in txt:
        if line.startswith("############ ARM="):
            arm = ARMLAB[line.split("=", 1)[1].strip()]
        elif line.startswith("# pbmc_10k_v3"):
            f = line.replace("#", "").split()
            nq[arm] = {"donor1": int(f[1].split("=")[1]), "donor2": int(f[3].split("=")[1])}
        elif line and not line.startswith(("direction", "#")):
            f = line.split("\t")
            # row ids are "<a>-><b>_REAL" or "<a>-><b>_null_s<seed>"
            base, _, kind = f[0].partition("_")
            key = (arm, DIRLAB[base])
            if kind == "REAL":
                real[key] = dict(n_query=int(f[1]), n25=int(f[2]), n100=int(f[3]),
                                 f25=float(f[4]), f100=float(f[5]))
            else:
                nulls.setdefault(key, {"25": [], "100": []})
                nulls[key]["25"].append(float(f[4])); nulls[key]["100"].append(float(f[5]))
    for (a, d), r in real.items():
        # the arithmetic ceiling of a direction is the OTHER set's size / this set's size
        other = nq[a]["donor2"] if d.startswith("donor1") else nq[a]["donor1"]
        for w in (25, 100):
            rows.append(dict(tool=a, direction=d, window_bp=w, n_query=r["n_query"],
                             n_agree=r[f"n{w}"], concordance=r[f"f{w}"],
                             chance_mean=float(np.mean(nulls[(a, d)][str(w)])),
                             jaccard_strand_aware=np.nan, jaccard_intersect_bp=np.nan,
                             jaccard_union_bp=np.nan,
                             comparison="PBMC 10k v3 (donor 1) vs pbmc4k (donor 2)",
                             arithmetic_ceiling=min(1.0, other / r["n_query"]),
                             provenance="results/benchmark_tools/pbmc4k_donor2/concordance.txt "
                                        "(pre-registered protocol of manuscript/26 §5, verifier FIXED)"))
    return pd.DataFrame(rows)


rep_donor = donor_table()
repro = pd.concat([rep_v2, rep_v1, rep_donor], ignore_index=True)
repro.insert(0, "panel", "c")

# EXPECT asserts on the plotted cross-donor values (manuscript/26 §5, verifier FIXED)
def _dn(direction, w, col="concordance"):
    q = repro[(repro.tool == DEFAULT_ARM) & (repro.direction == direction) & (repro.window_bp == w)]
    assert len(q) == 1, (direction, w)
    return round(float(q[col].iloc[0]), 4)


assert (_dn("donor2->donor1", 100), _dn("donor2->donor1", 25)) == (0.8423, 0.8132)
assert (_dn("donor1->donor2", 100), _dn("donor1->donor2", 25)) == (0.3993, 0.3613)
assert _dn("donor1->donor2", 100, "arithmetic_ceiling") == 0.4443    # 20,672 / 46,524 (26 §5)
# SIDECAR-GENERATION CHECK: the chance level is a Legend sentence and a dotted
# reference line, not a labelled series -- its bound is asserted here.
assert float(repro.chance_mean.max()) <= 0.0095, float(repro.chance_mean.max())

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
# NOT PLOTTED since 2026-09-02 (directive 4: no v1 point, no v1->v2 arrow -- that is
# Fig S8 / Fig S12's story).  The row stays so the improvement number keeps a source.
crows.append(dict(tool="PeakATail v1", label="PeakATail v1", wall_h=w_v1, wall_str=el_v1,
                  peak_rss_gb=g_v1, kind="not_plotted_history",
                  note="final Stage-2 clip-seeded arm WITHOUT --ip-filter, code 4efeb125, --threads 16; "
                       "paired with the v2 no-IP arm below so the Stage-1d comparison is like for like. "
                       "NOT PLOTTED since 2026-09-02 -- see Fig S8 / Fig S12",
                  source="results/benchmark_tools/pbmc_10k_v3/peakatail_clipseeded_final/runtime_mem.txt"))
crows.append(dict(tool="PeakATail", label="PeakATail", wall_h=w_v2, wall_str=el_v2,
                  peak_rss_gb=g_v2, kind="peakatail",
                  note="code 9dfdefb, clip-seeded arm WITHOUT --ip-filter, --threads 16, measured with all "
                       "four arms running concurrently; this is the 34:37 / 12.53 GB pair quoted in 19 sec.3. "
                       "The IP-filtered default arm used in panels a-c was cheaper on both codes (audit rows below)",
                  source="results/benchmark_tools/pbmc_10k_v3/peakatail_clipseeded_final_v2/runtime_mem.txt"))
crows.append(dict(tool="PeakATail v2 uncontended", label="", wall_h=27 / 60 + 43 / 3600, wall_str="27:43",
                  peak_rss_gb=g_v2, kind="not_plotted_caveat",
                  note="single uncontended run quoted in 19 sec.3 (peak RSS not re-measured; v2 value reused). "
                       "NOT PLOTTED since 2026-09-02 -- the concurrency caveat is a Legend sentence",
                  source="manuscript/19_final_gate_v2.md section 3"))
crows.append(dict(tool="PeakATail shipped", label="PeakATail shipped (coverage-only, v1)",
                  wall_h=res["PeakATail"]["runtime_s"], wall_str="3:27:05",
                  peak_rss_gb=res["PeakATail"]["rss_gib"] * 1048576 / 1e6, kind="not_plotted_history",
                  note="pre-clip-seeding caller, v1-era benchmark run. NOT PLOTTED since 2026-09-02",
                  source="results/figures/manuscript/benchmark_headtohead.tsv (verified, manuscript/05)"))
# the IP-filtered arms -- the call sets panels a-c actually use.  NOT plotted (panel d shows the
# no-IP arm, the pair 19 sec.3 / 15 sec.5 quote), but recorded so the arm identity of every
# compute number is auditable.
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
comp.insert(1, "plotted", comp.kind.isin(["peakatail", "competitor"]))

# ===========================================================================
# FIGURE
# ===========================================================================
FIG_W, FIG_H = 7.09, 6.60
fig = plt.figure(figsize=(FIG_W, FIG_H))
gs_top = fig.add_gridspec(1, 2, left=0.086, right=0.988, top=0.9550, bottom=0.6600, wspace=0.235)
gs_bot = fig.add_gridspec(1, 2, left=0.178, right=0.988, top=0.4800, bottom=0.0980,
                          wspace=0.360, width_ratios=[1.0, 0.88])
axA = fig.add_subplot(gs_top[0, 0])
axB = fig.add_subplot(gs_top[0, 1])
axC = fig.add_subplot(gs_bot[0, 0])
axD = fig.add_subplot(gs_bot[0, 1])
callouts = {"a": 0, "b": 0, "c": 0, "d": 0}    # design law: <= 3 short callouts per panel
OBSTACLES = {}


def style(ax, grid_axis="both"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED, length=2.4, width=0.7, labelsize=TYPE["tick"])
    ax.grid(True, axis=grid_axis, color=GRID, lw=0.5, alpha=0.85)
    ax.set_axisbelow(True)


def panel_tag(ax, letter, title):
    ax.text(-0.005, 1.035, letter, transform=ax.transAxes, fontsize=TYPE["panel_letter"],
            fontweight="bold", va="bottom", ha="left", color=INK)
    ax.text(0.062, 1.038, sentence_case(title), transform=ax.transAxes,
            fontsize=TYPE["panel_title"], va="bottom", ha="left", color=INK)


def occupy(key, x, y, pad=4.0):
    OBSTACLES.setdefault(key, []).append((float(x), float(y), pad))


SLOTS = [(-7, -7, "right", "top"), (8, 8, "left", "bottom"), (8, -7, "left", "top"),
         (-7, 8, "right", "bottom"), (0, -10, "center", "top"), (0, 10, "center", "bottom"),
         (11, 0, "left", "center"), (-11, 0, "right", "center")]


def place(ax, key, s, xy, color=INK, weight="normal", slots=SLOTS, literal=False, count=True):
    """Auto-place a label into the first slot that touches no mark and no other text.

    Bounding-box discipline (DESIGN_DIRECTIVES.md item 1).  `literal` labels open
    with a relational operator, so directive 6's initial capital does not apply.
    """
    if count:
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
    txt = s if literal else sentence_case(s)
    ax_bb = ax.get_window_extent(renderer=rend)
    for dx, dy, ha, va in slots:
        a = ax.annotate(txt, xy=xy, xytext=(dx, dy), textcoords="offset points",
                        fontsize=TYPE["annotation"], color=color, ha=ha, va=va,
                        fontweight=weight, zorder=9)
        bb = a.get_window_extent(renderer=rend).expanded(1.08, 1.16)
        if any(bb.overlaps(b) for b in boxes) or not ax_bb.fully_contains(bb.x0, bb.y0) \
                or not ax_bb.fully_contains(bb.x1, bb.y1):
            a.remove()
            continue
        return a
    raise AssertionError(f"panel {key}: no free slot for label {s!r}")


# ---- panels a / b: matching-window sensitivity ----------------------------
plotted = cut[(cut.dataset == "pbmc") & (cut.cutoff_bp.isin(CUTOFFS))]
# draw order: the field first, our two arms last (they own the message)
ORDER = ["scUTRquant*", "SCAPTURE", "polyApipe", "Sierra", "scAPAtrap", SENS_ARM, DEFAULT_ARM]
TOOL_OF = {"scUTRquant*": "scUTRquant", "SCAPTURE": "SCAPTURE", "polyApipe": "polyApipe",
           "Sierra": "Sierra", "scAPAtrap": "scAPAtrap", SENS_ARM: "PeakATail",
           DEFAULT_ARM: "PeakATail"}
for ax, col, key, ylab in ((axA, "atlas_agreement_precision", "a", "atlas-agreement precision"),
                           (axB, "recall_detected_genes", "b", "detected-gene recall R_det")):
    # plain "R_det", not mathtext: a mathtext subscript renders at 0.7x the label
    # size (7.5 -> 5.25 pt), under the 6 pt floor the type guard cannot measure
    ax.set_label(key)
    for tool in ORDER:
        s = plotted[plotted.tool == tool].sort_values("cutoff_bp")
        st = TOOL_STYLE[TOOL_OF[tool]]
        ours = tool in (DEFAULT_ARM, SENS_ARM)
        dflt = tool == DEFAULT_ARM
        ax.plot(s.cutoff_bp.to_numpy(), s[col].to_numpy(),
                ls=(0, (3, 1.6)) if tool == "scUTRquant*" else "-", color=st["color"],
                lw=1.7 if dflt else (1.2 if ours else 1.0),
                marker="D" if ours else st["marker"],
                markersize=3.6 if ours else 2.8,
                markerfacecolor=(st["color"] if dflt else
                                 ("white" if ours else st.get("mfc", st["color"]))),
                markeredgecolor=st["color"] if ours else "white",
                markeredgewidth=1.0 if ours else 0.5,
                zorder=7 if dflt else (6 if ours else 3))
        for x, y in zip(s.cutoff_bp, s[col]):
            occupy(key, x, y, pad=3.4)
    ax.set_xscale("log")
    ax.set_xticks(CUTOFFS)
    ax.set_xticklabels([str(c) for c in CUTOFFS])
    ax.minorticks_off()
    ax.set_xlim(8.6, 168)
    ax.set_xlabel(sentence_case("matching window (bp, strand-matched)"))
    ax.set_ylabel(sentence_case(ylab))
    panel_tag(ax, key, "PBMC 10k v3")
    style(ax)
axA.set_ylim(0.0, 0.92)
axA.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
axB.set_ylim(0.0, 0.32)
axB.set_yticks([0.0, 0.1, 0.2, 0.3])
# direct labels on our two arms only (<= 4 direct labels, per the mark spec);
# the shared key below carries every other identity
for ax, col, key in ((axA, "atlas_agreement_precision", "a"),
                     (axB, "recall_detected_genes", "b")):
    for tool, lab, weight in ((DEFAULT_ARM, "≥2 mol", "bold"), (SENS_ARM, "≥1 mol", "normal")):
        r = plotted[(plotted.tool == tool) & (plotted.cutoff_bp == 100)].iloc[0]
        place(ax, key, lab, (100.0, float(r[col])), color=PA, weight=weight, literal=True,
              slots=[(9, 0, "left", "center"), (9, 7, "left", "bottom"), (9, -7, "left", "top")])

# ---- panel c: replicate agreement -----------------------------------------
axC.set_label("c")
MOUSE_CMP = "testis mouse 1 vs mouse 2"
DONOR_CMP = "PBMC 10k v3 (donor 1) vs pbmc4k (donor 2)"
MOUSE_ROWS = [("scAPAtrap (v1)", "scAPAtrap", "scAPAtrap"),
              ("Sierra (v1)", "Sierra", "Sierra"),
              (DEFAULT_ARM, "PeakATail", "PeakATail ≥2 mol"),
              (SENS_ARM, "PeakATail", "PeakATail ≥1 mol"),
              ("polyApipe (v1)", "polyApipe", "polyApipe")]
# One donor row: the shipping default.  The >=1-molecule arm's donor values are in
# the Legend and the audit TSV -- a second row would halve the bar pitch and force
# the tip labels to collide, and Fig S4 is the home of the full donor protocol.
DONOR_ROWS = [(DEFAULT_ARM, "PeakATail", "PeakATail ≥2 mol")]
# drawn bottom-up, so the donor block is listed first and reads UNDER the mice
GROUPS = [(DONOR_CMP, DONOR_ROWS, ("donor2->donor1", "donor1->donor2"), "PBMC donors"),
          (MOUSE_CMP, MOUSE_ROWS, ("mouse1->mouse2", "mouse2->mouse1"), "Testis mice")]
# bar pitch inside a row: wide enough that the two tip labels clear each other at
# 6 pt (the bounding-box check below fails the build if it is not)
bh = 0.40
y = 0.0
yticks, ylabels, group_marks = [], [], []
row_y = {}
# a direction whose query set is larger than the other replicate's cannot reach 1.0:
# its call-count ceiling is drawn as a tick, and the tip label clears it
CEIL = {(r.tool, r.direction): float(r.arithmetic_ceiling)
        for r in repro[(repro.comparison == DONOR_CMP) & (repro.window_bp == 100)].itertuples()
        if float(r.arithmetic_ceiling) < 0.999}
for gi, (cmpname, rowspec, dirs, glab) in enumerate(GROUPS):
    ybot = y
    for tool, style_key, lab in reversed(rowspec):
        st = TOOL_STYLE[style_key]
        s = repro[(repro.tool == tool) & (repro.comparison == cmpname) & (repro.window_bp == 100)]
        for j, dr in enumerate(dirs):
            v = float(s[s.direction == dr].concordance.iloc[0])
            yb = y + (0.5 - j) * bh
            axC.barh(yb, v, height=bh * 0.88, color=st["color"],
                     alpha=1.0 if j == 0 else 0.48, edgecolor="none", zorder=3)
            ceil = CEIL.get((tool, dr))
            if ceil is not None:
                axC.plot([ceil, ceil], [yb - 0.5 * bh, yb + 0.5 * bh], color=MUTED, lw=1.1,
                         zorder=7)
            axC.text(max(v, ceil or 0.0) + 0.016, yb, f"{v:.3f}", va="center", ha="left",
                     fontsize=TYPE["annotation_min"], color=INK, zorder=4)
        q = repro[(repro.tool == tool) & (repro.comparison == cmpname) & (repro.window_bp == 25)]
        for j, dr in enumerate(dirs):
            if len(q):
                axC.plot([float(q[q.direction == dr].concordance.iloc[0])], [y + (0.5 - j) * bh],
                         marker="D", markersize=3.0, color=INK, markeredgecolor="white",
                         markeredgewidth=0.6, zorder=6)
        yticks.append(y); ylabels.append(lab); row_y[(cmpname, tool)] = y
        y += 1.0
    group_marks.append((ybot - 0.50, y - 0.50, glab))
    y += 0.70 if gi == 0 else 0.0
# (the reverse donor direction is capped by call-count arithmetic, not by disagreement;
#  its ceiling tick is drawn with the bar above and keyed below -- no inline text, so
#  nothing can collide with it)
_ceil = CEIL[(DEFAULT_ARM, "donor1->donor2")]
axC.set_yticks(yticks)
axC.set_yticklabels([l if l.startswith("≥") else sentence_case(l) for l in ylabels],
                    fontsize=TYPE["tick"])
axC.set_xlim(0, 1.13)
axC.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
KEY_BAND = 1.45          # empty band under the bars that holds this panel's key
axC.set_ylim(-0.50 - KEY_BAND, y - 0.50)
axC.set_xlabel(sentence_case("replicate agreement (≤100 bp, strand-matched)"),
               fontsize=TYPE["tick"])
axC.plot([0.0095, 0.0095], [-0.46, y - 0.72], color=MUTED, lw=0.8, ls=(0, (1, 2)), zorder=2)
for y0, y1, glab in group_marks:
    axC.text(-0.395, (y0 + y1) / 2, sentence_case(glab), transform=axC.get_yaxis_transform(),
             rotation=90, ha="center", va="center", fontsize=TYPE["annotation"], color=MUTED)
    axC.plot([-0.345, -0.345], [y0, y1], transform=axC.get_yaxis_transform(),
             color=GRID, lw=1.1, clip_on=False, zorder=1)
axC.legend(handles=[Patch(facecolor=PAL["neutral"], edgecolor="none", label="First → second"),
                    Patch(facecolor=PAL["neutral"], alpha=0.48, edgecolor="none", label="Second → first"),
                    Line2D([], [], marker="D", ls="none", color=INK, markersize=3.0,
                           label="Same pair, ≤25 bp"),
                    Line2D([], [], ls=(0, (1, 2)), color=MUTED, lw=0.8, label="Chance ≤0.010"),
                    Line2D([], [], ls="-", color=MUTED, lw=1.1, label="Call-count ceiling")],
           loc="lower left", bbox_to_anchor=(-0.012, -0.012), frameon=False, ncol=2,
           handlelength=1.4, handletextpad=0.45, labelspacing=0.35, columnspacing=0.8,
           borderpad=0.1, fontsize=TYPE["annotation_min"])
panel_tag(axC, "c", "biological replicates")
style(axC, grid_axis="x")

# ---- panel d: compute, current caller vs competitors ----------------------
axD.set_label("d")
for r in comp[comp.plotted].itertuples():
    tool = "PeakATail" if r.kind == "peakatail" else ("scUTRquant" if r.label == "scUTRquant*" else r.label)
    st = TOOL_STYLE[tool]
    ours = r.kind == "peakatail"
    axD.plot(r.wall_h, r.peak_rss_gb, marker=st["marker"], ms=7.5 if ours else 6.0,
             mfc="none", mec="white", mew=2.4, ls="none", zorder=5)
    axD.plot(r.wall_h, r.peak_rss_gb, marker=st["marker"], ms=7.5 if ours else 6.0,
             mfc=st.get("mfc", st["color"]), mec=st["color"], mew=1.3, ls="none", zorder=6)
    occupy("d", np.log10(r.wall_h), np.log10(r.peak_rss_gb), pad=6.0)
axD.set_xscale("log"); axD.set_yscale("log")
axD.set_xlim(0.33, 20)
axD.set_ylim(3.3, 46)
axD.set_xticks([0.5, 1, 2, 4, 8, 16])
axD.set_xticklabels(["0.5", "1", "2", "4", "8", "16"])
axD.set_yticks([4, 10, 30])
axD.set_yticklabels(["4", "10", "30"])
axD.minorticks_off()
axD.set_xlabel(sentence_case("wall time (h, log)"))
axD.set_ylabel(sentence_case("peak RSS (GB, log)"))
panel_tag(axD, "d", "compute, PBMC 10k v3")
style(axD)
# a log-log scatter of six named points reads best with direct labels and no key
for r in comp[comp.plotted].itertuples():
    tool = "PeakATail" if r.kind == "peakatail" else ("scUTRquant" if r.label == "scUTRquant*" else r.label)
    st = TOOL_STYLE[tool]
    lab = "scUTRquant (catalog)" if tool == "scUTRquant" else tool
    place(axD, "d", lab, (r.wall_h, r.peak_rss_gb), color=st["color"],
          weight="bold" if r.kind == "peakatail" else "normal", count=False)

# ---- one shared key for panels a/b ----------------------------------------
row1 = [
    Line2D([], [], color=PA, lw=1.7, marker="D", ms=4.4, mfc=PA, mec=PA, mew=1.0,
           label="PeakATail ≥2 molecules — precision default"),
    Line2D([], [], color=PA, lw=1.2, marker="D", ms=4.0, mfc="white", mec=PA, mew=1.0,
           label="PeakATail ≥1 molecule — sensitivity arm"),
]
row2 = []
for tool in ["polyApipe", "scAPAtrap", "SCAPTURE", "Sierra", "scUTRquant"]:
    st = TOOL_STYLE[tool]
    row2.append(Line2D([], [], color=st["color"], lw=1.0,
                       ls=(0, (3, 1.6)) if tool == "scUTRquant" else "-",
                       marker=st["marker"], ms=3.6, mfc=st.get("mfc", st["color"]),
                       mec="white", mew=0.5,
                       label="scUTRquant (catalog)" if tool == "scUTRquant" else tool))
lg1 = fig.legend(handles=row1, loc="center", bbox_to_anchor=(0.537, 0.590), ncol=2, frameon=False,
                 handletextpad=0.5, columnspacing=2.0, borderpad=0.0, handlelength=1.8,
                 fontsize=TYPE["annotation"])
# NB: fig.legend() already registers the legend on the figure; the fig.add_artist(lg1)
# that used to stand here registered the SAME artist a second time, so row 1 was drawn
# twice -- invisible in the raster (identical pixels) but two stacked text/marker objects
# in the vector PDF (design-judge finding 2, 2026-09-03).
fig.legend(handles=row2, loc="center", bbox_to_anchor=(0.537, 0.5585), ncol=5, frameon=False,
           handletextpad=0.5, columnspacing=1.7, borderpad=0.0, handlelength=1.8,
           fontsize=TYPE["annotation"])

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

# SIDECAR-GENERATION CHECK (2026-09-03 no-loss audit): panel c's Legend quotes Fig 5's
# v2 spermatogenesis replication numbers.  Pin them to Fig 5's own audit TSV whenever it
# has been rendered, so the two main figures can never drift apart again.
_f5 = OUTDIR / "fig5_spermatogenesis_replication.tsv"
if _f5.exists():
    _r5 = pd.read_csv(_f5, sep="\t").set_index("pair")
    _got = tuple(int(_r5.loc[k, "replicated_same_dir"]) for k in ("RS_vs_SPC", "ES_vs_RS", "ES_vs_SPC"))
    assert _got == (11219, 6070, 9480), ("Fig 5's v2 replication counts moved; Fig 3's Legend "
                                         "quotes them and must be updated together", _got)
    assert int(_r5.n_null_pairings.sum()) == 15 and int(_r5.null_replicated_same_dir_total.sum()) == 0
    print("cross-figure check: Fig 3's Legend matches fig5_spermatogenesis_replication.tsv "
          "(11,219 / 6,070 / 9,480; 0 in all 15 null pairings)")


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
fig.savefig(OUTDIR / f"{NAME}.png", dpi=600)
print("wrote", OUTDIR / f"{NAME}.png")

# ---------------------------------------------------------------------------
# audit TSVs -- every value, plotted or kept for the record
# ---------------------------------------------------------------------------
p = OUTDIR / f"{NAME}.tsv"
pd.concat([sweep, cut], ignore_index=True).to_csv(p, sep="\t", index=False, float_format="%.6f")
print("wrote", p)
p = OUTDIR / f"{NAME}_reproducibility.tsv"
# row_y holds exactly the (comparison, tool) pairs panel c drew; the rest -- the shipped
# v1 caller and the second donor's >=1-molecule arm -- are kept for the record only
repro.insert(1, "plotted", [(r.comparison, r.tool) in row_y for r in repro.itertuples()])
repro.to_csv(p, sep="\t", index=False, float_format="%.6f")
print("wrote", p)
p = OUTDIR / f"{NAME}_compute.tsv"
comp.to_csv(p, sep="\t", index=False, float_format="%.6f")
print("wrote", p)

# ---------------------------------------------------------------------------
# sidecar caption -- the single source of the journal legend
# ---------------------------------------------------------------------------
_p10 = {t: float(plotted[(plotted.tool == t) & (plotted.cutoff_bp == 10)].atlas_agreement_precision.iloc[0])
        / float(plotted[(plotted.tool == t) & (plotted.cutoff_bp == 100)].atlas_agreement_precision.iloc[0])
        for t in ORDER}
_mr = repro[(repro.comparison == MOUSE_CMP) & (repro.window_bp == 100)]
rr = _mr[_mr.tool == DEFAULT_ARM].concordance
r1 = _mr[_mr.tool == SENS_ARM].concordance
_d12 = (float(_mr[(_mr.tool == DEFAULT_ARM) & (_mr.direction == "mouse1->mouse2")].concordance.iloc[0])
        - float(_mr[(_mr.tool == SENS_ARM) & (_mr.direction == "mouse1->mouse2")].concordance.iloc[0]))
_d21 = (float(_mr[(_mr.tool == DEFAULT_ARM) & (_mr.direction == "mouse2->mouse1")].concordance.iloc[0])
        - float(_mr[(_mr.tool == SENS_ARM) & (_mr.direction == "mouse2->mouse1")].concordance.iloc[0]))
_ch = _mr[_mr.tool.isin([t for t, _s, _l in MOUSE_ROWS])].chance_mean
_ch_lo, _ch_hi = float(_ch.min()), float(_ch.max())
_ch_pt = float(_mr[_mr.tool == DEFAULT_ARM].chance_mean.max())
_ch_sa = float(_mr[_mr.tool == "scAPAtrap (v1)"].chance_mean.max())
jd = repro[(repro.tool == DEFAULT_ARM) & (repro.comparison == MOUSE_CMP)]
j25 = float(jd[jd.window_bp == 25].jaccard_strand_aware.iloc[0])
j100 = float(jd[jd.window_bp == 100].jaccard_strand_aware.iloc[0])
_25 = repro[(repro.comparison == MOUSE_CMP) & (repro.window_bp == 25) & (repro.tool == DEFAULT_ARM)].concordance
dfl = {ds: sweep[(sweep.dataset == ds) & (sweep.min_molecules == 2)].iloc[0]
       for ds in ("pbmc", "mouse1", "mouse2")}
k10 = {ds: sweep[(sweep.dataset == ds) & (sweep.min_molecules == 10)].iloc[0]
       for ds in ("pbmc", "mouse1", "mouse2")}
v1 = comp[comp.tool == "PeakATail v1"].iloc[0]
v2 = comp[comp.tool == "PeakATail"].iloc[0]
un = comp[comp.tool == "PeakATail v2 uncontended"].iloc[0]
shp = comp[comp.tool == "PeakATail shipped"].iloc[0]
_v1i = comp[comp.tool == "PeakATail v1 IP arm"].iloc[0]
_v2i = comp[comp.tool == "PeakATail v2 IP arm"].iloc[0]
_dd = {(d, w): float(repro[(repro.tool == DEFAULT_ARM) & (repro.direction == d)
                           & (repro.window_bp == w)].concordance.iloc[0])
       for d in ("donor2->donor1", "donor1->donor2") for w in (25, 100)}
_ds = {(d, w): float(repro[(repro.tool == SENS_ARM) & (repro.direction == d)
                           & (repro.window_bp == w)].concordance.iloc[0])
       for d in ("donor2->donor1", "donor1->donor2") for w in (25, 100)}
_dnull = float(repro[repro.comparison == DONOR_CMP].chance_mean.max())

footer = (
    f"Code {CODE} (PRs #93/#96/#97 merged); every PeakATail value is from the v2 run recorded in "
    "manuscript/19_final_gate_v2.md and results/benchmark_tools/final_v2_verify/VERIFIED_v2.md (verifier verdict FIXED). "
    "DEFINITIONS — atlas-agreement precision @W = fraction of called sites whose nearest same-strand curated PolyASite 2.0 "
    "representative site lies ≤W bp away (bedtools closest -s -d -t first, point mode); detected-gene recall R_det@W = fraction of "
    "the detected-gene-restricted atlas (285,136 human / 126,686 mouse sites) with a called site ≤W bp on the same strand; "
    "full-atlas recall (569,005 / 301,006 sites) is in the audit TSV. Molecule support = distinct UMI-deduplicated poly(A)-clip "
    "molecules at a site (pas.bed column 5); the ≥2-molecule default and its P ≥ 0.50 gate were pre-registered (13 §1, committed "
    "01:18:59 on 2026-08-21, before the arms started at 02:59:36) and the operating curve those thresholds sweep is **Fig 2**. "
    "REPLICATION (panel c) — the cross-mouse testis evidence is 18. Its v2 (merged-code) section, verified SOUND 2026-09-02, "
    "is what Fig 5 draws: per-gene effects ρ = 0.641 (n = 923) and 11,219 / 6,070 / 9,480 same-direction replicated PAS for "
    "SPC→RS / RS→ES / SPC→ES, 0 in all 15 null pairings. The superseded v1-code (4efeb125) record of the same quantities is "
    "ρ = 0.655 (n = 917; null 0.001 ± 0.032) and 6,049 / 9,469 / 11,263; 18 §v1→v2 tabulates both sides. "
    "COMPETITORS — competitor precision, recall, reproducibility and compute are from their own verified runs (15 §3, 09, 05 and "
    "results/benchmark_tools/gse104556/depth_concordance.tsv); they were not re-run for v2, so panels c and d compare the v2 caller "
    "against v1-era competitor measurements, and no competitor was run on the second donor. scUTRquant* is catalog-based: shown, not "
    "ranked with the de novo tools. scTail is absent (unrunnable on this BAM, R1 = 28 bp); SCAPTURE: mouse 1 only (its mouse-2 "
    "site-level run completed and scored 0.672 — 15 §3; only the per-cell PASquant step failed, so it has no cross-mouse "
    "reproducibility row in c). "
    "SCOPE — PBMC 10k v3 is one donor and one CellRanger BAM; the second donor (pbmc4k) is a different chemistry and depth; the two "
    "mice are one study and one chemistry. Peak RSS is GNU time -v maximum resident set size as kbytes/1e6, the convention of 19 §3 "
    "and 15 §5; the plotted wall time was measured with all four arms running concurrently. "
    "Every plotted value: results/figures/manuscript/fig3_tradeoff{,_reproducibility,_compute}.tsv."
)

removed_md = f"""**What this figure deliberately no longer shows** (2026-09-02 publication design pass,
`manuscript/figures/DESIGN_DIRECTIVES.md` items 4 and 5; every number below is still computed, still asserted
and still written to this figure's audit TSVs).
*The trade surface* — the molecule-support sweep that was panel a is now **Fig 2 panels a/b**, so the paper
holds exactly one precision/recall plane. It is still recomputed by this script through the shared module
`scripts/manuscript_figures/_molsweep.py` and written to `fig3_tradeoff.tsv` with `plotted = False`, and the
assert that its ≥1 and ≥2 points reproduce `19_final_gate_v2.md` §1 to 4 dp on all three arms still runs here.
For the record: raising the threshold from the ≥2-molecule default to ≥10 molecules moves atlas-agreement precision @100 bp
{dfl['pbmc'].atlas_agreement_precision:.3f} → {k10['pbmc'].atlas_agreement_precision:.3f} (PBMC),
{dfl['mouse1'].atlas_agreement_precision:.3f} → {k10['mouse1'].atlas_agreement_precision:.3f} and
{dfl['mouse2'].atlas_agreement_precision:.3f} → {k10['mouse2'].atlas_agreement_precision:.3f} (mice) while R_det falls
{dfl['pbmc'].recall_detected_genes:.3f} → {k10['pbmc'].recall_detected_genes:.3f} /
{dfl['mouse1'].recall_detected_genes:.3f} → {k10['mouse1'].recall_detected_genes:.3f} /
{dfl['mouse2'].recall_detected_genes:.3f} → {k10['mouse2'].recall_detected_genes:.3f} — the shape of the surface, never a proposal.
The pre-registered ≥2-molecule defaults themselves are PBMC {dfl['pbmc'].atlas_agreement_precision:.4f} / {dfl['pbmc'].recall_detected_genes:.4f}
(n {int(dfl['pbmc'].n):,}); mouse 1 {dfl['mouse1'].atlas_agreement_precision:.4f} / {dfl['mouse1'].recall_detected_genes:.4f}; mouse 2
{dfl['mouse2'].atlas_agreement_precision:.4f} / {dfl['mouse2'].recall_detected_genes:.4f}.
The old panel also encoded the call count in marker area (n {int(sweep.n.min()):,}–{int(sweep.n.max()):,}); that axis is now Fig 2 panels c/d, and
every n is a column of `fig3_tradeoff.tsv`. Its 3-seed gene-body-shuffled null P@100 at the default was
{float(dfl['pbmc'].null_P_mean):.4f} (PBMC) and {float(dfl['mouse1'].null_P_mean):.4f} / {float(dfl['mouse2'].null_P_mean):.4f} (mice), against the plotted
{dfl['pbmc'].atlas_agreement_precision:.3f} / {dfl['mouse1'].atlas_agreement_precision:.3f} / {dfl['mouse2'].atlas_agreement_precision:.3f}.
*Version history* — panel d used to carry a v1 point and a v1 → v2 arrow, and panel c a shipped-caller row. End
users choose a caller, not a release: the improvement is **Fig S8** (compute) and **Fig S12** (versions), with one
Results sentence pointing there. The values keep their rows: v1 {v1.wall_str} / {v1.peak_rss_gb:.1f} GB → current
{v2.wall_str} / {v2.peak_rss_gb:.2f} GB, i.e. {v1.peak_rss_gb / v2.peak_rss_gb:.1f}× less peak RSS on the same
no-IP arm (`fig3_tradeoff_compute.tsv`, `kind = not_plotted_history`), so no ≥300 GB node is needed — and on the
IP-filtered arm the same pair is {_v1i.wall_str} / {_v1i.peak_rss_gb:.1f} GB → {_v2i.wall_str} / {_v2i.peak_rss_gb:.2f} GB
(`kind = audit_only`, matching 15 §5); the shipped
coverage-only caller reproduced {float(repro[(repro.tool == 'PeakATail shipped (v1)') & (repro.window_bp == 100)].concordance.min()):.3f}–{float(repro[(repro.tool == 'PeakATail shipped (v1)') & (repro.window_bp == 100)].concordance.max()):.3f} across the two mice
(`fig3_tradeoff_reproducibility.tsv`).
*The concurrency circle* — panel d's open "uncontended" marker is gone; the caveat it carried is the sentence in
the panel-d paragraph above, and the measurement ({un.wall_str}) is the `not_plotted_caveat` row of the compute TSV.
*The gate lines* — P ≥ 0.50 and the superseded P ≥ 0.38 belong to the precision/recall plane and are drawn once,
in Fig 2."""

cap_md = f"""# Fig 3 — `fig3_tradeoff` caption (generated by `scripts/manuscript_figures/fig3_tradeoff.py`)

## Legend

Figure 3 | Robustness: matching-window resolution, replicate reproducibility and compute. (The file stem
`fig3_tradeoff` is retained from the figure-rename pass for file identity; the trade curve itself moved to
Fig 2 in the 2026-09-02 design pass — see "What this figure deliberately no longer shows".)

**Panels a and b — matching-window sensitivity.** Atlas-agreement precision (a) and detected-gene recall (b) at
10 / 25 / 50 / 100 bp for the ≥2-molecule precision default, the ≥1-molecule sensitivity arm and the five
competitors on PBMC 10k v3, from the verified per-tool score TSVs (values at 200 bp and the mouse arms are in the
audit TSV, marked `plotted = False`). The check the panels exist for: the de novo precision ordering does not
change when the window is tightened to 10 bp, and the ratio P@10 / P@100 is {_p10[DEFAULT_ARM]:.2f} for the
default and {_p10['polyApipe']:.2f} for polyApipe against {_p10['SCAPTURE']:.2f} (SCAPTURE), {_p10['Sierra']:.2f}
(Sierra) and {_p10['scAPAtrap']:.2f} (scAPAtrap) — the precision lead is not an artefact of a loose window.
Catalog-based scUTRquant* is above the default at every window ({_p10['scUTRquant*']:.2f}) and is shown but not
ranked with the de novo tools.

**Panel c — replicate agreement.** The fraction of one replicate's call set with a strand-matched call in the
other within 100 bp (`bedtools closest -s -d -t first`), both directions, with the ≤25 bp value as a diamond.
*Biological replicates (top block):* the two GSE104556 testis mice, computed here for the default
({rr.min():.3f}–{rr.max():.3f} at 100 bp, {float(_25.min()):.3f}–{float(_25.max()):.3f} at 25 bp; strand-aware Jaccard of the ±window site
sets {j25:.3f} at ±25 bp and {j100:.3f} at ±100 bp) and for the ≥1-molecule arm ({r1.min():.3f}–{r1.max():.3f}) —
the ≥2-molecule threshold buys {_d12 * 100:.1f} / {_d21 * 100:.1f} points of replicate agreement over the
≥1-molecule arm. *Human donors (bottom block):* the pre-registered second-donor comparison (26 §5, verifier
FIXED) — {_dd[('donor2->donor1', 100)]*100:.1f}% of donor-2 default calls reproduce within 100 bp in donor 1
({_dd[('donor2->donor1', 25)]*100:.1f}% within 25 bp) against a genic-shuffle null of ≤{_dnull:.4f}; the
≥1-molecule arm gives {_ds[('donor2->donor1', 100)]*100:.1f}% / {_ds[('donor1->donor2', 100)]*100:.1f}%.
**The reverse direction is capped by call-count arithmetic, not by disagreement**: donor 1 has {int(repro[(repro.tool == DEFAULT_ARM) & (repro.direction == 'donor1->donor2') & (repro.window_bp == 100)].n_query.iloc[0]):,}
default calls against donor 2's {int(repro[(repro.tool == DEFAULT_ARM) & (repro.direction == 'donor2->donor1') & (repro.window_bp == 100)].n_query.iloc[0]):,}, so its ceiling is {_ceil:.4f} and the observed
{_dd[('donor1->donor2', 100)]:.4f} is {_dd[('donor1->donor2', 100)] / _ceil * 100:.0f}% of it — that is what the "ceiling" tick marks. No competitor
was run on the second donor, so the human block holds PeakATail arms only; **Fig S4 carries the full
pre-registered second-donor protocol** (gate, per-direction nulls, ceiling arithmetic and the matched-call-count
control that shows donor 2's higher headline precision is an operating-point effect, not an improvement).
Competitor and shipped-caller mouse rows are the v1-era head-to-head run
(`results/benchmark_tools/gse104556/depth_concordance.tsv`, verified in `05_figure_index.md`), so the mouse block
compares a v2 PeakATail arm against v1-era competitor measurements. Chance ≤0.010 (each query set against the
other replicate's three gene-body-shuffled nulls); the dotted line sits at 0.0095, just above the **highest**
per-tool chance level, and the six mouse values span {_ch_lo:.4f}–{_ch_hi:.4f} in the audit TSV, because a denser call set has a
higher chance of a nearby match — the default's own chance level is {_ch_pt:.4f} against scAPAtrap's {_ch_sa:.4f}, so the
raw ranking is the conservative reading of this panel. **Honest scope: the default is above polyApipe and above
PeakATail's own shipped caller, but below Sierra (0.790 / 0.834) and scAPAtrap (0.884 / 0.793) in both
directions — and scAPAtrap's lead is measured on a set its own `reducePeaks(min.cells=10, min.count=10)` step
has already depth-cleaned.** Supporting replication evidence is `18_spermatogenesis_final.md`, whose **v2 (merged-code) section was
re-run after #96/#97 and verified SOUND on 2026-09-02** — that is the record **Fig 5** draws: cross-mouse per-gene
switch effects ρ = 0.641 (n = 923) and 11,219 / 6,070 / 9,480 same-direction replicated PAS for
SPC→RS / RS→ES / SPC→ES, 0 in all 15 null pairings. The superseded v1-code (`4efeb125`) values of the same
quantities — ρ = 0.655 (n = 917; null 0.001 ± 0.032) and 6,049 / 9,469 / 11,263 — are kept here because earlier
drafts of this legend quoted them; 18 §"v1 → v2" tabulates both sides digit for digit. (Corrected 2026-09-03: this
legend previously said the v2 re-run "has not happened", which contradicted Fig 5 and the Fig 5 manifest row.)

**Panel d — compute.** Wall time against peak RSS on the PBMC 10k v3 BAM, log–log, for the **current** caller
(code {CODE}: {v2.wall_str}, {v2.peak_rss_gb:.2f} GB) beside the five competitors' own verified runs on the same
box (`results/figures/manuscript/benchmark_headtohead.tsv`, verified in `05_figure_index.md`). The plotted
PeakATail point is the clip-seeded arm run *without* `--ip-filter` (the pair `19` §3 quotes); the IP-filtered
default arm whose call sets panels a–c show was cheaper ({_v2i.wall_str} / {_v2i.peak_rss_gb:.2f} GB, read here
from that arm's own `runtime_mem.txt`) and is an audit row in `fig3_tradeoff_compute.tsv`, not a plotted point.
That wall time was measured with all four arms running concurrently, so peak RSS is the production number and the
wall time must be quoted with the concurrency disclosed (~28–35 min; the uncontended single run was {un.wall_str},
`19` §3). Each competitor row carries its run's retry / skipped-stage caveat in the compute TSV's `note` column —
scAPAtrap's {float(comp[comp.tool == 'scAPAtrap'].wall_h.iloc[0]):.2f} h is a resumed run that skipped three completed stages and so understates a
from-scratch run (12:58:50 of total machine time across both attempts). Peak RSS is GNU `time -v` maximum
resident set size reported as kbytes/1e6, the convention of `19` §3 and `15` §5.

{removed_md}

**Definitions, pre-registration and scope** (the former on-figure footer, verbatim — it travels with the
legend): {footer}

## Provenance

**Index paragraph.** `fig3_tradeoff` — Fig 3: the robustness of the shipping caller. The de novo precision
ordering is unchanged at a 10 bp matching window (P@10 / P@100 {_p10[DEFAULT_ARM]:.2f} for the default against
{_p10['scAPAtrap']:.2f} for scAPAtrap), the two testis mice agree on {rr.min():.3f}–{rr.max():.3f} of default sites at 100 bp and the two
human donors on {_dd[('donor2->donor1', 100)]:.3f} (donor 2 → donor 1) against a ≤0.010 chance level, and the caller runs the PBMC BAM in
{v2.wall_str} at {v2.peak_rss_gb:.2f} GB peak RSS — the cheapest wall time of any de novo tool benchmarked here.

**Not drawn — what Fig 3 still lacks.** The long-read re-ranking panel (Spearman ρ between atlas and Kinnex
precision orderings, 0.90 at 7 of 8 arms) is quarantined and was not regenerated for v2, so it is not in this
figure; the fairness table (≥95% gene-proximal for every tool) is a table, not a panel; and no competitor was
re-run on the v2 code or on the second donor, so no v2-vs-v2 reproducibility or compute comparison exists.

Sources: `manuscript/19_final_gate_v2.md` + `results/benchmark_tools/final_v2_verify/VERIFIED_v2.md`
(FIXED), `manuscript/26_second_donor_preregistration.md` + `results/benchmark_tools/pbmc4k_donor2/concordance.txt`,
`manuscript/18_spermatogenesis_final.md`, `manuscript/15_final_gate.md` §3/§5,
`manuscript/09_headtohead_results.md`, `manuscript/05_figure_index.md`, the shared sweep module
`scripts/manuscript_figures/_molsweep.py`, and the per-file `source` columns of
`results/figures/manuscript/fig3_tradeoff.tsv`, `fig3_tradeoff_reproducibility.tsv` and
`fig3_tradeoff_compute.tsv`.
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

print(f"code {CODE}: cross-mouse agreement@100 {rr.min():.3f}-{rr.max():.3f}; "
      f"cross-donor {_dd[('donor2->donor1', 100)]:.3f} / {_dd[('donor1->donor2', 100)]:.3f} "
      f"(ceiling {_ceil:.3f}); compute {v2.wall_str} / {v2.peak_rss_gb:.2f} GB; "
      f"callouts per panel {callouts}")

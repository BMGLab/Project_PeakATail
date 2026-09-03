"""Molecule-support sweep of the v2 IP-filtered arms -- ONE implementation, two callers.

DESIGN_DIRECTIVES.md item 5 folds the operating (trade) curve out of Fig 3 panel a
and into Fig 2's precision/recall plane, so exactly one P/R plane exists in the paper.
The numbers must not drift between the two scripts, so the computation that used to
live in `fig3_tradeoff.py::sweep_scores()` moved here verbatim (2026-09-02):

  * `fig2_accuracy.py` PLOTS the curve (panels a/b) and writes it to
    `results/figures/manuscript/fig2_accuracy_sweep.tsv`.
  * `fig3_tradeoff.py` still WRITES its rows into `fig3_tradeoff.tsv` for the record
    (marked `plotted = False`), so that audit TSV keeps every value it had.

Every point is recomputed by `scripts/benchmark_tools/score_tool.py` with the reference
arguments of `scripts/benchmark_tools/stage2_final_launch.sh` -- nothing is interpolated.
The >=1 and >=2 points are asserted against manuscript/19 section 1 (verifier FIXED).

Heavy steps cache into $FIG3_WORKDIR (name kept from the original home so the recorded
reproduction commands still work); delete the directory to force a recompute.
"""
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
BT = WD / "results/benchmark_tools"
WORK = Path(os.environ.get("FIG3_WORKDIR", "/mnt/ssd0/emaout/fig3_trade_repro_work"))

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
PREREG_K = 2          # the ONLY pre-registered point on the curve (13 section 1)
SENSITIVITY_K = 1     # the user-choosable sensitivity arm

DSLAB = {"pbmc": "PBMC 10k v3", "mouse1": "testis mouse 1", "mouse2": "testis mouse 2"}

# manuscript/19 section 1 (verified FIXED): (dataset, k) -> (n, P@100, R_det@100).
# If these fail the pipeline drifted -- both callers inherit the guard.
VERIFIED_19 = {
    ("pbmc", 2):   (46524, 0.7062, 0.1754), ("pbmc", 1):   (167565, 0.3520, 0.2685),
    ("mouse1", 2): (26255, 0.7450, 0.2048), ("mouse1", 1): (52792, 0.5686, 0.2806),
    ("mouse2", 2): (26526, 0.7572, 0.2080), ("mouse2", 1): (51842, 0.5880, 0.2826),
}


def _sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, executable="/bin/bash",
                          env=dict(os.environ, LC_ALL="C"), capture_output=True, text=True).stdout


def molecule_sweep(verbose=True):
    """DataFrame of the >=1/2/3/5/10 molecule-support points on all three v2 IP arms."""
    WORK.mkdir(parents=True, exist_ok=True)
    rows = []
    for ds, (sp, arm, _pfx) in ARM.items():
        out = WORK / ds
        out.mkdir(parents=True, exist_ok=True)
        for k in THRESHOLDS:
            lab = f"{ds}_ge{k}mol"
            tsv = out / f"score_{lab}.tsv"
            if not tsv.exists():
                bed = out / f"pas_ge{k}mol.bed"
                _sh(f"awk -F'\\t' -v K={k} '$5>=K' {arm/'pas.bed'} > {bed}")
                refs = HUMAN_REFS if sp == "human" else MOUSE_REFS
                subprocess.run([PY, SCORE, str(bed), lab, "--outdir", str(out)] + refs,
                               check=True, env=dict(os.environ, LC_ALL="C"),
                               stdout=subprocess.DEVNULL)
                print(f"  computed {tsv}")
            d = pd.read_csv(tsv, sep="\t")

            def pick(panel, ref, series="real"):
                return d[(d.panel == panel) & (d.series == series) & (d.reference == ref)
                         & (d.cutoff_bp == 100.0)]

            pr = pick("precision", "atlas_full").iloc[0]
            rd = pick("recall", "atlas_detected").iloc[0]
            rf = pick("recall", "atlas_full").iloc[0]
            f1 = pick("f1", "atlas_detected").iloc[0]
            nulls = pick("precision", "atlas_full", "null_genic").sort_values("replicate")["value"].tolist()
            rows.append(dict(dataset=ds, arm=f"v2 IP arm, >={k} molecules",
                             min_molecules=k, pre_registered=(k == PREREG_K),
                             cutoff_bp=100, n=int(pr.n_query), n_matched=int(pr.n_matched),
                             atlas_agreement_precision=float(pr.value),
                             recall_detected_genes=float(rd.value),
                             recall_detected_denominator=int(rd.n_query),
                             recall_full_atlas=float(rf.value),
                             recall_full_denominator=int(rf.n_query),
                             F1_detected_genes=float(f1.value),
                             null_P_seed1=nulls[0], null_P_seed2=nulls[1], null_P_seed3=nulls[2],
                             null_P_mean=float(np.mean(nulls)),
                             source=str(tsv.relative_to(WORK.parent)) + " (recomputed by score_tool.py)"))
    sweep = pd.DataFrame(rows)

    for (ds, k), (n, p, r) in VERIFIED_19.items():
        q = sweep[(sweep.dataset == ds) & (sweep.min_molecules == k)].iloc[0]
        assert int(q.n) == n, (ds, k, q.n, n)
        assert abs(q.atlas_agreement_precision - p) < 5e-5, (ds, k, q.atlas_agreement_precision, p)
        assert abs(q.recall_detected_genes - r) < 5e-5, (ds, k, q.recall_detected_genes, r)
    if verbose:
        print("molecule sweep: >=1 and >=2 points reproduce manuscript/19 section 1 "
              "to 4 dp on all three arms")
    return sweep

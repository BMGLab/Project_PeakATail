#!/usr/bin/env python3
"""
fig5_spermatogenesis.py -- the verified testis 3'-UTR control
(manuscript/18_spermatogenesis_final.md; both GSE104556 mice).

VERSION SWITCH -- env var SPERMATOGENESIS_VERSION picks which record is drawn:

  v2  (default)                       merged caller, frozen tree 9dfdefb3 (#96 + #97);
                                      results/stage3_spermatogenesis_v2/; the record of
                                      18's v2 section (verified SOUND, 2026-09-02)
  v1  (SPERMATOGENESIS_VERSION=v1)    final caller 4efeb125; results/stage3_spermatogenesis_final/;
                                      18's v1 body (verified FIXED, 2026-08-21), kept
                                      reproducible as the labelled comparison

      export LC_ALL=C; SPERMATOGENESIS_VERSION=v1 python3 .../fig5_spermatogenesis.py

The claim the figure carries, in the exact form that survives (18, v2 section):

    Per-gene 3'UTR usage shifts progressively proximal along spermatocyte -> round ->
    elongating spermatid in both mice: 31.4% and 30.5% of depth-guarded genes shorten
    monotonically across the three stages against label-shuffle nulls of 17.4% and 21.0%
    (z = 10.5 and 6.4, 20 shuffles), exceeding monotone lengthening (486 vs 318 and
    366 vs 291 genes), and the composition-controlled per-cell distal-usage residual falls
    at every step (Cliff's delta SPC vs ES = 0.55 and 0.63).  It is a per-gene, equal-weight
    statement: the per-gene medians are not monotone and the across-gene UMI-weighted distal
    index reverses at RS->ES.

Nothing is typed in that is not read from the sources below (the literal p-values in the
footer are quoted from manuscript/18 and asserted to be present in that file):

  <SRC>/mouseN/summary/pdui_summary.json        panels a,b,d1
  <SRC>/mouseN/summary/pdui_null_shuffles.tsv   panel a (20-shuffle null)
  <SRC>/mouseN/summary/pdui_per_cell.tsv        panel b (per-cell residual)
  <SRC>/summary/switch_cross_mouse_replication.tsv  panel c1 + switch TRUE hits
  <SRC>/summary/switch_null_summary.tsv         panel c1 (null calibration)
  <SRC>/summary/gene_cross_mouse_delta.tsv      panel c2 (per-gene scatter)
  <SRC>/summary/gene_replication_summary.json   panel c2 (rho + null) , panel d2
  <SRC>/summary/gene_literature_panel.tsv       panel d2 (16 genes)
  <SRC>/REPORT_PROVISIONAL.md                   provenance asserts
  manuscript/18_spermatogenesis_final.md        claim + cautions

where <SRC> = results/stage3_spermatogenesis_v2 (v2) or results/stage3_spermatogenesis_final (v1).

Panels d1/d2 are the negative findings and are drawn inside a dashed frame: they are NOT
claim carriers.

Outputs
-------
manuscript/figures/fig5_spermatogenesis.{png,pdf}      600 dpi / fonttype 42 (no Type 3)
manuscript/figures/fig5_spermatogenesis.caption.md     '## Legend' (journal legend, incl. the
                                                       binding cautions) + '## Provenance'
results/figures/manuscript/fig5_spermatogenesis.png
results/figures/manuscript/fig5_spermatogenesis.tsv                  panel a (monotone fractions vs null)
results/figures/manuscript/fig5_spermatogenesis_null_shuffles.tsv    panel a (all 20 null draws per mouse)
results/figures/manuscript/fig5_spermatogenesis_percell_resid.tsv    panel b (per-stage summary + Cliff deltas)
results/figures/manuscript/fig5_spermatogenesis_percell_values.tsv   panel b (every plotted cell)
results/figures/manuscript/fig5_spermatogenesis_replication.tsv      panel c1 (per stage pair, TRUE + null)
results/figures/manuscript/fig5_spermatogenesis_gene_scatter.tsv     panel c2 (every plotted gene)
results/figures/manuscript/fig5_spermatogenesis_negatives.tsv        panel d1 + d2
results/figures/manuscript/fig5_spermatogenesis_reference_lines.tsv  every number in title/footer + source
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")          # <= 4 threads; matplotlib/numpy stay single-threaded

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, FancyBboxPatch
import numpy as np
import pandas as pd

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

WD = "/mnt/ssd1/Projects/PeakATail_wd"
FIGDIR = f"{WD}/manuscript/figures"
TSVDIR = f"{WD}/results/figures/manuscript"
os.makedirs(FIGDIR, exist_ok=True); os.makedirs(TSVDIR, exist_ok=True)
NAME = "fig5_spermatogenesis"

# ---------------------------------------------------------------------------
# Which record: v2 (default; merged caller 9dfdefb3, 18's v2 section, verified
# SOUND 2026-09-02) or v1 (final caller 4efeb125, 18's v1 body, verified FIXED
# 2026-08-21, kept reproducible as the labelled comparison).
# ---------------------------------------------------------------------------
VERSION = os.environ.get("SPERMATOGENESIS_VERSION", "v2").lower()
VERSIONS = {
    "v1": dict(
        src="results/stage3_spermatogenesis_final",
        code="4efeb125", commit_full="4efeb1252e6d7c7d51b252b443b5be5947ae000f",
        caller_phrase="PeakATail final caller",
        report_commit_assert="commit `4efeb1252e6d7c7d51b252b443b5be5947ae000f` (`4efeb125`)",
        report_flag_assert="PRs #96",       # v1 report flags the v2 re-run as pending
        # the literal p-values quoted in the footer live in manuscript/18 (v1 body) --
        # asserted below so the figure can never drift from the write-up
        binom_p={"mouse1": "9.9e-9", "mouse2": "0.0125"},
        claim_asserts=("binomial\np 9.9e-9",
                       "p 0.0125; 1 of 40 null draws reached the same excess in absolute value"),
        expect_m2_null_draws_ge_excess=1,   # 18 (v1 body): "1 of 40"
        genes_with_pair_approx="~10,400",
        prm_share_text="~26% of the guarded ES UMIs",
        record_note="manuscript/18 v1 body (verdict FIXED, 2026-08-21); superseded by the v2 section",
        footer_status=("PROVISIONAL: code is the frozen tree 4efeb125, and a v2 re-run after PRs #96 "
                       "(IP-filter/strand) and #97 will redo every number here"),
        v2_pending_ref=("v2_rerun_pending", "PRs #96 (IP-filter/strand) and #97",
                        "REPORT_PROVISIONAL.md status line; manuscript/18 header"),
    ),
    "v2": dict(
        src="results/stage3_spermatogenesis_v2",
        code="9dfdefb3", commit_full="9dfdefb3eb353b0817ef79c4eb9ace6d6c8aab53",
        caller_phrase="PeakATail merged caller",
        report_commit_assert="commit `9dfdefb3eb353b0817ef79c4eb9ace6d6c8aab53` (`9dfdefb3`)",
        report_flag_assert="PR #96 IP-filter/strand + PR #97 clip-memory",
        # quoted from manuscript/18's v2 section (verifier-recomputed from the raw outputs)
        binom_p={"mouse1": "3.4e-9", "mouse2": "0.0039"},
        claim_asserts=("binomial\np 3.4e-9",
                       "p 0.0039; 0 of 40 null draws reached the\nexcess in absolute value"),
        expect_m2_null_draws_ge_excess=0,   # 18 (v2 section): "0 of 40"
        genes_with_pair_approx="~10,500",
        prm_share_text="25.2% / 14.5% of guarded ES UMIs (m1/m2)",
        record_note="manuscript/18 v2 section (verdict SOUND, 2026-09-02); v1 retained as the labelled comparison",
        footer_status=("RECORD: code is the frozen merged tree 9dfdefb3 (#96 IP-filter/strand + #97), verified "
                       "SOUND 2026-09-02; the v1 render (4efeb125) stays reproducible via SPERMATOGENESIS_VERSION=v1"),
        v2_pending_ref=("v1_comparison_render", "SPERMATOGENESIS_VERSION=v1",
                        "manuscript/18 v2 section; the v1-code render is the labelled comparison, not the record"),
    ),
}
assert VERSION in VERSIONS, f"SPERMATOGENESIS_VERSION must be one of {sorted(VERSIONS)}, got {VERSION!r}"
V = VERSIONS[VERSION]
SRC = f"{WD}/{V['src']}"
CODE = V["code"]                  # frozen caller tree (18 / REPORT_PROVISIONAL §0)
STAGES = ["SPC", "RS", "ES"]      # SPG is excluded from the claim (manuscript/11)
MICE = ["mouse1", "mouse2"]
MOUSE_LABEL = {"mouse1": "mouse 1", "mouse2": "mouse 2"}
CHANCE = 1.0 / 6.0                # a uniformly random ordering of 3 stages is monotone 1/6 of the time

# Okabe-Ito, as mandated.  Colour follows the ROLE and is fixed across panels; every series also
# carries a marker or fill style so nothing depends on hue alone.
BLUE, SKY, GREEN, ORANGE, PINK, VERM, GREY = "#0072B2", "#56B4E9", "#009E73", "#E69F00", "#CC79A7", "#D55E00", "#999999"
INK, MUTED, GRID, SURFACE, FAINT = "#1B2429", "#5A6B73", "#D8E0E3", "#FFFFFF", "#C6D0D4"
STAGE_COL = {"SPC": SKY, "RS": GREEN, "ES": BLUE}
DIR_COL = {"shortening": BLUE, "lengthening": ORANGE}
MOUSE_MARK = {"mouse1": "o", "mouse2": "s"}
MOUSE_LS = {"mouse1": "-", "mouse2": (0, (4, 2))}

# ---------------------------------------------------------------------------
# 1. read the sources (+ provenance asserts)
# ---------------------------------------------------------------------------
report_txt = open(f"{SRC}/REPORT_PROVISIONAL.md").read()
claim_txt = open(f"{WD}/manuscript/18_spermatogenesis_final.md").read()
assert V["report_commit_assert"] in report_txt, "report does not record the frozen commit"
assert "PROVISIONAL" in report_txt and V["report_flag_assert"] in report_txt, \
    "report does not carry the expected version flag"
# the literal p-values quoted in the footer live in manuscript/18 and nowhere else on disk
# (v1 body and v2 section respectively; both guard strings survive the v2 supersession
# because the v1 claim block is retained as the labelled comparison)
for s in V["claim_asserts"]:
    assert s in claim_txt, f"manuscript/18 no longer states the quoted excess fact: {s!r}"
BINOM_P = V["binom_p"]                                   # quoted from manuscript/18 (asserted above)

S = {m: json.load(open(f"{SRC}/{m}/summary/pdui_summary.json")) for m in MICE}
NULLS = {m: pd.read_csv(f"{SRC}/{m}/summary/pdui_null_shuffles.tsv", sep="\t") for m in MICE}
CELLS = {m: pd.read_csv(f"{SRC}/{m}/summary/pdui_per_cell.tsv", sep="\t") for m in MICE}
rep = pd.read_csv(f"{SRC}/summary/switch_cross_mouse_replication.tsv", sep="\t")
swnull = pd.read_csv(f"{SRC}/summary/switch_null_summary.tsv", sep="\t")
gene = pd.read_csv(f"{SRC}/summary/gene_cross_mouse_delta.tsv", sep="\t")
grep_ = json.load(open(f"{SRC}/summary/gene_replication_summary.json"))
lit = pd.read_csv(f"{SRC}/summary/gene_literature_panel.tsv", sep="\t")

for m in MICE:
    assert S[m]["info"]["n_shuffles"] == len(NULLS[m]) == 20
    assert S[m]["info"]["min_umi_per_stage"] == 50
    assert set(NULLS[m]["n_genes3_fixed"]) == {S[m]["true_3stage_pseudobulk"]["n_genes"]}
N_SHUF = 20
GUARD_UMI = S["mouse1"]["info"]["min_umi_per_stage"]

# ---------------------------------------------------------------------------
# 2. panel a: monotone fraction TRUE vs the 20-shuffle label null
# ---------------------------------------------------------------------------
a_rows = []
for m in MICE:
    t = S[m]["true_3stage_pseudobulk"]
    for d, tcol, ncol_fixed, ncol_free, ncount in (
        ("shortening", "frac_monotone_shortening", "frac_mono_short3_fixed", "frac_mono_short3", "n_monotone_shortening"),
        ("lengthening", "frac_monotone_lengthening", "frac_mono_long3_fixed", "frac_mono_long3", "n_monotone_lengthening"),
    ):
        nf = NULLS[m][ncol_fixed].to_numpy(); nr = NULLS[m][ncol_free].to_numpy()
        a_rows.append(dict(
            mouse=m, direction=d, n_genes_guarded=int(t["n_genes"]), n_monotone=int(t[ncount]),
            true_frac=float(t[tcol]),
            null_fixed_mean=float(nf.mean()), null_fixed_sd=float(nf.std(ddof=1)),
            null_fixed_min=float(nf.min()), null_fixed_max=float(nf.max()),
            null_free_mean=float(nr.mean()), null_free_sd=float(nr.std(ddof=1)),
            z_vs_null_fixed=float((t[tcol] - nf.mean()) / nf.std(ddof=1)),  # sample sd (ddof=1), as in REPORT_PROVISIONAL §2.2
            ratio_true_over_null=float(t[tcol] / nf.mean()),
            n_null_shuffles=N_SHUF, chance_1_in_6=CHANCE))
A = pd.DataFrame(a_rows)
A.to_csv(f"{TSVDIR}/{NAME}.tsv", sep="\t", index=False)


def arow(m, d):
    return A[(A.mouse == m) & (A.direction == d)].iloc[0]


# the direction-specific excess and how often the null reaches it (two-sided, all 40 draws:
# 20 fixed-gene-set + 20 free-guard) -- this is the "1/40 shuffles" of manuscript/18
exc_rows = []
for m in MICE:
    te = float(arow(m, "shortening")["true_frac"] - arow(m, "lengthening")["true_frac"])
    fx = (NULLS[m]["frac_mono_short3_fixed"] - NULLS[m]["frac_mono_long3_fixed"]).to_numpy()
    fr = (NULLS[m]["frac_mono_short3"] - NULLS[m]["frac_mono_long3"]).to_numpy()
    both = np.concatenate([fx, fr])
    exc_rows.append(dict(mouse=m, true_excess_frac=te,
                         n_monotone_shortening=int(arow(m, "shortening")["n_monotone"]),
                         n_monotone_lengthening=int(arow(m, "lengthening")["n_monotone"]),
                         binomial_p_quoted_from_manuscript_18=BINOM_P[m],
                         null_excess_min=float(both.min()), null_excess_max=float(both.max()),
                         n_null_draws=int(both.size),
                         n_null_draws_abs_ge_true=int((np.abs(both) >= te).sum())))
EXC = pd.DataFrame(exc_rows).set_index("mouse")
assert int(EXC.loc["mouse2", "n_null_draws_abs_ge_true"]) == V["expect_m2_null_draws_ge_excess"], \
    f"manuscript/18's {V['expect_m2_null_draws_ge_excess']}/40 no longer reproduces"

null_long = []
for m in MICE:
    for _, r in NULLS[m].iterrows():
        for d, cf, cr in (("shortening", "frac_mono_short3_fixed", "frac_mono_short3"),
                          ("lengthening", "frac_mono_long3_fixed", "frac_mono_long3")):
            null_long.append(dict(mouse=m, seed=int(r["seed"]), direction=d,
                                  frac_fixed_gene_set=float(r[cf]), frac_free_guard=float(r[cr]),
                                  n_genes_fixed=int(r["n_genes3_fixed"]), n_genes_free=int(r["n_genes3"]),
                                  cell_cliffs_resid_SPC_ES=float(r["cell_cliffs_resid_SPC_ES"])))
pd.DataFrame(null_long).to_csv(f"{TSVDIR}/{NAME}_null_shuffles.tsv", sep="\t", index=False)

# ---------------------------------------------------------------------------
# 3. panel b: composition-controlled per-cell distal-usage residual by stage
# ---------------------------------------------------------------------------
PB = {}
b_rows, b_vals = [], []
for m in MICE:
    c = CELLS[m][CELLS[m].stage.isin(STAGES)].copy()
    PB[m] = {s: c[c.stage == s]["mean_resid_pdui"].to_numpy() for s in STAGES}
    for s in STAGES:
        v = PB[m][s]
        assert len(v) == S[m]["info"]["stage_counts"][s]
        b_rows.append(dict(mouse=m, stage=s, n_cells=int(len(v)), median=float(np.median(v)),
                           q1=float(np.percentile(v, 25)), q3=float(np.percentile(v, 75)),
                           p5=float(np.percentile(v, 5)), p95=float(np.percentile(v, 95)),
                           minimum=float(v.min()), maximum=float(v.max()),
                           median_n_genes_per_cell=float(S[m]["per_cell"]["median_n_genes_per_cell"][s])))
    b_vals.append(c.assign(mouse=m)[["mouse", "cell", "stage", "mean_resid_pdui", "median_pdui", "wdi", "n_genes", "umi_at_pairs"]])
B = pd.DataFrame(b_rows)
CLIFF_PAIRS = [("SPC", "RS"), ("RS", "ES"), ("SPC", "ES")]
b_cliff = []
for m in MICE:
    for a_, b_ in CLIFF_PAIRS:
        k = f"mean_resid_pdui_{a_}_vs_{b_}"
        nullc = NULLS[m]["cell_cliffs_resid_SPC_ES"].to_numpy()
        b_cliff.append(dict(mouse=m, comparison=f"{a_}_vs_{b_}", n_a=int(S[m]["per_cell"][k]["n_a"]),
                            n_b=int(S[m]["per_cell"][k]["n_b"]),
                            cliffs_delta=float(S[m]["per_cell"][k]["cliffs_delta_a_minus_b"]),
                            mannwhitney_p=float(S[m]["per_cell"][k]["mannwhitney_p"]),
                            median_diff_b_minus_a=float(S[m]["per_cell"][k]["median_diff_b_minus_a"]),
                            null_cliff_SPC_vs_ES_min=float(nullc.min()) if (a_, b_) == ("SPC", "ES") else np.nan,
                            null_cliff_SPC_vs_ES_max=float(nullc.max()) if (a_, b_) == ("SPC", "ES") else np.nan,
                            n_null_shuffles=N_SHUF if (a_, b_) == ("SPC", "ES") else np.nan))
BC = pd.DataFrame(b_cliff)
pd.concat([B.assign(row="per_stage_summary"), BC.assign(row="cliffs_delta")], ignore_index=True) \
    .to_csv(f"{TSVDIR}/{NAME}_percell_resid.tsv", sep="\t", index=False)
pd.concat(b_vals, ignore_index=True).to_csv(f"{TSVDIR}/{NAME}_percell_values.tsv", sep="\t", index=False)


def cliff(m, a_, b_):
    return float(BC[(BC.mouse == m) & (BC.comparison == f"{a_}_vs_{b_}")]["cliffs_delta"].iloc[0])


NULL_CLIFF = {m: (float(NULLS[m]["cell_cliffs_resid_SPC_ES"].min()), float(NULLS[m]["cell_cliffs_resid_SPC_ES"].max())) for m in MICE}

# ---------------------------------------------------------------------------
# 4. panel c1: cross-mouse replication of switch-test hits (arm B0)
# ---------------------------------------------------------------------------
PAIRS = [("RS_vs_SPC", "SPC → RS"), ("ES_vs_RS", "RS → ES"), ("ES_vs_SPC", "SPC → ES")]
true_rep = rep[rep.pairing == "true"].set_index("pair")
null_rep = rep[rep.pairing != "true"]
c_rows = []
for pid, lab in PAIRS:
    t = true_rep.loc[pid]
    n = null_rep[null_rep.pair == pid]
    c_rows.append(dict(pair=pid, label=lab, n_pas_matched=int(t["n_matched"]),
                       hits_m1=int(t["n_hits_m1"]), hits_m2=int(t["n_hits_m2"]),
                       hits_m1_matched=int(t["n_hits_m1_matched"]), hits_m2_matched=int(t["n_hits_m2_matched"]),
                       replicated_same_dir=int(t["n_replicated_same_dir"]),
                       both_hit_opposite_dir=int(t["n_both_hit_opposite_dir"]),
                       sign_agreement=float(t["sign_agreement_among_both_hit"]),
                       frac_m1_matched_hits_replicated=float(t["frac_m1_matched_hits_replicated"]),
                       frac_m2_matched_hits_replicated=float(t["frac_m2_matched_hits_replicated"]),
                       expected_same_dir_if_independent=float(t["exp_same_dir_independent"]),
                       pas_shuffle_same_dir_mean=float(t["shuf_same_dir_mean"]),
                       pas_shuffle_same_dir_max=int(t["shuf_same_dir_max"]),
                       enrichment_over_independent=float(t["enrichment_same_dir_over_independent"]),
                       n_null_pairings=int(len(n)),
                       null_replicated_same_dir_total=int(n["n_replicated_same_dir"].sum()),
                       null_replicated_same_dir_max=int(n["n_replicated_same_dir"].max()),
                       null_hits_m1_total=int(n["n_hits_m1"].sum()), null_hits_m2_total=int(n["n_hits_m2"].sum())))
C = pd.DataFrame(c_rows)
C.to_csv(f"{TSVDIR}/{NAME}_replication.tsv", sep="\t", index=False)
N_NULL_PAIRINGS = int(C["n_null_pairings"].sum())              # 15 = 5 shuffles x 3 stage pairs
assert int(C["null_replicated_same_dir_total"].sum()) == 0
SW_TRUE = {"mouse1": int(true_rep["n_hits_m1"].sum()), "mouse2": int(true_rep["n_hits_m2"].sum())}
SW_NULLP05 = {m: float(swnull[swnull.mouse == m]["null_frac_p05_mean"].mean()) for m in MICE}
SW_NULL_FAM = int(len(swnull) * swnull["n_runs"].iloc[0])      # 30 BH families (3 pairs x 5 shuffles x 2 mice)
assert int(swnull["null_q05_total"].sum()) == 0 and int(swnull["null_runs_with_hit"].sum()) == 0

# ---------------------------------------------------------------------------
# 5. panel c2: per-gene effect, mouse 1 vs mouse 2
# ---------------------------------------------------------------------------
g = gene.copy()
g["category"] = np.where((g.m1_delta == 0) | (g.m2_delta == 0), "delta 0 in one mouse",
                np.where((g.m1_delta < 0) & (g.m2_delta < 0), "shorten in both",
                np.where((g.m1_delta > 0) & (g.m2_delta > 0), "lengthen in both", "discordant")))
counts = g["category"].value_counts().to_dict()
assert counts["shorten in both"] == grep_["both_shorten"] and counts["lengthen in both"] == grep_["both_lengthen"]
assert counts["discordant"] == grep_["strictly_discordant"] and counts["delta 0 in one mouse"] == grep_["zero_delta_one_mouse"]
assert len(g) == grep_["n_common"]
g[["gene_id", "symbol", "m1_delta", "m2_delta", "m1_mono_short", "m2_mono_short", "category"]] \
    .to_csv(f"{TSVDIR}/{NAME}_gene_scatter.tsv", sep="\t", index=False)
RHO, RHO_N = float(grep_["replicate_spearman_rho"]), int(grep_["n_common"])
NULL_RHO_M, NULL_RHO_SD, NULL_RHO_MAX = float(grep_["null_rho_mean"]), float(grep_["null_rho_sd"]), float(grep_["null_rho_max_abs"])
CAT_COL = {"shorten in both": BLUE, "lengthen in both": ORANGE, "discordant": GREY, "delta 0 in one mouse": FAINT}

# ---------------------------------------------------------------------------
# 6. panels d1 / d2: the negative findings
# ---------------------------------------------------------------------------
d_rows = []
for m in MICE:
    for s in STAGES:
        d_rows.append(dict(panel="d1_across_gene_umi_weighted_distal_index", mouse=m, stage=s,
                           median_wdi_across_cells=float(S[m]["per_cell"]["wdi_median_by_stage"][s]),
                           n_cells=int(S[m]["info"]["stage_counts"][s])))
    for a_, b_ in CLIFF_PAIRS:
        k = f"wdi_{a_}_vs_{b_}"
        d_rows.append(dict(panel="d1_across_gene_umi_weighted_distal_index", mouse=m, stage=f"cliff_{a_}_vs_{b_}",
                           cliffs_delta=float(S[m]["per_cell"][k]["cliffs_delta_a_minus_b"]),
                           mannwhitney_p=float(S[m]["per_cell"][k]["mannwhitney_p"])))
LIT_CATS = [("shorten in both", BLUE), ("lengthen in both", ORANGE), ("discordant", GREY), ("delta 0 in one mouse", FAINT)]
lit["category"] = np.where((lit.m1_delta == 0) | (lit.m2_delta == 0), "delta 0 in one mouse",
                  np.where((lit.m1_delta < 0) & (lit.m2_delta < 0), "shorten in both",
                  np.where((lit.m1_delta > 0) & (lit.m2_delta > 0), "lengthen in both", "discordant")))
LP = grep_["literature_panel"]
lc = lit["category"].value_counts().to_dict()
assert lc["shorten in both"] == LP["n_shorten_both"] and lc["lengthen in both"] == LP["n_lengthen_both"]
assert lc["discordant"] == LP["n_strictly_discordant"] and lc["delta 0 in one mouse"] == LP["n_zero_delta_one_mouse"]
assert len(lit) == LP["n_measurable_in_both"]
for _, r in lit.iterrows():
    d_rows.append(dict(panel="d2_literature_panel", mouse="both", stage=r["symbol"], gene_id=r["gene_id"],
                       m1_delta=float(r["m1_delta"]), m2_delta=float(r["m2_delta"]),
                       mean_delta=float(r["mean_delta"]), category=r["category"]))
pd.DataFrame(d_rows).to_csv(f"{TSVDIR}/{NAME}_negatives.tsv", sep="\t", index=False)
LIT_SYMS = {c: sorted(lit[lit.category == c]["symbol"].tolist()) for c, _ in LIT_CATS}

# ---------------------------------------------------------------------------
# 7. figure
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.size": 7.4, "axes.labelsize": 7.0, "axes.titlesize": 8.0,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.0,
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": MUTED,
    "xtick.color": MUTED, "ytick.color": MUTED, "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.titlelocation": "left", "axes.titlepad": 5, "axes.linewidth": 0.7,
    "font.family": "DejaVu Sans",
})
# Publication layout: the figure-level title, subtitle and cautions footer live in the
# .caption.md Legend now, so the canvas keeps only the panel band (height 11.4 -> 8.4 in;
# width 7.5 -> 7.1 in, the double-column norm, which the panels tolerate).
LEFT = 0.075
fig = plt.figure(figsize=(7.1, 8.4))
gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 0.80], width_ratios=[1.0, 1.0],
                      hspace=0.72, wspace=0.30, left=LEFT, right=0.980, top=0.968, bottom=0.070)
ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[0, 1])
ax_c1 = fig.add_subplot(gs[1, 0])
gs_c2 = gs[1, 1].subgridspec(2, 1, height_ratios=[1.0, 0.10], hspace=0.62)
ax_c2 = fig.add_subplot(gs_c2[0, 0])
ax_c2n = fig.add_subplot(gs_c2[1, 0])
ax_d1 = fig.add_subplot(gs[2, 0])
ax_d2 = fig.add_subplot(gs[2, 1])


def style(ax, axis="y"):
    ax.set_axisbelow(True); ax.grid(True, axis=axis, color=GRID, lw=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# --- a: monotone fraction vs the label-shuffle null ---------------------------
XA = {("mouse1", "shortening"): 0.0, ("mouse1", "lengthening"): 1.0,
      ("mouse2", "shortening"): 2.4, ("mouse2", "lengthening"): 3.4}
rng = np.random.default_rng(0)
ax_a.axhline(CHANCE, color=MUTED, lw=0.8, ls=(0, (1, 1.6)), zorder=2)
ax_a.text(3.78, CHANCE + 0.004, "chance 1/6 for a\nrandom ordering\nof 3 stages", fontsize=5.2, color=MUTED,
          ha="left", va="bottom", linespacing=1.25)
for m in MICE:
    for d in ("shortening", "lengthening"):
        x = XA[(m, d)]; r = arow(m, d)
        nv = NULLS[m]["frac_mono_short3_fixed" if d == "shortening" else "frac_mono_long3_fixed"].to_numpy()
        ax_a.scatter(x + rng.uniform(-0.14, 0.14, nv.size), nv, s=6.0, c=GREY, alpha=0.75,
                     linewidths=0, zorder=3)
        ax_a.plot([x - 0.22, x + 0.22], [nv.mean()] * 2, color=MUTED, lw=1.2, solid_capstyle="butt", zorder=4)
        ax_a.plot([x, x], [nv.mean(), r["true_frac"]], color=DIR_COL[d], lw=0.7, alpha=0.55, zorder=3)
        ax_a.plot([x], [r["true_frac"]], marker=MOUSE_MARK[m], ms=6.4, mfc=DIR_COL[d], mec=SURFACE, mew=0.8, zorder=6)
        ax_a.text(x, r["true_frac"] + 0.016, f"{r['true_frac']:.3f}", ha="center", va="bottom",
                  fontsize=6.2, color=DIR_COL[d], fontweight="bold", zorder=6)
        ax_a.text(x, r["true_frac"] + 0.040, f"z {r['z_vs_null_fixed']:.1f}", ha="center", va="bottom",
                  fontsize=5.6, color=MUTED, zorder=6)
for m in MICE:
    xl, xr = XA[(m, "shortening")], XA[(m, "lengthening")]
    yb = 0.400
    ax_a.plot([xl, xl, xr, xr], [yb - 0.012, yb, yb, yb - 0.012], color=INK, lw=0.6, zorder=5)
    ax_a.text((xl + xr) / 2, yb + 0.006, f"excess: {int(EXC.loc[m, 'n_monotone_shortening']):,} vs "
              f"{int(EXC.loc[m, 'n_monotone_lengthening']):,} genes\nbinomial p {BINOM_P[m]}",
              ha="center", va="bottom", fontsize=5.4, color=INK, linespacing=1.25, zorder=6)
    ax_a.text((xl + xr) / 2, -0.175, f"{MOUSE_LABEL[m]}\n{int(arow(m, 'shortening')['n_genes_guarded']):,} guarded genes",
              transform=ax_a.get_xaxis_transform(), ha="center", va="top", fontsize=6.4, color=INK,
              fontweight="bold", linespacing=1.2)
ax_a.set_xticks(list(XA.values()))
ax_a.set_xticklabels([f"{d}\n{int(arow(m, d)['n_monotone']):,}" for (m, d) in XA], fontsize=6.0, linespacing=1.2)
for t, (m, d) in zip(ax_a.get_xticklabels(), XA):
    t.set_color(DIR_COL[d])
ax_a.set_xlim(-0.55, 4.75); ax_a.set_ylim(0.10, 0.475)
ax_a.set_yticks([0.10, 0.15, 0.20, 0.25, 0.30, 0.35])
ax_a.set_ylabel("fraction of depth-guarded genes strictly\nmonotone across SPC → RS → ES", labelpad=3)
style(ax_a)
ax_a.set_title("a  Claim carrier: per-gene monotone shortening vs its null")
ax_a.legend(handles=[
    Line2D([], [], ls="none", marker="o", ms=5.0, mfc=BLUE, mec=SURFACE, mew=0.7, label="TRUE labels, shortening (mouse 1 ● / mouse 2 ■)"),
    Line2D([], [], ls="none", marker="o", ms=5.0, mfc=ORANGE, mec=SURFACE, mew=0.7, label="TRUE labels, lengthening"),
    Line2D([], [], ls="none", marker="o", ms=2.8, mfc=GREY, mec=GREY, label=f"{N_SHUF} label shuffles (fixed gene set); bar = null mean"),
], loc="upper center", bbox_to_anchor=(0.5, -0.365), frameon=False, ncol=1, handlelength=1.2,
    labelspacing=0.25, borderaxespad=0.0, handletextpad=0.5)

# --- b: composition-controlled per-cell residual -------------------------------
XB = {("mouse1", "SPC"): 0.0, ("mouse1", "RS"): 0.85, ("mouse1", "ES"): 1.70,
      ("mouse2", "SPC"): 3.05, ("mouse2", "RS"): 3.90, ("mouse2", "ES"): 4.75}
ax_b.axhline(0.0, color=MUTED, lw=0.7, ls=(0, (1, 1.6)), zorder=2)
for (m, s), x in XB.items():
    v = PB[m][s]
    parts = ax_b.violinplot([v], positions=[x], widths=0.66, showextrema=False, showmedians=False)
    for pc in parts["bodies"]:
        pc.set_facecolor(STAGE_COL[s]); pc.set_alpha(0.30); pc.set_edgecolor(STAGE_COL[s]); pc.set_linewidth(0.7); pc.set_zorder(3)
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    ax_b.plot([x, x], [np.percentile(v, 5), np.percentile(v, 95)], color=STAGE_COL[s], lw=0.8, zorder=4)
    ax_b.plot([x, x], [q1, q3], color=STAGE_COL[s], lw=3.0, solid_capstyle="butt", zorder=5)
    ax_b.plot([x], [med], marker=MOUSE_MARK[m], ms=4.4, mfc=SURFACE, mec=STAGE_COL[s], mew=1.1, zorder=6)
for m in MICE:
    xs = [XB[(m, s)] for s in STAGES]
    meds = [np.median(PB[m][s]) for s in STAGES]
    ax_b.plot(xs, meds, color=INK, lw=0.8, ls=MOUSE_LS[m], alpha=0.8, zorder=5)
    for (a_, b_), yb in zip([("SPC", "RS"), ("RS", "ES")], [0.0655, 0.0655]):
        xl, xr = XB[(m, a_)], XB[(m, b_)]
        ax_b.plot([xl, xl, xr, xr], [yb - 0.004, yb, yb, yb - 0.004], color=MUTED, lw=0.6, zorder=5)
        ax_b.text((xl + xr) / 2, yb + 0.0015, f"δ {cliff(m, a_, b_):.2f}", ha="center", va="bottom",
                  fontsize=5.6, color=INK, zorder=6)
    xl, xr = XB[(m, "SPC")], XB[(m, "ES")]
    ax_b.plot([xl, xl, xr, xr], [0.0805 - 0.004, 0.0805, 0.0805, 0.0805 - 0.004], color=INK, lw=0.7, zorder=5)
    ax_b.text((xl + xr) / 2, 0.0820, f"SPC vs ES δ {cliff(m, 'SPC', 'ES'):.2f}\n{N_SHUF}-shuffle null "
              f"[{NULL_CLIFF[m][0]:+.2f}, {NULL_CLIFF[m][1]:+.2f}]", ha="center", va="bottom",
              fontsize=5.4, color=INK, linespacing=1.2, zorder=6)
    ax_b.text((XB[(m, 'SPC')] + XB[(m, 'ES')]) / 2, -0.175, f"{MOUSE_LABEL[m]}\n{sum(S[m]['info']['stage_counts'][s] for s in STAGES):,} labelled cells",
              transform=ax_b.get_xaxis_transform(), ha="center", va="top", fontsize=6.4, color=INK,
              fontweight="bold", linespacing=1.2)
ax_b.set_xticks(list(XB.values()))
ax_b.set_xticklabels([f"{s}\n{S[m]['info']['stage_counts'][s]:,}" for (m, s) in XB], fontsize=6.0, linespacing=1.2)
for t, (m, s) in zip(ax_b.get_xticklabels(), XB):
    t.set_color(STAGE_COL[s])
ax_b.set_xlim(-0.60, 5.35); ax_b.set_ylim(-0.052, 0.108)
ax_b.set_yticks([-0.04, -0.02, 0.0, 0.02, 0.04])
ax_b.set_ylabel("per-cell distal-usage residual", labelpad=3)
style(ax_b)
ax_b.set_title("b  Composition-controlled per-cell residual")
# (the residual definition and the violin/IQR/whisker glyph key moved to the caption Legend)

# --- c1: cross-mouse replication of switch hits ---------------------------------
xc = np.arange(len(C)); w = 0.34
ax_c1.bar(xc - w / 2, C["replicated_same_dir"], width=w * 0.94, color=BLUE, zorder=3, label="replicated in both mice, same direction")
ax_c1.bar(xc + w / 2, C["expected_same_dir_if_independent"], width=w * 0.94, facecolor="none",
          edgecolor=MUTED, lw=0.9, zorder=3, label="expected if the two mice were independent")
for i, r in C.iterrows():
    ax_c1.text(i - w / 2, r["replicated_same_dir"] + 320, f"{int(r['replicated_same_dir']):,}",
               ha="center", va="bottom", fontsize=6.2, color=BLUE, fontweight="bold")
    ax_c1.text(i - w / 2, r["replicated_same_dir"] + 1180, f"{r['sign_agreement']:.1%} sign\nagreement",
               ha="center", va="bottom", fontsize=5.3, color=INK, linespacing=1.15)
    ax_c1.text(i + w / 2, r["expected_same_dir_if_independent"] + 320, f"{r['enrichment_over_independent']:.1f}×",
               ha="center", va="bottom", fontsize=5.8, color=MUTED)
    ax_c1.plot([i], [0], marker="v", ms=4.2, mfc=VERM, mec=VERM, zorder=6, clip_on=False)
ax_c1.set_xticks(xc)
ax_c1.set_xticklabels([f"{r['label']}\n{int(r['n_pas_matched']):,} PAS matched\n"
                       f"({int(r['hits_m1_matched']):,} / {int(r['hits_m2_matched']):,} hit)" for _, r in C.iterrows()],
                      fontsize=5.9, linespacing=1.25)
ax_c1.set_ylim(0, 17500); ax_c1.set_yticks([0, 4000, 8000, 12000, 16000])
ax_c1.set_yticklabels(["0", "4k", "8k", "12k", "16k"])
ax_c1.spines["left"].set_bounds(0, 16000)
ax_c1.set_ylabel("PAS with q < 0.05 in both mice,\nsame sign of Δ proportion", labelpad=3)
style(ax_c1)
ax_c1.set_title("c  Cross-mouse replication: PAS hits (left) and per-gene effects (right)")
# (the arm-B0 switch-test description and PAS-matching rule moved to the caption Legend)
ax_c1.legend(handles=[Patch(facecolor=BLUE, label="replicated in both mice, same direction"),
                      Patch(facecolor="none", edgecolor=MUTED, lw=0.9, label="expected if the mice were independent (analytic; PAS-shuffle agrees)"),
                      Line2D([], [], ls="none", marker="v", ms=4.2, mfc=VERM, mec=VERM,
                             label=f"0 replicated in all {N_NULL_PAIRINGS} null pairings (5 shuffles × 3 pairs)")],
             loc="upper center", bbox_to_anchor=(0.5, -0.30), frameon=False, ncol=1, handlelength=1.2,
             labelspacing=0.25, borderaxespad=0.0, handletextpad=0.5)

# --- c2: per-gene effect scatter -------------------------------------------------
for cat in ("delta 0 in one mouse", "discordant", "lengthen in both", "shorten in both"):
    sub = g[g.category == cat]
    ax_c2.scatter(sub["m1_delta"], sub["m2_delta"], s=5.0, c=CAT_COL[cat], alpha=0.75, linewidths=0,
                  zorder=(3 if cat in ("delta 0 in one mouse", "discordant") else 4))
ax_c2.axhline(0, color=MUTED, lw=0.6, zorder=2); ax_c2.axvline(0, color=MUTED, lw=0.6, zorder=2)
ax_c2.plot([-1, 1], [-1, 1], color=GRID, lw=0.8, ls=(0, (3, 2)), zorder=2)
ax_c2.set_xlim(-1.04, 1.04); ax_c2.set_ylim(-1.04, 1.04)
ax_c2.set_xticks([-1, -0.5, 0, 0.5, 1]); ax_c2.set_yticks([-1, -0.5, 0, 0.5, 1])
ax_c2.set_xlabel("mouse 1  Δ PDUI (ES − SPC)", labelpad=2)
ax_c2.set_ylabel("mouse 2  Δ PDUI (ES − SPC)", labelpad=2)
TBOX = dict(boxstyle="round,pad=0.18", fc=SURFACE, ec="none", alpha=0.82)
ax_c2.text(-0.98, -0.98, f"shorten in both\n{counts['shorten in both']:,} ({counts['shorten in both'] / len(g):.1%})",
           ha="left", va="bottom", fontsize=5.6, color=BLUE, linespacing=1.2, zorder=6, bbox=TBOX)
ax_c2.text(0.98, 0.98, f"lengthen in both\n{counts['lengthen in both']:,} ({counts['lengthen in both'] / len(g):.1%})",
           ha="right", va="top", fontsize=5.6, color=ORANGE, linespacing=1.2, zorder=6, bbox=TBOX)
ax_c2.text(0.98, -0.98, f"discordant {counts['discordant']:,}\nΔ = 0 in one mouse {counts['delta 0 in one mouse']:,}",
           ha="right", va="bottom", fontsize=5.6, color=MUTED, linespacing=1.2, zorder=6, bbox=TBOX)
ax_c2.text(-0.98, 0.98, f"Spearman ρ = {RHO:.3f}\nn = {RHO_N:,} genes guarded in both\nPearson r = "
           f"{float(grep_['replicate_pearson_r']):.3f}", ha="left", va="top", fontsize=6.0, color=INK,
           fontweight="bold", linespacing=1.3, zorder=6, bbox=TBOX)
style(ax_c2, axis="both")

# null band for that rho, as its own strip under the scatter
ax_c2n.add_patch(plt.Rectangle((-NULL_RHO_MAX, -1), 2 * NULL_RHO_MAX, 2, color=GREY, alpha=0.30, lw=0, zorder=2))
ax_c2n.add_patch(plt.Rectangle((NULL_RHO_M - NULL_RHO_SD, -1), 2 * NULL_RHO_SD, 2, color=GREY, alpha=0.80, lw=0, zorder=3))
ax_c2n.plot([RHO], [0], marker="D", ms=5.2, mfc=BLUE, mec=SURFACE, mew=0.8, zorder=5, clip_on=False)
ax_c2n.text(RHO, 1.35, f"{RHO:.3f}", ha="center", va="bottom", fontsize=5.6, color=BLUE, fontweight="bold")
ax_c2n.text(NULL_RHO_M, 1.35, f"null {NULL_RHO_M:.3f} ± {NULL_RHO_SD:.3f}", ha="center", va="bottom", fontsize=5.2, color=MUTED)
ax_c2n.set_xlim(-0.22, 0.80); ax_c2n.set_ylim(-1, 1)
ax_c2n.set_yticks([]); ax_c2n.set_xticks([-0.2, 0, 0.2, 0.4, 0.6, 0.8])
ax_c2n.tick_params(axis="x", labelsize=5.4, pad=1.5, length=2)
for sp in ("top", "right", "left"):
    ax_c2n.spines[sp].set_visible(False)
ax_c2n.set_xlabel(f"cross-mouse Spearman ρ of the per-gene effect, against\n"
                  f"the gene-correspondence null (200 permutations; dark band\n"
                  f"± 1 sd, light band the full range, max |ρ| {NULL_RHO_MAX:.3f})",
                  fontsize=5.4, labelpad=3)

# --- d1: the across-gene UMI-weighted index (NEGATIVE) ---------------------------
xd = np.arange(3)
for m, col, dy in zip(MICE, (VERM, PINK), (-1, 1)):
    ys = [S[m]["per_cell"]["wdi_median_by_stage"][s] for s in STAGES]
    ax_d1.plot(xd, ys, color=col, lw=1.4, ls=MOUSE_LS[m], marker=MOUSE_MARK[m], ms=5.0,
               mfc=col, mec=SURFACE, mew=0.7, zorder=4, label=MOUSE_LABEL[m])
    for x, y in zip(xd, ys):
        ax_d1.text(x, y + dy * 0.0055, f"{y:.3f}", ha="center", va=("top" if dy < 0 else "bottom"),
                   fontsize=5.5, color=col, zorder=5)
ax_d1.text(1.42, 0.694, f"reversal at RS → ES\nCliff δ {S['mouse1']['per_cell']['wdi_RS_vs_ES']['cliffs_delta_a_minus_b']:.2f}"
           f" / {S['mouse2']['per_cell']['wdi_RS_vs_ES']['cliffs_delta_a_minus_b']:.2f}",
           ha="center", va="top", fontsize=5.6, color=VERM, fontweight="bold", linespacing=1.25)
ax_d1.set_xticks(xd); ax_d1.set_xticklabels(STAGES, fontsize=6.5)
ax_d1.set_xlim(-0.35, 2.45); ax_d1.set_ylim(0.555, 0.700)
ax_d1.set_yticks([0.56, 0.60, 0.64, 0.68])
ax_d1.set_ylabel("median per-cell distal index,\nUMIs pooled ACROSS genes", labelpad=3)
style(ax_d1)
ax_d1.set_title("d  Negative findings — not claim carriers", color=VERM)
ax_d1.legend(loc="upper left", frameon=False, handlelength=2.0, labelspacing=0.25, borderaxespad=0.2)
# (the protamine-ceiling explanation and the wdi name-collision caution moved to the caption Legend)

# --- d2: the 16-gene literature panel (NEGATIVE) ----------------------------------
yl = np.arange(len(LIT_CATS))
for j, (cat, col) in enumerate(LIT_CATS):
    n = len(LIT_SYMS[cat])
    ax_d2.barh(j, n, color=col, height=0.56, zorder=3)
    ax_d2.text(n + 0.16, j, f"{n}   " + ", ".join(LIT_SYMS[cat]), va="center", ha="left", fontsize=5.5, color=INK)
ax_d2.set_yticks(yl); ax_d2.set_yticklabels([c.replace("delta 0", "Δ = 0") for c, _ in LIT_CATS], fontsize=6.0)
for t, (_, col) in zip(ax_d2.get_yticklabels(), LIT_CATS):
    t.set_color(col if col != FAINT else MUTED)
ax_d2.invert_yaxis(); ax_d2.set_xlim(0, 10.8); ax_d2.set_xticks([0, 2, 4, 6])
ax_d2.set_xlabel(f"genes, of the {LP['n_measurable_in_both']} of the {LP['n_in_panel']}-gene curated panel\n"
                 f"measurable in both mice", labelpad=2)
style(ax_d2, axis="x")
ax_d2.set_title("   the 16-gene literature panel does not reproduce", color=VERM)
# (the quarantined-6/6 context and the Prm/Tnp ceiling-PDUI caution moved to the caption Legend)

# dashed frame around row 3 so the negatives cannot be read as part of the claim
fig.canvas.draw()
rend = fig.canvas.get_renderer()
from matplotlib.transforms import Bbox
bb = Bbox.union([ax_d1.get_tightbbox(rend), ax_d2.get_tightbbox(rend)]) \
     .transformed(fig.transFigure.inverted())
pad_x, pad_y = 0.010, 0.010
frame = FancyBboxPatch((bb.x0 - pad_x, bb.y0 - pad_y), bb.width + 2 * pad_x, bb.height + 2 * pad_y,
                       boxstyle="round,pad=0.002,rounding_size=0.006", transform=fig.transFigure, fill=False,
                       edgecolor=VERM, linewidth=0.9, linestyle=(0, (4, 2)), zorder=0, clip_on=False)
fig.add_artist(frame)

# --- no on-figure title / subtitle / footer: journal style ------------------------
# The figure-level title, the subtitle paragraph and the numbered cautions footer all
# live in the caption sidecar's '## Legend' now (surgery pass, 2026-09-02).  m1s/m2s
# stay: the VERDICT print and the reference-lines TSV still use them.
m1s, m2s = arow("mouse1", "shortening"), arow("mouse2", "shortening")

# --- save ---------------------------------------------------------------------------
fig.savefig(f"{TSVDIR}/{NAME}.png", dpi=600, bbox_inches="tight", pad_inches=0.07)
for ext in ("png", "pdf"):
    p = f"{FIGDIR}/{NAME}.{ext}"
    fig.savefig(p, bbox_inches="tight", pad_inches=0.07, **({"dpi": 600} if ext == "png" else {}))
    print("wrote", p)

# --- reference lines: every number that appears only in the caption Legend ------------
# (formerly the on-figure title/footer; the TSV is unchanged so the audit trail holds)
ref = [
    dict(item="claim_sentence_source", value="manuscript/18_spermatogenesis_final.md", source=f"{V['record_note']}; the figure states this claim and no stronger one"),
    dict(item="code_frozen_tree", value=CODE, source=f"REPORT_PROVISIONAL.md §0 (commit {V['commit_full']})"),
    dict(item="binomial_p_excess_mouse1", value=BINOM_P["mouse1"],
         source=f"manuscript/18 (shortening-vs-lengthening excess, {int(EXC.loc['mouse1', 'n_monotone_shortening'])} vs {int(EXC.loc['mouse1', 'n_monotone_lengthening'])})"),
    dict(item="binomial_p_excess_mouse2", value=BINOM_P["mouse2"],
         source=f"manuscript/18 ({int(EXC.loc['mouse2', 'n_monotone_shortening'])} vs {int(EXC.loc['mouse2', 'n_monotone_lengthening'])})"),
    dict(item="null_draws_reaching_mouse2_excess", value=f"{int(EXC.loc['mouse2', 'n_null_draws_abs_ge_true'])}/{int(EXC.loc['mouse2', 'n_null_draws'])}",
         source=f"computed here from mouse2/summary/pdui_null_shuffles.tsv (|short-long| over 20 fixed + 20 free-guard draws); reproduces 18's '{V['expect_m2_null_draws_ge_excess']}/40'"),
    dict(item="median_pdui_mouse1_SPC_RS_ES", value="/".join(f"{S['mouse1']['true_3stage_pseudobulk']['median_pdui'][s]:.4f}" for s in STAGES),
         source="mouse1/summary/pdui_summary.json (NOT monotone: RS > SPC)"),
    dict(item="median_pdui_mouse2_SPC_RS_ES", value="/".join(f"{S['mouse2']['true_3stage_pseudobulk']['median_pdui'][s]:.4f}" for s in STAGES),
         source="mouse2/summary/pdui_summary.json (NOT monotone: RS > SPC)"),
    dict(item="mean_pdui_mouse1_SPC_RS_ES", value="/".join(f"{S['mouse1']['true_3stage_pseudobulk']['mean_pdui'][s]:.4f}" for s in STAGES),
         source="mouse1/summary/pdui_summary.json (monotone decreasing; not plotted)"),
    dict(item="mean_pdui_mouse2_SPC_RS_ES", value="/".join(f"{S['mouse2']['true_3stage_pseudobulk']['mean_pdui'][s]:.4f}" for s in STAGES),
         source="mouse2/summary/pdui_summary.json (monotone decreasing; not plotted)"),
    dict(item="spg_cells_mouse1_mouse2", value=f"{S['mouse1']['info']['stage_counts']['SPG']}/{S['mouse2']['info']['stage_counts']['SPG']}",
         source="pdui_summary.json info.stage_counts (SPG excluded from the claim per manuscript/11)"),
    dict(item="switch_true_q05_hits_mouse1", value=SW_TRUE["mouse1"], source="summary/switch_cross_mouse_replication.tsv, sum of n_hits_m1 over the 3 TRUE rows"),
    dict(item="switch_true_q05_hits_mouse2", value=SW_TRUE["mouse2"], source="summary/switch_cross_mouse_replication.tsv, sum of n_hits_m2 over the 3 TRUE rows"),
    dict(item="switch_null_frac_p05_mouse1", value=round(SW_NULLP05["mouse1"], 5), source="summary/switch_null_summary.tsv, mean over the 3 pairs (5 shuffles each)"),
    dict(item="switch_null_frac_p05_mouse2", value=round(SW_NULLP05["mouse2"], 5), source="summary/switch_null_summary.tsv, mean over the 3 pairs (5 shuffles each)"),
    dict(item="switch_null_bh_families_with_a_hit", value=f"0/{SW_NULL_FAM}", source="summary/switch_null_summary.tsv (null_q05_total and null_runs_with_hit both 0)"),
    dict(item="replicate_pearson_r", value=round(float(grep_["replicate_pearson_r"]), 4), source="summary/gene_replication_summary.json"),
    dict(item="binom_p_both_shorten_vs_both_lengthen", value=round(float(grep_["binom_p_shorten_vs_lengthen"]), 5),
         source=f"summary/gene_replication_summary.json ({grep_['both_shorten']} vs {grep_['both_lengthen']}; per-gene reproducibility and net shortening are different claims)"),
    dict(item="genes_monotone_shortening_in_both_mice", value=f"{grep_['n_monotone_shortening_in_both']} ({grep_['frac_monotone_shortening_in_both']:.3f})",
         source=f"summary/gene_replication_summary.json (expected under independence {float(grep_['expected_frac_mono_both_if_independent']):.3f})"),
    dict(item="label_agreement_gex_vs_pas_clusters", value="0.9878 / 0.9877", source="REPORT_PROVISIONAL.md §0 (3-stage set)"),
    dict(item="protamine_share_of_guarded_ES_umis", value=V["prm_share_text"], source="manuscript/18 and REPORT_PROVISIONAL.md caveat 8 (Prm2 + Prm1)"),
    dict(item=V["v2_pending_ref"][0], value=V["v2_pending_ref"][1], source=V["v2_pending_ref"][2]),
]
pd.DataFrame(ref).to_csv(f"{TSVDIR}/{NAME}_reference_lines.tsv", sep="\t", index=False)

# --- edge check: nothing may be clipped; the outer 8 px must be blank -------------------
try:
    from PIL import Image
    im = np.asarray(Image.open(f"{FIGDIR}/{NAME}.png").convert("L"))
    edge = 8
    border = np.concatenate([im[:edge, :].ravel(), im[-edge:, :].ravel(), im[:, :edge].ravel(), im[:, -edge:].ravel()])
    ink = int((border < 250).sum())
    print(f"edge check: {im.shape[1]}x{im.shape[0]} px, ink pixels in the outer {edge} px = {ink}")
    assert ink == 0, "INK IN THE OUTER 8 PX -- something is clipped"
except ImportError:
    print("edge check skipped (no PIL)")

print(A.to_string(index=False))
print(C[["label", "replicated_same_dir", "sign_agreement", "enrichment_over_independent", "null_replicated_same_dir_total"]].to_string(index=False))
print(f"VERDICT: monotone shortening {m1s['true_frac']:.3f} (z {m1s['z_vs_null_fixed']:.1f}) / {m2s['true_frac']:.3f} "
      f"(z {m2s['z_vs_null_fixed']:.1f}); per-cell residual Cliff SPC vs ES {cliff('mouse1', 'SPC', 'ES'):.3f} / "
      f"{cliff('mouse2', 'SPC', 'ES'):.3f}; cross-mouse rho {RHO:.3f} (n {RHO_N})")

# ---------------------------------------------------------------------------
# sidecar caption -- script-written (the figure convention requires the
# script, never a hand edit, to own the caption).  Version-selected: the v1
# text is the verified 2026-08-21 caption (content unchanged), the v2 text the
# verified 2026-09-02 record (manuscript/18 v2 section).  At the 2026-09-02
# legend surgery both were restructured into '## Legend' (the journal legend:
# 'Figure 5 | title', panel-by-panel description, and the binding cautions
# that used to be the on-figure footer) + '## Provenance' (sources, commits,
# audit TSVs); every sentence formerly printed on the image lives here now.
# ---------------------------------------------------------------------------
CAPTION_V1 = r"""# Fig 5 — `fig5_spermatogenesis` caption (generated by `scripts/manuscript_figures/fig5_spermatogenesis.py`; renamed 2026-09-02, see FIGURE_MAP.tsv)

## Legend

Figure 5 | **Per-gene 3′UTR usage shifts progressively proximal along spermatocyte → round → elongating spermatid, in both mice.** GSE104556 testis, both mice, PeakATail final caller (frozen tree `4efeb125`); stage labels are marker-panel argmax over leiden clusters of the STARsolo **gene-expression** matrix, joined to the PAS matrix by exact cell-barcode match (1,294 / 1,364 labelled cells; 98.8% agreement with the independent PAS-clustering labels). PDUI = distal / (proximal + distal), higher = longer 3′UTR; proximal/distal = first/last PAS of the gene in transcription order, strand taken from `pasbed.bed` and re-asserted per row (this is the strand/counts re-validation of 13 §2, which the `per_gene` path **passes**). The per-gene statistic is the **pseudobulk PDUI per stage** over genes with ≥ 50 UMI at the PAS pair in **every** stage (1,553 / 1,194 genes of ~10,400 with a pair). **31.5% and 30.2%** of guarded genes shorten monotonically across SPC > RS > ES (**489 / 1,553** and **361 / 1,194**) against label-shuffle nulls of **17.5% and 21.0%** (z = **10.5** and **5.8**, 20 shuffles, one label per cell, composition preserved), exceeding monotone lengthening (489 vs 325 and 361 vs 296 genes), and the composition-controlled per-cell distal-usage residual falls at every step (Cliff's δ SPC vs ES = **0.56** and **0.61**, outside the entire 20-shuffle null range in both mice). The claim is a **per-gene, equal-weight** statement: the per-gene medians are not monotone and the across-gene UMI-weighted distal index reverses at RS → ES (panel d).

**Panel a (the claim carrier)** — fraction of depth-guarded genes strictly monotone across the three
stages, per mouse and per direction, with the 20 label-shuffle draws as grey dots (fixed TRUE gene set;
bar = null mean) and the TRUE value as a filled marker (● mouse 1, ■ mouse 2): shortening 0.315 vs
0.175 ± 0.013 (z 10.5, 1.80×) and 0.302 vs 0.210 ± 0.016 (z 5.8, 1.44×); lengthening 0.209 vs
0.178 ± 0.010 (z 3.1, 1.17×) and 0.248 vs 0.204 ± 0.021 (z 2.1, 1.21×). Dotted line = 1/6, the chance
rate for a random ordering of three stages. The bracket carries the direction-specific evidence: the
**excess** of shortening over lengthening (489 vs 325, binomial p 9.9e-9; 361 vs 296, p 0.0125).
**Panel b** — the composition-controlled per-cell index: for each cell, the mean over genes of
(per-cell PDUI − that gene's 3-stage pseudobulk PDUI), which removes the gene-composition confound
(the median cell carries 296 / 326.5 / 116.5 genes over SPC / RS / ES in mouse 1, so raw per-cell
medians are not comparable across stages). Medians fall at every step in both mice
(+0.0017 → −0.0036 → −0.0115 and +0.0028 → −0.0010 → −0.0084); Cliff's δ SPC vs RS 0.28 / 0.28,
RS vs ES 0.41 / 0.46, SPC vs ES 0.56 / 0.61 against 20-shuffle null ranges of [−0.11, +0.10] and
[−0.07, +0.07]. Violin = all cells, thick bar = IQR, whisker = 5–95th percentile, open marker =
median; δ = Cliff's delta (earlier stage more distal). **Panel c** — reliability payload. Left: PAS called switching at q < 0.05 in **both**
mice with the same sign of Δ proportion (arm B0: `fisher --count-mode cells --marker-top-n 0`), per
stage pair — 11,263 / 6,049 / 9,469 of 36,826 / 35,644 / 35,372 coordinate-matched PAS (gene + strand,
≤ 100 bp, greedy 1:1), sign agreement 99.7–99.8%, 3.2–5.1× the independence expectation (hollow bars;
the 10-draw PAS-shuffle mean agrees with the analytic value to within 1.3%), and **0 replicated in all 15 null
pairings** (▼). TRUE switch hits over the three pairs: 50,415 / 52,043; label-shuffle null 3.2% / 3.2%
at p < 0.05 with **0 q < 0.05 hits in all 30 null BH families**. Right: the per-gene effect
Δ PDUI (ES − SPC) in mouse 1 vs mouse 2 over the 917 genes guarded in both — Spearman ρ = **0.655**
(Pearson r = 0.853) against a gene-correspondence null of 0.001 ± 0.032 (200 permutations, max |ρ|
0.110; strip below the scatter). 338 genes (36.9%) shorten in both mice, 278 (30.3%) lengthen in both,
244 are strictly discordant, 57 have Δ exactly 0 in one mouse. **Panel d (boxed — NEGATIVE findings,
not claim carriers)** — left: the **across-gene UMI-weighted** per-cell distal index rises at RS → ES
in both mice (0.589 → 0.578 → 0.638 and 0.641 → 0.636 → 0.672; Cliff δ RS vs ES −0.95 / −0.86), because
protamine transcripts alone carry ~26% of the guarded ES UMIs at ceiling PDUI. Right: of the 16 genes
of the 27-gene curated literature panel measurable in both mice, only **4 shorten in both**
(Prm3, Ybx2, Ppp1cc, Nsun7), **5 lengthen in both** (Prm1, Tnp2, Odf1, Spata19, Smcp), 5 are
discordant and 2 uninformative.

**Cautions, binding on any use of this figure (they are part of the legend and must travel with it).**
- **Monotone LENGTHENING is also above its null** (z 3.1 / 2.1). The label shuffle destroys ordered
  structure of *either* sign, so "above null" alone is not directional evidence; the direction-specific
  claim rests on the **excess** of shortening over lengthening — strong in mouse 1 (binomial p 9.9e-9),
  **modest in mouse 2** (p 0.0125; 1 of 40 null draws — 20 fixed-gene-set + 20 free-guard — reached the
  same excess in absolute value). Never write "shortening is above null" without the lengthening column.
- **The per-gene MEDIANS are not monotone**: RS sits marginally above SPC in both mice
  (0.7867 → 0.7944 → 0.7320 and 0.8780 → 0.8863 → 0.8787), and in mouse 2 the SPC and ES medians are
  indistinguishable. The means *are* monotone (0.593/0.576/0.556 and 0.635/0.620/0.602). The ordering is
  carried by the monotone fraction (panel a), the paired per-gene tests and the per-cell residual
  (panel b) — never quote the medians as the gradient.
- **The statistic is the depth-weighted pseudobulk per-gene PDUI**, and that choice must be stated
  wherever these numbers are used. The cell-weighted arm (per-gene median of per-cell PDUI) is
  **unusable at this depth**: the guard leaves 214 / 197 genes and every median delta is 0.0 to 5 dp
  (REPORT_PROVISIONAL §2.3b).
- **Wilcoxon p-values are quoted only next to their shuffle-null column** and are therefore not on the
  figure. On the fixed TRUE gene set the best (smallest) of the 20 shuffle p-values reaches 8.2e-4
  (mouse 1 SPC→ES) and 2.4e-5 (mouse 2 SPC→ES), so a bare Wilcoxon p from this design is not
  interpretable on its own.
- **SPG is excluded from the claim** (manuscript/11): SPG → SPC *lengthens* under exactly the
  marker-argmax labelling used here, and rests on 64 / 68 cells. All 4-stage numbers are supplementary.
- **Gene-level and UMI-weighted summaries disagree in direction** (panel d, left) and both framings must
  appear wherever these data are used. This `wdi` is **not** the `wdi` of manuscript/11 (documented name
  collision): 11's is equal-weight per gene over all PAS, this one pools UMIs *across* genes on the
  proximal/distal endpoint pair only. The discrepancy is reproduced here, **not resolved**. Fig 5 must
  not use across-gene UMI weighting.
- **The literature-gene leg is dropped or shown with these numbers.** The quarantined "6/6 literature
  genes in the expected direction" was computed on a different index over a different gene set, so it is
  not a straight contradiction — but this panel is not a positive control, and Prm1 / Prm2 / Tnp1 / Tnp2
  sit at |Δ| < 0.003 because their PDUI is at ceiling in the stages where they are expressed.
- **Per-gene reproducibility and net shortening are different claims** and must not be conflated: the
  cross-mouse effect correlation is strong (ρ 0.66 / r 0.85, null 0.00 ± 0.03) and monotone shortening
  co-occurs in both mice in 167 genes (18.2%) against 9.5% expected under independence, while the
  aggregate sign-level bias is weak (338 both-shorten vs 278 both-lengthen, binomial p 0.017).
- **"Replication" here is two mice of one study** (GSE104556), one protocol, one chemistry: it controls
  animal-level and sampling noise, not study- or protocol-level artefacts. Labels are argmax GEX labels,
  not ground truth — cluster-level misassignment moves whole blocks of cells.
- **The switch test's null is conservative, not exact** (manuscript/14; Fisher discreteness puts a large
  fraction of p-values at exactly 1), so the TRUE hit counts are **not** to be read as "N true switches";
  only 5 label shuffles per mouse there (20 for the PDUI arm). Cross-mouse PAS matching is
  coordinate-based, so only the matched subset (~60% of tested PAS) can replicate — the denominator is
  on the figure.
- **The ≥ 50-UMI depth guard** keeps 1,553 / 1,194 of ~10,400 genes with a PAS pair; the surviving set is
  biased toward highly expressed genes. The shuffle null is computed on the *same* guarded set, so the
  comparison is internally valid, but the gene set is not a random sample of the transcriptome.
- **PROVISIONAL.** Code is the frozen tree `4efeb125`, not tool HEAD, and a **v2 re-run after PRs #96
  (IP-filter / strand) and #97** will regenerate every number here. Not yet final for the manuscript.
- Not covered by this run and therefore not on the figure: the `wul` vs `wdi` reconciliation, the
  feature-class decomposition, a 200-permutation null for the PDUI arm, and the 3′UTR-scoped PDUI arm —
  the last cannot be produced with this caller (`switch length --isoform-agg per_isoform` emits
  degenerate proximal == distal pairs; issue draft `github/issue11_per_isoform_degenerate.md`).

## Provenance

**Index entry (applied to `manuscript/05_figure_index.md` at the 2026-09-02 rename pass).**

> ## fig5_spermatogenesis — Fig 5, testis 3′-UTR control on the final caller (2026-08-21, verified FIXED)
>
> (a) the claim carrier: 31.5% / 30.2% of depth-guarded genes shorten monotonically across SPC > RS > ES
> against 20-shuffle label nulls of 17.5% / 21.0% (z 10.5 / 5.8), with monotone lengthening shown
> alongside (0.209 / 0.248, z 3.1 / 2.1) because the shuffle destroys ordered structure of either sign —
> the direction-specific evidence is the **excess** (489 vs 325, p 9.9e-9; 361 vs 296, p 0.0125);
> (b) the composition-controlled per-cell distal-usage residual falls at every step in both mice
> (Cliff's δ SPC vs ES 0.56 / 0.61, outside the whole 20-shuffle null range); (c) reliability payload —
> 11,263 / 6,049 / 9,469 PAS replicate across mice with the same direction (99.7–99.8% sign agreement,
> 3.2–5.1× over independence, **0 in all 15 null pairings**; arm B0 TRUE hits 50,415 / 52,043, 0 hits in
> 30 null BH families) and the per-gene effect correlates across mice at ρ 0.655 (n 917, null
> 0.001 ± 0.032); (d) boxed **negatives** — the across-gene UMI-weighted distal index *reverses* at
> RS → ES (Cliff −0.95 / −0.86; protamine ceiling-PDUI UMIs dominate ES) and the 16-gene literature panel
> does not reproduce (4 shorten / 5 lengthen / 5 discordant / 2 uninformative). **Caveats that travel with
> it:** lengthening is also above null, so quote the excess, not "above null"; mouse-2 direction excess is
> modest (p 0.0125, 1/40 shuffles); per-gene medians are **not** monotone and must never be quoted as the
> gradient; the statistic is depth-weighted pseudobulk PDUI (the cell-weighted arm is unusable at this
> depth); SPG excluded (it lengthens under these labels); two mice of one study with argmax GEX labels;
> switch null conservative, not exact; PROVISIONAL on frozen code `4efeb125`, v2 re-run after PRs #96/#97
> pending. Caption: `figures/fig5_spermatogenesis.caption.md`; write-up
> [18_spermatogenesis_final.md](18_spermatogenesis_final.md).

**`figures/README.md` table row (applied at the 2026-09-02 rename pass).**

> \| `fig5_spermatogenesis` \| Fig 5a–c: per-gene 3′UTR shortening along SPC → RS → ES in both testis mice — 31.5% / 30.2% of depth-guarded genes monotone vs 17.5% / 21.0% label-shuffle null (z 10.5 / 5.8); cross-mouse ρ 0.655; UMI-weighted index and literature panel reported as negatives \|

Sources: `manuscript/18_spermatogenesis_final.md` (verdict **FIXED**, 2026-08-21 — the figure states that
claim and no stronger one) and `results/stage3_spermatogenesis_final/` —
`REPORT_PROVISIONAL.md` §0–§6, `mouse{1,2}/summary/pdui_summary.json`,
`mouse{1,2}/summary/pdui_null_shuffles.tsv`, `mouse{1,2}/summary/pdui_per_cell.tsv`,
`summary/switch_cross_mouse_replication.tsv`, `summary/switch_null_summary.tsv`,
`summary/gene_cross_mouse_delta.tsv`, `summary/gene_replication_summary.json`,
`summary/gene_literature_panel.tsv`. Audit TSVs (every plotted value):
`results/figures/manuscript/fig5_spermatogenesis{,_null_shuffles,_percell_resid,_percell_values,_replication,_gene_scatter,_negatives,_reference_lines}.tsv`.
Palette: Okabe-Ito (#0072B2 #009E73 #D55E00 #CC79A7 #E69F00 #56B4E9 #999999); every series also carries a
marker or fill style, so nothing depends on hue alone. PNG 600 dpi, PDF vector with subsetted TrueType
(fonttype 42, no Type 3). The on-figure title, subtitle and cautions footer of the working-phase render
moved into this Legend at the 2026-09-02 surgery pass; the image keeps panel letters, short titles,
axis labels and data annotations only.

**Note for whoever quotes the residual test.** An earlier draft of `manuscript/18` gave the per-cell
residual Cliff's δ a "permutation p < 0.005"; that wording was corrected in the same commit as this
figure and `manuscript/18` now states the null-range form used here. The on-disk null for that
statistic is the **20-shuffle** column
`cell_cliffs_resid_SPC_ES`, whose empirical p floor is 1/21 ≈ 0.048; the observed δ (0.562 / 0.608) lies
outside the entire null range ([−0.114, +0.097] and [−0.066, +0.071]), which is what the figure states.
The Mann–Whitney p for the same comparison is 3.0e-39 / 2.6e-34, but it treats cells as independent.
"""

CAPTION_V2 = r"""# Fig 5 — `fig5_spermatogenesis` caption (generated by `scripts/manuscript_figures/fig5_spermatogenesis.py`; v2 record, verified SOUND 2026-09-02 — manuscript/18 v2 section)

## Legend

Figure 5 | **Per-gene 3′UTR usage shifts progressively proximal along spermatocyte → round → elongating spermatid, in both mice.** GSE104556 testis, both mice, PeakATail merged caller (frozen tree `9dfdefb3`, PR #96 IP-filter/strand + PR #97; the v1-code render on `4efeb125` stays reproducible via `SPERMATOGENESIS_VERSION=v1` and is retained in 18 as the labelled comparison); stage labels are marker-panel argmax over leiden clusters of the STARsolo **gene-expression** matrix, byte-identical to the v1 run, joined to the PAS matrix by exact cell-barcode match (1,294 / 1,364 labelled cells; 98.8% agreement with the independent PAS-clustering labels). PDUI = distal / (proximal + distal), higher = longer 3′UTR; proximal/distal = first/last PAS of the gene in transcription order, strand taken from `pasbed.bed` and re-asserted per row (the strand/counts re-validation of 13 §2, which the `per_gene` path **passes on the merged code**). The per-gene statistic is the **pseudobulk PDUI per stage** over genes with ≥ 50 UMI at the PAS pair in **every** stage (1,548 / 1,200 genes of ~10,500 with a pair). **31.4% and 30.5%** of guarded genes shorten monotonically across SPC > RS > ES (**486 / 1,548** and **366 / 1,200**) against label-shuffle nulls of **17.4% and 21.0%** (z = **10.5** and **6.4**, 20 shuffles, one label per cell, composition preserved), exceeding monotone lengthening (486 vs 318 and 366 vs 291 genes), and the composition-controlled per-cell distal-usage residual falls at every step (Cliff's δ SPC vs ES = **0.55** and **0.63**, outside the entire 20-shuffle null range in both mice). The claim is a **per-gene, equal-weight** statement: the per-gene medians are not monotone and the across-gene UMI-weighted distal index reverses at RS → ES (panel d).

**Panel a (the claim carrier)** — fraction of depth-guarded genes strictly monotone across the three
stages, per mouse and per direction, with the 20 label-shuffle draws as grey dots (fixed TRUE gene set;
bar = null mean) and the TRUE value as a filled marker (● mouse 1, ■ mouse 2): shortening 0.314 vs
0.174 ± 0.013 (z 10.5, 1.80×) and 0.305 vs 0.210 ± 0.015 (z 6.4, 1.45×); lengthening 0.205 vs
0.177 ± 0.010 (z 2.8, 1.16×) and 0.242 vs 0.204 ± 0.019 (z 2.1, 1.19×). Dotted line = 1/6, the chance
rate for a random ordering of three stages. The bracket carries the direction-specific evidence: the
**excess** of shortening over lengthening (486 vs 318, binomial p 3.4e-9; 366 vs 291, p 0.0039).
**Panel b** — the composition-controlled per-cell index: for each cell, the mean over genes of
(per-cell PDUI − that gene's 3-stage pseudobulk PDUI), which removes the gene-composition confound
(the median cell carries 292.5 / 326 / 119 genes over SPC / RS / ES in mouse 1, so raw per-cell
medians are not comparable across stages). Medians fall at every step in both mice
(+0.0018 → −0.0036 → −0.0114 and +0.0031 → −0.0014 → −0.0081); Cliff's δ SPC vs RS 0.26 / 0.32,
RS vs ES 0.41 / 0.46, SPC vs ES 0.55 / 0.63 against 20-shuffle null ranges of [−0.10, +0.11] and
[−0.08, +0.07]. Violin = all cells, thick bar = IQR, whisker = 5–95th percentile, open marker =
median; δ = Cliff's delta (earlier stage more distal). **Panel c** — reliability payload. Left: PAS called switching at q < 0.05 in **both**
mice with the same sign of Δ proportion (arm B0: `fisher --count-mode cells --marker-top-n 0`), per
stage pair — 11,219 / 6,070 / 9,480 of 36,840 / 35,668 / 35,382 coordinate-matched PAS (gene + strand,
≤ 100 bp, greedy 1:1), sign agreement 99.7–99.8%, 3.25–5.05× the independence expectation (hollow bars;
the 10-draw PAS-shuffle mean agrees with the analytic value to within 0.8%), and **0 replicated in all 15 null
pairings** (▼). TRUE switch hits over the three pairs: 50,354 / 52,111; label-shuffle null 3.2% / 3.2%
at p < 0.05 with **0 q < 0.05 hits in all 30 null BH families**. Right: the per-gene effect
Δ PDUI (ES − SPC) in mouse 1 vs mouse 2 over the 923 genes guarded in both — Spearman ρ = **0.641**
(Pearson r = 0.842) against a gene-correspondence null of 0.004 ± 0.033 (200 permutations, max |ρ|
0.088; strip below the scatter). 339 genes (36.7%) shorten in both mice, 273 (29.6%) lengthen in both,
253 are strictly discordant, 58 have Δ exactly 0 in one mouse. **Panel d (boxed — NEGATIVE findings,
not claim carriers)** — left: the **across-gene UMI-weighted** per-cell distal index rises at RS → ES
in both mice (0.587 → 0.572 → 0.627 and 0.640 → 0.636 → 0.672; Cliff δ RS vs ES −0.93 / −0.87), because
protamine transcripts alone carry 25.2% (mouse 1) / 14.5% (mouse 2) of the guarded ES UMIs at ceiling
PDUI. Right: of the 16 genes of the 27-gene curated literature panel measurable in both mice, only
**4 shorten in both** (Prm3, Ybx2, Ppp1cc, Nsun7), **5 lengthen in both** (Prm1, Tnp2, Odf1, Spata19,
Smcp), 5 are discordant and 2 uninformative — the identical classification to v1.

**Cautions, binding on any use of this figure (they are part of the legend and must travel with it).**
- **Monotone LENGTHENING is also above its null** (z 2.8 / 2.1). The label shuffle destroys ordered
  structure of *either* sign, so "above null" alone is not directional evidence; the direction-specific
  claim rests on the **excess** of shortening over lengthening — strong in mouse 1 (binomial p 3.4e-9),
  weaker in mouse 2 though clearer than on v1 (p 0.0039; 0 of 40 null draws — 20 fixed-gene-set +
  20 free-guard — reached the excess in absolute value, where v1 had 1 of 40). Never write "shortening
  is above null" without the lengthening column.
- **The per-gene MEDIANS are not monotone**: RS sits marginally above SPC in both mice
  (0.7953 → 0.8014 → 0.7364 and 0.8787 → 0.8874 → 0.8813), and in mouse 2 the ES median sits marginally
  above SPC. The means *are* monotone (0.595/0.578/0.558 and 0.634/0.619/0.602). The ordering is
  carried by the monotone fraction (panel a), the paired per-gene tests and the per-cell residual
  (panel b) — never quote the medians as the gradient.
- **The statistic is the depth-weighted pseudobulk per-gene PDUI**, and that choice must be stated
  wherever these numbers are used. The cell-weighted arm (per-gene median of per-cell PDUI) is
  **unusable at this depth**: the guard leaves 215 / 194 genes and every median delta is 0.0 to 5 dp
  (v2 REPORT_PROVISIONAL §2.3b; its SPC→RS rank-biserial even flips sign vs v1 at trivial magnitude,
  disclosed in 18's v2 section).
- **Wilcoxon p-values are quoted only next to their shuffle-null column** and are therefore not on the
  figure. On the fixed TRUE gene set the best (smallest) of the 20 shuffle p-values reaches 2.6e-4
  (mouse 1 SPC→ES) and 1.4e-4 (mouse 2 SPC→ES), so a bare Wilcoxon p from this design is not
  interpretable on its own.
- **SPG is excluded from the claim** (manuscript/11): SPG → SPC *lengthens* under exactly the
  marker-argmax labelling used here, and rests on 64 / 68 cells. All 4-stage numbers are supplementary
  (mouse 1's supplementary SPG→SPC Wilcoxon crosses 0.05 on v2 — 0.0509 vs v1's 0.0354 — disclosed in
  18's v2 section; it touches nothing in the claim).
- **Gene-level and UMI-weighted summaries disagree in direction** (panel d, left) and both framings must
  appear wherever these data are used. This `wdi` is **not** the `wdi` of manuscript/11 (documented name
  collision): 11's is equal-weight per gene over all PAS, this one pools UMIs *across* genes on the
  proximal/distal endpoint pair only. The discrepancy is reproduced on the merged code, **not resolved**.
  Fig 5 must not use across-gene UMI weighting.
- **The literature-gene leg is dropped or shown with these numbers.** The quarantined "6/6 literature
  genes in the expected direction" was computed on a different index over a different gene set, so it is
  not a straight contradiction — but this panel is not a positive control, and Prm1 / Prm2 / Tnp1 / Tnp2
  sit at |Δ| < 0.003 because their PDUI is at ceiling in the stages where they are expressed.
- **Per-gene reproducibility and net shortening are different claims** and must not be conflated: the
  cross-mouse effect correlation is strong (ρ 0.64 / r 0.84, null 0.00 ± 0.03) and monotone shortening
  co-occurs in both mice in 165 genes (17.9%) against 9.5% expected under independence, while the
  aggregate sign-level bias is weak (339 both-shorten vs 273 both-lengthen, binomial p 0.0085).
- **"Replication" here is two mice of one study** (GSE104556), one protocol, one chemistry: it controls
  animal-level and sampling noise, not study- or protocol-level artefacts. Labels are argmax GEX labels,
  not ground truth — cluster-level misassignment moves whole blocks of cells.
- **The switch test's null is conservative, not exact** (manuscript/14; Fisher discreteness puts a large
  fraction of p-values at exactly 1), so the TRUE hit counts are **not** to be read as "N true switches";
  only 5 label shuffles per mouse there (20 for the PDUI arm; pooled null p<0.05 3.169% / 3.156%).
  Cross-mouse PAS matching is coordinate-based, so only the matched subset (~60% of tested PAS) can
  replicate — the denominator is on the figure.
- **The ≥ 50-UMI depth guard** keeps 1,548 / 1,200 of ~10,500 genes with a PAS pair; the surviving set is
  biased toward highly expressed genes. The shuffle null is computed on the *same* guarded set, so the
  comparison is internally valid, but the gene set is not a random sample of the transcriptome.
- **This is the record.** Code is the frozen merged tree `9dfdefb3`; every headline was recounted to the
  digit by an adversarial verifier from the raw v2 outputs (18's v2 section, verdict SOUND, 2026-09-02).
  The only v1→v2 input change is 100% minus-strand pasbed churn (m1 −1,505/+1,248, m2 −1,302/+1,279,
  zero plus-strand changes) — the #96 footprint — and no headline moved beyond the small,
  churn-consistent tier.
- Not covered by this run and therefore not on the figure: the `wul` vs `wdi` reconciliation, the
  feature-class decomposition, a 200-permutation null for the PDUI arm, and the 3′UTR-scoped PDUI arm —
  the last still cannot be produced (issue #98 reproduces bit-for-bit on the merged code: `per_isoform`
  emits degenerate proximal == distal pairs on both mice, same offending-row counts as v1).

## Provenance

**Index entry (for `manuscript/05_figure_index.md`).**

> ## fig5_spermatogenesis — Fig 5, testis 3′-UTR control on the merged caller (v2 record, verified SOUND 2026-09-02)
>
> (a) the claim carrier: 31.4% / 30.5% of depth-guarded genes shorten monotonically across SPC > RS > ES
> against 20-shuffle label nulls of 17.4% / 21.0% (z 10.5 / 6.4), with monotone lengthening shown
> alongside (0.205 / 0.242, z 2.8 / 2.1) because the shuffle destroys ordered structure of either sign —
> the direction-specific evidence is the **excess** (486 vs 318, p 3.4e-9; 366 vs 291, p 0.0039, 0/40
> null draws); (b) the composition-controlled per-cell distal-usage residual falls at every step in both
> mice (Cliff's δ SPC vs ES 0.55 / 0.63, outside the whole 20-shuffle null range); (c) reliability
> payload — 11,219 / 6,070 / 9,480 PAS replicate across mice with the same direction (99.7–99.8% sign
> agreement, 3.25–5.05× over independence, **0 in all 15 null pairings**; arm B0 TRUE hits
> 50,354 / 52,111, 0 hits in 30 null BH families) and the per-gene effect correlates across mice at
> ρ 0.641 (n 923, null 0.004 ± 0.033); (d) boxed **negatives** — the across-gene UMI-weighted distal
> index *reverses* at RS → ES (Cliff −0.93 / −0.87; protamine ceiling-PDUI UMIs dominate ES: 25.2% /
> 14.5%) and the 16-gene literature panel does not reproduce (4 shorten / 5 lengthen / 5 discordant /
> 2 uninformative — identical to v1). **Caveats that travel with it:** lengthening is also above null,
> so quote the excess, not "above null"; per-gene medians are **not** monotone and must never be quoted
> as the gradient; the statistic is depth-weighted pseudobulk PDUI (the cell-weighted arm is unusable at
> this depth); SPG excluded (it lengthens under these labels); two mice of one study with argmax GEX
> labels; switch null conservative, not exact; issue #98 still blocks the 3′UTR-scoped arm. **v2 record
> on frozen merged code `9dfdefb3`, verified SOUND 2026-09-02** (18's v2 section; v1 render 4efeb125
> reproducible via `SPERMATOGENESIS_VERSION=v1`). Caption: `figures/fig5_spermatogenesis.caption.md`;
> write-up [18_spermatogenesis_final.md](18_spermatogenesis_final.md).

**`figures/README.md` table row.**

> \| `fig5_spermatogenesis` \| Fig 5a–c: per-gene 3′UTR shortening along SPC → RS → ES in both testis mice — 31.4% / 30.5% of depth-guarded genes monotone vs 17.4% / 21.0% label-shuffle null (z 10.5 / 6.4); cross-mouse ρ 0.641; UMI-weighted index and literature panel reported as negatives (v2 record, code 9dfdefb3) \|

Sources: `manuscript/18_spermatogenesis_final.md` **v2 section** (verdict **SOUND**, 2026-09-02 — the
figure states that claim and no stronger one; the v1 body is the labelled comparison) and
`results/stage3_spermatogenesis_v2/` — `REPORT_PROVISIONAL.md` §0–§8,
`mouse{1,2}/summary/pdui_summary.json`, `mouse{1,2}/summary/pdui_null_shuffles.tsv`,
`mouse{1,2}/summary/pdui_per_cell.tsv`, `summary/switch_cross_mouse_replication.tsv`,
`summary/switch_null_summary.tsv`, `summary/gene_cross_mouse_delta.tsv`,
`summary/gene_replication_summary.json`, `summary/gene_literature_panel.tsv`. Audit TSVs (every
plotted value):
`results/figures/manuscript/fig5_spermatogenesis{,_null_shuffles,_percell_resid,_percell_values,_replication,_gene_scatter,_negatives,_reference_lines}.tsv`.
Palette: Okabe-Ito (#0072B2 #009E73 #D55E00 #CC79A7 #E69F00 #56B4E9 #999999); every series also carries a
marker or fill style, so nothing depends on hue alone. PNG 600 dpi, PDF vector with subsetted TrueType
(fonttype 42, no Type 3). The on-figure title, subtitle and cautions footer of the working-phase render
moved into this Legend at the 2026-09-02 surgery pass; the image keeps panel letters, short titles,
axis labels and data annotations only.

**Note for whoever quotes the residual test.** The on-disk null for the per-cell residual Cliff's δ is
the **20-shuffle** column `cell_cliffs_resid_SPC_ES`, whose empirical p floor is 1/21 ≈ 0.048; the
observed δ (0.554 / 0.632) lies outside the entire null range ([−0.100, +0.108] and [−0.078, +0.073]),
which is what the figure states. The Mann–Whitney p for the same comparison is 3.3e-38 / 5.6e-37, but
it treats cells as independent.
"""

CAPTION_MD = CAPTION_V1 if VERSION == "v1" else CAPTION_V2
with open(f"{FIGDIR}/{NAME}.caption.md", "w") as fh:
    fh.write(CAPTION_MD)
print("wrote", f"{FIGDIR}/{NAME}.caption.md", f"({VERSION})")

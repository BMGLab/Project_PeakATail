#!/usr/bin/env python3
"""
laughney_switches.py -- manuscript Fig 6 ("laughney_switches"): reliable
cell-type APA switches in a tumour cohort (Laughney lung adenocarcinoma),
Stage-3 patient-level replication.

VERSION SWITCH -- env var LAUGHNEY_SWITCHES_VERSION picks which Stage-3 chain
is plotted.  Same registry pattern as scripts/manuscript_figures/final_benchmark.py.
  v1_code  (default; THE CURRENT RECORD)  PeakATail code 4efeb125 (pre-Stage-1d),
      tree /mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v2/.
      Source of truth: manuscript/20_stage3_replication.md (verifier verdict FIXED).
  v2_code  (LAUGHNEY_SWITCHES_VERSION=v2_code)  PeakATail code 9dfdefb
      (#93 + #96 IP-filter minus-strand fix + #97), tree
      /mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v3/.
      *** PENDING: that chain is still running as this script is written.  The
      paths are wired but nothing under stage3_laughney_v3 has been read.  When
      the chain lands and a verified write-up exists, re-run with the env var
      set -- no code edit -- and then fill in EXPECT/doc strings for v2_code. ***

SOURCE OF TRUTH -- every plotted value is READ or RECOMPUTED from the primary
replication tree (.../replication/primary_noMetBone/), the per-GSM switch
summaries (.../switch/<GSM>/summary.json), the verifier's null recount
(.../replication/verify/), and, for the genomic-context panel, the same Ensembl
GTF the cohort run used.  Nothing numeric is typed by hand except facts that
exist only in the verified write-up (code commit, the pre-registration, the
cohort label counts, the 33%-of-top-30 gene-name audit) -- those are cited to 20
on the figure.  Consistency asserts below check the headline counts against the
verified table in 20.

PANELS
  a  Replication funnel, real vs 10 patient-wise label-shuffle nulls on ONE log
     axis: tested (pair, PAS) hypotheses -> called BH q<0.05 in >=1 patient ->
     called in >=2 patients (any direction) -> replicated (>=2 patients, same
     direction, any opposite-direction patient vetoes) -> + effect floor
     |dprop| >= 0.1.  The null replicates nothing at any stage past the first.
  b  Per cell-type pair (top 16 of the 47 pairs tested in >=2 patients):
     replicated (pair, PAS) at K=2 with the K=3 subset overlaid; patients tested
     annotated on each row.
  c  Distribution of the consensus |delta proportion| of the 14,480 replicated
     switches, with the pre-registered 0.1 effect floor marked.
  d  Patient support: how many patients back each replicated switch.  The K=3
     sensitivity set is exactly the replication_count >= 3 tail (4,937).
  e  Honesty panel: genomic context of the 5,395 distinct replicated PAS
     (recomputed here from the GTF, exclusive partition) + the gene-name caveat.

CRITICAL (20 disclosure 1) -- NO ranked list of named top genes is drawn.  33%
of the top-30 gene-level rows name a spanning/readthrough model rather than the
gene whose 3' UTR holds the PAS (tool issue #99); the named list must not be
published until re-assignment.  The re-assignment-safe subset (PAS verified
here to lie in their OWN assigned gene's annotated 3' UTR) is written to an
audit TSV only, and is not plotted.

OUTPUTS
  manuscript/figures/laughney_switches.{png,pdf}      (PNG 300 dpi, PDF fonttype 42)
  manuscript/figures/laughney_switches.caption.md     (sidecar caption + index paragraph)
  results/figures/manuscript/laughney_switches_funnel.tsv          (panel a)
  results/figures/manuscript/laughney_switches_per_pair.tsv        (panel b, all 59 pairs)
  results/figures/manuscript/laughney_switches_effects.tsv         (panel c histogram)
  results/figures/manuscript/laughney_switches_support.tsv         (panel d)
  results/figures/manuscript/laughney_switches_context.tsv         (panel e, per-PAS)
  results/figures/manuscript/laughney_switches_context_summary.tsv (panel e, plotted values)
  results/figures/manuscript/laughney_switches_cohort.tsv          (cohort facts on the figure)
  results/figures/manuscript/laughney_switches_own3utr_genes.tsv   (AUDIT ONLY, not plotted)

Run:  export LC_ALL=C; python3 scripts/manuscript_figures/laughney_switches.py
      export LC_ALL=C; LAUGHNEY_SWITCHES_VERSION=v2_code python3 .../laughney_switches.py
"""
import collections
import json
import os
import textwrap
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "DejaVu Sans"
matplotlib.rcParams["font.size"] = 7.5
matplotlib.rcParams["axes.titlesize"] = 8
matplotlib.rcParams["axes.labelsize"] = 7.5
matplotlib.rcParams["legend.fontsize"] = 6.3

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
OUTDIR = WD / "results/figures/manuscript"
FIGDIR = WD / "manuscript/figures"
NAME = "laughney_switches"
OUTDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)

INK, MUTED, GRID = "#1B2429", "#5A6B73", "#D8E0E3"
# Okabe-Ito colourblind-safe palette (manuscript/figures/README.md convention).
BLUE, GREEN, VERM = "#0072B2", "#009E73", "#D55E00"
PINK, ORANGE, SKY, GREY = "#CC79A7", "#E69F00", "#56B4E9", "#999999"

FDR = 0.05
FLOOR = 0.10
K_PRIMARY, K_SENS = 2, 3
N_NULL_COMBOS = 10

# ---------------------------------------------------------------------------
# Which Stage-3 chain: v1_code (default, the current record) or v2_code (wired,
# pending).  Only the tree root and the hand-cited write-up facts differ.
# ---------------------------------------------------------------------------
VERSION = os.environ.get("LAUGHNEY_SWITCHES_VERSION", "v1_code").lower()
VERSIONS = {
    "v1_code": dict(
        tree=Path("/mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v2"),
        commit="4efeb125", commit_note="pre-Stage-1d",
        doc="20", doc_path="manuscript/20_stage3_replication.md",
        verdict="verifier verdict FIXED",
        record_note=("v1-code record; the Stage-3 v2 chain on code 9dfdefb "
                     "(stage3_laughney_v3) is running and will supersede these numbers "
                     "(re-run with LAUGHNEY_SWITCHES_VERSION=v2_code)"),
        # headline counts from the verified table in 20 -- asserted below
        expect=dict(n_patients=12, n_gsm=15, n_pairs=59, n_features=1977134,
                    n_replicated=14480, n_replicated_floor=13826,
                    n_replicated_K3=4937, n_replicated_K3_floor=4776,
                    n_gene_K2=9579, n_gene_K2_floor=9303,
                    n_distinct_pas=5395, n_distinct_genes=2687, n_pairs_with_hits=47,
                    universe=74954, sens_replicated=14541, sens_floor=13888,
                    sens_patients=13),
        # facts that exist only in the write-up (20) -- cited, never recomputed here
        cohort=dict(curated_cells=29063, confirmed_cells=18651, confirmed_frac=0.642,
                    n_gsm_cohort=17, n_patients_cohort=14, pas_space=505197,
                    top30_bad_gene_frac=0.33),
    ),
    "v2_code": dict(
        # PENDING -- paths wired, nothing under this tree has been read.
        tree=Path("/mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v3"),
        commit="9dfdefb", commit_note="#93 + #96 + #97",
        doc="(pending)", doc_path="(pending -- no verified Stage-3 v2 write-up yet)",
        verdict="NOT YET VERIFIED",
        record_note=("v2-code chain (9dfdefb); PROVISIONAL until a verified Stage-3 v2 "
                     "write-up exists -- do not quote"),
        expect=None,          # no verified table to check against yet
        cohort=None,          # fill in from the v2 write-up when it lands
    ),
}
assert VERSION in VERSIONS, f"LAUGHNEY_SWITCHES_VERSION must be one of {sorted(VERSIONS)}, got {VERSION!r}"
V = VERSIONS[VERSION]
TREE = V["tree"]
REPL = TREE / "replication" / "primary_noMetBone"      # pre-registered primary (MetBone excluded)
SENS = TREE / "replication"                            # sensitivity tree, one level up
SWITCH = TREE / "switch"                               # per-GSM switch summaries
VERIFY = TREE / "replication" / "verify"               # verifier's independent recount
GTF = Path("/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.99.gtf")

# Refuse to read a half-finished chain: replication_filter.py drops DONE.ok into a
# config directory only after that config is complete.  This is a stat, not a read.
CONFIGS = ("pas_K2", "pas_K3", "gene_K2")
if not REPL.exists():
    raise SystemExit(
        f"{VERSION}: {REPL} does not exist.\n"
        "  v2_code is pre-wired but PENDING -- the Stage-3 v2 chain must finish first.")
_missing = [c for c in CONFIGS if not (REPL / c / "DONE.ok").exists()]
if _missing:
    raise SystemExit(
        f"{VERSION}: {REPL} is INCOMPLETE -- no DONE.ok in {', '.join(_missing)}.\n"
        "  Refusing to plot a still-running chain.  Wait for the chain to finish.")
EXPECT = V["expect"]

def cite(v1_text, v2_text="(pending)"):
    """Facts that live only in a write-up: available for v1_code, pending for v2_code."""
    return v1_text if VERSION == "v1_code" else v2_text


# ---------------------------------------------------------------------------
# 1. summary.json of the pre-registered primary config (pas, K=2)
# ---------------------------------------------------------------------------
def load_summary(cfg):
    with open(REPL / cfg / "summary.json") as fh:
        return json.load(fh)


S_K2 = load_summary("pas_K2")
S_K3 = load_summary("pas_K3")
S_G2 = load_summary("gene_K2")

PARAMS = S_K2["params"]
assert PARAMS["level"] == "pas" and PARAMS["strategy"] == "fisher"
assert PARAMS["min_samples"] == K_PRIMARY and S_K3["params"]["min_samples"] == K_SENS
assert abs(PARAMS["fdr"] - FDR) < 1e-12 and abs(PARAMS["effect_floor"] - FLOOR) < 1e-12
assert PARAMS["unit_of_replication"] == "group" and PARAMS["discordant_policy"] == "exclude"

SAMPLE_GROUPS = PARAMS["sample_groups"]                 # GSM -> patient
GSMS = sorted(SAMPLE_GROUPS)
PATIENTS = sorted(set(SAMPLE_GROUPS.values()))
N_GSM, N_PATIENTS = len(GSMS), len(PATIENTS)
_libs_per_patient = collections.Counter(SAMPLE_GROUPS.values())
N_MULTILIB_PATIENTS = sum(1 for v in _libs_per_patient.values() if v > 1)
MAX_LIBS_PER_PATIENT = max(_libs_per_patient.values())
PAIRS = S_K2["pairs_seen"]
N_PAIRS = len(PAIRS)
N_FEATURES = S_K2["real"]["n_features"]
N_REPL = S_K2["real"]["n_replicated"]
N_REPL_FLOOR = S_K2["real"]["n_replicated_floor"]
N_REPL_K3 = S_K3["real"]["n_replicated"]
N_REPL_K3_FLOOR = S_K3["real"]["n_replicated_floor"]
EMP_P = S_K2["null"]["replicated"]["empirical_p_ge_observed"]
NULL_PER_COMBO = S_K2["null"]["replicated"]["null_per_combo"]
assert S_K2["null"]["n_combos"] == N_NULL_COMBOS and len(NULL_PER_COMBO) == N_NULL_COMBOS
assert max(NULL_PER_COMBO) == 0, NULL_PER_COMBO
assert max(S_K3["null"]["replicated"]["null_per_combo"]) == 0
assert max(S_G2["null"]["replicated"]["null_per_combo"]) == 0

if EXPECT:
    assert (N_PATIENTS, N_GSM, N_PAIRS) == (EXPECT["n_patients"], EXPECT["n_gsm"], EXPECT["n_pairs"])
    assert N_FEATURES == EXPECT["n_features"], N_FEATURES
    assert (N_REPL, N_REPL_FLOOR) == (EXPECT["n_replicated"], EXPECT["n_replicated_floor"])
    assert (N_REPL_K3, N_REPL_K3_FLOOR) == (EXPECT["n_replicated_K3"], EXPECT["n_replicated_K3_floor"])
    assert (S_G2["real"]["n_replicated"], S_G2["real"]["n_replicated_floor"]) == \
        (EXPECT["n_gene_K2"], EXPECT["n_gene_K2_floor"])
    assert abs(EMP_P - 1.0 / (N_NULL_COMBOS + 1)) < 1e-9, EMP_P
else:
    print(f"WARNING {VERSION}: no verified table to check against -- numbers are PROVISIONAL.")

# ---------------------------------------------------------------------------
# 2. panel a -- the funnel, recomputed row by row from all_features.tsv
#    columns used: n_groups_called, replication_count, direction_consistent,
#    passes_replication, passes_replication_floor.
# ---------------------------------------------------------------------------
FUNNEL_COLS = ["n_groups_called", "replication_count", "direction_consistent",
               "passes_replication", "passes_replication_floor"]


def funnel_counts(cfg, K):
    af = REPL / cfg / "all_features.tsv"
    n_rows = 0
    c_ge1 = c_geK = rc_geK = rep = repf = 0
    with open(af) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {k: hdr.index(k) for k in FUNNEL_COLS}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            n_rows += 1
            if int(f[ix["n_groups_called"]] or 0) >= 1:
                c_ge1 += 1
            if int(f[ix["n_groups_called"]] or 0) >= K:
                c_geK += 1
            if int(f[ix["replication_count"]] or 0) >= K:
                rc_geK += 1
            if f[ix["passes_replication"]] == "True":
                rep += 1
            if f[ix["passes_replication_floor"]] == "True":
                repf += 1
    return dict(tested=n_rows, called_ge1=c_ge1, called_geK=c_geK,
                same_dir_geK=rc_geK, replicated=rep, replicated_floor=repf)


FUN = funnel_counts("pas_K2", K_PRIMARY)
assert FUN["tested"] == N_FEATURES and FUN["replicated"] == N_REPL \
    and FUN["replicated_floor"] == N_REPL_FLOOR, FUN
# The 14,643 called in >=2 patients lose 163 features before the replicated set:
#   N_SPLIT_DIRECTION -- called in >=2 patients but no single direction reached K
#                        (e.g. one patient up, one patient down)
#   N_DISCORDANT_VETO -- reached K in one direction but an opposite-direction
#                        patient vetoed the call (discordant_policy="exclude")
N_SPLIT_DIRECTION = FUN["called_geK"] - FUN["same_dir_geK"]
N_DISCORDANT_VETO = FUN["same_dir_geK"] - FUN["replicated"]
assert N_SPLIT_DIRECTION + N_DISCORDANT_VETO == FUN["called_geK"] - FUN["replicated"]

# null side of the same funnel ------------------------------------------------
null_ctrl = pd.read_csv(REPL / "pas_K2" / "null_control.tsv", sep="\t")
assert len(null_ctrl) == N_NULL_COMBOS
assert null_ctrl.n_replicated.max() == 0 and null_ctrl.n_replicated_floor.max() == 0
NULL_TESTED_MEAN = float(null_ctrl.n_features.mean())

# the verifier's independent recount of every nominal q<FDR hit in every null
# run (one row per GSM x permutation x pair x PAS) -- restricted to the 15 GSMs
# of the pre-registered primary.
DISCLOSED_GSM = "GSM3516672-StageIB"        # 20: trips the per-GSM null rule via BH discreteness
qhits = pd.read_csv(VERIFY / "null_qhits_raw.tsv", sep="\t", header=None,
                    names=["gsm", "perm", "pair", "pas_id", "q", "effect"])
qhits = qhits[qhits.gsm.isin(GSMS)]
ncounts = pd.read_csv(VERIFY / "null_counts_per_gsm_perm.tsv", sep="\t", header=None,
                      names=["gsm", "perm", "n_tests", "n_qhits"])
ncounts = ncounts[ncounts.gsm.isin(GSMS)]
NULL_TESTS_TOTAL = int(ncounts.n_tests.sum())
NULL_QHITS_TOTAL = int(ncounts.n_qhits.sum())
assert NULL_QHITS_TOTAL == len(qhits), (NULL_QHITS_TOTAL, len(qhits))
_dg = ncounts[ncounts.gsm == DISCLOSED_GSM]
assert len(_dg) == N_NULL_COMBOS, DISCLOSED_GSM
DISCLOSED_NULL_TESTS = int(_dg.n_tests.sum())
DISCLOSED_NULL_QHITS = int(_dg.n_qhits.sum())
PERMS = sorted(ncounts.perm.unique())
assert len(PERMS) == N_NULL_COMBOS
# per null combination: distinct (pair, PAS) called in >=1 patient
_per_combo = qhits.drop_duplicates(["perm", "pair", "pas_id"]).groupby("perm").size()
NULL_CALLED_GE1 = [int(_per_combo.get(p, 0)) for p in PERMS]
NULL_CALLED_GE1_MEAN = float(np.mean(NULL_CALLED_GE1))
# nothing in any null was called in two patients, hence 0 at every later stage
_two_patient_null = qhits.groupby(["perm", "pair", "pas_id"]).gsm.nunique()
assert (_two_patient_null <= 1).all(), "a null (pair,PAS) was called in >1 library"

# real single-patient call totals, straight from the per-GSM switch summaries
per_sample_called = S_K2["per_sample_called"]
REAL_CALLS_TOTAL = int(sum(per_sample_called.values()))
sw = {}
for g in GSMS:
    with open(SWITCH / g / "summary.json") as fh:
        sw[g] = json.load(fh)
REAL_TESTS_TOTAL = int(sum(d["true"]["n_tests"] for d in sw.values()))
CELLS_TESTED = int(sum(d["n_cells_final"] for d in sw.values()))
CELLS_CONFIRMED_MATCHED = int(sum(d["n_cells_matched"] for d in sw.values()))
UNIVERSE = sw[GSMS[0]]["n_universe_cohort"]
assert all(d["n_universe_cohort"] == UNIVERSE for d in sw.values())
assert all(d["tier_filter"] == "tier1_ge2" for d in sw.values())
if EXPECT:
    assert UNIVERSE == EXPECT["universe"], UNIVERSE
    assert REAL_CALLS_TOTAL == sum(per_sample_called.values())

FUNNEL_STAGES = [
    ("tested (pair, PAS) hypotheses", FUN["tested"], NULL_TESTED_MEAN, "mean of 10 null combinations"),
    (f"called BH q<{FDR:g} in ≥1 patient", FUN["called_ge1"], NULL_CALLED_GE1_MEAN,
     f"mean of 10; per-combination {min(NULL_CALLED_GE1)}–{max(NULL_CALLED_GE1)}"),
    (f"called in ≥{K_PRIMARY} patients (any direction)", FUN["called_geK"], 0.0, "0 in all 10"),
    (f"replicated: ≥{K_PRIMARY} patients, same direction,\nno opposite-direction patient",
     FUN["replicated"], 0.0, "0 in all 10"),
    (f"+ effect floor |Δproportion| ≥ {FLOOR:g}", FUN["replicated_floor"], 0.0, "0 in all 10"),
]

# ---------------------------------------------------------------------------
# 3. panel b -- per cell-type pair, K=2 with the K=3 subset
# ---------------------------------------------------------------------------
def per_pair(cfg):
    df = pd.read_csv(REPL / f"per_pair_{cfg}.tsv", sep="\t")
    return df[df.pair != "ALL"].copy()


pp2, pp3 = per_pair("pas_K2"), per_pair("pas_K3")
assert len(pp2) == N_PAIRS and len(pp3) == N_PAIRS
pairs = pp2.merge(pp3[["pair", "replicated_q", "replicated_q_floor"]], on="pair",
                  suffixes=("_K2", "_K3"))
assert pairs.replicated_q_K2.sum() == N_REPL and pairs.replicated_q_K3.sum() == N_REPL_K3
pairs["multi_patient"] = pairs.n_patients_tested >= K_PRIMARY
MULTI = pairs[pairs.multi_patient].sort_values("replicated_q_K2", ascending=False)
N_MULTI_PAIRS = len(MULTI)
N_PAIRS_WITH_HITS = int((pairs.replicated_q_K2 > 0).sum())
N_SINGLE_PATIENT_PAIRS = N_PAIRS - N_MULTI_PAIRS
if EXPECT:
    assert N_PAIRS_WITH_HITS == EXPECT["n_pairs_with_hits"], N_PAIRS_WITH_HITS
N_SHOW = 16
TOPB = MULTI.head(N_SHOW).iloc[::-1]                    # bottom-to-top for barh
REST = MULTI.iloc[N_SHOW:]
REST_N, REST_SUM = len(REST), int(REST.replicated_q_K2.sum())
REST_MIN, REST_MAX = int(REST.replicated_q_K2.min()), int(REST.replicated_q_K2.max())
pairs["plotted_panel_b"] = pairs.pair.isin(MULTI.head(N_SHOW).pair)


def pretty_pair(p):
    a, b = p.split("_vs_")
    return f"{a.replace('_', ' ')} – {b.replace('_', ' ')}"


# ---------------------------------------------------------------------------
# 4. panels c/d -- effect size and patient support of the replicated set
# ---------------------------------------------------------------------------
rep = pd.read_csv(REPL / "pas_K2" / "replicated.tsv", sep="\t",
                  usecols=["c1", "c2", "pas_id", "gene_id", "chrom", "start", "end", "strand",
                           "n_groups_tested", "replication_count", "consensus_direction",
                           "min_q_called", "mean_effect_consensus", "passes_replication_floor"])
assert len(rep) == N_REPL
rep["abs_effect"] = rep.mean_effect_consensus.abs()
BINW = 0.02
bins = np.arange(0.0, 1.0 + BINW, BINW)
hist_all, _ = np.histogram(rep.abs_effect, bins=bins)
MEDIAN_EFFECT = float(rep.abs_effect.median())
MEAN_EFFECT = float(rep.abs_effect.mean())
N_CONSENSUS_BELOW_FLOOR = int((rep.abs_effect < FLOOR).sum())
# the pre-registered floor is applied PER PATIENT and the patients are re-counted,
# so it removes more than the consensus-below-floor rows -- state both.
N_REMOVED_BY_FLOOR = N_REPL - N_REPL_FLOOR

support = rep.replication_count.value_counts().sort_index()
SUPPORT_K = support.index.to_numpy()
SUPPORT_N = support.to_numpy()
# the K=3 sensitivity set is exactly the replication_count >= 3 tail
assert int(support[support.index >= K_SENS].sum()) == N_REPL_K3, support
N_DISTINCT_PAS = rep.pas_id.nunique()
N_DISTINCT_GENES = rep.gene_id.nunique()
N_PAIRS_REP = rep.groupby(["c1", "c2"]).ngroups
if EXPECT:
    assert (N_DISTINCT_PAS, N_DISTINCT_GENES, N_PAIRS_REP) == \
        (EXPECT["n_distinct_pas"], EXPECT["n_distinct_genes"], EXPECT["n_pairs_with_hits"])

# ---------------------------------------------------------------------------
# 5. panel e -- genomic context of the distinct replicated PAS, recomputed here
#    from the same Ensembl GTF the cohort run used (gtf_cache/manifest.json).
#    One streaming pass; only gene / exon / three_prime_utr features that
#    overlap one of the PAS intervals are kept.  Strand must match.
# ---------------------------------------------------------------------------
gtf_manifest = json.loads((TREE / "cohort_run/run/gtf_cache/manifest.json").read_text())
assert Path(gtf_manifest["gtf_path"]) == GTF, (gtf_manifest["gtf_path"], GTF)

pas_int = (rep.drop_duplicates("pas_id")[["pas_id", "chrom", "start", "end", "strand", "gene_id"]]
           .set_index("pas_id"))
assert len(pas_int) == N_DISTINCT_PAS

BUCKET = 1 << 17
buckets = collections.defaultdict(list)
for pid, r in pas_int.iterrows():
    for b in range(int(r.start) // BUCKET, int(r.end) // BUCKET + 1):
        buckets[(r.chrom, b)].append(pid)

gene_bodies = collections.defaultdict(set)     # pas_id -> {gene_id} same-strand gene bodies
utr3_genes = collections.defaultdict(set)      # pas_id -> {gene_id} same-strand 3' UTRs
exonic = set()                                 # pas_id overlapping any same-strand exon
n_gtf = 0
with open(GTF) as fh:
    for line in fh:
        if line[0] == "#":
            continue
        n_gtf += 1
        t = line.split("\t", 8)
        feat = t[2]
        if feat not in ("gene", "exon", "three_prime_utr"):
            continue
        chrom, s, e, strand = t[0], int(t[3]) - 1, int(t[4]), t[6]
        hit = set()
        for b in range(s // BUCKET, e // BUCKET + 1):
            hit.update(buckets.get((chrom, b), ()))
        if not hit:
            continue
        a = t[8]
        i = a.find('gene_id "')
        gid = a[i + 9:a.find('"', i + 9)] if i >= 0 else ""
        for pid in hit:
            r = pas_int.loc[pid]
            if r.strand != strand or int(r.start) >= e or int(r.end) <= s:
                continue
            if feat == "gene":
                gene_bodies[pid].add(gid)
            elif feat == "exon":
                exonic.add(pid)
            else:
                utr3_genes[pid].add(gid)
print(f"GTF scan: {n_gtf:,} feature lines, {N_DISTINCT_PAS:,} PAS classified")


def classify(pid):
    own = pas_int.loc[pid, "gene_id"]
    u = utr3_genes.get(pid, set())
    if own in u:
        return "own-gene 3' UTR"
    if u:
        return "other-gene 3' UTR"
    if pid in exonic:
        return "exonic, not 3' UTR"
    if gene_bodies.get(pid):
        return "intronic (same-strand gene body)"
    return "outside any same-strand gene"


CONTEXT_ORDER = ["own-gene 3' UTR", "other-gene 3' UTR", "exonic, not 3' UTR",
                 "intronic (same-strand gene body)", "outside any same-strand gene"]
CONTEXT_COLOR = dict(zip(CONTEXT_ORDER, [GREEN, PINK, SKY, ORANGE, GREY]))
ctx = pd.DataFrame(dict(
    pas_id=pas_int.index,
    chrom=pas_int.chrom.values, start=pas_int.start.values, end=pas_int.end.values,
    strand=pas_int.strand.values, assigned_gene_id=pas_int.gene_id.values))
ctx["context"] = [classify(p) for p in ctx.pas_id]
ctx["n_samestrand_gene_bodies"] = [len(gene_bodies.get(p, ())) for p in ctx.pas_id]
ctx["utr3_gene_ids"] = [",".join(sorted(utr3_genes.get(p, ()))) for p in ctx.pas_id]
ctx["in_own_gene_body"] = [pas_int.loc[p, "gene_id"] in gene_bodies.get(p, ()) for p in ctx.pas_id]
ctx["exonic_any_gene"] = ctx.pas_id.isin(exonic)

CTX_N = ctx.context.value_counts().reindex(CONTEXT_ORDER).fillna(0).astype(int)
assert int(CTX_N.sum()) == N_DISTINCT_PAS
C_OWN3, C_OTH3, C_EX, C_INTRON, C_OUT = CONTEXT_ORDER
N_OWN3 = int(CTX_N[C_OWN3])
N_OTH3 = int(CTX_N[C_OTH3])
N_EXNOT3 = int(CTX_N[C_EX])
N_INTRON = int(CTX_N[C_INTRON])
N_OUTSIDE = int(CTX_N[C_OUT])
N_ANY3 = N_OWN3 + N_OTH3
FRAC_ANY_3UTR = N_ANY3 / N_DISTINCT_PAS
FRAC_EXONIC = float(ctx.exonic_any_gene.mean())
FRAC_GENEBODY = float((ctx.n_samestrand_gene_bodies > 0).mean())
FRAC_INTRONIC = N_INTRON / N_DISTINCT_PAS
FRAC_OUTSIDE = N_OUTSIDE / N_DISTINCT_PAS
FRAC_OWN_3UTR = N_OWN3 / N_DISTINCT_PAS
N_MULTI_GENE = int((ctx.n_samestrand_gene_bodies >= 2).sum())
FRAC_MULTI_GENE = N_MULTI_GENE / N_DISTINCT_PAS
N_3UTR_MISASSIGNED = N_OTH3
FRAC_3UTR_MISASSIGNED = N_3UTR_MISASSIGNED / N_ANY3

# Regression check (v1_code only; the recomputed values are the ones drawn).
# NOTE ON PROVENANCE: 20 disclosure 1 quotes only two of these -- 62.5% ("own-gene
# 3' UTRs") and 96.8% (same-strand gene bodies).  The exonic / intronic / outside
# splits appear in NO write-up: they are recomputed here from the GTF, and the
# figure and caption say so.  The three unsourced targets below are this script's
# own first-run values, kept only so a silent drift would fail loudly.
if VERSION == "v1_code":
    for got, doc, what, src in (
            (FRAC_ANY_3UTR, 0.625, "3' UTR", "20 disclosure 1"),
            (FRAC_GENEBODY, 0.968, "gene body", "20 disclosure 1"),
            (FRAC_EXONIC, 0.764, "exonic", "recomputed here, not in any write-up"),
            (FRAC_INTRONIC, 0.20, "intronic", "recomputed here, not in any write-up"),
            (FRAC_OUTSIDE, 0.032, "outside", "recomputed here, not in any write-up")):
        assert abs(got - doc) < 0.01, f"{what}: recomputed {got:.3f} vs {doc:.3f} ({src})"

# AUDIT ONLY -- the re-assignment-safe subset (verified own-gene 3' UTR).
# NOT plotted and NOT a ranked list: 20 disclosure 1 forbids publishing named
# top-switch genes until the PAS->gene assignment is fixed (tool issue #99).
gene_names = pd.read_csv(REPL.parent / "gene_names.GRCh38.99.tsv", sep="\t", header=None,
                         names=["gene_id", "gene_name", "gene_biotype"])
safe = ctx[ctx.context == C_OWN3][["pas_id", "chrom", "start", "end", "strand",
                                              "assigned_gene_id", "n_samestrand_gene_bodies"]]
safe = safe.merge(gene_names, left_on="assigned_gene_id", right_on="gene_id", how="left") \
           .drop(columns=["gene_id"])
safe = safe.merge(rep.groupby("pas_id").agg(n_pairs_replicated=("c1", "size"),
                                            max_patient_support=("replication_count", "max"),
                                            max_abs_effect=("abs_effect", "max")).reset_index(),
                  on="pas_id", how="left")
safe = safe.sort_values(["chrom", "start"])            # genomic order, NOT ranked by effect
safe["disclosure"] = "own-gene 3'UTR verified vs GTF; audit only, not a published ranked list (20 disclosure 1, issue #99)"

# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(8.4, 11.6))
FW = 604.8  # figure width in points, for text-width budgeting


def style(ax, grid_axis="both"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED, length=2.5, labelsize=6.6)
    ax.grid(True, axis=grid_axis, color=GRID, lw=0.5, alpha=0.7)
    ax.set_axisbelow(True)


def panel_title(ax, letter, text, pad=6):
    ax.set_title(f"{letter}   {text}", loc="left", fontweight="bold", pad=pad)


# ---- a: funnel ------------------------------------------------------------
axA = fig.add_axes([0.275, 0.800, 0.700, 0.122])
XMIN, XMAX = 0.55, 1.3e7
nstage = len(FUNNEL_STAGES)
ys = np.arange(nstage)[::-1].astype(float)
bh = 0.34
for y, (lab, real, null, note) in zip(ys, FUNNEL_STAGES):
    yr, yn = y + bh / 2 + 0.015, y - bh / 2 - 0.015
    axA.barh(yr, real - XMIN, left=XMIN, height=bh, color=BLUE, edgecolor="none", zorder=3)
    axA.text(real * 1.30, yr, f"{real:,.0f}", va="center", ha="left",
             fontsize=6.6, color=INK, zorder=5)
    if null >= 1.0:
        axA.barh(yn, null - XMIN, left=XMIN, height=bh, color=GREY, edgecolor="none", zorder=3)
        axA.text(null * 1.30, yn, f"{null:,.0f}" if null >= 10 else f"{null:.1f}",
                 va="center", ha="left", fontsize=6.2, color=MUTED, zorder=5)
    else:
        axA.plot([XMIN * 1.16], [yn], marker="o", ms=3.0, mfc="white", mec=GREY,
                 mew=0.9, zorder=4)
        axA.text(XMIN * 1.75, yn, "0", va="center", ha="left", fontsize=6.2,
                 color=MUTED, zorder=5)
axA.set_yticks(ys)
axA.set_yticklabels([s[0] for s in FUNNEL_STAGES], fontsize=6.5, color=INK)
axA.set_xscale("log")
axA.set_xlim(XMIN, XMAX)
axA.set_ylim(-0.62, nstage - 0.38)
axA.set_xticks([1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7])
axA.set_xticklabels(["1", "10", "100", "1k", "10k", "100k", "1M", "10M"])
axA.set_xlabel("number of (pair, PAS) hypotheses  (log scale)")
style(axA, grid_axis="x")
axA.legend([Patch(facecolor=BLUE), Patch(facecolor=GREY)],
           ["real labels", f"label-shuffle null (mean of {N_NULL_COMBOS})"],
           loc="lower right", bbox_to_anchor=(1.006, -0.045), frameon=False,
           handlelength=1.1, handleheight=0.85, borderpad=0.2, labelspacing=0.3)
panel_title(axA, "a", "Replication funnel, real vs patient-wise label-shuffle null (one axis)")
fig.text(0.275, 0.772, textwrap.fill(
         f"Null = {N_NULL_COMBOS} patient-wise label-shuffle combinations through the identical pipeline; "
         f"bars are the mean. Single-patient calls per null combination range "
         f"{min(NULL_CALLED_GE1)}–{max(NULL_CALLED_GE1)}; nothing was called in two patients in any "
         f"combination, so every later stage is exactly 0 in all {N_NULL_COMBOS}.", 140),
         fontsize=5.9, color=MUTED, va="top", ha="left", linespacing=1.36)

# ---- b: per cell-type pair ------------------------------------------------
axB = fig.add_axes([0.275, 0.505, 0.700, 0.215])
yb = np.arange(len(TOPB))
axB.barh(yb, TOPB.replicated_q_K2, height=0.66, color=SKY, edgecolor="none", zorder=3,
         label=f"K={K_PRIMARY} (pre-registered): ≥2 patients, same direction")
axB.barh(yb, TOPB.replicated_q_K3, height=0.34, color=BLUE, edgecolor="none", zorder=4,
         label=f"K={K_SENS} (sensitivity): ≥3 patients")
for y, r in zip(yb, TOPB.itertuples()):
    axB.text(r.replicated_q_K2 + N_REPL * 0.004, y, f"{int(r.replicated_q_K2):,}",
             va="center", ha="left", fontsize=6.0, color=INK, zorder=5)
axB.set_yticks(yb)
axB.set_yticklabels([f"{pretty_pair(r.pair)}  ({int(r.n_patients_tested)} pt)"
                     for r in TOPB.itertuples()], fontsize=6.4, color=INK)
axB.set_xlim(0, float(TOPB.replicated_q_K2.max()) * 1.16)
axB.set_ylim(-0.75, len(TOPB) - 0.35)
axB.set_xlabel("replicated (pair, PAS) switches")
style(axB, grid_axis="x")
axB.legend(loc="lower right", frameon=False, handlelength=1.1, handleheight=0.85,
           borderpad=0.2, labelspacing=0.3)
panel_title(axB, "b", f"Top {N_SHOW} of the {N_MULTI_PAIRS} cell-type pairs tested in ≥{K_PRIMARY} "
                      f"patients   ('n pt' = patients tested)")
_note_b = (f"The other {REST_N} multi-patient pairs hold {REST_SUM:,} more ({REST_MIN}–{REST_MAX} each); "
           f"{N_SINGLE_PATIENT_PAIRS} of the {N_PAIRS} pairs were tested in one patient only and cannot replicate.")
assert len(_note_b) <= 140, f"panel-b note is {len(_note_b)} chars; >140 overflows the panel width"
fig.text(0.975, 0.4735, _note_b, fontsize=5.9, color=MUTED, va="top", ha="right")

# ---- c: effect sizes ------------------------------------------------------
axC = fig.add_axes([0.085, 0.345, 0.385, 0.094])
centers = bins[:-1] + BINW / 2
cols = [VERM if c < FLOOR else BLUE for c in centers]
axC.bar(centers, hist_all, width=BINW * 0.92, color=cols, edgecolor="none", zorder=3)
CMAX = float(hist_all.max())
axC.set_ylim(0, CMAX * 1.45)
axC.axvline(FLOOR, color=INK, lw=0.8, ls=(0, (4, 2)), zorder=4)
axC.text(FLOOR + 0.018, CMAX * 1.30, f"pre-registered effect floor |Δprop| ≥ {FLOOR:g}",
         fontsize=5.9, color=INK, va="center", ha="left", zorder=6)
axC.text(0.985, 0.985, f"median {MEDIAN_EFFECT:.2f}   mean {MEAN_EFFECT:.2f}",
         transform=axC.transAxes, fontsize=6.0, color=INK, va="top", ha="right")
axC.set_xlim(0, 1.0)
axC.set_xlabel("consensus |Δ proportion| of the replicated switch")
axC.set_ylabel("replicated (pair, PAS)")
style(axC, grid_axis="y")
panel_title(axC, "c", f"Effect size of the {N_REPL:,} replicated switches")

# ---- d: patient support ---------------------------------------------------
axD = fig.add_axes([0.590, 0.345, 0.385, 0.094])
dcols = [SKY if k < K_SENS else BLUE for k in SUPPORT_K]
axD.bar(SUPPORT_K, SUPPORT_N, width=0.72, color=dcols, edgecolor="none", zorder=3)
for k, n in zip(SUPPORT_K, SUPPORT_N):
    axD.text(k, n + SUPPORT_N.max() * 0.025, f"{n:,}", ha="center", va="bottom",
             fontsize=5.9, color=INK, zorder=5)
axD.set_xticks(SUPPORT_K)
axD.set_xlim(SUPPORT_K.min() - 0.75, SUPPORT_K.max() + 0.75)
axD.set_ylim(0, SUPPORT_N.max() * 1.34)
axD.set_xlabel("patients supporting the switch (same direction)")
axD.set_ylabel("replicated (pair, PAS)")
style(axD, grid_axis="y")
axD.legend([Patch(facecolor=BLUE)],
           [f"≥{K_SENS} patients = the K={K_SENS} sensitivity set ({N_REPL_K3:,})"],
           loc="upper right", frameon=False, handlelength=1.1, handleheight=0.85,
           borderpad=0.2)
panel_title(axD, "d", "Patient support of the replicated set")

# ---- e: honesty panel -----------------------------------------------------
fig.text(0.085, 0.3030,
         f"e   Genomic context of the {N_DISTINCT_PAS:,} replicated PAS "
         f"(recomputed here from the Ensembl GRCh38.99 GTF)",
         fontsize=8, fontweight="bold", color=INK, va="top", ha="left")
axE = fig.add_axes([0.085, 0.258, 0.890, 0.026])
left = 0.0
for cat in CONTEXT_ORDER:
    frac = CTX_N[cat] / N_DISTINCT_PAS
    axE.barh(0, frac, left=left, height=1.0, color=CONTEXT_COLOR[cat],
             edgecolor="white", linewidth=0.9, zorder=3)
    if frac > 0.055:
        axE.text(left + frac / 2, 0, f"{frac*100:.1f}%", ha="center", va="center",
                 fontsize=6.4, color="white", fontweight="bold", zorder=5)
    left += frac
axE.set_xlim(0, 1.0)
axE.set_ylim(-0.6, 0.6)
axE.set_yticks([])
axE.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
axE.set_xticklabels(["0", "25%", "50%", "75%", "100%"])
for s in ("top", "right", "left"):
    axE.spines[s].set_visible(False)
axE.spines["bottom"].set_color(GRID)
axE.tick_params(colors=MUTED, length=2.5, labelsize=6.6)
fig.legend([Patch(facecolor=CONTEXT_COLOR[c]) for c in CONTEXT_ORDER],
           [f"{c} — {int(CTX_N[c]):,} ({CTX_N[c]/N_DISTINCT_PAS*100:.1f}%)" for c in CONTEXT_ORDER],
           loc="upper left", bbox_to_anchor=(0.085, 0.2410), frameon=False, ncol=3,
           handlelength=1.1, handleheight=0.85, borderpad=0.2, columnspacing=1.4,
           labelspacing=0.32, fontsize=6.2)

honesty = (
    "Exclusive, strand-matched partition of the distinct replicated PAS. "
    f"Roll-ups: {FRAC_GENEBODY*100:.1f}% of the replicated PAS sit inside a same-strand gene body, "
    f"{FRAC_EXONIC*100:.1f}% are exonic, {FRAC_ANY_3UTR*100:.1f}% are in some gene's 3' UTR, "
    f"{FRAC_INTRONIC*100:.1f}% are intronic and {FRAC_OUTSIDE*100:.1f}% fall outside any same-strand gene "
    + cite("(of these, only 62.5% and 96.8% appear in 20 disclosure 1, both reproduced here to "
           "within 0.1 pp; the rest are recomputed for this figure). ",
           "(no verified v2 write-up to compare against yet). ")
    + f"BUT only {FRAC_OWN_3UTR*100:.1f}% lie in the 3' UTR of the gene the caller ASSIGNED them to: "
    f"{N_3UTR_MISASSIGNED:,} of the {N_ANY3:,} 3'-UTR PAS ({FRAC_3UTR_MISASSIGNED*100:.1f}%) lie in a "
    f"DIFFERENT gene's 3' UTR, and {N_MULTI_GENE:,} ({FRAC_MULTI_GENE*100:.1f}%) lie inside ≥2 overlapping "
    "same-strand gene models. NO ranked list of named top-switch genes is shown here: "
    + cite("33% of the top-30 gene-level rows in 20 name a spanning/readthrough model rather than the gene "
           "whose 3' UTR holds the PAS, and seven of them have zero assigned PAS cohort-wide (tool issue "
           "#99). ", "the gene-name audit must be repeated on the v2 chain. ")
    + "Replication statistics are unaffected — they are computed on PAS identifiers, never on gene names."
)
_honesty_wrapped = textwrap.fill(honesty, 160)
_n_honesty_lines = _honesty_wrapped.count("\n") + 1
assert _n_honesty_lines <= 7, f"honesty block is {_n_honesty_lines} lines; >7 collides with the footer"
fig.text(0.085, 0.2105, _honesty_wrapped, fontsize=5.9, color=INK, va="top", ha="left",
         linespacing=1.36)

# ---- title ----------------------------------------------------------------
fig.text(0.085, 0.9920,
         "Fig 6 | Reliable cell-type APA switches in a tumour cohort",
         fontsize=11.5, fontweight="bold", color=INK, va="top", ha="left")
fig.text(0.085, 0.9725, textwrap.fill(
         f"Laughney lung adenocarcinoma, {N_PATIENTS} patients ({N_GSM} scRNA-seq libraries): "
         f"{N_REPL:,} of {N_FEATURES:,} tested (pair, PAS) hypotheses replicate in ≥{K_PRIMARY} "
         f"patients in the same direction — and none replicate in any of the {N_NULL_COMBOS} "
         f"patient-wise label-shuffle nulls.", 128),
         fontsize=7.6, color=MUTED, va="top", ha="left", linespacing=1.4)

# ---- footer ---------------------------------------------------------------
footer = (
    f"Design: one cohort run so all {N_GSM} libraries share a single PAS identifier space; "
    f"pre-registered precision-first universe (IP-filtered ∧ tier-1 ∧ cohort clip-molecule sum ≥2) "
    f"= {UNIVERSE:,} PAS; "
    + cite(f"{V['cohort']['curated_cells']:,} curated cells → {V['cohort']['confirmed_cells']:,} "
           f"confirmed ({V['cohort']['confirmed_frac']*100:.1f}%) by the pre-registered label policy (20), "
           if V["cohort"] else "", "")
    + f"of which {CELLS_TESTED:,} entered the tests in the {N_PATIENTS} primary patients; "
    f"{N_PAIRS} cell-type pairs; test = Fisher on cells with no marker pre-selection (14); unit of "
    f"replication = patient, so the {N_MULTILIB_PATIENTS} patients with "
    f"{MAX_LIBS_PER_PATIENT} libraries each count once. "
    f"RECORD: PeakATail code {V['commit']} ({V['commit_note']}) — {V['record_note']}; "
    f"source of truth {V['doc_path']} ({V['verdict']}). "
    f"NULL: not one feature replicates in ANY of the {N_NULL_COMBOS} patient-wise "
    f"label-shuffle combinations run through the identical pipeline: 0 replicated in every one of the 10, over a "
    f"comparable {NULL_TESTED_MEAN:,.0f} (pair, PAS) hypotheses tested per combination "
    f"({NULL_TESTS_TOTAL/1e6:.1f}M null tests produced {NULL_QHITS_TOTAL} nominal q<{FDR:g} calls, none "
    f"co-occurring in two patients). Ten permutations resolve only to an empirical "
    f"p ≤ {EMP_P:.3f}; report this as 'none in {N_NULL_COMBOS} nulls', NEVER as an FDR estimate. "
    f"EXCLUSION: GSM3516664-MetBone is excluded by pre-registration (13 addendum item 8), but its stated "
    f"premise was WRONG — the 0.015% clip rate came from the caller sampling only the first 200,000 CB "
    f"reads of a coordinate-sorted BAM (head of chr1); MetBone in fact has 306k clip molecules "
    f"(tool issue #99). Never repeat the 'no clip evidence' justification. "
    + cite(f"The {EXPECT['sens_patients']}-patient sensitivity run that includes it gives "
           f"{EXPECT['sens_replicated']:,} / {EXPECT['sens_floor']:,}, a "
           f"{abs(EXPECT['sens_replicated']-N_REPL)/EXPECT['sens_replicated']*100:.1f}% difference confined "
           f"to one pair. " if EXPECT else "", "")
    + f"GENE NAMES: no named top-gene list may be published until PAS→gene re-assignment (panel e; "
    f"20 disclosure 1, tool issue #99) — the replication statistics themselves are computed on PAS ids. "
    f"CAVEATS: one cohort, one chemistry, one caller — no orthogonal 3'-end assay confirms these sites; "
    f"pairs tested in more patients replicate more, so panel b tracks cohort composition as much as "
    f"biology; one library ({DISCLOSED_GSM.split(chr(45))[0]}) trips the per-GSM null rule only through BH "
    f"discreteness ({DISCLOSED_NULL_QHITS} hits in {DISCLOSED_NULL_TESTS/1e6:.1f}M null tests) and is kept "
    f"and disclosed (20). "
    f"Panel c plots the consensus effect, but the pre-registered floor is applied PER PATIENT and the "
    f"patients re-counted, so it removes {N_REMOVED_BY_FLOOR} switches ({N_REPL:,} → {N_REPL_FLOOR:,}) — "
    f"more than the {N_CONSENSUS_BELOW_FLOOR} whose consensus alone falls below {FLOOR:g}. "
    f"Every plotted value: results/figures/manuscript/{NAME}*.tsv."
)
_footer_wrapped = textwrap.fill(footer, 172)
_n_footer_lines = _footer_wrapped.count("\n") + 1
assert _n_footer_lines <= 15, f"footer is {_n_footer_lines} wrapped lines; >15 collides with panel e"
fig.text(0.030, 0.010, _footer_wrapped, fontsize=5.8, color=INK, va="bottom", ha="left",
         linespacing=1.36)

for ext, kw in (("png", dict(dpi=300)), ("pdf", {})):
    p = FIGDIR / f"{NAME}.{ext}"
    fig.savefig(p, **kw)
    print("wrote", p)
p = OUTDIR / f"{NAME}.png"
fig.savefig(p, dpi=300)
print("wrote", p)

# ---------------------------------------------------------------------------
# TSVs: every plotted value
# ---------------------------------------------------------------------------
SRC_REPL = str(REPL)
funnel_rows = []
for stage, (lab, real, null, note) in enumerate(FUNNEL_STAGES, start=1):
    funnel_rows.append(dict(panel="a", stage=stage, stage_label=lab.replace("\n", " "),
                            series="real", n=real, note="",
                            source=f"{SRC_REPL}/pas_K2/all_features.tsv"))
    funnel_rows.append(dict(panel="a", stage=stage, stage_label=lab.replace("\n", " "),
                            series=f"null_mean_of_{N_NULL_COMBOS}", n=null, note=note,
                            source=(f"{SRC_REPL}/pas_K2/null_control.tsv"
                                    if stage in (1, 4, 5) else f"{VERIFY}/null_qhits_raw.tsv")))
fun_df = pd.DataFrame(funnel_rows)
fun_df.loc[len(fun_df)] = dict(panel="a", stage=2, stage_label="null per-combination detail",
                               series="null_called_ge1_per_combo",
                               n=float("nan"), note=";".join(map(str, NULL_CALLED_GE1)),
                               source=f"{VERIFY}/null_qhits_raw.tsv")
fun_df.loc[len(fun_df)] = dict(
    panel="a", stage=4, stage_label="dropped: called in >=2 patients but no single direction reached K",
    series="real", n=N_SPLIT_DIRECTION,
    note=f"stage 3 -> stage 4 loses {FUN['called_geK'] - FUN['replicated']} = "
         f"{N_SPLIT_DIRECTION} split-direction + {N_DISCORDANT_VETO} vetoed",
    source=f"{SRC_REPL}/pas_K2/all_features.tsv")
fun_df.loc[len(fun_df)] = dict(
    panel="a", stage=4, stage_label="dropped: reached K in one direction but an opposite-direction patient vetoed it",
    series="real", n=N_DISCORDANT_VETO, note="discordant_policy=exclude",
    source=f"{SRC_REPL}/pas_K2/all_features.tsv")
p = OUTDIR / f"{NAME}_funnel.tsv"
fun_df.to_csv(p, sep="\t", index=False, float_format="%.4f")
print("wrote", p)

pairs_out = pairs.rename(columns={"replicated_q_K2": "replicated_K2",
                                  "replicated_q_floor_K2": "replicated_K2_floor",
                                  "replicated_q_K3": "replicated_K3",
                                  "replicated_q_floor_K3": "replicated_K3_floor"})
pairs_out["pair_label"] = pairs_out.pair.map(pretty_pair)
pairs_out["panel"] = "b"
pairs_out["source"] = f"{SRC_REPL}/per_pair_pas_K2.tsv + per_pair_pas_K3.tsv"
p = OUTDIR / f"{NAME}_per_pair.tsv"
pairs_out[["panel", "pair", "pair_label", "n_patients_tested", "n_gsm_tested", "n_features",
           "replicated_K2", "replicated_K2_floor", "replicated_K3", "replicated_K3_floor",
           "empirical_p_q", "multi_patient", "plotted_panel_b", "source"]].to_csv(
    p, sep="\t", index=False)
print("wrote", p)

eff = pd.DataFrame(dict(panel="c", bin_left=bins[:-1], bin_right=bins[1:],
                        n_replicated=hist_all,
                        below_effect_floor=bins[:-1] + BINW / 2 < FLOOR))
eff["source"] = f"{SRC_REPL}/pas_K2/replicated.tsv (mean_effect_consensus, abs)"
p = OUTDIR / f"{NAME}_effects.tsv"
eff.to_csv(p, sep="\t", index=False, float_format="%.4f")
print("wrote", p)

sup = pd.DataFrame(dict(panel="d", replication_count=SUPPORT_K, n_replicated=SUPPORT_N))
sup["in_K3_sensitivity_set"] = sup.replication_count >= K_SENS
sup["source"] = f"{SRC_REPL}/pas_K2/replicated.tsv (replication_count)"
p = OUTDIR / f"{NAME}_support.tsv"
sup.to_csv(p, sep="\t", index=False)
print("wrote", p)

ctx_out = ctx.copy()
ctx_out["panel"] = "e"
ctx_out["gtf"] = str(GTF)
p = OUTDIR / f"{NAME}_context.tsv"
ctx_out.to_csv(p, sep="\t", index=False)
print("wrote", p)

ctx_sum = pd.DataFrame([dict(panel="e", category=c, n=int(CTX_N[c]),
                             fraction=float(CTX_N[c] / N_DISTINCT_PAS)) for c in CONTEXT_ORDER])
for k, n, f_ in (("roll-up: any same-strand gene body", int((ctx.n_samestrand_gene_bodies > 0).sum()), FRAC_GENEBODY),
                 ("roll-up: exonic (any same-strand gene)", int(ctx.exonic_any_gene.sum()), FRAC_EXONIC),
                 ("roll-up: any gene's 3' UTR", N_ANY3, FRAC_ANY_3UTR),
                 ("caveat: inside >=2 overlapping same-strand genes", N_MULTI_GENE, FRAC_MULTI_GENE),
                 ("caveat: 3'UTR PAS assigned to a different gene", N_3UTR_MISASSIGNED, FRAC_3UTR_MISASSIGNED)):
    ctx_sum.loc[len(ctx_sum)] = dict(panel="e", category=k, n=n, fraction=f_)
ctx_sum["denominator"] = N_DISTINCT_PAS
ctx_sum["source"] = f"recomputed: {SRC_REPL}/pas_K2/replicated.tsv vs {GTF}"
p = OUTDIR / f"{NAME}_context_summary.tsv"
ctx_sum.to_csv(p, sep="\t", index=False, float_format="%.5f")
print("wrote", p)

cohort_rows = [
    ("version", VERSION, ""), ("peakatail_code", V["commit"], V["commit_note"]),
    ("n_patients_with_multiple_libraries", N_MULTILIB_PATIENTS,
     "libraries per patient: " + ",".join(f"{k}={v}" for k, v in sorted(_libs_per_patient.items()))),
    ("tree", str(TREE), "read-only"),
    ("source_of_truth", V["doc_path"], V["verdict"]),
    ("n_patients", N_PATIENTS, ",".join(PATIENTS)),
    ("n_libraries_gsm", N_GSM, ""), ("n_celltype_pairs", N_PAIRS, ""),
    ("n_pairs_multi_patient", N_MULTI_PAIRS, f">={K_PRIMARY} patients"),
    ("n_pairs_with_replicated", N_PAIRS_WITH_HITS, ""),
    ("pas_universe", UNIVERSE, "IP-pass AND tier-1 AND cohort clip-molecule sum >=2"),
    ("cells_entering_tests", CELLS_TESTED, "sum of n_cells_final over the 15 primary libraries"),
    ("cells_confirmed_label_matched", CELLS_CONFIRMED_MATCHED, "sum of n_cells_matched"),
    ("real_tests_total", REAL_TESTS_TOTAL, "sum of true.n_tests over the 15 libraries"),
    ("real_calls_total_q<0.05", REAL_CALLS_TOTAL, "sum of per_sample_called"),
    ("null_tests_total", NULL_TESTS_TOTAL, f"15 libraries x {N_NULL_COMBOS} permutations"),
    ("null_nominal_q_hits_total", NULL_QHITS_TOTAL, "none co-occurring in two patients"),
    ("n_tested_pair_pas", N_FEATURES, ""),
    ("n_replicated_K2", N_REPL, ""), ("n_replicated_K2_floor", N_REPL_FLOOR, ""),
    ("n_replicated_K3", N_REPL_K3, ""), ("n_replicated_K3_floor", N_REPL_K3_FLOOR, ""),
    ("n_gene_level_K2", S_G2["real"]["n_replicated"], ""),
    ("n_distinct_pas", N_DISTINCT_PAS, ""), ("n_distinct_genes", N_DISTINCT_GENES, ""),
    ("split_direction_dropped", N_SPLIT_DIRECTION,
     "called in >=2 patients but no single direction reached K"),
    ("discordance_vetoed", N_DISCORDANT_VETO,
     "reached K in one direction but an opposite-direction patient vetoed it"),
    ("empirical_p_floor", EMP_P, f"1/({N_NULL_COMBOS}+1); NOT an FDR"),
    # ---- panel c / footer values printed on the figure ----
    ("median_abs_consensus_effect", round(MEDIAN_EFFECT, 6), "panel c annotation; replicated.tsv"),
    ("mean_abs_consensus_effect", round(MEAN_EFFECT, 6), "panel c annotation; replicated.tsv"),
    ("n_consensus_below_effect_floor", N_CONSENSUS_BELOW_FLOOR,
     f"|mean_effect_consensus| < {FLOOR:g}; footer"),
    ("n_removed_by_effect_floor", N_REMOVED_BY_FLOOR,
     "floor applied per patient then patients re-counted; footer"),
    ("null_tested_mean_per_combo", round(NULL_TESTED_MEAN, 1), "panel a grey bar; null_control.tsv"),
    ("null_called_ge1_mean_per_combo", NULL_CALLED_GE1_MEAN, "panel a grey bar; null_qhits_raw.tsv"),
    (f"null_tests_{DISCLOSED_GSM}", DISCLOSED_NULL_TESTS, "footer BH-discreteness disclosure"),
    (f"null_qhits_{DISCLOSED_GSM}", DISCLOSED_NULL_QHITS, "footer BH-discreteness disclosure"),
] + ([
    # ---- printed on the figure but NOT derived here: cited to the write-up ----
    ("curated_cells", V["cohort"]["curated_cells"], "CITED 20 (verified vs labels/manifest.json)"),
    ("confirmed_cells", V["cohort"]["confirmed_cells"], "CITED 20 (verified vs labels/manifest.json)"),
    ("confirmed_fraction", V["cohort"]["confirmed_frac"], "CITED 20"),
    ("top30_gene_rows_misnamed_frac", V["cohort"]["top30_bad_gene_frac"],
     "CITED 20 disclosure 1; NOT recomputed here (would require the forbidden ranked list)"),
    ("sens_replicated_13_patients", EXPECT["sens_replicated"], "CITED 20 / PRIMARY_vs_SENSITIVITY.txt"),
    ("sens_replicated_floor_13_patients", EXPECT["sens_floor"], "CITED 20 / PRIMARY_vs_SENSITIVITY.txt"),
    ("metbone_clip_rate_quoted", "0.015%",
     "CITED 20 disclosure 2; a first-200,000-CB-read sampling artefact, premise WRONG"),
    ("metbone_actual_clip_molecules", "306k", "CITED 20 disclosure 2 / tool issue #99"),
] if (V["cohort"] and EXPECT) else [])
p = OUTDIR / f"{NAME}_cohort.tsv"
pd.DataFrame(cohort_rows, columns=["key", "value", "note"]).to_csv(p, sep="\t", index=False)
print("wrote", p)

p = OUTDIR / f"{NAME}_own3utr_genes.tsv"
safe.to_csv(p, sep="\t", index=False, float_format="%.5f")
print("wrote", p, f"({len(safe):,} re-assignment-safe PAS; AUDIT ONLY, not plotted)")

# ---------------------------------------------------------------------------
# sidecar caption (+ the index paragraph for 05_figure_index.md, which is
# owned by another process and is NOT edited here)
# ---------------------------------------------------------------------------
cap_md = f"""# Fig 6 — `laughney_switches` caption (generated by `scripts/manuscript_figures/laughney_switches.py`)

**Fig 6 | Reliable cell-type APA switches in a tumour cohort.** Across {N_PATIENTS} patients
({N_GSM} scRNA-seq libraries) of the Laughney lung-adenocarcinoma cohort, processed in a single cohort run so that
all libraries share one PAS identifier space, PeakATail called cell-type APA switches in {N_PAIRS} cell-type pairs
over a pre-registered precision-first PAS universe (clip-supported, internal-priming-filtered, ≥2 clip molecules;
{UNIVERSE:,} sites). Requiring a switch to be called (Fisher on cells, BH q<{FDR:g}) in the same direction in at least
{K_PRIMARY} independent patients, with any opposite-direction patient vetoing the call, **{N_REPL:,} of
{N_FEATURES:,} tested (pair, PAS) hypotheses replicate ({N_REPL/N_FEATURES*100:.2f}%)**, of which {N_REPL_FLOOR:,}
also exceed an effect floor of |Δproportion| ≥ {FLOOR:g}; these are {N_DISTINCT_PAS:,} distinct PAS in
{N_DISTINCT_GENES:,} genes across {N_PAIRS_WITH_HITS} pairs. **Under ten patient-wise label-shuffle nulls run
through the identical pipeline, nothing replicated in any combination.** Requiring three patients retains
{N_REPL_K3:,}.

**Panel a** — the funnel on one log axis, real (blue) against the mean of the {N_NULL_COMBOS} label-shuffle nulls
(grey): {FUN['tested']:,} tested (pair, PAS) hypotheses → {FUN['called_ge1']:,} called at BH q<{FDR:g} in ≥1 patient
→ {FUN['called_geK']:,} called in ≥{K_PRIMARY} patients in any direction → {FUN['replicated']:,} replicated in the
same direction with no opposite-direction patient (the {FUN['called_geK'] - FUN['replicated']} lost here split as
{N_SPLIT_DIRECTION} in which no single direction reached {K_PRIMARY} patients and {N_DISCORDANT_VETO} vetoed because a
patient went the other way) →
{FUN['replicated_floor']:,} also over the |Δproportion| ≥ {FLOOR:g} floor. The null tests a comparable
{NULL_TESTED_MEAN:,.0f} hypotheses per combination and reaches a mean of {NULL_CALLED_GE1_MEAN:.1f} single-patient
calls (range {min(NULL_CALLED_GE1)}–{max(NULL_CALLED_GE1)}) and **zero** at every later stage.
**Panel b** — the top {N_SHOW} of the {N_MULTI_PAIRS} cell-type pairs tested in ≥{K_PRIMARY} patients, K={K_PRIMARY}
(light) with the K={K_SENS} subset overlaid (dark); "n pt" is the number of patients in which the pair was testable.
The other {REST_N} multi-patient pairs hold {REST_SUM:,} more switches ({REST_MIN}–{REST_MAX} each);
{N_SINGLE_PATIENT_PAIRS} pairs were tested in a single patient and can never replicate.
**Panel c** — consensus |Δproportion| of the {N_REPL:,} replicated switches (median {MEDIAN_EFFECT:.2f}, mean
{MEAN_EFFECT:.2f}); the dashed line is the pre-registered floor. **Panel d** — patient support; the K={K_SENS}
sensitivity set is exactly the ≥3-patient tail ({N_REPL_K3:,}).
**Panel e (honesty panel)** — genomic context of the {N_DISTINCT_PAS:,} distinct replicated PAS, recomputed for this
figure from `{GTF}` as an exclusive, strand-matched partition:
{FRAC_OWN_3UTR*100:.1f}% own-gene 3′ UTR, {N_OTH3/N_DISTINCT_PAS*100:.1f}% another gene's 3′ UTR,
{N_EXNOT3/N_DISTINCT_PAS*100:.1f}% exonic but not 3′ UTR, {FRAC_INTRONIC*100:.1f}% intronic,
{FRAC_OUTSIDE*100:.1f}% outside any same-strand gene — i.e. {FRAC_GENEBODY*100:.1f}% inside a same-strand gene body,
{FRAC_EXONIC*100:.1f}% exonic, {FRAC_ANY_3UTR*100:.1f}% in some gene's 3′ UTR. **Only two of these five roll-ups exist in
20**: disclosure 1 quotes 62.5% and 96.8%, and both are reproduced here to within 0.1 pp; the exonic, intronic and
outside-any-gene figures appear in no write-up and are recomputed for this figure from the GTF. The 62.5% additionally
needs **a wording correction that matters**: 20 describes it as PAS "in own-gene 3′ UTRs", whereas the recomputation shows that
{FRAC_ANY_3UTR*100:.1f}% lie in *some* gene's 3′ UTR and only {FRAC_OWN_3UTR*100:.1f}% in the 3′ UTR of the gene the
caller actually assigned them to. Use the own-gene figure ({FRAC_OWN_3UTR*100:.1f}%) whenever the claim is about
gene-level interpretation.

**Gene names (20 disclosure 1, tool issue #99).** No ranked list of named top-switch genes is shown. 33% of the
top-30 gene-level rows in 20 name a spanning or readthrough model rather than the gene whose 3′ UTR holds the PAS
(CD68 labelled SENP3-EIF4A1, PTPRCAP under CORO1B, FKBP11 under AC073610.2; seven of the named genes have zero
assigned PAS cohort-wide), and an ambient-immunoglobulin signature (IGLL5 in non-plasma pairs) is present.
Recomputed here: {N_3UTR_MISASSIGNED:,} of the {N_ANY3:,}
3′-UTR PAS ({FRAC_3UTR_MISASSIGNED*100:.1f}%) lie in a *different* gene's 3′ UTR from the one the caller assigned,
and {N_MULTI_GENE:,} ({FRAC_MULTI_GENE*100:.1f}%) lie inside ≥2 overlapping same-strand gene models. Replication
statistics are computed on PAS ids and are unaffected. The re-assignment-safe subset — PAS verified against the GTF
to lie in their own assigned gene's 3′ UTR — is written to
`results/figures/manuscript/{NAME}_own3utr_genes.tsv` in genomic order for audit only; it is not a ranked list and
must not be published as one until the assignment is fixed.

**Record and controls.** PeakATail code `{V['commit']}` ({V['commit_note']}); {V['record_note']}. Source of truth
`{V['doc_path']}` ({V['verdict']}). The null resolves only to an empirical p ≤ {EMP_P:.3f} (the
{N_NULL_COMBOS}-permutation floor): report it as "none in {N_NULL_COMBOS} nulls", **never** as an FDR estimate.
{NULL_TESTS_TOTAL/1e6:.1f}M null tests across the {N_GSM} libraries produced {NULL_QHITS_TOTAL} nominal q<{FDR:g}
calls, none co-occurring in two patients. GSM3516664-MetBone is excluded by pre-registration (13 addendum item 8),
but **its stated premise was wrong**: the 0.015% clip rate came from the caller sampling only the first 200,000 CB
reads of a coordinate-sorted BAM (head of chr1) — MetBone in fact has 306k clip molecules (tool issue #99). The
"no clip evidence" justification must never be repeated; the {EXPECT['sens_patients'] if EXPECT else '13'}-patient
sensitivity run that includes it gives {format(EXPECT['sens_replicated'], ',') if EXPECT else 'n/a'} /
{format(EXPECT['sens_floor'], ',') if EXPECT else 'n/a'}, a 0.4% difference confined to one pair.

Sources: `{V['doc_path']}` ({V['verdict']}) and, read directly, `{SRC_REPL}/` (`pas_K2`, `pas_K3`, `gene_K2`:
`all_features.tsv`, `replicated.tsv`, `null_control.tsv`, `summary.json`; `per_pair_pas_K{{2,3}}.tsv`), the per-GSM
switch summaries `{SWITCH}/<GSM>/summary.json`, the verifier's null recount `{VERIFY}/`, and `{GTF}` for panel e.
Every plotted value: `results/figures/manuscript/{NAME}_{{funnel,per_pair,effects,support,context,context_summary,cohort}}.tsv`.

---

## Index paragraph (for `manuscript/05_figure_index.md` — paste there; this script does not edit that file)

**Fig 6 — `laughney_switches` — reliable cell-type APA switches in a tumour cohort.** Stage-3 patient-level
replication on the Laughney lung-adenocarcinoma cohort, PeakATail code `{V['commit']}`, source of truth
`{V['doc_path']}` ({V['verdict']}). Headline: across {N_PATIENTS} patients / {N_GSM} libraries and {N_PAIRS}
cell-type pairs over a {UNIVERSE:,}-PAS precision-first universe, **{N_REPL:,} of {N_FEATURES:,} tested (pair, PAS)
hypotheses replicate in ≥{K_PRIMARY} patients in the same direction ({N_REPL/N_FEATURES*100:.2f}%),
{N_REPL_FLOOR:,} of them over the |Δproportion| ≥ {FLOOR:g} floor, covering {N_DISTINCT_PAS:,} PAS in
{N_DISTINCT_GENES:,} genes across {N_PAIRS_WITH_HITS} pairs; requiring three patients retains {N_REPL_K3:,};
and nothing replicates in any of {N_NULL_COMBOS} patient-wise label-shuffle nulls** ({NULL_TESTS_TOTAL/1e6:.1f}M
null tests → {NULL_QHITS_TOTAL} nominal q<{FDR:g} calls, none in two patients). Caveats that must travel with the
figure: (i) the null resolves only to empirical p ≤ {EMP_P:.3f}, so say "none in {N_NULL_COMBOS} nulls" and never
quote an FDR; (ii) these are v1-code numbers and the Stage-3 v2 chain on `9dfdefb` will supersede them
(`LAUGHNEY_SWITCHES_VERSION=v2_code`); (iii) the MetBone exclusion is pre-registered but its stated premise was a
clip-rate sampling artefact (issue #99), and the 13-patient sensitivity run differs by 0.4%; (iv) **no named
top-gene list may be published** until PAS→gene re-assignment — 33% of the top-30 gene rows name a
spanning/readthrough model (issue #99), and {FRAC_3UTR_MISASSIGNED*100:.1f}% of 3′-UTR PAS are assigned to a
different gene than the one whose 3′ UTR they occupy; (v) pairs tested in more patients replicate more, so panel b
tracks cohort composition as much as biology.
"""
p = FIGDIR / f"{NAME}.caption.md"
p.write_text(cap_md)
print("wrote", p)

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

print(f"VERSION={VERSION} code={V['commit']}: {N_REPL:,} replicated (pair, PAS) in {N_PATIENTS} patients, "
      f"{N_REPL_FLOOR:,} over the effect floor, {N_REPL_K3:,} at K=3, "
      f"0 in all {N_NULL_COMBOS} nulls (empirical p <= {EMP_P:.3f})")

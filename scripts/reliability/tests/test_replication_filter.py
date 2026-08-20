#!/usr/bin/env python3
"""
Unit tests for scripts/reliability/replication_filter.py on SYNTHETIC inputs.

Run:
    /mnt/ssd1/Projects/PeakATail_wd/tools/PeakATail/.venv/bin/python -m pytest \
        scripts/reliability/tests -v
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
SCRIPT = HERE.parent / "replication_filter.py"


def _load():
    spec = importlib.util.spec_from_file_location("_rf_under_test", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["_rf_under_test"] = mod
    spec.loader.exec_module(mod)
    return mod


rf = _load()

DIFF_COLS = ["pas_id", "gene_id", "chrom", "start", "end", "strand", "cluster1", "cluster2",
             "pvalue", "qvalue", "n_cells", "n_cells_cluster1", "n_cells_cluster2",
             "n_cells_expr_cluster1", "n_cells_expr_cluster2", "n_reads_pas_cluster1",
             "n_reads_pas_cluster2", "n_reads_gene_cluster1", "n_reads_gene_cluster2",
             "odds_ratio", "delta_proportion", "log2fc"]

# synthetic PAS universe: pas_id -> (gene, chrom, start, end, strand)
PAS = {
    "1": ("GENE_A", "1", 1000, 1020, "+"),   # proximal of GENE_A
    "2": ("GENE_A", "1", 2000, 2020, "+"),   # distal of GENE_A
    "3": ("GENE_B", "2", 5000, 5020, "-"),   # distal of GENE_B (minus strand: smallest start)
    "4": ("GENE_B", "2", 6000, 6020, "-"),   # proximal of GENE_B
    "5": ("GENE_C", "3", 100, 120, "+"),
    "6": ("GENE_C", "3", 300, 320, "+"),
    "7": ("GENE_D", "4", 100, 120, "+"),
    "8": ("GENE_D", "4", 300, 320, "+"),
}


def make_pair_tsv(path: Path, c1: str, c2: str, rows: dict, strategy="fisher"):
    """rows: pas_id -> (qvalue, delta_proportion) in the (c1, c2) orientation."""
    recs = []
    for pid, (q, d) in rows.items():
        g, ch, s, e, st = PAS[pid]
        recs.append(dict(pas_id=pid, gene_id=g, chrom=ch, start=s, end=e, strand=st,
                         cluster1=c1, cluster2=c2, pvalue=q / 2, qvalue=q, n_cells=100,
                         n_cells_cluster1=50, n_cells_cluster2=50, n_cells_expr_cluster1=20,
                         n_cells_expr_cluster2=20, n_reads_pas_cluster1=30, n_reads_pas_cluster2=30,
                         n_reads_gene_cluster1=60, n_reads_gene_cluster2=60, odds_ratio=1.5,
                         delta_proportion=d, log2fc=-d * 3))   # log2fc has the OPPOSITE sign
    path.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(recs, columns=DIFF_COLS).to_csv(path / f"{strategy}_{c1}_vs_{c2}.tsv",
                                                 sep="\t", index=False)


@pytest.fixture
def three_samples(tmp_path):
    """Three samples, pair (T, N). Canonical orientation is (N, T) because 'N' < 'T',
    so s1 and s2 (written as T_vs_N) are the FLIPPED ones and s3 (N_vs_T) is not.
    Effects below are quoted in each file's own orientation; canonical effects
    are the negatives for s1/s2.

    PAS 1: called +0.3 in all three          -> replicated, count 3
    PAS 2: called +0.3 in s1, -0.3 in s2     -> discordant, count 1
    PAS 3: called -0.3 in s1 only            -> count 1
    PAS 4: called +0.3 in s1, +0.05 in s2    -> count 2 but floor count 1
    PAS 5: q=0.5 everywhere                  -> never called
    PAS 6: only present in s1, s2 (not s3)   -> n_samples_tested 2, called both +
    PAS 7: called +0.3 in s1 and s3 (s3 file is flipped: written as -0.3 in N_vs_T)
    PAS 8: effect NaN in s1 (q small), +0.3 in s2 -> count 1 (NaN never called)
    """
    s1 = tmp_path / "s1" / "differential"
    s2 = tmp_path / "s2" / "differential"
    s3 = tmp_path / "s3" / "differential"
    make_pair_tsv(s1, "T", "N", {"1": (1e-5, 0.3), "2": (1e-4, 0.3), "3": (1e-3, -0.3),
                                 "4": (1e-3, 0.3), "5": (0.5, 0.2), "6": (1e-3, 0.25),
                                 "7": (1e-3, 0.3), "8": (1e-3, float("nan"))})
    make_pair_tsv(s2, "T", "N", {"1": (1e-4, 0.3), "2": (1e-4, -0.3), "3": (0.9, -0.3),
                                 "4": (1e-3, 0.05), "5": (0.5, 0.2), "6": (1e-3, 0.25),
                                 "7": (0.7, 0.3), "8": (1e-3, 0.3)})
    # sample 3 wrote the pair the other way round: effects are in (N, T) orientation
    make_pair_tsv(s3, "N", "T", {"1": (1e-3, -0.3), "2": (0.5, 0.0), "3": (0.9, 0.3),
                                 "4": (0.9, 0.0), "5": (0.5, -0.2),
                                 "7": (1e-3, -0.3), "8": (0.9, 0.0)})
    return tmp_path


def _run_cli(args):
    return rf.main([str(a) for a in args])


# ---------------------------------------------------------------- parsing
def test_parse_pair_filename_handles_underscores():
    assert rf.parse_pair_filename("fisher_Normal_vs_StageIA.tsv", None) == ("fisher", "Normal", "StageIA")
    assert rf.parse_pair_filename("nb_pairwise_T_cell_vs_B_cell.tsv", None) == ("nb_pairwise", "T_cell", "B_cell")
    assert rf.parse_pair_filename("fisher_A_vs_B.tsv", "nb_pairwise") is None
    assert rf.parse_pair_filename("nb_multi_omnibus.tsv", None) is None
    assert rf.parse_pair_filename("switch_diff_long.tsv", None) is None


def test_canonical_pair_flips():
    assert rf.canonical_pair("T", "N") == ("N", "T", True)
    assert rf.canonical_pair("N", "T") == ("N", "T", False)


def test_pick_effect_col_prefers_delta_proportion():
    assert rf.pick_effect_col(["qvalue", "log2fc", "delta_proportion"]) == "delta_proportion"
    assert rf.pick_effect_col(["qvalue", "log2fc"]) == "log2fc"
    with pytest.raises(KeyError):
        rf.pick_effect_col(["qvalue"])
    with pytest.raises(KeyError):
        rf.pick_effect_col(["qvalue", "log2fc"], "dpdui")


def test_bh_qvalues_monotone_and_nan_safe():
    from scipy.stats import false_discovery_control
    p = np.array([0.01, 0.04, 0.03, np.nan, 0.5])
    q = rf.bh_qvalues(p)
    assert np.isnan(q[3])
    ok = np.isfinite(p)
    assert q[ok] == pytest.approx(false_discovery_control(p[ok], method="bh"))
    assert q[0] == pytest.approx(0.04)            # 0.01*4/1
    assert q[1] == pytest.approx(0.04 * 4 / 3)    # 0.0533 (no smaller value to its right)
    assert q[4] == pytest.approx(0.5)
    assert np.all(np.diff(np.sort(q[np.isfinite(q)])) >= 0)


def test_load_diff_sample_canonicalises_orientation(three_samples):
    ldf, strat, eff = rf.load_diff_sample("s3", three_samples / "s3")
    assert strat == "fisher" and eff == "delta_proportion"
    assert set(zip(ldf["c1"], ldf["c2"])) == {("N", "T")}
    assert not ldf["flipped"].any()                         # N_vs_T is already canonical
    r = ldf.set_index("feature_id")
    assert r.loc["1", "effect"] == pytest.approx(-0.3)
    ldf1, _, _ = rf.load_diff_sample("s1", three_samples / "s1")
    assert set(zip(ldf1["c1"], ldf1["c2"])) == {("N", "T")}
    assert ldf1["flipped"].all()                            # T_vs_N file -> flipped
    r1 = ldf1.set_index("feature_id")
    assert r1.loc["1", "effect"] == pytest.approx(-0.3)     # +0.3 in T_vs_N -> -0.3 canonical
    assert r1.loc["7", "effect"] == pytest.approx(-0.3)
    assert r1.loc["3", "effect"] == pytest.approx(0.3)


def test_load_diff_sample_explicit_log2fc_column(three_samples):
    ldf, _, eff = rf.load_diff_sample("s1", three_samples / "s1", effect_col="log2fc")
    assert eff == "log2fc"
    # file log2fc = -0.9 (T_vs_N); flipped to canonical (N, T) -> +0.9
    assert ldf.set_index("feature_id").loc["1", "effect"] == pytest.approx(0.9)


# --------------------------------------------------------- replication logic
def test_replication_counts_and_flags(three_samples):
    frames = [rf.load_diff_sample(s, three_samples / s)[0] for s in ("s1", "s2", "s3")]
    long_df = pd.concat(frames, ignore_index=True)
    prm = rf.Params(min_samples=2, fdr=0.05, effect_floor=0.1)
    t = rf.replication_table(long_df, prm, ["s1", "s2", "s3"]).set_index("feature_id")

    assert t.loc["1", "replication_count"] == 3
    assert t.loc["1", "consensus_direction"] == "-"          # canonical (N, T): higher in T
    assert bool(t.loc["1", "direction_consistent"]) is True
    assert bool(t.loc["1", "passes_replication"]) is True
    assert bool(t.loc["1", "passes_replication_floor"]) is True
    assert t.loc["1", "n_samples_tested"] == 3
    assert t.loc["1", "utr_direction"] == "higher_in_c2"          # c2 == T after canonicalisation

    assert t.loc["2", "n_pos"] == 1 and t.loc["2", "n_neg"] == 1
    assert t.loc["2", "replication_count"] == 1
    assert bool(t.loc["2", "direction_consistent"]) is False
    assert t.loc["2", "consensus_direction"] == "tie"
    assert bool(t.loc["2", "passes_replication"]) is False

    assert t.loc["3", "replication_count"] == 1
    assert bool(t.loc["3", "passes_replication"]) is False

    assert t.loc["4", "replication_count"] == 2
    assert bool(t.loc["4", "passes_replication"]) is True
    assert t.loc["4", "replication_count_floor"] == 1
    assert bool(t.loc["4", "effect_floor_pass"]) is False
    assert bool(t.loc["4", "passes_replication_floor"]) is False

    assert t.loc["5", "n_samples_called"] == 0
    assert t.loc["5", "consensus_direction"] == ""
    assert bool(t.loc["5", "passes_replication"]) is False

    assert t.loc["6", "n_samples_tested"] == 2
    assert t.loc["6", "replication_count"] == 2
    assert pd.isna(t.loc["6", "q__s3"])

    assert t.loc["7", "replication_count"] == 2           # s1 + flipped s3
    assert bool(t.loc["7", "passes_replication"]) is True
    assert t.loc["7", "effect__s3"] == pytest.approx(-0.3)
    assert t.loc["7", "effect__s1"] == pytest.approx(-0.3)

    assert t.loc["8", "replication_count"] == 1            # NaN effect in s1 is not a call
    assert bool(t.loc["8", "called__s1"]) is False

    # per-sample columns carry the raw q values
    assert t.loc["1", "q__s1"] == pytest.approx(1e-5)
    assert t.loc["1", "min_q_called"] == pytest.approx(1e-5)
    assert t.loc["1", "mean_effect_consensus"] == pytest.approx(-0.3)


def test_discordant_policy_allow_and_K(three_samples):
    frames = [rf.load_diff_sample(s, three_samples / s)[0] for s in ("s1", "s2", "s3")]
    long_df = pd.concat(frames, ignore_index=True)
    # K=3: only PAS 1 survives
    t3 = rf.replication_table(long_df, rf.Params(min_samples=3), ["s1", "s2", "s3"]).set_index("feature_id")
    assert set(t3.index[t3["passes_replication"]]) == {"1"}
    # allow discordant + K=1: PAS 2 passes
    ta = rf.replication_table(long_df, rf.Params(min_samples=1, discordant_policy="allow"),
                              ["s1", "s2", "s3"]).set_index("feature_id")
    assert bool(ta.loc["2", "passes_replication"]) is True
    te = rf.replication_table(long_df, rf.Params(min_samples=1, discordant_policy="exclude"),
                              ["s1", "s2", "s3"]).set_index("feature_id")
    assert bool(te.loc["2", "passes_replication"]) is False


def test_gene_level_collapse(three_samples):
    frames = [rf.load_diff_sample(s, three_samples / s)[0] for s in ("s1", "s2", "s3")]
    long_df = pd.concat(frames, ignore_index=True)
    pas_t = rf.replication_table(long_df, rf.Params(), ["s1", "s2", "s3"])
    g = rf.collapse_to_gene(pas_t).set_index("gene_id")
    assert g.loc["GENE_A", "representative_pas_id"] == "1"
    assert g.loc["GENE_A", "n_pas_tested"] == 2
    assert g.loc["GENE_A", "n_pas_replicated"] == 1
    assert bool(g.loc["GENE_A", "passes_replication"]) is True
    # distal PAS of GENE_A (+ strand) is pas 2, which is discordant -> not replicated
    assert bool(g.loc["GENE_A", "distal_pas_replicated"]) is False
    assert g.loc["GENE_A", "distal_utr_direction"] == ""
    # GENE_B (- strand): distal = smallest start = pas 3 (count 1) -> gene replicates via pas 4
    assert g.loc["GENE_B", "representative_pas_id"] == "4"
    assert bool(g.loc["GENE_B", "distal_pas_replicated"]) is False
    # GENE_D: pas 7 replicated, pas 8 not; distal = pas 8 -> not replicated
    assert g.loc["GENE_D", "representative_pas_id"] == "7"
    assert bool(g.loc["GENE_D", "passes_replication"]) is True


def test_gene_level_distal_direction(tmp_path):
    """Distal PAS replicated with + effect (higher in c1) => c1_longer."""
    for s in ("a", "b"):
        make_pair_tsv(tmp_path / s / "differential", "N", "T",
                      {"1": (1e-3, -0.3), "2": (1e-3, 0.3)})
    frames = [rf.load_diff_sample(s, tmp_path / s)[0] for s in ("a", "b")]
    t = rf.replication_table(pd.concat(frames), rf.Params(), ["a", "b"])
    g = rf.collapse_to_gene(t).set_index("gene_id")
    assert bool(g.loc["GENE_A", "distal_pas_replicated"]) is True
    assert g.loc["GENE_A", "distal_utr_direction"] == "c1_longer"
    # with log2fc (opposite sign convention) the same biology reads the same way
    frames = [rf.load_diff_sample(s, tmp_path / s, effect_col="log2fc")[0] for s in ("a", "b")]
    t = rf.replication_table(pd.concat(frames), rf.Params(effect_sign=-1, effect_floor=0.5), ["a", "b"])
    g = rf.collapse_to_gene(t).set_index("gene_id")
    assert g.loc["GENE_A", "distal_utr_direction"] == "c1_longer"


# ------------------------------------------------------------------ nulls
def test_expand_null_specs_forms(tmp_path):
    for s in ("s1", "s2"):
        for k in range(3):
            (tmp_path / "null" / s / f"perm_{k:02d}" / "differential").mkdir(parents=True)
    m = rf.expand_null_specs([str(tmp_path / "null" / "{sample}" / "perm_*")], ["s1", "s2"])
    assert len(m["s1"]) == 3 and len(m["s2"]) == 3
    m = rf.expand_null_specs([f"s2={tmp_path}/null/s2/perm_0[01]"], ["s1", "s2"])
    assert len(m["s1"]) == 0 and len(m["s2"]) == 2
    m = rf.expand_null_specs([f"{tmp_path}/null/s1/perm_*", f"{tmp_path}/null/s2/perm_00"], ["s1", "s2"])
    assert len(m["s1"]) == 3 and len(m["s2"]) == 1
    with pytest.raises(KeyError):
        rf.expand_null_specs(["zz=/nowhere/*"], ["s1"])
    with pytest.raises(ValueError):
        rf.expand_null_specs(["/a/*", "/b/*"], ["s1"])


def test_null_combinations():
    rng = np.random.default_rng(0)
    c = rf.null_combinations({"a": 3, "b": 5, "c": 0}, None, rng)
    assert len(c) == 3 and all(set(x) == {"a", "b"} for x in c) and c[2] == {"a": 2, "b": 2}
    c = rf.null_combinations({"a": 3, "b": 5}, 10, rng)
    assert len(c) == 10 and all(0 <= x["a"] < 3 and 0 <= x["b"] < 5 for x in c)
    assert rf.null_combinations({"a": 0}, None, rng) == []


def _make_null_dirs(root: Path, n_perms: int, hot: bool):
    """Null perms: nothing significant (hot=False) or PAS 1 replicating everywhere."""
    for s in ("s1", "s2", "s3"):
        for k in range(n_perms):
            d = root / s / f"perm_{k:02d}" / "differential"
            q1 = 1e-3 if hot else 0.6
            make_pair_tsv(d, "T", "N", {"1": (q1, 0.3), "2": (0.7, 0.1), "3": (0.8, -0.2)})


def test_cli_end_to_end_with_quiet_null(three_samples, tmp_path):
    nroot = tmp_path / "nulls"
    _make_null_dirs(nroot, 4, hot=False)
    out = tmp_path / "out"
    rc = _run_cli(["--sample", f"s1={three_samples / 's1'}", "--sample", f"s2={three_samples / 's2'}",
                   "--sample", f"s3={three_samples / 's3'}", "--out", out,
                   "--null-dirs", f"{nroot}/{{sample}}/perm_*", "--note", "synthetic"])
    assert rc == 0
    summ = json.loads((out / "summary.json").read_text())
    assert summ["params"]["effect_col"] == "delta_proportion"
    assert summ["params"]["strategy"] == "fisher"
    assert summ["n_samples"] == 3
    assert summ["pairs_seen"] == ["N_vs_T"]
    assert summ["real"]["n_replicated"] == 4          # PAS 1, 4, 6, 7
    assert summ["real"]["n_replicated_floor"] == 3    # PAS 4 fails the floor
    assert summ["real"]["per_pair"]["N_vs_T"]["n_replicated"] == 4
    assert summ["null"]["enabled"] and summ["null"]["n_combos"] == 4
    assert summ["null"]["replicated"]["null_mean"] == 0.0
    assert summ["null"]["replicated"]["empirical_p_ge_observed"] == pytest.approx(1 / 5)
    assert summ["null"]["replicated"]["expected_false_replicated_fraction"] == 0.0
    assert summ["note"] == "synthetic"
    rep = pd.read_csv(out / "replicated.tsv", sep="\t", dtype={"feature_id": str})
    assert sorted(rep["feature_id"]) == ["1", "4", "6", "7"]
    allf = pd.read_csv(out / "all_features.tsv", sep="\t", dtype={"feature_id": str})
    assert len(allf) == 8
    nulltsv = pd.read_csv(out / "null_control.tsv", sep="\t")
    assert list(nulltsv["n_replicated"]) == [0, 0, 0, 0]
    assert {"perm__s1", "perm__s2", "perm__s3"} <= set(nulltsv.columns)


def test_cli_hot_null_and_random_combos(three_samples, tmp_path):
    nroot = tmp_path / "nulls"
    _make_null_dirs(nroot, 2, hot=True)
    out = tmp_path / "out"
    rc = _run_cli(["--sample", three_samples / "s1", "--sample", three_samples / "s2",
                   "--sample", three_samples / "s3", "--out", out,
                   "--null-dirs", f"{nroot}/{{sample}}/perm_*", "--null-combos", 5, "--seed", 1])
    assert rc == 0
    summ = json.loads((out / "summary.json").read_text())
    assert [s["name"] for s in summ["samples"]] == ["s1", "s2", "s3"]   # names from dir basenames
    assert summ["null"]["n_combos"] == 5
    assert summ["null"]["replicated"]["null_per_combo"] == [1] * 5
    assert summ["null"]["replicated"]["expected_false_replicated_fraction"] == pytest.approx(1 / 4)


def test_cli_gene_level_and_pairs_restriction(three_samples, tmp_path):
    # add a second pair to s1 only; restrict to N,T -> second pair ignored
    make_pair_tsv(three_samples / "s1" / "differential", "T", "X", {"1": (1e-3, 0.3)})
    out = tmp_path / "out"
    rc = _run_cli(["--sample", three_samples / "s1", "--sample", three_samples / "s2",
                   "--sample", three_samples / "s3", "--out", out, "--level", "gene",
                   "--pairs", "T,N", "--min-samples", 2])
    assert rc == 0
    summ = json.loads((out / "summary.json")
                      .read_text())
    assert summ["pairs_seen"] == ["N_vs_T"]
    g = pd.read_csv(out / "all_features.tsv", sep="\t")
    assert set(g["gene_id"]) == {"GENE_A", "GENE_B", "GENE_C", "GENE_D"}
    # GENE_A (pas 1), GENE_B (pas 4), GENE_C (pas 6), GENE_D (pas 7)
    assert summ["real"]["n_replicated"] == 4
    assert summ["real"]["n_replicated_floor"] == 3           # GENE_B's pas 4 fails the floor


def test_cli_rejects_nb_multi(tmp_path):
    d = tmp_path / "s1" / "differential"
    d.mkdir(parents=True)
    pd.DataFrame({"pas_id": ["1"], "pvalue": [0.1], "qvalue": [0.2]}).to_csv(
        d / "nb_multi_omnibus.tsv", sep="\t", index=False)
    with pytest.raises((SystemExit, FileNotFoundError)):
        _run_cli(["--sample", tmp_path / "s1", "--sample", tmp_path / "s1", "--out", tmp_path / "o"])


def test_cli_rejects_duplicate_sample_names(three_samples, tmp_path):
    with pytest.raises(SystemExit):
        _run_cli(["--sample", f"x={three_samples / 's1'}", "--sample", f"x={three_samples / 's2'}",
                  "--out", tmp_path / "o"])


# ------------------------------------------------------------------ PDUI
PDUI_COLS = ["gene_id", "transcript_id", "proximal_pas_id", "distal_pas_id", "cell", "pdui",
             "proximal_reads", "distal_reads", "total_reads", "cluster", "proximal_chrom",
             "proximal_start", "proximal_end", "proximal_strand", "distal_chrom", "distal_start",
             "distal_end", "distal_strand"]


def make_pdui_table(path: Path, rng, shift: dict, n_cells=40, genes=("G1", "G2", "G3")):
    """Per-cell pdui_classic.tsv with clusters A, B; shift[gene] = mean PDUI(A) - mean PDUI(B)."""
    recs = []
    for g in genes:
        for cl in ("A", "B"):
            base = 0.5 + (shift.get(g, 0.0) / 2 if cl == "A" else -shift.get(g, 0.0) / 2)
            vals = np.clip(rng.normal(base, 0.05, n_cells), 0, 1)
            for i, v in enumerate(vals):
                recs.append(dict(gene_id=g, transcript_id="_gene_", proximal_pas_id=1, distal_pas_id=2,
                                 cell=f"{cl}_{i}", pdui=v, proximal_reads=5, distal_reads=5,
                                 total_reads=10, cluster=cl, proximal_chrom="1", proximal_start=100,
                                 proximal_end=120, proximal_strand="+", distal_chrom="1",
                                 distal_start=900, distal_end=920, distal_strand="+"))
    # a few cells with NaN pdui (zero reads) must be ignored
    recs.append(dict(gene_id="G1", transcript_id="_gene_", proximal_pas_id=1, distal_pas_id=2,
                     cell="A_nan", pdui=np.nan, proximal_reads=0, distal_reads=0, total_reads=0,
                     cluster="A", proximal_chrom="1", proximal_start=100, proximal_end=120,
                     proximal_strand="+", distal_chrom="1", distal_start=900, distal_end=920,
                     distal_strand="+"))
    path.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(recs, columns=PDUI_COLS).to_csv(path / "pdui_classic.tsv", sep="\t", index=False)


def test_summarise_pdui_cells_detects_shift(tmp_path):
    rng = np.random.default_rng(0)
    make_pdui_table(tmp_path / "p1", rng, {"G1": 0.3, "G2": -0.3, "G3": 0.0})
    raw = pd.read_csv(tmp_path / "p1" / "pdui_classic.tsv", sep="\t")
    s = rf.summarise_pdui_cells(raw, min_cells=10).set_index("feature_id")
    assert set(zip(s["c1"], s["c2"])) == {("A", "B")}
    assert s.loc["G1", "dpdui"] == pytest.approx(0.3, abs=0.05)
    assert s.loc["G2", "dpdui"] == pytest.approx(-0.3, abs=0.05)
    assert s.loc["G1", "qvalue"] < 1e-6 and s.loc["G2", "qvalue"] < 1e-6
    assert s.loc["G3", "qvalue"] > 0.05
    assert s.loc["G1", "n_cells_c1"] == 40          # NaN cell dropped
    assert s.loc["G1", "chrom"] == "1" and s.loc["G1", "strand"] == "+"


def test_summarise_pdui_min_cells_and_shuffle_null(tmp_path):
    rng = np.random.default_rng(1)
    make_pdui_table(tmp_path / "p1", rng, {"G1": 0.4}, n_cells=12)
    raw = pd.read_csv(tmp_path / "p1" / "pdui_classic.tsv", sep="\t")
    assert rf.summarise_pdui_cells(raw, min_cells=20).empty
    real = rf.summarise_pdui_cells(raw, min_cells=10).set_index("feature_id")
    assert real.loc["G1", "qvalue"] < 0.01
    null = rf.summarise_pdui_cells(raw, min_cells=10, rng=np.random.default_rng(5)).set_index("feature_id")
    assert abs(null.loc["G1", "dpdui"]) < abs(real.loc["G1", "dpdui"])


def test_pdui_cli_with_internal_null(tmp_path):
    rng = np.random.default_rng(2)
    # G1 lengthened in A in both samples; G2 opposite directions; G3 null
    make_pdui_table(tmp_path / "p1", rng, {"G1": 0.3, "G2": 0.3})
    make_pdui_table(tmp_path / "p2", rng, {"G1": 0.3, "G2": -0.3})
    out = tmp_path / "out"
    rc = _run_cli(["--sample", tmp_path / "p1", "--sample", tmp_path / "p2", "--out", out,
                   "--null-perms", 3, "--seed", 3])
    assert rc == 0
    summ = json.loads((out / "summary.json").read_text())
    assert summ["params"]["input_kind"] == "pdui"
    assert summ["params"]["effect_col"] == "dpdui"
    assert summ["null"]["mode"] == "internal_label_shuffle"
    assert summ["null"]["n_combos"] == 3
    t = pd.read_csv(out / "all_features.tsv", sep="\t").set_index("feature_id")
    assert bool(t.loc["G1", "passes_replication_floor"]) is True
    assert t.loc["G1", "utr_direction"] == "higher_in_c1"
    assert bool(t.loc["G2", "direction_consistent"]) is False
    assert bool(t.loc["G2", "passes_replication"]) is False
    assert bool(t.loc["G3", "passes_replication"]) is False
    assert summ["real"]["n_replicated"] == 1
    assert summ["null"]["replicated"]["null_max"] <= 1


def test_pdui_presummarised_table(tmp_path):
    rows = pd.DataFrame({"gene_id": ["G1", "G1", "G2"], "cluster1": ["B", "A", "A"],
                         "cluster2": ["A", "B", "B"], "dpdui": [-0.25, 0.25, 0.02],
                         "qvalue": [1e-3, 1e-3, 0.5]})
    (tmp_path / "p1").mkdir()
    (tmp_path / "p2").mkdir()
    rows.iloc[[0, 2]].to_csv(tmp_path / "p1" / "pdui_summary.tsv", sep="\t", index=False)
    rows.iloc[[1, 2]].to_csv(tmp_path / "p2" / "pdui_summary.tsv", sep="\t", index=False)
    out = tmp_path / "out"
    rc = _run_cli(["--sample", f"p1={tmp_path / 'p1' / 'pdui_summary.tsv'}",
                   "--sample", f"p2={tmp_path / 'p2' / 'pdui_summary.tsv'}", "--out", out])
    assert rc == 0
    t = pd.read_csv(out / "all_features.tsv", sep="\t").set_index("feature_id")
    # p1 wrote B_vs_A with -0.25 -> canonical (A,B) +0.25 ; replicates with p2
    assert t.loc["G1", "replication_count"] == 2
    assert t.loc["G1", "effect__p1"] == pytest.approx(0.25)
    assert bool(t.loc["G1", "passes_replication_floor"]) is True
    assert bool(t.loc["G2", "passes_replication"]) is False
    with pytest.raises((SystemExit, ValueError)):
        _run_cli(["--sample", f"p1={tmp_path / 'p1' / 'pdui_summary.tsv'}",
                  "--sample", f"p2={tmp_path / 'p2' / 'pdui_summary.tsv'}", "--out", out,
                  "--null-perms", 2])


def test_detect_input_kind(three_samples, tmp_path):
    assert rf.detect_input_kind(three_samples / "s1") == "diff"
    assert rf.detect_input_kind(three_samples / "s1" / "differential") == "diff"
    make_pdui_table(tmp_path / "p", np.random.default_rng(0), {})
    assert rf.detect_input_kind(tmp_path / "p") == "pdui"
    assert rf.detect_input_kind(tmp_path / "p" / "pdui_classic.tsv") == "pdui"
    with pytest.raises(FileNotFoundError):
        rf.detect_input_kind(tmp_path / "nope")


def test_wide_columns_when_a_sample_has_no_rows_in_a_pair(tmp_path):
    """Regression: a sample present in --sample but absent from the long table
    (no rows for any feature, e.g. a null perm that tested nothing) used to
    break the per-sample column assignment (length-1 scalar vs n rows)."""
    make_pair_tsv(tmp_path / "a" / "differential", "N", "T", {"1": (1e-3, 0.3), "2": (0.5, 0.1)})
    ldf, _, _ = rf.load_diff_sample("a", tmp_path / "a")
    t = rf.replication_table(ldf, rf.Params(), ["a", "ghost"]).set_index("feature_id")
    assert len(t) == 2
    assert t["q__ghost"].isna().all() and t["effect__ghost"].isna().all()
    assert t["called__ghost"].isna().all()
    assert bool(t.loc["1", "called__a"]) is True and bool(t.loc["2", "called__a"]) is False
    # a bool column must survive for a sample with only one tested feature
    make_pair_tsv(tmp_path / "b" / "differential", "N", "T", {"1": (1e-3, 0.3)})
    ldf2, _, _ = rf.load_diff_sample("b", tmp_path / "b")
    t2 = rf.replication_table(pd.concat([ldf, ldf2]), rf.Params(), ["a", "b"]).set_index("feature_id")
    assert bool(t2.loc["1", "called__b"]) is True
    assert pd.isna(t2.loc["2", "called__b"])
    assert t2.loc["1", "replication_count"] == 2


def test_shuffle_cell_labels_is_cell_level(tmp_path):
    """Each cell keeps ONE label across all its genes; the label multiset over
    cells is preserved; the per-gene split is not artificially balanced."""
    rng = np.random.default_rng(0)
    make_pdui_table(tmp_path / "p", rng, {"G1": 0.3}, n_cells=30)
    raw = pd.read_csv(tmp_path / "p" / "pdui_classic.tsv", sep="\t")
    raw = raw[raw["pdui"].notna()]                 # drop the fixture's NaN-PDUI cell
    # make the label split unbalanced: drop most B cells of every gene
    b_cells = sorted(raw.loc[raw["cluster"] == "B", "cell"].unique())[:25]
    raw = raw[~raw["cell"].isin(b_cells)]
    sh = rf.shuffle_cell_labels(raw, "cluster", np.random.default_rng(1))
    per_cell = sh.groupby("cell")["cluster"].nunique()
    assert (per_cell == 1).all()
    orig = raw.drop_duplicates("cell")["cluster"].value_counts().to_dict()
    new = sh.drop_duplicates("cell")["cluster"].value_counts().to_dict()
    assert orig == new == {"A": 30, "B": 5}
    assert (sh.drop_duplicates("cell").set_index("cell")["cluster"]
            != raw.drop_duplicates("cell").set_index("cell")["cluster"]).any()
    # per gene the split stays 30/5 (row-wise shuffling would have given ~17/18)
    g1 = sh[sh["gene_id"] == "G1"].drop_duplicates("cell")["cluster"].value_counts().to_dict()
    assert g1 == {"A": 30, "B": 5}
    with pytest.raises(ValueError):
        rf.shuffle_cell_labels(raw.drop(columns=["cell"]), "cluster", rng)


def test_null_universe_restriction(three_samples, tmp_path):
    """Null perms that replicate a feature ABSENT from the real run count only
    with --null-universe all; the default 'real' universe ignores them."""
    nroot = tmp_path / "nulls"
    for s in ("s1", "s2", "s3"):
        for k in range(2):
            d = nroot / s / f"perm_{k:02d}" / "differential"
            # pas 1 is in the real universe; GENE_C's pas 5 is too; but a
            # foreign pas id "99" is not -> inject it via a custom row
            make_pair_tsv(d, "T", "N", {"1": (0.6, 0.3), "5": (0.6, 0.2)})
            f = d / "fisher_T_vs_N.tsv"
            df = pd.read_csv(f, sep="\t", dtype=str)
            extra = df.iloc[[0]].copy()
            extra["pas_id"] = "99"; extra["qvalue"] = "0.001"; extra["delta_proportion"] = "0.5"
            pd.concat([df, extra]).to_csv(f, sep="\t", index=False)
    base = ["--sample", three_samples / "s1", "--sample", three_samples / "s2",
            "--sample", three_samples / "s3", "--null-dirs", f"{nroot}/{{sample}}/perm_*"]
    out_real = tmp_path / "o_real"
    assert _run_cli(base + ["--out", out_real]) == 0
    sr = json.loads((out_real / "summary.json").read_text())
    assert sr["null"]["universe"] == "real"
    assert sr["null"]["replicated"]["null_per_combo"] == [0, 0]
    assert sr["null"]["n_features_null_mean"] <= sr["null"]["n_features_real"]
    out_all = tmp_path / "o_all"
    assert _run_cli(base + ["--out", out_all, "--null-universe", "all"]) == 0
    sa = json.loads((out_all / "summary.json").read_text())
    assert sa["null"]["replicated"]["null_per_combo"] == [1, 1]
    assert sa["null"]["n_features_null_mean"] == 3.0
    # restrict_to_universe unit behaviour
    ldf, _, _ = rf.load_diff_sample("s1", three_samples / "s1")
    assert len(rf.restrict_to_universe(ldf, {("N", "T", "1")})) == 1
    assert rf.restrict_to_universe(ldf, set()).empty

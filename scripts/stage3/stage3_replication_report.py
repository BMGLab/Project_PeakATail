#!/usr/bin/env python3
"""
stage3_replication_report.py -- human-readable report over replication_filter.py
outputs for the Stage-3 Laughney cohort (manuscript/13 addendum item 2;
unit of replication = PATIENT).

Reads <rep_root>/{pas,gene}_{K2,K3}/{summary.json, all_features.tsv,
replicated.tsv, null_control.tsv} and writes
    <rep_root>/REPORT.md                         the report (PROVISIONAL banner)
    <rep_root>/per_pair_<cfg>.tsv                per cell-type pair table per config
    <rep_root>/top_gene_switches_<cfg>.tsv       ranked replicated gene-level switches
    <rep_root>/known_apa_hits_<cfg>.tsv          replicated genes on the APA watch-list
    <rep_root>/report_summary.json               machine-readable totals

Median delta-proportion per replicated feature = median of effect__<sample>
over the samples called in the consensus direction (for a two-GSM patient
both GSMs contribute; the count of PATIENTS is replication_count).

Every number is PROVISIONAL until the verifier passes.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

CONFIGS = ("pas_K2", "gene_K2", "pas_K3", "gene_K3")

# Watch-list of genes with well-documented APA / 3'UTR-isoform regulation
# (literature heuristic, NOT a validation set): Sandberg 2008 Science
# (proliferation shortening), Mayr & Bartel 2009 Cell (cancer shortening),
# Masamha 2014 Nature (CFIm25 targets), Gruber & Zavolan 2019 Nat Rev Genet,
# Takagaki 1996 Cell (IgM secreted/membrane), Tian & Manley 2017, Lianoglou
# 2013 Genes Dev (ubiquitous vs tissue-specific), Xia 2014 Nat Commun (DaPars),
# Singh 2018 (CD47, immune APA), Mitschka & Mayr 2022.
KNOWN_APA_GENES = {
    "IGHM": "IgM secreted vs membrane isoform (B cell -> plasma cell); Takagaki 1996",
    "IGHG1": "Ig heavy chain secreted/membrane poly(A) choice in plasma cells",
    "IGHA1": "Ig heavy chain secreted/membrane poly(A) choice in plasma cells",
    "CSTF2": "CstF-64, drives proximal PAS choice; Takagaki 1996",
    "NUDT21": "CFIm25, global 3'UTR length regulator; Masamha 2014",
    "CPSF6": "CFIm68, 3'UTR length regulator; Gruber 2012",
    "PABPN1": "suppresses proximal PAS; Jenal 2012",
    "CD47": "3'UTR isoform-dependent localisation; Berkovits & Mayr 2015",
    "CCND1": "cyclin D1 3'UTR shortening in cancer; Mayr & Bartel 2009",
    "CCND2": "3'UTR shortening target; Sandberg 2008",
    "DICER1": "3'UTR shortening; Mayr & Bartel 2009",
    "IGF2BP1": "IMP-1 3'UTR shortening in cancer; Mayr & Bartel 2009",
    "HMGA2": "let-7 escape via 3'UTR shortening; Mayr 2007",
    "TIMP2": "CFIm25 target; Masamha 2014",
    "PTEN": "APA-regulated tumour suppressor; Xia 2014",
    "MDM2": "3'UTR isoform; Xia 2014",
    "E2F1": "proliferation-linked 3'UTR shortening; Sandberg 2008",
    "CDC6": "proliferation-linked shortening; Sandberg 2008",
    "ELAVL1": "HuR, 3'UTR autoregulation; Dai 2012",
    "BCL2": "3'UTR isoform regulation; Lianoglou 2013",
    "TP53": "3'UTR isoform; Xia 2014",
    "RAC1": "3'UTR shortening in cancer; Xia 2014",
    "VMA21": "DaPars example; Xia 2014",
    "TNFRSF10B": "APA in immune signalling",
    "IL7R": "alternative 3' end usage in T cells",
    "PTPRC": "CD45 isoforms (splicing, reported APA-linked)",
    "CD5": "APA in T cells (Singh 2018 immune APA atlas)",
    "CD44": "3'UTR isoforms; immune APA",
    "MAP3K7": "immune APA; Singh 2018",
    "IRF7": "immune APA",
    "STAT3": "3'UTR shortening; Xia 2014",
    "POLR2K": "DaPars example",
    "SPATA13": "DaPars example",
    "FUS": "APA autoregulation",
    "CPEB1": "polyadenylation regulator",
    "CELF1": "APA regulator",
    "CPSF1": "core CPSF",
    "FIP1L1": "Fip1, stem-cell APA; Lackford 2014",
    "TRA2B": "immune APA; Singh 2018",
    "CHMP2A": "APA example",
}


def load_cfg(root: Path, cfg: str) -> dict | None:
    d = root / cfg
    if not (d / "summary.json").is_file():
        return None
    out = {"dir": d, "summary": json.loads((d / "summary.json").read_text())}
    for f in ("all_features", "replicated", "null_control"):
        p = d / f"{f}.tsv"
        out[f] = pd.read_csv(p, sep="\t", low_memory=False) if p.is_file() else None
    return out


def per_pair_table(cfg: dict) -> pd.DataFrame:
    s = cfg["summary"]
    rows = []
    nullpp = (s.get("null") or {}).get("per_pair", {}) or {}
    for pair, v in sorted(s["real"]["per_pair"].items()):
        nb = nullpp.get(pair, {})
        r = {"pair": pair, "n_patients_tested": v.get("n_groups_tested"),
             "n_gsm_tested": v.get("n_samples_tested"), "n_features": v["n_features"],
             "replicated_q": v["n_replicated"], "replicated_q_floor": v["n_replicated_floor"]}
        for key, tag in (("replicated", "q"), ("replicated_floor", "q_floor")):
            b = nb.get(key, {})
            r[f"null_mean_{tag}"] = b.get("null_mean")
            r[f"null_sd_{tag}"] = b.get("null_sd")
            r[f"null_max_{tag}"] = b.get("null_max")
            r[f"empirical_p_{tag}"] = b.get("empirical_p_ge_observed")
            r[f"exp_false_frac_{tag}"] = b.get("expected_false_replicated_fraction")
        rows.append(r)
    t = pd.DataFrame(rows)
    # totals row
    nul = s.get("null") or {}
    tot = {"pair": "ALL", "n_patients_tested": s["params"].get("n_groups"),
           "n_gsm_tested": s["n_samples"], "n_features": s["real"]["n_features"],
           "replicated_q": s["real"]["n_replicated"], "replicated_q_floor": s["real"]["n_replicated_floor"]}
    for key, tag in (("replicated", "q"), ("replicated_floor", "q_floor")):
        b = nul.get(key, {}) if nul.get("enabled") else {}
        tot[f"null_mean_{tag}"] = b.get("null_mean"); tot[f"null_sd_{tag}"] = b.get("null_sd")
        tot[f"null_max_{tag}"] = b.get("null_max")
        tot[f"empirical_p_{tag}"] = b.get("empirical_p_ge_observed")
        tot[f"exp_false_frac_{tag}"] = b.get("expected_false_replicated_fraction")
    return pd.concat([t, pd.DataFrame([tot])], ignore_index=True)


def median_consensus_effect(row: pd.Series, samples: list[str]) -> float:
    cons = row.get("consensus_direction")
    vals = []
    for smp in samples:
        c = row.get(f"called__{smp}")
        e = row.get(f"effect__{smp}")
        if c is True or c == "True" or c == 1:
            if pd.notna(e) and ((cons == "+" and e > 0) or (cons == "-" and e < 0)):
                vals.append(float(e))
    return float(np.median(vals)) if vals else float("nan")


def direction_text(row: pd.Series) -> str:
    """'+' on delta_proportion = higher PAS usage in c1 (prop(c1) - prop(c2))."""
    ud = row.get("utr_direction", "")
    if ud == "higher_in_c1":
        return f"PAS usage higher in {row['c1']}"
    if ud == "higher_in_c2":
        return f"PAS usage higher in {row['c2']}"
    return ""


def top_genes(cfg: dict, names: dict, n: int = 15) -> pd.DataFrame:
    rep = cfg["replicated"]
    if rep is None or rep.empty:
        return pd.DataFrame()
    samples = [s["name"] for s in cfg["summary"]["samples"]]
    t = rep.copy()
    t["gene_name"] = t["gene_id"].map(lambda g: names.get(str(g), ("", ""))[0])
    t["gene_biotype"] = t["gene_id"].map(lambda g: names.get(str(g), ("", ""))[1])
    t["pair"] = t["c1"] + "_vs_" + t["c2"]
    t["direction"] = t.apply(direction_text, axis=1)
    t["median_dprop_consensus"] = t.apply(lambda r: median_consensus_effect(r, samples), axis=1)
    t["on_apa_watchlist"] = t["gene_name"].isin(KNOWN_APA_GENES)
    t = t.sort_values(["replication_count", "replication_count_floor", "min_q_called"],
                      ascending=[False, False, True])
    cols = ["gene_id", "gene_name", "gene_biotype", "pair", "direction", "replication_count",
            "replication_count_floor", "n_groups_tested", "n_samples_tested", "median_dprop_consensus",
            "min_q_called", "passes_replication_floor"]
    if "representative_pas_id" in t.columns:
        cols += ["representative_pas_id", "n_pas_tested", "n_pas_replicated", "distal_pas_replicated",
                 "distal_utr_direction"]
    else:
        cols += ["pas_id", "chrom", "start", "end", "strand"]
    cols += ["on_apa_watchlist"]
    return t[cols].reset_index(drop=True), t


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rep-root", required=True)
    ap.add_argument("--gene-names", required=True, help="TSV gene_id gene_name biotype (no header)")
    ap.add_argument("--top-n", type=int, default=15)
    ap.add_argument("--title", default="Stage-3 Laughney replication filter (unit = patient)")
    args = ap.parse_args()
    root = Path(args.rep_root)
    names = {}
    with open(args.gene_names) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                names[p[0]] = (p[1], p[2] if len(p) > 2 else "")
    cfgs = {c: load_cfg(root, c) for c in CONFIGS}
    md = [f"# {args.title}", "",
          "**PROVISIONAL -- every number below is provisional until the verifier passes.**",
          f"Generated {time.strftime('%Y-%m-%d %H:%M:%S')} by stage3_replication_report.py from `{root}`.", ""]
    js = {"generated": time.strftime("%Y-%m-%dT%H:%M:%S"), "rep_root": str(root), "status": "PROVISIONAL",
          "configs": {}}
    for c, cfg in cfgs.items():
        if cfg is None:
            md += [f"## {c}: MISSING (no summary.json)", ""]
            js["configs"][c] = None
            continue
        s = cfg["summary"]
        prm = s["params"]
        nul = s.get("null") or {}
        md += [f"## {c}  (level={prm['level']}, K={prm['min_samples']} patients, fdr={prm['fdr']}, "
               f"floor={prm['effect_floor']}, discordant={prm['discordant_policy']}, unit={prm.get('unit_of_replication')})", ""]
        md += [f"- samples (GSMs): {s['n_samples']}; groups (patients): {prm.get('n_groups')}; "
               f"pairs seen: {len(s['pairs_seen'])}; (pair, feature) rows: {s['real']['n_features']}"]
        md += [f"- replicated q-only: **{s['real']['n_replicated']}**; q + |dprop|>=floor: **{s['real']['n_replicated_floor']}**"]
        if nul.get("enabled"):
            for key, tag in (("replicated", "q-only"), ("replicated_floor", "q+floor")):
                b = nul[key]
                md += [f"- null ({tag}): combos {nul['n_combos']}, per-combo {b['null_per_combo']}, "
                       f"mean {b['null_mean']:.2f} sd {b['null_sd']:.2f} max {b['null_max']}, "
                       f"empirical p(null >= obs) {b['empirical_p_ge_observed']:.3f}, "
                       f"expected_false_replicated_fraction {b['expected_false_replicated_fraction']}"]
            if nul.get("samples_without_null"):
                md += [f"- WARNING samples without null: {nul['samples_without_null']}"]
        else:
            md += ["- null: not run"]
        pp = per_pair_table(cfg)
        pp.to_csv(root / f"per_pair_{c}.tsv", sep="\t", index=False)
        md += ["", f"### Per cell-type pair ({c})", "",
               "| pair | patients tested | GSMs | features | repl. q | repl. q+floor | null mean q | emp. p q | null mean q+floor | emp. p q+floor |",
               "|---|---|---|---|---|---|---|---|---|---|"]
        for _, r in pp.iterrows():
            f = lambda v, d=2: ("" if v is None or (isinstance(v, float) and np.isnan(v)) else (f"{v:.{d}f}" if isinstance(v, float) else str(v)))
            md += [f"| {r['pair']} | {r['n_patients_tested']} | {r['n_gsm_tested']} | {r['n_features']} | {r['replicated_q']} | "
                   f"{r['replicated_q_floor']} | {f(r['null_mean_q'])} | {f(r['empirical_p_q'],3)} | "
                   f"{f(r['null_mean_q_floor'])} | {f(r['empirical_p_q_floor'],3)} |"]
        md += [""]
        js["configs"][c] = {"params": prm, "n_samples": s["n_samples"], "real": s["real"],
                            "null_totals": {k: nul.get(k) for k in ("n_combos", "replicated", "replicated_floor")} if nul.get("enabled") else None,
                            "per_sample_called": s.get("per_sample_called")}
        if prm["level"] == "gene":
            res = top_genes(cfg, names, args.top_n)
            if isinstance(res, tuple):
                top, full = res
                full_cols = top.columns
                full[full_cols].to_csv(root / f"top_gene_switches_{c}.tsv", sep="\t", index=False)
                hits = full[full["on_apa_watchlist"]][full_cols]
                hits.to_csv(root / f"known_apa_hits_{c}.tsv", sep="\t", index=False)
                md += [f"### Top {args.top_n} replicated gene-level switches ({c}; ranked by patients, then floor count, then min q)", "",
                       "| gene | pair | direction (representative PAS) | n patients | n patients (floor) | median dprop | min q | distal PAS replicated | distal UTR direction | watch-list |",
                       "|---|---|---|---|---|---|---|---|---|---|"]
                for _, r in top.head(args.top_n).iterrows():
                    gname = r["gene_name"] or r["gene_id"]
                    md += [f"| {gname} ({r['gene_id']}) | {r['pair']} | {r['direction']} | {r['replication_count']} | "
                           f"{r['replication_count_floor']} | {r['median_dprop_consensus']:.3f} | {r['min_q_called']:.2e} | "
                           f"{r['distal_pas_replicated']} | {r['distal_utr_direction'] or '-'} | {'yes' if r['on_apa_watchlist'] else ''} |"]
                md += ["", f"APA watch-list genes among ALL replicated gene-level switches ({c}): "
                       + (", ".join(sorted(set(f"{g} [{p}]" for g, p in zip(hits['gene_name'], hits['pair'])))) if len(hits) else "none"), ""]
                js["configs"][c]["top_genes"] = top.head(args.top_n).to_dict(orient="records")
                js["configs"][c]["apa_watchlist_hits"] = hits[["gene_name", "pair", "direction", "replication_count"]].to_dict(orient="records")
            else:
                md += ["(no replicated gene-level switches)", ""]
    (root / "REPORT.md").write_text("\n".join(md) + "\n")
    (root / "report_summary.json").write_text(json.dumps(js, indent=2, default=lambda o: None if (isinstance(o, float) and np.isnan(o)) else (o.item() if hasattr(o, "item") else str(o))))
    print("\n".join(md))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

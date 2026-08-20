#!/usr/bin/env python3
"""
Unit tests for scripts/reliability/trusted_novel_pas.py on a SYNTHETIC genome.

The genome background is C/G only, so no poly(A) hexamer (all contain >= 3 A) and no A-run
can occur by accident on either strand; every signal is implanted explicitly.  Minus-strand
implants are written as the reverse complement on the forward strand, which is exactly the
case that an asymmetric getfasta window gets wrong (motif_validation.py).

Run:
    /mnt/ssd1/Projects/PeakATail_wd/tools/PeakATail/.venv/bin/python -m pytest \
        /mnt/ssd1/Projects/PeakATail_wd/scripts/reliability/tests -v
"""
import csv
import importlib.util
import json
import random
import shutil
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
SCRIPT = HERE.parent / "trusted_novel_pas.py"

pytestmark = pytest.mark.skipif(shutil.which("bedtools") is None, reason="bedtools not on PATH")


def _load():
    spec = importlib.util.spec_from_file_location("_tnp_under_test", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["_tnp_under_test"] = mod
    spec.loader.exec_module(mod)
    return mod


tnp = _load()

RC = {"A": "T", "C": "G", "G": "C", "T": "A"}


def revcomp(s):
    return "".join(RC[c] for c in reversed(s))


def write_fasta_with_fai(path: Path, contigs):
    """Write a 60-column FASTA and a hand-built .fai (name, len, offset, linebases, linewidth)."""
    lines = []
    fai = []
    offset = 0
    for name, seq in contigs:
        header = f">{name}\n"
        offset += len(header)
        fai.append(f"{name}\t{len(seq)}\t{offset}\t60\t61")
        wrapped = [seq[i:i + 60] + "\n" for i in range(0, len(seq), 60)]
        lines.append(header)
        lines.extend(wrapped)
        offset += sum(len(w) for w in wrapped)
    path.write_text("".join(lines))
    Path(str(path) + ".fai").write_text("\n".join(fai) + "\n")


def implant(seq: list, start: int, motif: str):
    for i, ch in enumerate(motif):
        seq[start + i] = ch


@pytest.fixture(scope="module")
def synth(tmp_path_factory):
    """Synthetic genome + PAS/clip/atlas/truth BEDs with known expected outcomes."""
    d = tmp_path_factory.mktemp("synth")
    rng = random.Random(7)
    c1 = [rng.choice("CG") for _ in range(3000)]
    c2 = [rng.choice("CG") for _ in range(500)]
    # --- '+' sites: hexamer occupying r = -21..-16  -> forward [c-21, c-15)
    implant(c1, 500 - 21, "AATAAA")               # P1 strong (atlas site at 450 '+', 50 bp)
    implant(c1, 700 - 30, "ATTAAA")               # P2 strong hexamer ...
    implant(c1, 700 + 12, "AAAAAA")               # ... but A-run at r=+12..+17 -> IP
    implant(c1, 1300 - 21, "AGTAAA")              # W1 weak-tier hexamer only
    #                                              N1 (1500): nothing implanted
    implant(c1, 1700 - 21, "AATAAA")              # S1 strong but support 1
    implant(c1, 1900 - 21, "AATAAA")              # C1 strong but not clip-supported
    implant(c2, 250 - 21, "AATAAA")               # O1 strong; atlas site at c2:240 on '-'
    # --- '-' sites: r = c - g; hexamer r = -21..-16 -> g = c+16..c+21 -> forward [c+16, c+22)
    implant(c1, 900 + 16, revcomp("AATAAA"))      # M1 strong   (forward TTTATT)
    implant(c1, 1100 + 16, revcomp("AATAAA"))     # M2 strong hexamer ...
    implant(c1, 1100 - 25, "T" * 7)               # ... T-run forward [1075,1082) = A-run at r=+19..+25 -> IP
    fasta = d / "genome.fa"
    write_fasta_with_fai(fasta, [("c1", "".join(c1)), ("c2", "".join(c2))])

    pas = [  # chrom, cleavage, name, support, strand
        ("c1", 500, "P1", 5, "+"), ("c1", 700, "P2", 5, "+"), ("c1", 900, "M1", 5, "-"),
        ("c1", 1100, "M2", 5, "-"), ("c1", 1300, "W1", 5, "+"), ("c1", 1500, "N1", 5, "+"),
        ("c1", 1700, "S1", 1, "+"), ("c1", 1900, "C1", 5, "+"), ("c1", 10, "E1", 5, "+"),
        ("c9", 100, "X1", 5, "+"), ("c2", 250, "O1", 5, "+"), ("c1", 900, "PM", 5, "+"),
    ]
    pas_bed = d / "pas.bed"
    pas_bed.write_text("".join(f"{c}\t{p}\t{p + 1}\t{n}\t{s}\t{st}\n" for c, p, n, s, st in pas))
    clip_bed = d / "clip.bed"
    clip_bed.write_text("".join(f"{c}\t{p}\t{p + 1}\t{n}\t{s}\t{st}\n" for c, p, n, s, st in pas if n != "C1"))
    atlas_bed = d / "atlas.bed6"
    atlas_bed.write_text("c1\t450\t451\ta1\t1\t+\nc1\t2500\t2501\ta3\t1\t+\nc2\t240\t241\ta2\t1\t-\n")
    sizes = d / "chrom.sizes"
    sizes.write_text("c1\t3000\nc2\t500\n")
    truth = d / "truth.bed"   # M1 exact, W1 at +10, O1 at +30 (outside 25)
    truth.write_text("c1\t900\t901\tt1\t50\t-\nc1\t1310\t1311\tt2\t50\t+\nc2\t280\t281\tt3\t50\t+\n")
    truth_flip = d / "truth_flip.bed"   # M1 position but on '+': strand-aware miss
    truth_flip.write_text("c1\t900\t901\tt1\t50\t+\n")
    return dict(dir=d, fasta=fasta, pas=pas_bed, clip=clip_bed, atlas=atlas_bed, sizes=sizes,
                truth=truth, truth_flip=truth_flip)


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ----------------------------------------------------------------------------- pure functions
def test_hexamer_window_boundaries_fully_inside():
    W = 40
    base = ["C"] * (2 * W + 1)

    def with_hex_at(r_start, hexamer="AATAAA"):
        s = base[:]
        for i, ch in enumerate(hexamer):
            s[r_start + W + i] = ch
        return "".join(s)

    assert tnp.hexamer_call(with_hex_at(-21), W)[0] is True
    assert tnp.hexamer_call(with_hex_at(-40), W)[0] is True          # starts exactly at -40
    assert tnp.hexamer_call(with_hex_at(-10), W)[0] is True          # ends exactly at -5
    assert tnp.hexamer_call(with_hex_at(-41), W)[0] is False         # one base too far upstream
    assert tnp.hexamer_call(with_hex_at(-9), W)[0] is False          # ends at -4: outside
    assert tnp.hexamer_call(with_hex_at(+5), W)[0] is False          # downstream never counts
    strong, any12, hits = tnp.hexamer_call(with_hex_at(-21, "AGTAAA"), W)
    assert (strong, any12, hits) == (False, True, ["AGTAAA"])
    assert tnp.hexamer_call(with_hex_at(-21, "ATTAAA"), W)[0] is True


def test_ip_call_run_and_fraction():
    W = 40
    base = ["C"] * (2 * W + 1)

    def with_downstream(s_rel, motif):
        s = base[:]
        for i, ch in enumerate(motif):
            s[s_rel + W + i] = ch
        return "".join(s)

    assert tnp.ip_call(with_downstream(12, "A" * 6), W) is True
    assert tnp.ip_call(with_downstream(12, "A" * 5), W) is False
    assert tnp.ip_call(with_downstream(25, "A" * 6), W) is True      # ends exactly at +30
    assert tnp.ip_call(with_downstream(26, "A" * 6), W) is False     # spills past +30
    assert tnp.ip_call(with_downstream(-5, "A" * 6), W) is False     # upstream A-run is not IP
    # 70 % rule on the 21-nt window: 15/21 = 0.714 flags, 14/21 = 0.667 does not
    assert tnp.ip_call(with_downstream(10, "A" * 5 + "C" + "A" * 5 + "C" + "A" * 5 + "CCCC"), W) is True
    assert tnp.ip_call(with_downstream(10, "A" * 5 + "C" + "A" * 4 + "C" + "A" * 5 + "CCCCC"), W) is False


def test_read_bed6_rejects_bad_strand(tmp_path):
    p = tmp_path / "bad.bed"
    p.write_text("c1\t10\t11\tx\t1\t.\n")
    with pytest.raises(ValueError):
        tnp.read_bed6(p)


def test_wilson_interval():
    lo, hi = tnp.wilson(70, 100)
    assert 0.60 < lo < 0.70 < hi < 0.79
    assert tnp.wilson(0, 0)[0] != tnp.wilson(0, 0)[0]  # nan


# ----------------------------------------------------------------------------- bedtools seams
def test_fetch_windows_strand_mapping_and_edges(synth, tmp_path):
    sites = tnp.read_bed6(synth["pas"])
    lens = tnp.read_fai(synth["fasta"])
    W = 40
    seqs = tnp.fetch_windows(sites, synth["fasta"], lens, W, tmp_path)
    by = {s.name: s for s in sites}
    # minus-strand site: reading the transcript strand, AATAAA sits at r=-21..-16 -> i=19..24
    assert seqs[by["M1"].uid][19:25] == "AATAAA"
    # plus-strand site at the SAME base sees the forward strand: TTTATT at r=+16..+21
    assert seqs[by["PM"].uid][56:62] == "TTTATT"
    assert "AATAAA" not in tnp.rel_slice(seqs[by["PM"].uid], W, -40, -5)
    # plus-strand hexamer at r=-21
    assert seqs[by["P1"].uid][19:25] == "AATAAA"
    # M2: forward T-run becomes an A-run at r=+19..+25 on the transcript strand
    assert tnp.rel_slice(seqs[by["M2"].uid], W, 19, 25) == "A" * 7
    # edge window and contig absent from the FASTA are not returned
    assert by["E1"].uid not in seqs
    assert by["X1"].uid not in seqs
    assert all(len(s) == 2 * W + 1 for s in seqs.values())


def test_closest_distance_point_semantics(tmp_path):
    q = tmp_path / "q.bed"
    q.write_text("c1\t10\t11\tq\t0\t+\n")
    r = tmp_path / "r.bed"
    r.write_text("c1\t10\t11\tx\t0\t-\nc1\t20\t21\ty\t0\t+\n")
    d = tnp.closest_distances(q, r, tmp_path, "t", strand_aware=True)
    assert d["q"] == 10                                   # |10-20| on the same strand
    d2 = tnp.closest_distances(q, r, tmp_path, "t2", strand_aware=False)
    assert d2["q"] == 0                                   # opposite-strand site at the same base
    r2 = tmp_path / "r2.bed"
    r2.write_text("c1\t10\t11\tx\t0\t-\n")
    assert tnp.closest_distances(q, r2, tmp_path, "t3", strand_aware=True)["q"] is None


# ----------------------------------------------------------------------------- CLI: call
def _call(synth, out, extra=()):
    args = ["call", "--pas-bed", str(synth["pas"]), "--clip-bed", str(synth["clip"]),
            "--support-unit", "molecules", "--genome", str(synth["fasta"]),
            "--atlas", str(synth["atlas"]), "--outdir", str(out), *extra]
    assert tnp.main(args) == 0
    return {r["stage"]: r for r in read_tsv(out / "funnel.tsv")}


def test_call_funnel_counts_and_membership(synth, tmp_path):
    out = tmp_path / "call"
    funnel = _call(synth, out)
    assert int(funnel["input"]["n_pass"]) == 12
    assert int(funnel["clip_supported"]["n_pass"]) == 11            # C1
    assert int(funnel["support_ge2_molecules"]["n_pass"]) == 10     # S1
    assert int(funnel["sequence_testable"]["n_pass"]) == 8          # E1 edge, X1 contig missing
    assert int(funnel["not_internal_priming"]["n_pass"]) == 6       # P2 (+), M2 (-)
    assert int(funnel["hexamer_any12_-40..-5"]["n_pass"]) == 4      # N1, PM
    assert int(funnel["atlas_novel_ge100bp"]["n_pass"]) == 3        # P1 (50 bp from atlas, same strand)
    assert int(funnel["trusted_novel"]["n_pass"]) == 3
    assert int(funnel["trusted_novel_strong"]["n_pass"]) == 2       # W1 is weak tier
    trusted = [l.split("\t") for l in (out / "trusted_novel.bed").read_text().splitlines() if not l.startswith("#")]
    assert sorted(r[3] for r in trusted) == ["M1", "O1", "W1"]
    tiers = {r[3]: r[8] for r in trusted}
    assert tiers == {"M1": "strong", "O1": "strong", "W1": "weak"}
    strong = [l.split("\t")[3] for l in (out / "trusted_novel.strong.bed").read_text().splitlines() if not l.startswith("#")]
    assert sorted(strong) == ["M1", "O1"]
    ann = {r["name"]: r for r in read_tsv(out / "site_annotations.tsv")}
    assert ann["P2"]["fail_stage"] == "not_internal_priming" and ann["P2"]["ip_flag"] == "1"
    assert ann["M2"]["fail_stage"] == "not_internal_priming" and ann["M2"]["ip_flag"] == "1"
    assert ann["PM"]["fail_stage"] == "hexamer_any12_-40..-5"
    assert ann["P1"]["atlas_dist_bp"] == "50" and ann["P1"]["fail_stage"] == "atlas_novel_ge100bp"
    assert ann["O1"]["atlas_dist_bp"] == "" and ann["O1"]["atlas_novel"] == "1"   # opposite strand ignored
    assert ann["E1"]["ip_source"] == "untestable" and ann["X1"]["ip_source"] == "untestable"
    assert ann["M1"]["hex_hits"] == "AATAAA" and ann["W1"]["hex_hits"] == "AGTAAA"
    stage_files = sorted(p.name for p in (out / "stages").glob("*.bed"))
    assert stage_files[0] == "00_input.bed" and stage_files[-1].endswith("trusted_novel_strong.bed")
    assert len(stage_files) == 9
    m = json.loads((out / "run_manifest.json").read_text())
    assert m["n_trusted_novel"] == 3 and m["stale_input_dry_run"] is False
    # stage BEDs are cleavage-point BED6 carrying the original names
    tn = [p for p in (out / "stages").glob("*_trusted_novel.bed")]
    assert len(tn) == 1 and tn[0].name == "07_trusted_novel.bed"   # 6 criteria stages + input, no contig stage
    line = tn[0].read_text().splitlines()[0].split("\t")
    assert len(line) == 6 and int(line[2]) - int(line[1]) == 1


def test_call_atlas_ignore_strand_drops_o1(synth, tmp_path):
    funnel = _call(synth, tmp_path / "call_is", ["--atlas-ignore-strand"])
    assert int(funnel["atlas_novel_ge100bp"]["n_pass"]) == 2
    assert int(funnel["trusted_novel_strong"]["n_pass"]) == 1


def test_call_restrict_chroms_stage(synth, tmp_path):
    funnel = _call(synth, tmp_path / "call_rc", ["--restrict-chroms", str(synth["sizes"])])
    assert int(funnel["on_listed_contigs"]["n_pass"]) == 11        # X1 on c9 dropped here
    assert int(funnel["sequence_testable"]["n_pass"]) == 8
    assert int(funnel["trusted_novel"]["n_pass"]) == 3


def test_call_support_unit_reads_warns_and_marks_dry_run(synth, tmp_path, capsys):
    out = tmp_path / "call_reads"
    args = ["call", "--pas-bed", str(synth["pas"]), "--clip-bed", str(synth["clip"]),
            "--support-unit", "reads", "--genome", str(synth["fasta"]),
            "--atlas", str(synth["atlas"]), "--outdir", str(out)]
    assert tnp.main(args) == 0
    err = capsys.readouterr().err
    assert "READS, not distinct molecules" in err
    m = json.loads((out / "run_manifest.json").read_text())
    assert m["stale_input_dry_run"] is True and any("READS" in w for w in m["warnings"])


def test_call_refuses_hg19_atlas(synth, tmp_path):
    bad = tmp_path / "polyadb.hg19.bed6"
    bad.write_text(synth["atlas"].read_text())
    with pytest.raises(SystemExit):
        _call(synth, tmp_path / "call_hg19", ["--atlas", str(bad)])
    _call(synth, tmp_path / "call_hg19_ok", ["--atlas", str(bad), "--allow-unsafe-atlas"])


def test_call_ip_from_annot_overrides_and_falls_back(synth, tmp_path):
    annot = tmp_path / "annotatedpas.bed"
    rows = []
    for line in synth["pas"].read_text().splitlines():
        f = line.split("\t")
        flag = {"M1": "True", "P2": "False"}.get(f[3], "")   # M1 forced IP, P2 forced clean, rest blank
        rows.append("\t".join(f + ["gene", "", "", flag, ""]))
    annot.write_text("\n".join(rows) + "\n")
    out = tmp_path / "call_annot"
    funnel = _call(synth, out, ["--ip-from-annot", f"{annot}:10"])
    ann = {r["name"]: r for r in read_tsv(out / "site_annotations.tsv")}
    assert ann["M1"]["ip_source"] == "annot" and ann["M1"]["fail_stage"] == "not_internal_priming"
    assert ann["P2"]["ip_source"] == "annot" and ann["P2"]["ip_flag"] == "0"
    assert ann["M2"]["ip_source"] == "recomputed" and ann["M2"]["ip_flag"] == "1"
    # M1 is lost (annot says IP); P2 is forced clean, has ATTAAA and no atlas site within
    # 100 bp, so it joins: the set is O1, P2, W1 (still 3, but a different 3)
    assert int(funnel["trusted_novel"]["n_pass"]) == 3
    trusted = sorted(l.split("\t")[3] for l in (out / "trusted_novel.bed").read_text().splitlines() if not l.startswith("#"))
    assert trusted == ["O1", "P2", "W1"]


# ----------------------------------------------------------------------------- CLI: validate
def test_validate_observed_fraction_and_null(synth, tmp_path):
    call_out = tmp_path / "call_v"
    _call(synth, call_out)
    out = tmp_path / "validate"
    args = ["validate", "--query", f"trusted={call_out / 'trusted_novel.bed'}",
            "--query", f"all={call_out / 'stages' / '00_input.bed'}",
            "--truth", f"t={synth['truth']}", "--truth", f"flip={synth['truth_flip']}",
            "--chrom-sizes", str(synth["sizes"]), "--seeds", "4", "--outdir", str(out)]
    assert tnp.main(args) == 0
    rows = {(r["query"], r["truth"]): r for r in read_tsv(out / "validation.tsv")}
    r = rows[("trusted", "t")]
    assert int(r["n_query"]) == 3 and int(r["n_hit"]) == 2            # M1 exact, W1 +10, O1 +30 misses
    assert abs(float(r["frac_hit"]) - 2 / 3) < 1e-6
    assert float(r["wilson95_lo"]) < 2 / 3 < float(r["wilson95_hi"])
    assert int(r["null_seeds"]) == 4
    assert float(r["null_mean"]) < float(r["frac_hit"])
    assert int(rows[("trusted", "flip")]["n_hit"]) == 0               # strand-aware miss
    # the 'all' query drops X1 (c9 not in chrom sizes): 11 scored
    assert int(rows[("all", "t")]["n_query"]) == 11
    seeds = read_tsv(out / "validation_null_seeds.tsv")
    assert len(seeds) == 2 * 2 * 4
    m = json.loads((out / "run_manifest.json").read_text())
    assert m["params"]["seeds"] == [1, 2, 3, 4]


def test_validate_ignore_strand_hits_flipped_truth(synth, tmp_path):
    call_out = tmp_path / "call_v2"
    _call(synth, call_out)
    out = tmp_path / "validate_is"
    args = ["validate", "--query", f"trusted={call_out / 'trusted_novel.bed'}",
            "--truth", f"flip={synth['truth_flip']}", "--chrom-sizes", str(synth["sizes"]),
            "--seeds", "2", "--ignore-strand", "--outdir", str(out)]
    assert tnp.main(args) == 0
    r = read_tsv(out / "validation.tsv")[0]
    assert int(r["n_hit"]) == 1 and r["strand_aware"] == "0"

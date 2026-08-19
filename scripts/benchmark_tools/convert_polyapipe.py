#!/usr/bin/env python3
"""
convert_polyapipe.py -- polyApipe *_polyA_peaks.gff -> standard PAS point BED6.

STANDARD FORM (shared by every benchmarked tool):
  BED6, one line per inferred PAS, the interval is the single 1-bp cleavage
  site ([point-1, point) in 0-based half-open), name = tool's peak id,
  score = tool's depth/support, strand = +/-.

GEOMETRY DERIVATION (verified in polyApipe source, polyApipe.py 0.1.0,
process_polyA_ends_to_peaks, lines ~331-457):
  * forward ('+') read: 3' end = read.reference_end-1 (0-based); region stored
    as [pos-region_size, pos]; at GFF write the END becomes pos+1 (1-based) and
    the peak name uses the GFF END  => PAS point = GFF end for strand '+'/'f'.
  * reverse ('-') read: 3' end = read.reference_start (0-based); region stored
    as [pos, pos+region_size]; at GFF write the START becomes pos+1 (1-based)
    and the peak name uses the GFF START => PAS point = GFF start for '-'/'r'.
  The 250 bp window is UPSTREAM (5') of the PAS on the mRNA in both cases; the
  peak name chrom_pos_f|r embeds the 1-based PAS coordinate, which we assert
  against the strand-derived value for every feature.

Outputs (same directory as the input GFF unless --outdir):
  pas_all.bed : every peak, including misprime="True"
  pas.bed     : DEFAULT for scoring -- excludes misprime="True"
Both carry the derivation as leading '#' comment lines.
"""
import argparse
import re
import sys
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gff")
    ap.add_argument("--outdir", default=None)
    a = ap.parse_args()
    gff = Path(a.gff)
    outdir = Path(a.outdir) if a.outdir else gff.parent

    rx = re.compile(r'peak="([^"]+)"; peakdepth="(\d+)"; misprime="(True|False)";')
    rows = []
    n_bad = 0
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, start, end, strand, attr = f[0], int(f[3]), int(f[4]), f[6], f[8]
            m = rx.search(attr)
            if not m:
                n_bad += 1
                continue
            name, depth, misprime = m.group(1), int(m.group(2)), m.group(3) == "True"
            # PAS point, 1-based: '+' -> GFF end ; '-' -> GFF start
            pas1 = end if strand == "+" else start
            # assert the name embeds the same coordinate
            want = f"{chrom}_{pas1}_{'f' if strand == '+' else 'r'}"
            assert name == want, f"name/geometry mismatch: {name} vs {want}"
            rows.append((chrom, pas1 - 1, pas1, name, depth, strand, misprime))
    assert n_bad == 0, f"{n_bad} unparseable feature lines"
    rows.sort(key=lambda r: (r[0], r[1]))

    hdr = [
        "# polyApipe inferred PAS points (BED6, 0-based half-open, 1-bp intervals)",
        f"# source: {gff}",
        "# derivation (verified in polyApipe.py 0.1.0 source, process_polyA_ends_to_peaks):",
        "#   GFF feature = 250 bp window UPSTREAM (5' on mRNA) of the polyA site.",
        "#   strand '+' ('f'): PAS = GFF end   (window end;   name embeds this coord)",
        "#   strand '-' ('r'): PAS = GFF start (window start; name embeds this coord)",
        "#   name chrom_pos_f|r re-derived and asserted for every feature.",
        "# score column = polyApipe peakdepth.",
    ]
    variants = {
        "pas_all.bed": (rows, "# variant: ALL peaks incl. misprime=True"),
        "pas.bed": ([r for r in rows if not r[6]],
                    "# variant: DEFAULT -- misprime=True EXCLUDED"),
    }
    for fn, (rs, note) in variants.items():
        p = outdir / fn
        with open(p, "w") as out:
            out.write("\n".join(hdr + [note]) + "\n")
            for r in rs:
                out.write("\t".join(map(str, r[:6])) + "\n")
        print(f"wrote {p}  n={len(rs)}")

if __name__ == "__main__":
    sys.exit(main())

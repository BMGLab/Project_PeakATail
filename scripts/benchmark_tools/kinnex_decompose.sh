#!/usr/bin/env bash
# kinnex_decompose.sh <tag> <kinnex_truth_dir> <pos_counts.tsv>
# Where does each tool's PAS call actually land, given long-read evidence?
#   genuine_pas   : within 100 bp (strand-matched) of a >=5-UMI NON-internal-priming terminus peak
#   ip_artifact   : not that, but within 100 bp of a >=5-UMI INTERNAL-PRIMING terminus peak
#   no_lr_support : neither  (no reproducible long-read 3' end nearby at all)
# Also emits the sequence-based truth calibration (canonical hexamer vs molecular support,
# with the internal-priming peaks as a depth-matched negative control).
set -euo pipefail
export LC_ALL=C
TAG=$1; G=$2; POSC=$3
WD=/mnt/ssd1/Projects/PeakATail_wd
GEN=$WD/data/references/chrom.sizes.nochr.filt
REF=/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.dna.primary_assembly.fa
PAS2=$WD/data/references/atlases/polyasite2.GRCh38.96.rep_sites.bed6
S=/mnt/ssd0/emaout/peakatail_benchmark/kinnex/scratch/${TAG}_dec; mkdir -p $S
awk -F'\t' 'NR==FNR{ok[$1]=1;next} ok[$1]' $GEN $G/${TAG}_truth_t5.point.bed | sort -k1,1 -k2,2n > $S/truth5.bed
awk -F'\t' 'NR==FNR{ok[$1]=1;next} ok[$1]' $GEN $G/${TAG}_decoy.point.bed    | sort -k1,1 -k2,2n > $S/decoy.bed
sort -k1,1 -k2,2n $PAS2 > $S/pas2.bed

OUT=$G/${TAG}_call_support_decomposition.tsv
{ echo -e "# Kinnex $TAG long-read support decomposition of pbmc_10k_v3 PAS calls (100 bp, strand-matched)"
  echo -e "# CAVEAT: different donor and 10x chemistry from pbmc_10k_v3 -- site-level truth, not cell-matched."
  echo -e "tool\tn_calls\tgenuine_pas\tip_artifact\tno_lr_support\tip_share_of_supported"; } > $OUT
for t in scutrquant scapture polyapipe sierra scapatrap peakatail; do
  B=$WD/results/benchmark_tools/pbmc_10k_v3/$t/pas.bed
  [ -s "$B" ] || continue
  grep -v '^#' $B | awk -F'\t' 'NR==FNR{ok[$1]=1;next} ok[$1]' $GEN - \
   | awk -F'\t' -v OFS='\t' '{if($6=="+"){s=$3-1;e=$3}else{s=$2;e=$2+1} print $1,s,e,NR,0,$6}' \
   | sort -k1,1 -k2,2n > $S/q.bed
  paste <(bedtools closest -s -d -t first -a $S/q.bed -b $S/truth5.bed 2>/dev/null | awk -F'\t' '{print $NF}') \
        <(bedtools closest -s -d -t first -a $S/q.bed -b $S/decoy.bed  2>/dev/null | awk -F'\t' '{print $NF}') \
   | awk -v T=$t -v OFS='\t' '{n++; g=($1>=0&&$1<=100); i=($2>=0&&$2<=100)
        if(g) ng++; else if(i) ni++; else nn++}
      END{printf "%s\t%d\t%.6f\t%.6f\t%.6f\t%.6f\n", T,n,ng/n,ni/n,nn/n,(ng+ni?ni/(ng+ni):0)}' >> $OUT
done
echo "wrote $OUT"; cat $OUT

QC=$G/${TAG}_truth_calibration.tsv
{ echo -e "# Sequence-only calibration of the long-read truth set (NO atlas used to choose thresholds)."
  echo -e "# hexamer = fraction of representative sites with AATAAA or ATTAAA in -40..-10."
  echo -e "# For scale: the SAME code scores PolyASite 2.0 rep sites at 0.4111 and protein-coding TES at 0.3799."
  echo -e "set\tumi_bin\tn_pas\thexamer_canonical"; } > $QC
python3 - "$REF" "$G/${TAG}_truth_pas.bed" "$G/${TAG}_decoy_pas.bed" >> $QC <<'PY'
import sys
gp, tb, db = sys.argv[1:4]
seqs={}; name=None; buf=[]
for line in open(gp,'rb'):
    if line[:1]==b'>':
        if name: seqs[name]=b''.join(buf).upper()
        name=line[1:].split()[0].decode(); buf=[]
    else: buf.append(line.rstrip())
if name: seqs[name]=b''.join(buf).upper()
BINS=[(5,9),(10,19),(20,49),(50,99),(100,199),(200,499),(500,999),(1000,4999),(5000,10**9)]
TR=bytes.maketrans(b'ACGT',b'TGCA')
for label, path in (("truth_non_internal_priming", tb), ("decoy_internal_priming", db)):
    agg={b:[0,0] for b in BINS}
    for line in open(path):
        f=line.rstrip('\n').split('\t')
        c=f[0]; s=f[5]; u=int(f[4]); p=int(f[6])
        g=seqs.get(c)
        if g is None: continue
        w = g[max(0,p-40):max(0,p-9)] if s=='+' else g[p+10:p+41].translate(TR)[::-1]
        can = 1 if (b'AATAAA' in w or b'ATTAAA' in w) else 0
        for b in BINS:
            if b[0]<=u<=b[1]: agg[b][0]+=1; agg[b][1]+=can; break
    for b in BINS:
        n,c=agg[b]
        if n: print(f"{label}\t{b[0]}-{'inf' if b[1]>10**8 else b[1]}\t{n}\t{c/n:.6f}")
PY
echo "wrote $QC"; cat $QC

TH=$G/${TAG}_truth_threshold_qc.tsv
{ echo -e "# Truth-set size and independent QC at each molecular-support threshold."
  echo -e "umi_threshold\tn_pas\tfrac_within25bp_polyasite2\tfrac_within100bp_polyasite2"; } > $TH
for T in 5 20 100 500; do
  B=$G/${TAG}_truth_t${T}.point.bed; N=$(wc -l < $B)
  bedtools closest -s -d -t first -a <(sort -k1,1 -k2,2n $B) -b $S/pas2.bed 2>/dev/null \
   | awk -F'\t' -v N=$N -v T=$T -v OFS='\t' '{d=$NF; if(d>=0&&d<=25)c++; if(d>=0&&d<=100)e++}
       END{printf "%d\t%d\t%.6f\t%.6f\n", T,N,c/N,e/N}' >> $TH
done
echo "wrote $TH"; cat $TH

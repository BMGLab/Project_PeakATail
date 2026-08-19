# Per-molecule 3' terminus extraction from a dedup FLNC alignment (SAM on stdin).
# Emits: chrom <TAB> terminus0 <TAB> strand <TAB> CB <TAB> UMI   (1 row = 1 UMI/molecule)
# Terminus = cleavage site = reference END for '+' alignments, reference START for '-'.
# (isoseq refine already trimmed the polyA tail, so the read 3' end IS the cleavage junction.)
# vars: WL=real-cell barcode list (one per line, whitelist orientation), CLIP=max 3' softclip,
#       STRIP=1 to convert chr-prefixed GENCODE names to Ensembl no-chr.
BEGIN{
  FS="\t"; OFS="\t"
  if(CLIP=="") CLIP=30
  while((getline ln < WL) > 0){ if(ln!="") wl[ln]=1; nwl++ }
  if(STRIP==1){ for(i=1;i<=22;i++) keep[i]=1; keep["X"]=1; keep["Y"]=1; keep["MT"]=1 }
}
{
  n_in++
  flag=$2+0
  chrom=$3; pos=$4+0; cig=$6
  if(STRIP==1){
    if(substr(chrom,1,3)!="chr") { n_offchrom++; next }
    chrom=substr(chrom,4); if(chrom=="M") chrom="MT"
    if(!(chrom in keep)) { n_offchrom++; next }
  }
  # --- CIGAR: reference span + terminal clips ---
  nf=split(cig, ops, /[0-9]+/, nums)   # ops[i+1] = operator for length nums[i]
  reflen=0; leadclip=0; trailclip=0
  for(i=1;i<nf;i++){
    L=nums[i]+0; O=ops[i+1]
    if(O=="M"||O=="D"||O=="N"||O=="="||O=="X") reflen+=L
    if((O=="S"||O=="H")){ if(i==1) leadclip+=L; if(i==nf-1) trailclip+=L }
  }
  # --- strand + terminus ---
  if(and(flag,16)){ strand="-"; term1=pos;              clip3=leadclip  }
  else            { strand="+"; term1=pos+reflen-1;     clip3=trailclip }
  # --- tags ---
  cb=""; umi=""; rc="NA"
  for(i=12;i<=NF;i++){
    t=substr($i,1,5)
    if(t=="CB:Z:") cb=substr($i,6)
    else if(t=="XM:Z:") umi=substr($i,6)
    else if(substr($i,1,5)=="rc:i:") rc=substr($i,6)
  }
  if(clip3>CLIP){ n_clip++; next }
  if(!(cb in wl)){ n_nocb++; if(rc=="1") n_rc1_notwl++; next }
  if(rc!="1") n_wl_rc0++
  n_out++
  print chrom, term1-1, strand, cb, umi
}
END{
  printf("[termini] whitelist_cells=%d records_in=%d off_chrom=%d dropped_3p_softclip_gt%d=%d dropped_cb_not_realcell=%d kept=%d\n",
         nwl, n_in, n_offchrom+0, CLIP, n_clip+0, n_nocb+0, n_out+0) > "/dev/stderr"
  printf("[termini] rc-tag concordance: kept_with_rc_not1=%d  dropped_cb_but_rc1=%d\n",
         n_wl_rc0+0, n_rc1_notwl+0) > "/dev/stderr"
}

# Make a STARsolo BAM CellRanger-like for patched kallisto 0.46.2sq `bus --bam`:
# kallisto segfaults unless every record carries raw-tag CR/CY/UR/UY (CellRanger invariant).
# STARsolo emits only corrected CB/UB (with "-" for unassigned). We copy CB->CR, UB->UR,
# synthesize all-F quals (CY/UY), and DROP records whose CB or UB is missing/"-"
# (~2-3% of records; these lack a corrected barcode/UMI and would be discarded by
# bustools correct/count anyway).
BEGIN { FS = OFS = "\t" }
/^@/ { print; next }
{
  cb = ""; ub = ""
  for (i = 12; i <= NF; i++) {
    if ($i ~ /^CB:Z:/) cb = substr($i, 6)
    else if ($i ~ /^UB:Z:/) ub = substr($i, 6)
  }
  if (cb == "" || cb == "-" || ub == "" || ub == "-") next
  q1 = sprintf("%*s", length(cb), ""); gsub(/ /, "F", q1)
  q2 = sprintf("%*s", length(ub), ""); gsub(/ /, "F", q2)
  print $0, "CR:Z:" cb, "CY:Z:" q1, "UR:Z:" ub, "UY:Z:" q2
}

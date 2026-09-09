#!/usr/bin/env bash
# Assemble submission/1_UPLOAD from the built sources. Idempotent: safe to re-run.
#
# 1_UPLOAD is a DERIVED folder. Never edit anything inside it — edit the draft in
# manuscript/00_draft, re-run submission/4_BUILD/build.sh, then re-run this script.
set -euo pipefail
export LC_ALL=C
WD=/mnt/ssd1/Projects/PeakATail_wd
S="$WD/submission"
U="$S/1_UPLOAD"
FIG="$WD/manuscript/figures"

rm -rf "$U"
mkdir -p "$U/04_Main_figures" "$U/05_Supplementary_figures"

# --- the two manuscript documents -----------------------------------------
cp "$S/4_BUILD/output/PeakATail_manuscript_submission.docx" "$U/01_Manuscript.docx"
cp "$S/4_BUILD/output/PeakATail_supplementary.docx"         "$U/02_Additional_file_1_Supplementary.docx"

# --- cover letter, as .docx and as plain text for pasting into a portal ----
pandoc -f markdown -t docx "$S/3_SUPPORTING/cover_letter.md" -o "$U/03_Cover_letter.docx"
cp "$S/3_SUPPORTING/cover_letter.md" "$U/03_Cover_letter.md"

# --- main figures, renamed to what a journal portal expects ----------------
declare -A MAIN=(
  [1]=fig1_overview [2]=fig2_accuracy [3]=fig3_tradeoff
  [4]=fig4_calibration [5]=fig5_spermatogenesis [6]=fig6_cohort
)
for n in 1 2 3 4 5 6; do
  src="${MAIN[$n]}"
  cp "$FIG/$src.pdf" "$U/04_Main_figures/Figure_$n.pdf"
  cp "$FIG/$src.png" "$U/04_Main_figures/Figure_$n.png"
done

# --- supplementary figures -------------------------------------------------
for f in "$FIG"/figS*.pdf; do
  b=$(basename "$f" .pdf)                 # figS10_clustering
  n=$(echo "$b" | sed -E 's/^figS([0-9]+)_.*/\1/')
  cp "$f" "$U/05_Supplementary_figures/Figure_S$n.pdf"
  [ -f "$FIG/$b.png" ] && cp "$FIG/$b.png" "$U/05_Supplementary_figures/Figure_S$n.png"
done

# --- manifest, so what was shipped is recorded -----------------------------
{
  echo "# 1_UPLOAD manifest"
  echo
  echo "Assembled by submission/assemble_upload.sh. Do not edit files here."
  echo
  printf '%-46s %10s  %s\n' FILE BYTES SHA256
  find "$U" -type f ! -name MANIFEST.md | sort | while read -r f; do
    printf '%-46s %10s  %s\n' "${f#$U/}" "$(stat -c%s "$f")" "$(sha256sum "$f" | cut -c1-16)"
  done
} > "$U/MANIFEST.md"

echo "assembled $U"
find "$U" -type f | wc -l | xargs echo "  files:"
du -sh "$U" | cut -f1 | xargs echo "  size: "

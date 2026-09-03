#!/usr/bin/env bash
# Build the three PeakATail submission Word documents from the verified draft.
#
#   bash submission/build/build.sh          # run from anywhere; WD is fixed below
#
# Inputs  : manuscript/00_draft/MAIN.md, manuscript/00_draft/SUPPLEMENTARY.md,
#           manuscript/figures/fig{1..6}_*.png
# Outputs : submission/PeakATail_manuscript_submission.docx
#           submission/PeakATail_manuscript_reading_copy.docx
#           submission/PeakATail_supplementary.docx
#
# The build is idempotent and derives everything from the two Markdown files;
# re-run it after any edit to the draft.  It never writes into manuscript/.
set -euo pipefail
export LC_ALL=C
export OMP_NUM_THREADS=1

WD=/mnt/ssd1/Projects/PeakATail_wd
B="$WD/submission/build"
OUT="$WD/submission"
PANDOC=/usr/bin/pandoc

# markdown reader:
#  -implicit_figures : a lone image stays an inline image in a styled paragraph
#                      instead of becoming a pandoc figure with an empty caption.
#  -smart            : REQUIRED. With smart on, pandoc rewrites '--' as an en dash
#                      and straight quotes as curly ones; that silently corrupted
#                      the bare command-line flag '--peak-workers' in Results.
#                      Off, the verified text passes through character for character.
FMT='markdown-implicit_figures-smart'

echo "== 0. pin pandoc's own default reference.docx as the style base =="
"$PANDOC" --print-default-data-file reference.docx > "$B/reference_base.docx"

echo "== 1. build one reference.docx per output profile =="
python3 "$B/make_reference.py" submission    "$B/reference_base.docx" "$B/reference_submission.docx"
python3 "$B/make_reference.py" reading       "$B/reference_base.docx" "$B/reference_reading.docx"
python3 "$B/make_reference.py" supplementary "$B/reference_base.docx" "$B/reference_supplementary.docx"

echo "== 2. preprocess the verified Markdown into pandoc build files =="
python3 "$B/prepare_md.py" "$WD"

echo "== 3. pandoc -> docx =="
"$PANDOC" -f "$FMT" -t docx \
  --reference-doc="$B/reference_submission.docx" \
  --resource-path="$WD/manuscript/figures:$WD" \
  -o "$OUT/PeakATail_manuscript_submission.docx" \
  "$B/main_submission.md"

"$PANDOC" -f "$FMT" -t docx \
  --reference-doc="$B/reference_reading.docx" \
  --resource-path="$WD/manuscript/figures:$WD" \
  -o "$OUT/PeakATail_manuscript_reading_copy.docx" \
  "$B/main_reading.md"

"$PANDOC" -f "$FMT" -t docx \
  --reference-doc="$B/reference_supplementary.docx" \
  --resource-path="$WD/manuscript/figures:$WD" \
  -o "$OUT/PeakATail_supplementary.docx" \
  "$B/supplementary.md"

echo "== 4. postprocess: continuous line numbers + page-number footer =="
python3 "$B/postprocess_docx.py" "$OUT/PeakATail_manuscript_submission.docx" \
  --title "PeakATail: precision-first poly(A)-site calling and calibrated alternative polyadenylation analysis in single-cell RNA-seq"
python3 "$B/postprocess_docx.py" "$OUT/PeakATail_manuscript_reading_copy.docx" \
  --no-line-numbers --title "PeakATail manuscript - reading copy (figures inline)"
python3 "$B/postprocess_docx.py" "$OUT/PeakATail_supplementary.docx" \
  --no-line-numbers --title "PeakATail - Additional file 1: Supplementary figure legends"

echo "== 5. verify =="
python3 "$B/verify_docx.py" "$WD"

echo
echo "Built:"
ls -la "$OUT"/PeakATail_manuscript_submission.docx \
       "$OUT"/PeakATail_manuscript_reading_copy.docx \
       "$OUT"/PeakATail_supplementary.docx

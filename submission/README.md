# Submission package — PeakATail methods paper

**Everything you upload to the journal is in `1_UPLOAD/`. Nothing else needs to be sent.**

Target: Genome Biology, Method track — but see the open cost question in `OPEN_ITEMS.md`.

---

## The folders

| Folder | What it is | Do you touch it? |
|---|---|---|
| **`1_UPLOAD/`** | Exactly the files the journal receives, named as a portal expects | **No — it is generated.** Edit the draft, rebuild, re-assemble |
| `2_REVIEW_COPY/` | Reading copy with the six main figures inline, for co-authors | Read it; never upload it |
| `3_SUPPORTING/` | Cover-letter source, references, licences, statements, checklists | Yes — these are sources |
| `4_BUILD/` | The machinery that turns the draft into Word files, plus its `output/` | Only if the build changes |
| **`OPEN_ITEMS.md`** | **What still blocks submission, and who owns each item** | **Read this first** |

## What goes where in the portal

| Portal field | File |
|---|---|
| Manuscript | `1_UPLOAD/01_Manuscript.docx` — line-numbered, no embedded images |
| Additional file 1 | `1_UPLOAD/02_Additional_file_1_Supplementary.docx` |
| Cover letter | `1_UPLOAD/03_Cover_letter.docx`, or paste `03_Cover_letter.md` |
| Figures 1–6 | `1_UPLOAD/04_Main_figures/Figure_N.pdf` (vector; `.png` at 600 dpi alongside) |
| Supplementary figures | `1_UPLOAD/05_Supplementary_figures/Figure_SN.pdf` |

`1_UPLOAD/MANIFEST.md` records every shipped file with its size and hash.

## Regenerating

`1_UPLOAD/` is derived. Never edit inside it — the next assemble overwrites everything.

```bash
cd /mnt/ssd1/Projects/PeakATail_wd
bash submission/4_BUILD/build.sh        # draft -> Word files in 4_BUILD/output/
bash submission/assemble_upload.sh      # Word + figures -> 1_UPLOAD/
```

The draft itself lives outside this folder, in `manuscript/00_draft/MAIN.md` and
`SUPPLEMENTARY.md`. That is the single source of truth for the text.

## What is verified

- Body 11 pt, headings bold 12 pt, legends 10 pt, body justified — checked by reading the
  built files back, not by trusting the builder (`4_BUILD/check_typography.py`).
- Submission copy carries 0 embedded images; reading copy carries the 6 mains inline.
- `build.sh` runs 38 structural checks on every build.
- Abstract 342 words against Genome Biology's 350.

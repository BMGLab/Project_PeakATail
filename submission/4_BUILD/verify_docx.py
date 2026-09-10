#!/usr/bin/env python3
"""Open the built .docx files back up and check them against the source.

    python3 verify_docx.py <WD>

Checks, per document:
  * heading sequence (section order)
  * embedded image count  (submission: 0, reading: 6, supplementary: 0)
  * author-block superscript runs and the corresponding-author ORCID
  * [[CITE]] / [[PLACEHOLDER]] token counts against the source Markdown
  * continuous line numbering and the page-number footer field
  * every source paragraph present in the document text (no-loss)
Exits non-zero on any failure.
"""
import os
import re
import sys
import zipfile

import docx
from docx.oxml.ns import qn

FAIL = []
OK = []


def check(cond, msg):
    (OK if cond else FAIL).append(msg)
    print(("  PASS  " if cond else "  FAIL  ") + msg)


def doc_text(d):
    """Paragraphs AND table cells.

    Markdown table rows become docx TABLE CELLS, not paragraphs, so a
    paragraph-only reading made every table invisible to the no-loss check --
    Supplementary Table T7 could have shipped empty and this would have passed.
    """
    parts = [p.text for p in d.paragraphs]
    for tb in d.tables:
        for row in tb.rows:
            parts.append(" ".join(c.text for c in row.cells))
    return "\n".join(parts)


def image_count(path):
    z = zipfile.ZipFile(path)
    media = [n for n in z.namelist() if n.startswith("word/media/")]
    body = z.read("word/document.xml").decode("utf-8")
    blips = re.findall(r"<a:blip[^>]*r:embed", body)
    z.close()
    return len(media), len(blips)


def emu_to_in(v):
    return v / 914400.0


def strip_md(s):
    s = re.sub(r"```\{=openxml\}.*?```", "", s, flags=re.S)
    s = re.sub(r"^:::.*$", "", s, flags=re.M)
    s = s.replace("**", "").replace("`", "")
    s = re.sub(r"\^([^^\s][^^]*?)\^", r"\1", s)      # ^1^ -> 1
    s = re.sub(r"!\[\]\([^)]*\)\{[^}]*\}", "", s)     # images
    s = s.replace("\\*", "*").replace("\\_", "_")
    s = re.sub(r"^\s*>\s?", "", s, flags=re.M)        # block quotes
    s = re.sub(r"^\s*[-*]\s+", "", s, flags=re.M)     # bullets
    s = re.sub(r"^\s*\d+\.\s+", "", s, flags=re.M)    # ordered list markers
    s = re.sub(r"^#+\s*", "", s, flags=re.M)
    # a markdown table row renders as a row of cells; drop the pipes so the probe
    # matches the text docx actually holds, and drop separator rows entirely
    s = re.sub(r"^[ \t]*\|[ \t:|-]+\|[ \t]*$", "", s, flags=re.M)
    s = re.sub(r"^[ \t]*\|[ \t]*(.*?)[ \t]*\|[ \t]*$", lambda m: m.group(1).replace("|", " "), s, flags=re.M)
    s = re.sub(r"\s+", " ", s)
    return s.strip()


def paragraph_noloss(build_md, d, label):
    body = re.sub(r"\s+", " ", doc_text(d))
    missing = []
    for blk in re.split(r"\n\s*\n", open(build_md, encoding="utf-8").read()):
        t = strip_md(blk)
        if len(t) < 25:
            continue
        probe = t if len(t) <= 160 else t[:160]
        if probe not in body:
            missing.append(probe[:90])
    check(not missing,
          "%s: every source paragraph present in the document (%d checked, "
          "%d missing)" % (label, len(re.split(r'\n\s*\n', open(build_md, encoding='utf-8').read())), len(missing)))
    for m in missing[:6]:
        print("        missing: %r" % m)


def tokens(s):
    return (len(re.findall(r"\[\[CITE", s)), len(re.findall(r"\[\[PLACEHOLDER", s)))


def sectpr_flags(path):
    z = zipfile.ZipFile(path)
    body = z.read("word/document.xml").decode("utf-8")
    has_ln = "w:lnNumType" in body
    footers = [n for n in z.namelist() if re.match(r"word/footer\d*\.xml", n)]
    foot = "".join(z.read(n).decode("utf-8") for n in footers)
    z.close()
    return has_ln, ("PAGE" in foot), ("NUMPAGES" in foot), len(footers)


def main(wd):
    B = os.path.join(wd, "submission", "4_BUILD")
    OUTD = os.path.join(wd, "submission", "4_BUILD", "output")

    jobs = [
        ("SUBMISSION", os.path.join(OUTD, "PeakATail_manuscript_submission.docx"),
         os.path.join(B, "main_submission.md"), 0, True),
        ("READING", os.path.join(OUTD, "PeakATail_manuscript_reading_copy.docx"),
         os.path.join(B, "main_reading.md"), 6, False),
        ("SUPPLEMENTARY", os.path.join(OUTD, "PeakATail_supplementary.docx"),
         os.path.join(B, "supplementary.md"), 0, False),
    ]

    for label, path, build_md, want_imgs, want_lines in jobs:
        print("\n=== %s : %s ===" % (label, os.path.basename(path)))
        d = docx.Document(path)
        heads = [(p.style.name, p.text) for p in d.paragraphs
                 if p.style.name.startswith("Heading") or p.style.name == "Title"]
        print("  section order (%d headings):" % len(heads))
        show = ("Title", "Heading 1", "Heading 2") if label == "SUPPLEMENTARY" \
            else ("Title", "Heading 1")
        for st, tx in heads:
            if st in show:
                print("      %-9s %s" % (st, tx))

        media, blips = image_count(path)
        check(media == want_imgs and blips == want_imgs,
              "%s: %d embedded image part(s), %d drawing reference(s) "
              "(expected %d)" % (label, media, blips, want_imgs))

        if want_imgs:
            widths = []
            for sh in d.inline_shapes:
                widths.append((emu_to_in(sh.width), emu_to_in(sh.height)))
            print("      image sizes (in): %s"
                  % ", ".join("%.2f x %.2f" % w for w in widths))
            sec = d.sections[0]
            col = emu_to_in(sec.page_width - sec.left_margin - sec.right_margin)
            usable_h = emu_to_in(sec.page_height - sec.top_margin - sec.bottom_margin)
            check(all(w <= col + 1e-6 for w, _ in widths),
                  "%s: every image width <= %.2f in text column" % (label, col))
            check(all(h <= usable_h + 1e-6 for _, h in widths),
                  "%s: every image height <= %.2f in usable page height"
                  % (label, usable_h))

        txt = doc_text(d)
        src = open(build_md, encoding="utf-8").read()
        c_src, p_src = tokens(src)
        c_doc, p_doc = tokens(txt)
        check(c_src == c_doc, "%s: [[CITE]] tokens %d in build md == %d in docx"
              % (label, c_src, c_doc))
        check(p_src == p_doc,
              "%s: [[PLACEHOLDER]] tokens %d in build md == %d in docx"
              % (label, p_src, p_doc))

        has_ln, has_page, has_np, nfoot = sectpr_flags(path)
        check(has_ln == want_lines,
              "%s: continuous line numbering %s"
              % (label, "present" if has_ln else "absent"))
        check(has_page and has_np and nfoot >= 1,
              "%s: footer carries PAGE and NUMPAGES fields (%d footer part(s))"
              % (label, nfoot))

        paragraph_noloss(build_md, d, label)

        if label != "SUPPLEMENTARY":
            # author block
            auth = [p for p in d.paragraphs if p.text.startswith("Amir Amiri Tabat")]
            check(len(auth) == 1, "%s: author line present exactly once" % label)
            if auth:
                sup = [r.text for r in auth[0].runs if r.font.superscript]
                check(sup == ["1", "2,*"],
                      "%s: author superscripts %r (expected ['1', '2,*'])"
                      % (label, sup))
            affs = [p for p in d.paragraphs
                    if re.match(r"^[12] Department of ", p.text)]
            check(len(affs) == 2, "%s: two affiliation lines, each its own "
                  "paragraph" % label)
            for a in affs:
                sup = [r.text for r in a.runs if r.font.superscript]
                check(len(sup) == 1 and sup[0] in ("1", "2"),
                      "%s: affiliation marker %r is superscript" % (label, sup))
            check("ORCID: 0000-0002-9725-7536" in txt,
                  "%s: corresponding-author ORCID present" % label)
            check("yasin.kaymaz@ege.edu.tr" in txt,
                  "%s: corresponding-author e-mail present" % label)

        # section order assertions
        h1 = [t for s, t in heads if s == "Heading 1"]
        if label == "SUBMISSION":
            want = ["Abstract", "Key points", "1 Background", "2 Results",
                    "3 Discussion", "4 Conclusions", "5 Methods", "Declarations",
                    "References", "Figure legends", "Supplementary information"]
            check(h1 == want, "%s: H1 order == %s" % (label, want))
            body = "\n".join(x.text for x in d.paragraphs)
            n_doi = body.count("doi:10.")
            check(n_doi >= 20,
                  "%s: reference list is populated (%d doi entries)" % (label, n_doi))
        elif label == "READING":
            want = ["Abstract", "Key points", "1 Background", "2 Results",
                    "3 Discussion", "4 Conclusions", "5 Methods", "Declarations",
                    "References", "Supplementary information"]
            check(h1 == want, "%s: H1 order == %s "
                  "(legends are inline, so no collected section)" % (label, want))
            legend_starts = [p.text[:12] for p in d.paragraphs
                             if p.style.name == "FigureLegend"]
            check(len(legend_starts) == 6,
                  "%s: 6 FigureLegend paragraphs, one per inline figure (%s)"
                  % (label, legend_starts))
        else:
            check(len(h1) >= 1, "%s: %d H1 section(s)" % (label, len(h1)))

    print("\n%d checks passed, %d failed." % (len(OK), len(FAIL)))
    if FAIL:
        for f in FAIL:
            print("  FAILED: " + f)
        sys.exit(1)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "/mnt/ssd1/Projects/PeakATail_wd")

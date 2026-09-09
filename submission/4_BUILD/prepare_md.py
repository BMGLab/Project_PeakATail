#!/usr/bin/env python3
"""Turn the verified manuscript Markdown into three pandoc-ready build files.

    python3 prepare_md.py <WD>

Reads   WD/manuscript/00_draft/MAIN.md
        WD/manuscript/00_draft/SUPPLEMENTARY.md
Writes  WD/submission/4_BUILD/main_submission.md
        WD/submission/4_BUILD/main_reading.md
        WD/submission/4_BUILD/supplementary.md

Rules (assembly only — no claim, number, citation or author is altered):
  * HTML comments are dropped (a submission file carries no working notes).
  * The H1 title becomes a Title-styled block; the author / affiliation /
    corresponding lines become TitlePageLine-styled blocks so their pandoc
    superscripts (^1^, ^2,\\*^) survive as real superscript runs.
  * '## Abstract' and '## Key points' are promoted to H1 so that H1 == major
    section and H2 == subsection throughout.
  * Plain section numbers are written into the heading text (pandoc 2.9's docx
    writer ignores --number-sections). Numbered: Background, Results,
    Discussion, Conclusions, Methods.  Unnumbered: Abstract, Key points,
    Declarations, Figure legends, Supplementary information.
  * Reading copy only: each of Fig 1-6 is inserted after the body paragraph that
    first cites it, followed by that figure's legend taken verbatim from the
    manuscript's own '# Figure legends' section; the collected legend section is
    then dropped from the reading copy because every legend now sits inline.
"""
import os
import re
import sys

NUMBERED = ["Background", "Results", "Discussion", "Conclusions", "Methods"]

FIGS = [
    (1, "fig1_overview.png"),
    (2, "fig2_accuracy.png"),
    (3, "fig3_tradeoff.png"),
    (4, "fig4_calibration.png"),
    (5, "fig5_spermatogenesis.png"),
    (6, "fig6_cohort.png"),
]

READING_IMG_WIDTH_IN = 6.4   # A4 with 2 cm margins gives a 6.69 in text column

PAGEBREAK = ("```{=openxml}\n"
             '<w:p><w:r><w:br w:type="page"/></w:r></w:p>\n'
             "```")


def strip_html_comments(text):
    found = re.findall(r"<!--.*?-->", text, flags=re.S)
    text = re.sub(r"<!--.*?-->\n?", "", text, flags=re.S)
    text = re.sub(r"\n{3,}", "\n\n", text)
    return text, found


def blocks(text):
    """Split into blank-line-separated blocks, preserving order."""
    return [b for b in re.split(r"\n\s*\n", text.strip("\n"))]


def number_headings(blks):
    """Write plain section numbers into heading text. Returns new block list."""
    out = []
    top = 0
    sub = 0
    in_numbered = False
    for b in blks:
        m = re.match(r"^(#{1,6})\s+(.*)$", b)
        if not m:
            out.append(b)
            continue
        hashes, title = m.group(1), m.group(2).rstrip()
        level = len(hashes)
        if level == 1:
            if title in NUMBERED:
                top += 1
                sub = 0
                in_numbered = True
                out.append("%s %d %s" % (hashes, top, title))
            else:
                in_numbered = False
                out.append(b)
        elif level == 2 and in_numbered:
            sub += 1
            out.append("%s %d.%d %s" % (hashes, top, sub, title))
        else:
            out.append(b)
    return out


def split_sections(blks):
    """-> list of (h1_title_or_None, [blocks]) in document order."""
    secs = []
    cur_title, cur = None, []
    for b in blks:
        m = re.match(r"^#\s+(.*)$", b)
        if m:
            secs.append((cur_title, cur))
            cur_title, cur = m.group(1).strip(), [b]
        else:
            cur.append(b)
    secs.append((cur_title, cur))
    return [s for s in secs if s[1]]


def main(wd):
    src = os.path.join(wd, "manuscript", "00_draft")
    build = os.path.join(wd, "submission", "4_BUILD")
    figdir = os.path.join(wd, "manuscript", "figures")

    raw = open(os.path.join(src, "MAIN.md"), encoding="utf-8").read()
    text, comments = strip_html_comments(raw)
    print("MAIN.md: stripped %d HTML comment block(s)" % len(comments))
    for c in comments:
        print("   [stripped] %s..." % c[:70].replace("\n", " "))

    blks = blocks(text)

    # --- title page ---
    assert blks[0].startswith("# "), blks[0][:60]
    title = blks[0][2:].strip()
    front = []
    i = 1
    while i < len(blks) and not blks[i].startswith("#"):
        front.append(blks[i])
        i += 1
    rest = blks[i:]

    # --- promote Abstract / Key points to H1 ---
    rest = [re.sub(r"^## (Abstract|Key points)\s*$", r"# \1", b) for b in rest]

    rest = number_headings(rest)

    # --- pull the legends out of the collected '# Figure legends' section ---
    secs = split_sections(rest)
    legends = {}
    for stitle, sblocks in secs:
        if stitle == "Figure legends":
            for b in sblocks:
                m = re.match(r"^\*\*Fig (\d+) \|", b)
                if m:
                    legends[int(m.group(1))] = b
    print("legends found in '# Figure legends': %s"
          % sorted(legends))

    # Each physical line of the front matter is its own paragraph: the two
    # affiliation lines must not be folded into one by a markdown soft break.
    title_page = ['::: {custom-style="Title"}', title, ":::", ""]
    n_front = 0
    for f in front:
        for line in f.split("\n"):
            if line.strip():
                title_page += ['::: {custom-style="TitlePageLine"}',
                               line.strip(), ":::", ""]
                n_front += 1
    title_page += [PAGEBREAK, ""]
    print("title page: 1 title + %d front-matter lines" % n_front)

    # ---------------- submission copy ----------------
    sub_rest = []
    for b in rest:
        if b.startswith("# Figure legends"):
            sub_rest.append(PAGEBREAK)   # legends start on their own page
        sub_rest.append(b)
    sub_md = "\n".join(title_page) + "\n" + "\n\n".join(sub_rest) + "\n"
    open(os.path.join(build, "main_submission.md"), "w", encoding="utf-8").write(sub_md)

    # ---------------- reading copy ----------------
    read_blocks = []
    placed = set()
    for stitle, sblocks in secs:
        if stitle == "Figure legends":
            continue            # every legend is inline in this copy
        for b in sblocks:
            read_blocks.append(b)
            if b.startswith("#"):
                continue
            for n, fname in FIGS:
                if n in placed or n not in legends:
                    continue
                if re.search(r"\bFig %d(?![0-9])" % n, b):
                    path = os.path.join(figdir, fname).replace("\\", "/")
                    read_blocks.append(
                        '::: {custom-style="FigureImage"}\n'
                        '![](%s){width=%.2fin}\n:::' % (path, READING_IMG_WIDTH_IN))
                    read_blocks.append(
                        '::: {custom-style="FigureLegend"}\n%s\n:::' % legends[n])
                    placed.add(n)
                    print("reading copy: Fig %d placed after a paragraph in "
                          "section %r" % (n, stitle))
                    break
    missing = [n for n, _ in FIGS if n not in placed]
    assert not missing, "no citation anchor found for Fig %s" % missing
    read_md = "\n".join(title_page) + "\n" + "\n\n".join(read_blocks) + "\n"
    open(os.path.join(build, "main_reading.md"), "w", encoding="utf-8").write(read_md)

    # ---------------- supplementary ----------------
    sraw = open(os.path.join(src, "SUPPLEMENTARY.md"), encoding="utf-8").read()
    stext, scomments = strip_html_comments(sraw)
    print("SUPPLEMENTARY.md: stripped %d HTML comment block(s)" % len(scomments))
    sblks = blocks(stext)
    assert sblks[0].startswith("# "), sblks[0][:60]
    stitle = sblks[0][2:].strip()
    supp_md = ('::: {custom-style="Title"}\n%s\n:::\n\n' % stitle
               + "\n\n".join(sblks[1:]) + "\n")
    open(os.path.join(build, "supplementary.md"), "w", encoding="utf-8").write(supp_md)

    print("\nwrote main_submission.md, main_reading.md, supplementary.md in %s" % build)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "/mnt/ssd1/Projects/PeakATail_wd")

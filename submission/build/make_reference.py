#!/usr/bin/env python3
"""Build a pandoc reference.docx (styles + page setup) for one output profile.

Usage:  python3 make_reference.py <profile> <base.docx> <out.docx>
        profile in {submission, reading, supplementary}

The base is pandoc's own default reference.docx (`pandoc --print-default-data-file
reference.docx`), so every style pandoc's docx writer expects survives; this script
only overrides the ones we care about and appends a few custom styles.

Nothing here touches manuscript content.
"""
import sys
import shutil
import zipfile
from lxml import etree

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
NS = {"w": W}


def q(tag):
    return "{%s}%s" % (W, tag)


PROFILES = {
    # body_pt, line (240ths of a line; 240=single, 276=1.15, 360=1.5, 480=double),
    # space_after (twips), h1_pt, h2_pt, h3_pt, title_pt, margin (twips)
    "submission": dict(body_pt=12, line=360, after=120, h1_pt=14, h2_pt=12, h3_pt=12,
                       title_pt=16, margin=1440, legend_pt=12),
    "reading": dict(body_pt=11, line=276, after=160, h1_pt=14, h2_pt=12, h3_pt=11,
                    title_pt=16, margin=1134, legend_pt=10),
    "supplementary": dict(body_pt=11, line=276, after=160, h1_pt=14, h2_pt=12, h3_pt=11,
                          title_pt=15, margin=1247, legend_pt=11),
}

SERIF = "Times New Roman"
MONO = "Courier New"


def frag(xml):
    return etree.fromstring(
        '<w:root xmlns:w="%s">%s</w:root>' % (W, xml)
    )


def half(pt):
    return str(int(round(pt * 2)))


def build(profile, base_path, out_path):
    cfg = PROFILES[profile]
    shutil.copyfile(base_path, out_path + ".tmp")

    zin = zipfile.ZipFile(base_path)
    parts = {n: zin.read(n) for n in zin.namelist()}
    zin.close()

    # ---------------- styles.xml ----------------
    st = etree.fromstring(parts["word/styles.xml"])

    # docDefaults: serif body font, size, line spacing
    dd = st.find(q("docDefaults"))
    rpr = dd.find(q("rPrDefault")).find(q("rPr"))
    rf = rpr.find(q("rFonts"))
    for a in list(rf.attrib):
        del rf.attrib[a]
    rf.set(q("ascii"), SERIF)
    rf.set(q("hAnsi"), SERIF)
    rf.set(q("eastAsia"), SERIF)
    rf.set(q("cs"), SERIF)
    rpr.find(q("sz")).set(q("val"), half(cfg["body_pt"]))
    rpr.find(q("szCs")).set(q("val"), half(cfg["body_pt"]))

    ppr = dd.find(q("pPrDefault")).find(q("pPr"))
    for ch in list(ppr):
        ppr.remove(ch)
    sp = etree.SubElement(ppr, q("spacing"))
    sp.set(q("before"), "0")
    sp.set(q("after"), str(cfg["after"]))
    sp.set(q("line"), str(cfg["line"]))
    sp.set(q("lineRule"), "auto")

    def style(sid):
        for s in st.findall(q("style")):
            if s.get(q("styleId")) == sid:
                return s
        return None

    def set_ppr(s, before=None, after=None, line=None, jc=None, keepnext=None,
                indent_left=None):
        p = s.find(q("pPr"))
        if p is None:
            p = etree.SubElement(s, q("pPr"))
            s.remove(p)
            # pPr must precede rPr in a w:style
            idx = 0
            for i, ch in enumerate(s):
                if ch.tag == q("rPr"):
                    idx = i
                    break
                idx = i + 1
            s.insert(idx, p)
        old = p.find(q("spacing"))
        if old is not None:
            p.remove(old)
        sp = etree.Element(q("spacing"))
        sp.set(q("before"), str(before if before is not None else 0))
        sp.set(q("after"), str(after if after is not None else cfg["after"]))
        sp.set(q("line"), str(line if line is not None else cfg["line"]))
        sp.set(q("lineRule"), "auto")
        p.insert(0, sp)
        if jc is not None:
            oldj = p.find(q("jc"))
            if oldj is not None:
                p.remove(oldj)
            j = etree.SubElement(p, q("jc"))
            j.set(q("val"), jc)
        if indent_left is not None:
            ind = etree.SubElement(p, q("ind"))
            ind.set(q("left"), str(indent_left))
        if keepnext:
            if p.find(q("keepNext")) is None:
                p.insert(0, etree.Element(q("keepNext")))
        return p

    def set_rpr(s, pt=None, bold=None, italic=None, font=SERIF, drop_color=True):
        r = s.find(q("rPr"))
        if r is None:
            r = etree.SubElement(s, q("rPr"))
        if drop_color:
            for c in r.findall(q("color")):
                r.remove(c)
        rfo = r.find(q("rFonts"))
        if rfo is None:
            rfo = etree.SubElement(r, q("rFonts"))
        for a in list(rfo.attrib):
            del rfo.attrib[a]
        rfo.set(q("ascii"), font)
        rfo.set(q("hAnsi"), font)
        rfo.set(q("eastAsia"), font)
        rfo.set(q("cs"), font)
        for tag, val in (("b", bold), ("i", italic)):
            for e in r.findall(q(tag)):
                r.remove(e)
            for e in r.findall(q(tag + "Cs")):
                r.remove(e)
            if val:
                etree.SubElement(r, q(tag))
                etree.SubElement(r, q(tag + "Cs"))
        if pt is not None:
            for tag in ("sz", "szCs"):
                for e in r.findall(q(tag)):
                    r.remove(e)
                e = etree.SubElement(r, q(tag))
                e.set(q("val"), half(pt))
        return r

    # Body styles
    for sid in ("BodyText", "FirstParagraph"):
        s = style(sid)
        if s is not None:
            set_ppr(s)
    c = style("Compact")
    if c is not None:
        set_ppr(c, after=int(cfg["after"] * 0.4))

    # Headings: plain black serif, no theme colour
    for sid, pt, bold, italic in (
        ("Heading1", cfg["h1_pt"], True, False),
        ("Heading2", cfg["h2_pt"], True, False),
        ("Heading3", cfg["h3_pt"], True, True),
        ("Heading4", cfg["h3_pt"], False, True),
        ("Heading5", cfg["h3_pt"], False, True),
        ("Heading6", cfg["h3_pt"], False, True),
    ):
        s = style(sid)
        if s is None:
            continue
        set_rpr(s, pt=pt, bold=bold, italic=italic)
        before = 320 if sid == "Heading1" else 240
        set_ppr(s, before=before, after=80, line=240, keepnext=True)

    t = style("Title")
    if t is not None:
        set_rpr(t, pt=cfg["title_pt"], bold=True)
        set_ppr(t, before=0, after=240, line=240, jc="left", keepnext=True)

    a = style("Author")
    if a is not None:
        set_rpr(a, pt=cfg["body_pt"])
        set_ppr(a, before=0, after=120, line=240, jc="left")

    for sid in ("VerbatimChar", "SourceCode"):
        s = style(sid)
        if s is not None:
            set_rpr(s, pt=cfg["body_pt"] - 1, font=MONO, drop_color=True)

    # Custom styles used by custom-style divs in the build markdown
    custom = """
    <w:style w:type="paragraph" w:customStyle="1" w:styleId="TitlePageLine">
      <w:name w:val="TitlePageLine"/><w:basedOn w:val="Normal"/><w:qFormat/>
      <w:pPr><w:spacing w:before="0" w:after="60" w:line="240" w:lineRule="auto"/></w:pPr>
      <w:rPr><w:rFonts w:ascii="{serif}" w:hAnsi="{serif}" w:eastAsia="{serif}" w:cs="{serif}"/>
        <w:sz w:val="{bodyhalf}"/><w:szCs w:val="{bodyhalf}"/></w:rPr>
    </w:style>
    <w:style w:type="paragraph" w:customStyle="1" w:styleId="FigureImage">
      <w:name w:val="FigureImage"/><w:basedOn w:val="Normal"/><w:qFormat/>
      <w:pPr><w:keepNext/><w:spacing w:before="240" w:after="80" w:line="240" w:lineRule="auto"/>
        <w:jc w:val="center"/></w:pPr>
    </w:style>
    <w:style w:type="paragraph" w:customStyle="1" w:styleId="FigureLegend">
      <w:name w:val="FigureLegend"/><w:basedOn w:val="Normal"/><w:qFormat/>
      <w:pPr><w:spacing w:before="0" w:after="280" w:line="240" w:lineRule="auto"/>
        <w:ind w:left="284" w:right="284"/></w:pPr>
      <w:rPr><w:rFonts w:ascii="{serif}" w:hAnsi="{serif}" w:eastAsia="{serif}" w:cs="{serif}"/>
        <w:sz w:val="{leghalf}"/><w:szCs w:val="{leghalf}"/></w:rPr>
    </w:style>
    """.format(serif=SERIF, bodyhalf=half(cfg["body_pt"]), leghalf=half(cfg["legend_pt"]))
    for node in frag(custom):
        st.append(node)

    parts["word/styles.xml"] = etree.tostring(
        st, xml_declaration=True, encoding="UTF-8", standalone=True)

    # ---------------- document.xml: page size / margins ----------------
    doc = etree.fromstring(parts["word/document.xml"])
    body = doc.find(q("body"))
    sect = body.find(q("sectPr"))
    if sect is None:
        sect = etree.SubElement(body, q("sectPr"))
    for ch in list(sect):
        sect.remove(ch)
    pgsz = etree.SubElement(sect, q("pgSz"))          # A4 portrait
    pgsz.set(q("w"), "11906")
    pgsz.set(q("h"), "16838")
    m = cfg["margin"]
    pgm = etree.SubElement(sect, q("pgMar"))
    for k, v in (("top", m), ("right", m), ("bottom", m), ("left", m),
                 ("header", 708), ("footer", 708), ("gutter", 0)):
        pgm.set(q(k), str(v))
    parts["word/document.xml"] = etree.tostring(
        doc, xml_declaration=True, encoding="UTF-8", standalone=True)

    zout = zipfile.ZipFile(out_path, "w", zipfile.ZIP_DEFLATED)
    for n, dat in parts.items():
        zout.writestr(n, dat)
    zout.close()
    import os
    os.remove(out_path + ".tmp")
    print("wrote %s (profile=%s, %d pt serif, line=%d/240, A4, margin=%d twips)"
          % (out_path, profile, cfg["body_pt"], cfg["line"], cfg["margin"]))


if __name__ == "__main__":
    build(sys.argv[1], sys.argv[2], sys.argv[3])

#!/usr/bin/env python3
"""Replace abbreviated author lists with the publisher-deposited full names
from Crossref. This is NOT expansion-from-memory: every given name below is
read out of the DOI's own Crossref deposit, fetched in crossref_records.json.

Safety rules:
  * only rewrite when the surname sequence in the .bib matches Crossref's
    surname sequence exactly (same order, same count);
  * only rewrite when Crossref strictly ADDS information (at least one given
    name longer than the .bib's initials);
  * never rewrite an entry whose .bib names are already full.
"""
import json
import re
import unicodedata

BIB = "/mnt/ssd1/Projects/PeakATail_wd/submission/references.bib"
REC = "/mnt/ssd1/Projects/PeakATail_wd/submission/.verify/crossref_records.json"

ACCENT = {"é": "{\\'e}", "í": "{\\'i}", "’": "'", " ": " "}


def to_bibtex(s):
    for k, v in ACCENT.items():
        s = s.replace(k, v)
    s = re.sub(r"\s+", " ", s).strip()
    assert all(ord(c) < 128 for c in s), "non-ascii survived: %r" % s
    return s


def surname(a):
    a = a.split(",")[0]
    for k, v in {"é": "e", "í": "i", "’": "'"}.items():
        a = a.replace(k, v)
    a = a.replace("{\\'e}", "e").replace("{\\'i}", "i")
    a = a.replace("{", "").replace("}", "").replace("\\", "")
    return re.sub(r"\s+", " ", a).strip().lower()


txt = open(BIB, encoding="utf-8").read()
recs = {r["key"]: r for r in json.load(open(REC))}
changed = []

for key, r in recs.items():
    cr = r.get("cr_authors")
    if not cr:
        continue
    bib_list = [a.strip() for a in r["bib_author"].split(" and ") if a.strip()]
    if len(bib_list) != len(cr):
        continue
    if [surname(a) for a in bib_list] != [surname(a) for a in cr]:
        continue
    gains = 0
    for b, c in zip(bib_list, cr):
        bg = b.split(",", 1)[1].strip() if "," in b else ""
        cg = c.split(",", 1)[1].strip() if "," in c else ""
        if len(cg.replace(".", "").replace(" ", "")) > len(bg.replace(".", "").replace(" ", "")):
            gains += 1
    if gains == 0:
        continue
    new = " and ".join(to_bibtex(a) for a in cr)
    # locate this entry's author field and replace it
    m = re.search(r"@\w+\{" + re.escape(key) + r",", txt)
    assert m, key
    seg_start = m.end()
    d = 1
    p = txt.index("{", m.start())
    d = 0
    q = p
    while q < len(txt):
        if txt[q] == "{":
            d += 1
        elif txt[q] == "}":
            d -= 1
            if d == 0:
                break
        q += 1
    body = txt[seg_start:q]
    am = re.search(r"(author\s*=\s*)\{", body)
    assert am, key
    k = am.end() - 1
    dd = 0
    pp = k
    while pp < len(body):
        if body[pp] == "{":
            dd += 1
        elif body[pp] == "}":
            dd -= 1
            if dd == 0:
                break
        pp += 1
    old = body[k + 1:pp]
    newbody = body[:k + 1] + new + body[pp:]
    txt = txt[:seg_start] + newbody + txt[q:]
    changed.append((key, len(cr), old[:60], new[:60]))

open(BIB, "w", encoding="utf-8").write(txt)
print("rewrote %d author lists from Crossref deposits:" % len(changed))
for k, n, o, nw in changed:
    print("  %-28s n=%-3d %s...  ->  %s..." % (k, n, o, nw))

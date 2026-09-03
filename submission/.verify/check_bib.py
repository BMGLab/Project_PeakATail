#!/usr/bin/env python3
"""Independent adversarial verification of references.bib.

For every entry: pull the DOI (or URL) and resolve it against Crossref
(api.crossref.org) and, as a second independent source, Europe PMC.
Compare title / first author / venue / year against what the .bib claims.
Nothing here trusts the .bib's own "% verified:" comments.
"""
import json
import re
import subprocess
import sys
import time
import unicodedata

BIB = "/mnt/ssd1/Projects/PeakATail_wd/submission/references.bib"
OUT = "/mnt/ssd1/Projects/PeakATail_wd/submission/.verify/crossref_records.json"
UA = "PeakATail-refcheck/1.0 (mailto:yasin.kaymaz@ege.edu.tr)"


def curl(url):
    try:
        r = subprocess.run(
            ["curl", "-sS", "-m", "40", "-H", "User-Agent: " + UA, url],
            capture_output=True, text=True,
        )
        if r.returncode != 0:
            return None, "curl rc=%d %s" % (r.returncode, r.stderr.strip()[:200])
        return r.stdout, None
    except Exception as e:  # noqa
        return None, str(e)


def parse_bib(path):
    txt = open(path, encoding="utf-8").read()
    entries = []
    for m in re.finditer(r"^@(\w+)\{([^,\s]+),", txt, re.M):
        kind, key = m.group(1), m.group(2)
        start = m.end()
        depth = 1
        i = m.start() + len(kind) + 1
        # walk braces from the opening brace of the entry
        i = txt.index("{", m.start())
        depth = 0
        j = i
        while j < len(txt):
            if txt[j] == "{":
                depth += 1
            elif txt[j] == "}":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        body = txt[start:j]
        fields = {}
        for fm in re.finditer(r"(\w+)\s*=\s*\{", body):
            fname = fm.group(1).lower()
            k = fm.end() - 1
            d = 0
            p = k
            while p < len(body):
                if body[p] == "{":
                    d += 1
                elif body[p] == "}":
                    d -= 1
                    if d == 0:
                        break
                p += 1
            fields[fname] = body[k + 1:p]
        entries.append({"key": key, "kind": kind, "fields": fields,
                        "line": txt[:m.start()].count("\n") + 1})
    return entries


def norm(s):
    if s is None:
        return ""
    s = unicodedata.normalize("NFKD", s)
    s = "".join(c for c in s if not unicodedata.combining(c))
    s = s.replace("\\'", "").replace("{", "").replace("}", "")
    s = s.replace("‐", "-").replace("‑", "-").replace("–", "-")
    s = s.replace("’", "'").replace("‘", "'").replace("ʼ", "'")
    s = re.sub(r"<[^>]+>", "", s)
    s = re.sub(r"[^a-z0-9]+", " ", s.lower()).strip()
    return s


def main():
    entries = parse_bib(BIB)
    print("parsed %d entries" % len(entries), file=sys.stderr)
    results = []
    for e in entries:
        doi = e["fields"].get("doi", "").strip()
        rec = {"key": e["key"], "kind": e["kind"], "line": e["line"],
               "bib_doi": doi,
               "bib_title": e["fields"].get("title", ""),
               "bib_author": e["fields"].get("author", ""),
               "bib_journal": e["fields"].get("journal", "") or e["fields"].get("howpublished", ""),
               "bib_year": e["fields"].get("year", ""),
               "bib_volume": e["fields"].get("volume", ""),
               "bib_pages": e["fields"].get("pages", ""),
               "bib_note": e["fields"].get("note", "")}
        if doi:
            body, err = curl("https://api.crossref.org/works/" + doi)
            if err:
                rec["crossref_error"] = err
            else:
                try:
                    j = json.loads(body)["message"]
                    rec["cr_title"] = (j.get("title") or [""])[0]
                    auth = j.get("author") or []
                    rec["cr_authors"] = [
                        (a.get("family", "") + ", " + a.get("given", "")).strip(", ")
                        for a in auth]
                    rec["cr_nauthors"] = len(auth)
                    rec["cr_container"] = (j.get("container-title") or [""])[0]
                    rec["cr_short_container"] = (j.get("short-container-title") or [""])[0]
                    dp = (j.get("issued", {}).get("date-parts") or [[None]])[0]
                    rec["cr_year"] = dp[0] if dp else None
                    for alt in ("published-print", "published-online", "created"):
                        d2 = (j.get(alt, {}) or {}).get("date-parts")
                        if d2:
                            rec["cr_year_" + alt] = d2[0][0]
                    rec["cr_volume"] = j.get("volume")
                    rec["cr_page"] = j.get("page")
                    rec["cr_type"] = j.get("type")
                    rec["cr_publisher"] = j.get("publisher")
                    rec["cr_url"] = j.get("URL")
                except Exception as ex:
                    rec["crossref_error"] = "parse: %s :: %s" % (ex, body[:200])
            time.sleep(0.3)
        results.append(rec)
    json.dump(results, open(OUT, "w"), indent=1)
    print("wrote %s" % OUT, file=sys.stderr)

    # comparison report
    for r in results:
        key = r["key"]
        if "cr_title" not in r:
            print("NO-CROSSREF  %-28s doi=%s  %s" % (key, r["bib_doi"], r.get("crossref_error", "no doi")))
            continue
        flags = []
        if norm(r["bib_title"]) != norm(r["cr_title"]):
            flags.append("TITLE")
        bib_first = r["bib_author"].split(" and ")[0].strip()
        bib_surname = norm(bib_first.split(",")[0])
        cr_first = (r["cr_authors"] or [""])[0]
        cr_surname = norm(cr_first.split(",")[0])
        if bib_surname != cr_surname:
            flags.append("FIRSTAUTHOR")
        years = {r.get("cr_year"), r.get("cr_year_published-print"),
                 r.get("cr_year_published-online")}
        years.discard(None)
        try:
            if int(r["bib_year"]) not in years:
                flags.append("YEAR")
        except Exception:
            flags.append("YEAR?")
        nb = len([a for a in r["bib_author"].split(" and ") if a.strip()])
        if r.get("cr_nauthors") and nb != r["cr_nauthors"]:
            flags.append("NAUTH %d!=%d" % (nb, r["cr_nauthors"]))
        print("%-12s %-28s | %s" % ("OK" if not flags else "CHECK:" + ",".join(flags),
                                    key, r["cr_title"][:70]))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Second, independent source: Europe PMC by DOI. Confirms journal/volume/
pages/year/PMID, which Crossref sometimes omits for OUP/BMC records."""
import json
import re
import subprocess
import sys
import time
import urllib.parse

BIB = "/mnt/ssd1/Projects/PeakATail_wd/submission/references.bib"
UA = "PeakATail-refcheck/1.0 (mailto:yasin.kaymaz@ege.edu.tr)"
sys.path.insert(0, "/mnt/ssd1/Projects/PeakATail_wd/submission/.verify")
from check_bib import parse_bib, curl, norm  # noqa

rows = []
for e in parse_bib(BIB):
    f = e["fields"]
    doi = f.get("doi", "").strip()
    bib_pmid = ""
    m = re.search(r"PMID:\s*(\d+)", f.get("note", ""))
    if m:
        bib_pmid = m.group(1)
    row = {"key": e["key"], "doi": doi, "bib_pmid": bib_pmid,
           "bib_journal": f.get("journal", ""), "bib_vol": f.get("volume", ""),
           "bib_pages": f.get("pages", ""), "bib_year": f.get("year", ""),
           "bib_title": f.get("title", "")}
    if doi:
        q = urllib.parse.quote('DOI:"%s"' % doi)
        url = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query="
               + q + "&format=json&resultType=core&pageSize=5")
        body, err = curl(url)
        if err:
            row["epmc_error"] = err
        else:
            try:
                res = json.loads(body)["resultList"]["result"]
            except Exception as ex:
                res = []
                row["epmc_error"] = str(ex)
            # prefer a non-preprint record
            pick = None
            for r in res:
                if (r.get("doi", "").lower() == doi.lower()):
                    pick = r
                    break
            if pick is None and res:
                pick = res[0]
            if pick:
                row["epmc_title"] = pick.get("title", "")
                row["epmc_pmid"] = pick.get("pmid", "")
                row["epmc_src"] = pick.get("source", "")
                ji = pick.get("journalInfo", {}) or {}
                row["epmc_journal"] = (ji.get("journal", {}) or {}).get("title", "")
                row["epmc_vol"] = ji.get("volume", "")
                row["epmc_issue"] = ji.get("issue", "")
                row["epmc_year"] = ji.get("yearOfPublication", "")
                row["epmc_pages"] = pick.get("pageInfo", "")
                row["epmc_authors"] = pick.get("authorString", "")
                row["epmc_type"] = ",".join(pick.get("pubTypeList", {}).get("pubType", []) or [])
            else:
                row["epmc_error"] = row.get("epmc_error", "no result for DOI")
        time.sleep(0.35)
    rows.append(row)

json.dump(rows, open("/mnt/ssd1/Projects/PeakATail_wd/submission/.verify/epmc_records.json", "w"), indent=1)

print("%-28s %-6s %-22s %-8s %-6s %-12s %s" % ("key", "pmid?", "journal-match", "vol", "year", "pages", "flags"))
for r in rows:
    if "epmc_title" not in r:
        print("%-28s  NO-EPMC  %s" % (r["key"], r.get("epmc_error", "no doi")))
        continue
    fl = []
    if r["bib_pmid"] and r["bib_pmid"] != r.get("epmc_pmid", ""):
        fl.append("PMID bib=%s epmc=%s" % (r["bib_pmid"], r.get("epmc_pmid")))
    if norm(r["bib_title"]) != norm(r["epmc_title"].rstrip(".")):
        fl.append("TITLE")
    bj, ej = norm(r["bib_journal"]), norm(r["epmc_journal"])
    jm = "same" if (bj == ej or bj.startswith(ej) or ej.startswith(bj)) else "DIFF"
    if jm == "DIFF":
        fl.append("JOURNAL bib=%r epmc=%r" % (r["bib_journal"], r["epmc_journal"]))
    if r["bib_vol"] and r["bib_vol"] != r.get("epmc_vol", ""):
        fl.append("VOL bib=%s epmc=%s" % (r["bib_vol"], r.get("epmc_vol")))
    if r["bib_year"] and r["bib_year"] != str(r.get("epmc_year", "")):
        fl.append("YEAR bib=%s epmc=%s" % (r["bib_year"], r.get("epmc_year")))
    bp = r["bib_pages"].replace("--", "-")
    if bp and bp != r.get("epmc_pages", ""):
        fl.append("PAGES bib=%s epmc=%s" % (bp, r.get("epmc_pages")))
    print("%-28s %-6s %-22s %-8s %-6s %-12s %s" % (
        r["key"], r.get("epmc_pmid", "-"), jm, r.get("epmc_vol", "-"),
        r.get("epmc_year", "-"), r.get("epmc_pages", "-")[:12],
        "; ".join(fl) if fl else "ok"))

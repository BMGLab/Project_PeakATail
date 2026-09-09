import json, sys, datetime
import yaml, jsonschema
from jsonschema import Draft7Validator, FormatChecker

cff_path, schema_path = sys.argv[1], sys.argv[2]

with open(cff_path, encoding="utf-8") as fh:
    doc = yaml.safe_load(fh)
print("YAML: parsed OK; top-level keys =", sorted(doc.keys()))

# CFF schema note: YAML date objects must be cast to str before validating.
def cast_dates(o):
    if isinstance(o, dict):  return {k: cast_dates(v) for k, v in o.items()}
    if isinstance(o, list):  return [cast_dates(v) for v in o]
    if isinstance(o, (datetime.date, datetime.datetime)): return o.isoformat()
    return o
doc = cast_dates(doc)

schema = json.load(open(schema_path, encoding="utf-8"))
v = Draft7Validator(schema, format_checker=FormatChecker())
errs = sorted(v.iter_errors(doc), key=lambda e: list(e.path))
if not errs:
    print("SCHEMA: VALID against CFF 1.2.0 (%d errors)" % 0)
else:
    print("SCHEMA: %d ERROR(S)" % len(errs))
    for e in errs:
        print("  path=%s : %s" % ("/".join(map(str, e.path)) or "<root>", e.message[:200]))

# explicit spot-checks
import re
print("\n-- spot checks --")
print("cff-version       :", doc.get("cff-version"))
print("type              :", doc.get("type"))
print("license (SPDX)    :", doc.get("license"), "in enum:",
      doc.get("license") in schema["definitions"]["license-enum"]["enum"])
print("version           :", repr(doc.get("version")))
print("'doi' key present :", "doi" in doc, "(must be False while placeholder)")
orc = schema["definitions"]["orcid"]["pattern"]
for who in ("authors", "contact"):
    for p in doc.get(who, []):
        if "orcid" in p:
            print("%-9s orcid ok  : %s -> %s" % (who, p["orcid"], bool(re.search(orc, p["orcid"]))))
pc = doc["preferred-citation"]
print("pref-cit type     :", pc["type"], "in enum:",
      pc["type"] in schema["definitions"]["reference"]["properties"]["type"]["enum"])
print("pref-cit status   :", pc["status"], "in enum:",
      pc["status"] in schema["definitions"]["reference"]["properties"]["status"]["enum"])
print("pref-cit required :", all(k in pc for k in ("authors", "title", "type")))
print("no extra top-level:", set(doc) <= set(schema["properties"]))
sys.exit(1 if errs else 0)

#!/usr/bin/env python3
"""Refresh the variant map on the homepage from any VCFART TSV (SNV rows).

usage:
  python docs/build_map.py BRCA1_ex17.tsv --from c.5075-12 --to c.5152+12 \
         --title "BRCA1 exon 17" [--html docs/index.html]

Rewrites the JSON between the MAP-DATA markers in docs/index.html.
"""

import argparse
import csv
import json
import re

POS_RE = re.compile(r"^c\.(.+?)([ACGT])>([ACGT])$")


def build(tsv, start, end, title):
    rows = []
    with open(tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["type"] != "SNV":
                continue
            tx, c = r["hgvs_c"].split(":", 1)
            m = POS_RE.match(c)
            if m:
                rows.append((int(r["pos"]), m.group(1), m.group(2), m.group(3), r, tx))
    if not rows:
        raise SystemExit("no SNV rows in the TSV")
    plain = sorted((g, int(p)) for g, p, *_ in rows if p.isdigit())
    minus = len(plain) > 1 and plain[0][1] > plain[-1][1]
    rows.sort(key=lambda x: -x[0] if minus else x[0])
    cols, index = [], {}
    for g, p, ref_t, alt_t, r, tx in rows:
        if p not in index:
            index[p] = len(cols)
            cols.append({"c": p, "ref": ref_t, "reg": r["region"], "alts": {}})
        csq = r["consequence"].split("&")[0]
        hp = "" if r["hgvs_p"] == "." else r["hgvs_p"]
        cols[index[p]]["alts"][alt_t] = [csq, hp, r["chrom"], g, r["ref"], r["alt"]]
    s, e = start.replace("c.", ""), end.replace("c.", "")
    if s not in index or e not in index:
        raise SystemExit(f"window {start}..{end} not found in the TSV")
    cols = cols[index[s]:index[e] + 1]
    return {"title": title, "tx": rows[0][5], "gene": rows[0][4]["gene"], "cols": cols}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tsv")
    ap.add_argument("--from", dest="start", required=True)
    ap.add_argument("--to", dest="end", required=True)
    ap.add_argument("--title", default="")
    ap.add_argument("--html", default="docs/index.html")
    a = ap.parse_args()
    data = build(a.tsv, a.start, a.end, a.title)
    html = open(a.html, encoding="utf-8").read()
    blob = json.dumps(data, separators=(",", ":"))
    new = re.sub(r"(/\*MAP-DATA\*/).*?(/\*END-MAP-DATA\*/)",
                 lambda m: m.group(1) + blob + m.group(2), html, flags=re.S)
    open(a.html, "w", encoding="utf-8").write(new)
    print(f"{len(data['cols'])} positions written to {a.html}")


if __name__ == "__main__":
    main()

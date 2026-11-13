#!/usr/bin/env python3
"""Cross-check VCFART HGVS against ANNOVAR refGene (exonic SNVs).

usage:
  perl table_annovar.pl X.vcf humandb/ -buildver hg38 -out X_annovar -remove \
       -protocol refGene -operation g -nastring . -vcfinput
  python compare_annovar.py X.tsv X_annovar.hg38_multianno.txt

Matches records through the VCF ID column (VCFART writes the c. description there),
converts VCFART HGVS to ANNOVAR style (c.541G>A -> c.G541A, p.Ala181Thr -> p.A181T,
Ter -> X, '=' -> same residue) and reports agreement for every coding SNV.
"""

import csv
import re
import sys

AA1 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
       "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
       "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "X"}


def c_to_annovar(c):
    m = re.match(r"c\.(\d+)([ACGT])>([ACGT])$", c)
    return f"c.{m.group(2)}{m.group(1)}{m.group(3)}" if m else None


def p_to_annovar(p):
    if p == "p.Met1?":                                   # start loss
        return "p.M1?"
    m = re.match(r"p\.Ter(\d+)([A-Z][a-z]{2})extTer(\d+|\?)$", p)
    if m:                                                # stop loss: ANNOVAR omits the extension
        return f"p.X{m.group(1)}{AA1[m.group(2)]}"
    m = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=)$", p)
    if not m:
        return None
    a, n, b = m.groups()
    return f"p.{AA1[a]}{n}{AA1[a] if b == '=' else AA1[b]}"


def main(tsv, multianno):
    va = {}
    tx_acc = None
    with open(tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tx, c = r["hgvs_c"].split(":", 1)
            tx_acc = tx.split(".")[0]
            va[c] = (r["hgvs_p"], r["consequence"])

    seen = compared = agree = 0
    skipped, diffs = [], []
    with open(multianno) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_aa = header.index("AAChange.refGene")
        for line in fh:
            row = line.rstrip("\n").split("\t")
            vid = row[-8] if row[-2] == "GT" else row[-6]
            if vid not in va:
                continue
            seen += 1
            if not re.match(r"c\.\d+[ACGT]>[ACGT]$", vid):
                continue                                   # not a coding SNV
            vp, csq = va[vid]
            want_c, want_p = c_to_annovar(vid), p_to_annovar(vp)
            ann = [x for x in row[i_aa].split(",") if f":{tx_acc}:" in x or f":{tx_acc}." in x]
            if want_p is None or not ann:
                skipped.append((vid, vp, ann[0] if ann else row[i_aa]))
                continue
            parts = ann[0].split(":")
            got_c = next((x for x in parts if x.startswith("c.")), "")
            got_p = next((x for x in parts if x.startswith("p.")), "")
            compared += 1
            if (got_c, got_p) == (want_c, want_p):
                agree += 1
            else:
                diffs.append((vid, vp, ann[0]))

    print(f"VCFART records found in ANNOVAR output : {seen} / {len(va)}")
    print(f"coding SNVs compared                   : {compared}")
    print(f"  identical c. and p.                  : {agree}")
    print(f"  different                            : {len(diffs)}")
    for d in diffs[:20]:
        print("    DIFF", *d, sep="\t")
    print(f"not comparable (MNV/indel p. etc.)        : {len(skipped)}")
    for s in skipped[:20]:
        print("    SKIP", *s, sep="\t")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])

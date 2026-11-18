# -*- coding: utf-8 -*-
"""What actually depends on 620 vs 602 in the discovery collection."""
import csv
import re
from collections import Counter

import numpy as np
import openpyxl

D = "data"
XLSX = D + "/Scapitis_620_cgMLST_2026.xlsx"
CLUST = D + "/wang2022_cgmlst_clusters_602.csv"
ANCHOR = "GCA_000712995.1"
OUTGROUP = "GCF_040739495.1"

# the six external genomes reported as profile-identical, and their discovery match
SIX = ["SRR7939128", "SRR8786603", "DRR015840", "SRR7939170",
       "GCA_015645265.1", "GCA_015645205.1"]


def canon(v):
    if v is None:
        return None
    if isinstance(v, float) and np.isnan(v):
        return None
    s = str(v).strip()
    if s.lower() in ("", "nan", "none", "-", "?"):
        return None
    m = re.fullmatch(r"(\d+)\.0+", s)
    return m.group(1) if m else s


vcanon = np.vectorize(canon, otypes=[object])

ws = openpyxl.load_workbook(XLSX, read_only=True).active
it = ws.iter_rows(values_only=True)
hdr = list(next(it))
i_id = hdr.index("Sample ID\n(Sample)")
i_st = hdr.index("ST\n(S. capitis MLST)")
cg = [i for i, h in enumerate(hdr) if h and h.endswith("(S. capitis cgMLST)")
      and not h.startswith("Perc.") and not h.startswith("Complex Type")]

ids, sts, prof = [], [], []
for row in it:
    if row[i_id] is None:
        continue
    ids.append(str(row[i_id]).strip())
    sts.append(canon(row[i_st]))
    prof.append([row[c] for c in cg])
P = vcanon(np.array(prof, dtype=object))
idx = {a: i for i, a in enumerate(ids)}

lab = {r["Accession"].strip(): r["Wang_Cluster"].strip()
       for r in csv.DictReader(open(CLUST, newline="", encoding="utf-8-sig"))}
labelled = [a for a in ids if a in lab]
extra = [a for a in ids if a not in lab]

print("=== composition ===")
print("export total        :", len(ids))
print("labelled (the 602)  :", len(labelled))
print("unlabelled (the 18) :", len(extra))
print("  outgroup present  :", OUTGROUP in extra)


def kind(a):
    if a.startswith("GCA"):
        return "GCA"
    if a.startswith("GCF"):
        return "GCF"
    return re.match(r"[A-Z]+", a).group(0)


print("\n=== accession types ===")
print("all 620 :", Counter(kind(a) for a in ids).most_common())
print("602     :", Counter(kind(a) for a in labelled).most_common())
print("18 extra:", Counter(kind(a) for a in extra).most_common())

a = P[idx[ANCHOR]]


def dist_to_anchor(j):
    b = P[j]
    both = (a != None) & (b != None)          # noqa: E711
    return int(np.sum(both & (a != b)))


print("\n=== the 18 unlabelled genomes: distance to LNZR-1 ===")
rows = sorted(((dist_to_anchor(idx[x]), x, sts[idx[x]]) for x in extra))
for d, x, s in rows[:6]:
    print("   %-18s ST=%-5s d=%d" % (x, s or "-", d))
print("   closest of the 18: %d alleles" % rows[0][0])
print("   any within 50 alleles: %s" % [x for d, x, s in rows if d <= 50])

print("\n=== nearest non-Cluster-L genome, 620 vs 602 ===")
L = {x for x, c in lab.items() if c.upper() == "L"}
cand = [(dist_to_anchor(idx[x]), x, sts[idx[x]], x in lab) for x in ids if x not in L]
cand.sort()
for d, x, s, inlab in cand[:4]:
    print("   %-18s ST=%-5s d=%-5d labelled=%s" % (x, s or "-", d, inlab))
print("   nearest within 620 :", cand[0][0])
print("   nearest within 602 :", next(d for d, x, s, i in cand if i))

print("\n=== profile-identical matches: is the discovery counterpart in the 602? ===")
for m in SIX:
    print("   %-18s in export=%-5s labelled=%s"
          % (m, m in idx, m in lab))

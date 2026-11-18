# -*- coding: utf-8 -*-
"""Check the within-Cluster-L nearest-neighbour claim and the 602/620 bookkeeping."""
import csv
import re

import numpy as np
import openpyxl

D = "data"
XLSX = D + "/Scapitis_620_cgMLST_2026.xlsx"
CLUST = D + "/wang2022_cgmlst_clusters_602.csv"


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

rows = list(csv.DictReader(open(CLUST, newline="", encoding="utf-8-sig")))
lab = {r["Accession"].strip(): r["Wang_Cluster"].strip() for r in rows}
L = [a for a, c in lab.items() if c.upper() == "L"]

print("export genomes            :", len(ids))
print("labelled in Wang table    :", len(lab))
print("  Cluster L               :", len(L))
print("  non-Cluster L           :", sum(1 for c in lab.values() if c.upper() != "L"))
print("unlabelled in export      :", sum(1 for a in ids if a not in lab))
print("ST6 in whole export       :", sum(1 for s in sts if s == "6"))
print("ST6 that are Cluster L    :", sum(1 for a in L if sts[idx[a]] == "6"))
print("Cluster L not ST6         :", [a for a in L if sts[idx[a]] != "6"])
print("ST6 not Cluster L         :", [ids[i] for i, s in enumerate(sts)
                                      if s == "6" and lab.get(ids[i], "").upper() != "L"])

Li = [idx[a] for a in L]
n = len(Li)
Dm = np.zeros((n, n), dtype=int)
for x in range(n):
    a = P[Li[x]]
    for y in range(x + 1, n):
        b = P[Li[y]]
        both = (a != None) & (b != None)          # noqa: E711
        d = int(np.sum(both & (a != b)))
        Dm[x, y] = Dm[y, x] = d

off = Dm + np.eye(n, dtype=int) * 10 ** 6
nn = off.min(axis=1)
print("\nwithin-Cluster-L pairwise distance: min %d  max %d  median %.0f"
      % (Dm[np.triu_indices(n, 1)].min(), Dm[np.triu_indices(n, 1)].max(),
         float(np.median(Dm[np.triu_indices(n, 1)]))))
print("nearest-neighbour distance : min %d  max %d  median %.0f"
      % (nn.min(), nn.max(), float(np.median(nn))))
print("nearest-neighbour values   :", sorted(nn.tolist()))
pairs = [(ids[Li[x]], ids[Li[y]]) for x in range(n) for y in range(x + 1, n)
         if Dm[x, y] == 0]
print("identical pairs (distance 0):", len(pairs))
for p in pairs:
    print("   ", p)

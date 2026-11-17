# -*- coding: utf-8 -*-
"""Seven-locus ST and cgMLST distance profile of the 19 discovery Cluster L genomes."""
import csv
import re
import sys
from collections import Counter

import numpy as np
import openpyxl

D = "data/capitis定序後重整之資料"
XLSX = D + "/ridom/Scapitis_620_cgMLST_2026.xlsx"
CLUST = D + "/投稿用/上傳到github/wang2022_cgmlst_clusters_602.csv"
ANCHOR = "GCA_000712995.1"


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

wb = openpyxl.load_workbook(XLSX, read_only=True)
ws = wb[wb.sheetnames[0]]
it = ws.iter_rows(values_only=True)
hdr = list(next(it))

i_id = hdr.index("Sample ID\n(Sample)")
i_st = hdr.index("ST\n(S. capitis MLST)")
i_cc = hdr.index("CC\n(S. capitis MLST)")
cg = [i for i, h in enumerate(hdr) if h and h.endswith("(S. capitis cgMLST)")
      and not h.startswith("Perc.") and not h.startswith("Complex Type")]
print("cgMLST target columns:", len(cg), file=sys.stderr)

ids, sts, ccs, prof = [], [], [], []
for row in it:
    if row[i_id] is None:
        continue
    ids.append(str(row[i_id]).strip())
    sts.append(canon(row[i_st]))
    ccs.append(canon(row[i_cc]))
    prof.append([row[c] for c in cg])

P = vcanon(np.array(prof, dtype=object))
print("genomes:", len(ids), file=sys.stderr)

L = [r["Accession"].strip() for r in
     csv.DictReader(open(CLUST, newline="", encoding="utf-8-sig"))
     if r["Wang_Cluster"].strip().upper() == "L"]

idx = {a: i for i, a in enumerate(ids)}
missing = [a for a in L if a not in idx]
print("Cluster L accessions:", len(L), " not found in export:", missing, file=sys.stderr)

rows = [(a, idx[a]) for a in L if a in idx]

print("\n== seven-locus MLST of the %d Cluster L genomes (Ridom Typer) ==" % len(rows))
print("ST:", Counter(sts[i] or "not assigned" for _, i in rows).most_common())
print("CC:", Counter(ccs[i] or "not assigned" for _, i in rows).most_common())

ai = idx[ANCHOR]
a = P[ai]


def dist(j):
    b = P[j]
    both = (a != None) & (b != None)          # noqa: E711
    return int(np.sum(both & (a != b))), int(np.sum(both))


d = [(acc, sts[i] or "-", ) + dist(i) for acc, i in rows]
print("\n== cgMLST allelic distance to LNZR-1 (%s) ==" % ANCHOR)
for acc, st, dd, sh in sorted(d, key=lambda x: x[2]):
    print("  %-18s ST=%-5s d=%3d  (shared targets %d)" % (acc, st, dd, sh))
vals = [x[2] for x in d]
print("range %d-%d, median %.0f" % (min(vals), max(vals), float(np.median(vals))))

# nearest non-Cluster-L genome
Lset = {idx[x] for x in L if x in idx}
out = [(dist(j)[0], ids[j], sts[j]) for j in range(len(ids)) if j not in Lset]
out.sort()
print("\n== nearest non-Cluster-L genomes ==")
for dd, acc, st in out[:5]:
    print("  %-18s ST=%-5s d=%d" % (acc, st or "-", dd))

# ST composition of the whole discovery set, for context
print("\n== ST composition of all %d discovery genomes ==" % len(ids))
print(Counter(s or "not assigned" for s in sts).most_common(12))

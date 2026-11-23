"""Build Microreact inputs for all 921 genomes.

Usage:
    python3 scripts/microreact.py <620_cgMLST.xlsx> <external_merged.csv> <clusters_602.csv>

Writes into results/microreact/:
    microreact_tree.nwk        UPGMA dendrogram of cgMLST allelic distances
    microreact_metadata.csv    one row per genome, Microreact column conventions
    microreact_matrix.csv      the 921 x 921 allelic distance matrix (optional upload)

The tree is a UPGMA (average-linkage) dendrogram of allelic distances, not a
phylogeny. Distances are counted across the 1,492 shared targets, ignoring any
target missing in either genome, exactly as dist_all.py does.
"""
import csv
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

XLSX, EXT, CLU = sys.argv[1], sys.argv[2], sys.argv[3]

BATCH = "Scapitis_ExtVal_304thirdBatch"
OUTGROUP_SID = "DSM20326"
OUTGROUP_ACC = "GCF_040739495.1"
ANCHOR = "GCA_000712995.1"          # LNZR-1
BORDERLINE = "GCF_039733025.1"      # ST46, 50 alleles from LNZR-1
MARKER_ALLELE = "4"                 # ftsZ cgMLST allele carried by LNZR-1
MISSING = {"", "nan", "none", "-", "?", "na", "n/a"}

OUT = "results/microreact"
os.makedirs(OUT, exist_ok=True)


def canon(v):
    if v is None:
        return None
    if isinstance(v, float) and v != v:
        return None
    s = str(v).strip()
    if s.lower() in MISSING:
        return None
    m = re.fullmatch(r"(\d+)\.0+", s)
    return m.group(1) if m else s


vcanon = np.vectorize(canon, otypes=[object])

# ---- load both allele matrices on the shared targets ----------------------
x = pd.read_excel(XLSX, dtype=str)
x = x.rename(columns={c: " ".join(str(c).split()) for c in x.columns})
S6 = " (S. capitis cgMLST)"
t620 = {c[:-len(S6)]: c for c in x.columns
        if c.endswith(S6) and not c.startswith(("Complex Type", "Perc."))}
sid620 = [c for c in x.columns if c.startswith("Sample ID")][0]
st620 = [c for c in x.columns if c.startswith("ST (")][0]

m = pd.read_csv(EXT, dtype=str, low_memory=False)
m = m[m.SourceBatch == BATCH].reset_index(drop=True)
m = m[m["Sample_ID_(Sample)"] != OUTGROUP_ACC].reset_index(drop=True)
SM = "_(S._capitis_cgMLST)"
tM = {c[:-len(SM)]: c for c in m.columns
      if c.endswith(SM) and not c.startswith(("Complex_Type", "Perc."))}

shared = sorted(set(t620) & set(tM))
assert len(shared) == 1492, len(shared)

D = vcanon(x[[t620[t] for t in shared]].to_numpy(dtype=object))
E = vcanon(m[[tM[t] for t in shared]].to_numpy(dtype=object))
ids = list(x[sid620].astype(str)) + list(m["Sample_ID_(Sample)"].astype(str))
coll = ["Discovery"] * len(D) + ["External"] * len(E)
sts = list(vcanon(x[st620].to_numpy(dtype=object))) + \
      list(vcanon(m["ST_(S._capitis_MLST)"].to_numpy(dtype=object)))
A = np.vstack([D, E])
n = len(ids)
print("genomes: %d (%d discovery, %d external) x %d targets"
      % (n, len(D), len(E), len(shared)))
assert n == 921, n

# ---- integer-encode each column so distances are cheap -------------------
M = np.full(A.shape, -1, dtype=np.int32)
for j in range(A.shape[1]):
    col = A[:, j]
    codes = {}
    for i, v in enumerate(col):
        if v is None:
            continue
        if v not in codes:
            codes[v] = len(codes)
        M[i, j] = codes[v]

# ---- pairwise allelic distance, ignoring targets missing in either -------
dist = np.zeros((n, n), dtype=np.int32)
comp = np.zeros((n, n), dtype=np.int32)
present = M >= 0
for i in range(n):
    both = present & present[i]
    diff = (M != M[i]) & both
    dist[i] = diff.sum(1)
    comp[i] = both.sum(1)
np.fill_diagonal(dist, 0)
assert (dist == dist.T).all(), "distance matrix not symmetric"
print("distance matrix done; off-diagonal range %d-%d"
      % (dist[~np.eye(n, dtype=bool)].min(), dist.max()))

# ---- UPGMA dendrogram -> Newick -----------------------------------------
Z = linkage(squareform(dist.astype(float), checks=False), method="average")


def to_newick(Z, labels):
    """Iterative Newick conversion, safe for deep trees."""
    n = len(labels)
    parts = {i: labels[i] for i in range(n)}
    height = {i: 0.0 for i in range(n)}
    for k, (a, b, d, _) in enumerate(Z):
        a, b = int(a), int(b)
        la = max(d / 2.0 - height[a], 0.0)
        lb = max(d / 2.0 - height[b], 0.0)
        parts[n + k] = "(%s:%.4f,%s:%.4f)" % (parts[a], la, parts[b], lb)
        height[n + k] = d / 2.0
        del parts[a], parts[b]
    return parts[n + len(Z) - 1] + ";"


nwk = to_newick(Z, ids)
open(os.path.join(OUT, "microreact_tree.nwk"), "w").write(nwk + "\n")
print("newick written: %d leaves, %d chars" % (nwk.count(":") - nwk.count("),") , len(nwk)))

# ---- metadata ------------------------------------------------------------
labels = {}
with open(CLU, newline="") as fh:
    for row in csv.DictReader(fh):
        labels[canon(row["Accession"])] = (canon(row["Wang_Cluster"]),
                                           canon(row["Country"]),
                                           canon(row["Year"]))
NO_LABEL = ("Not assigned", "Outgroup")

ftsz = {}
with open("results/tables/ftsz_allele_assignments.csv", newline="") as fh:
    for row in csv.DictReader(fh):
        ftsz[row["id"]] = canon(row["ftsZ"])

ai = ids.index(ANCHOR)
COLL_C = {"Discovery": "#0072B2", "External": "#D55E00"}
MARK_C = {"Carrier": "#CC79A7", "Non-carrier": "#BBBBBB", "Not typed": "#EEEEEE"}
LIN_C = {"Cluster L (discovery)": "#009E73", "ST6 (external)": "#009E73",
         "ST46, 50 alleles from LNZR-1": "#E69F00",
         "Outgroup": "#000000", "Non-lineage": "#BBBBBB"}

rows = []
for i, gid in enumerate(ids):
    cluster, country, year = labels.get(gid, (None, None, None))
    if cluster in NO_LABEL:
        cluster = None
    st = sts[i]
    if gid == OUTGROUP_SID:
        lineage = "Outgroup"
    elif cluster == "L":
        lineage = "Cluster L (discovery)"
    elif coll[i] == "External" and st == "6":
        lineage = "ST6 (external)"
    elif gid == BORDERLINE:
        lineage = "ST46, 50 alleles from LNZR-1"
    else:
        lineage = "Non-lineage"
    a = ftsz.get(gid)
    mark = "Not typed" if a is None else ("Carrier" if a == MARKER_ALLELE else "Non-carrier")
    rows.append({
        "id": gid,
        "Accession": OUTGROUP_ACC if gid == OUTGROUP_SID else gid,
        "Collection": coll[i],
        "Collection__colour": COLL_C[coll[i]],
        "Published_cgMLST_cluster": cluster or "",
        "Seven_locus_ST": st or "",
        "Lineage": lineage,
        "Lineage__colour": LIN_C[lineage],
        "ftsZ_marker_allele": mark,
        "ftsZ_marker_allele__colour": MARK_C[mark],
        "ftsZ_cgMLST_allele": a or "",
        "Distance_to_LNZR1": int(dist[i, ai]),
        "Targets_compared_with_LNZR1": int(comp[i, ai]),
        "Country": country or "",
        "Year": year or "",
    })

md = pd.DataFrame(rows)
md.to_csv(os.path.join(OUT, "microreact_metadata.csv"), index=False)

pd.DataFrame(dist, index=ids, columns=ids).to_csv(
    os.path.join(OUT, "microreact_matrix.csv"))

print("\nmetadata: %d rows x %d cols" % md.shape)
print(md["Collection"].value_counts().to_string())
print()
print(md["Lineage"].value_counts().to_string())
print()
print(md["ftsZ_marker_allele"].value_counts().to_string())
print("\ncarriers by collection:")
print(md[md.ftsZ_marker_allele == "Carrier"].groupby("Collection").size().to_string())

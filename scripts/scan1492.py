import sys, re, pandas as pd, numpy as np, warnings
warnings.filterwarnings("ignore")
x620_p, clus_p = sys.argv[1], sys.argv[2]

def canon(v):
    if v is None: return None
    if isinstance(v, float) and np.isnan(v): return None
    s = str(v).strip()
    if s.lower() in ("", "nan", "none", "-", "?"): return None
    m = re.fullmatch(r"(\d+)\.0+", s)
    return m.group(1) if m else s
vc = np.vectorize(canon, otypes=[object])

x = pd.read_excel(x620_p, dtype=str)
x = x.rename(columns={c: " ".join(str(c).split()) for c in x.columns})
SUF = " (S. capitis cgMLST)"
tcols = [c for c in x.columns if c.endswith(SUF) and not c.startswith(("Complex Type","Perc."))]
sid = [c for c in x.columns if c.startswith("Sample ID")][0]
ids = x[sid].astype(str).to_numpy()

cl = pd.read_csv(clus_p, dtype=str)
print("clusters csv cols:", list(cl.columns), "rows:", len(cl))
acc = [c for c in cl.columns if "acc" in c.lower() or "genome" in c.lower() or "id" in c.lower()][0]
clc = [c for c in cl.columns if "clust" in c.lower()][0]
cmap = dict(zip(cl[acc].astype(str).str.strip(), cl[clc].astype(str).str.strip()))
lab = np.array([cmap.get(i) for i in ids], dtype=object)
keep = np.array([v is not None and v not in ("Not assigned","Outgroup","nan") for v in lab])
isL  = np.array([v == "L" for v in lab])
print("labeled:", int(keep.sum()), " Cluster L:", int((isL & keep).sum()),
      " non-L:", int((keep & ~isL).sum()))

A = vc(x[tcols].to_numpy(dtype=object))
Lrows = np.where(isL & keep)[0]
res = []
for j, c in enumerate(tcols):
    col = A[:, j]
    lv = {col[i] for i in Lrows}
    if len(lv) != 1 or None in lv:
        continue                      # lineage not uniform, or missing in a lineage genome
    a = lv.pop()
    carr = np.array([v == a for v in col]) & keep
    TP = int((carr & isL).sum()); FP = int((carr & ~isL).sum())
    res.append((c[:-len(SUF)], a, TP+FP, FP))
r = pd.DataFrame(res, columns=["target","allele","carriers","false_pos"]).sort_values("false_pos")
nonL = int((keep & ~isL).sum()); nL = int((isL & keep).sum())
r["specificity_%"] = (100*(nonL - r.false_pos)/nonL).round(2)
r["PPV_%"] = (100*nL/r.carriers).round(1)
r.to_csv("results/tables/cgmlst_1492_screen.csv", index=False)
print("\ntargets with a 100%%-sensitivity candidate allele: %d of %d" % (len(r), len(tcols)))
print("\nTOP 15 by fewest false positives:")
print(r.head(15).to_string(index=False))
print("\nwhere does ftsZ rank?")
f = r[r.target.str.lower().str.contains("ftsz")]
print(f.to_string(index=False))
print("rank of ftsZ (1 = fewest FP):", (r.false_pos.values < f.false_pos.values[0]).sum() + 1 if len(f) else "n/a")
print("\nn targets with <= 11 false positives:", int((r.false_pos <= 11).sum()))

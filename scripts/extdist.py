import sys, re, pandas as pd, numpy as np, warnings
warnings.filterwarnings("ignore")
x620_p, merged_p = sys.argv[1], sys.argv[2]

MISSING = {"", "nan", "none", "-", "?", "na", "n/a"}

def canon(v):
    # Must handle float NaN explicitly: pandas Series.map() does not pass NaN
    # to the callable, so relying on str(v) == "nan" silently leaves missing
    # targets in place and inflates every pairwise distance. Kept identical to
    # scripts/dist_all.py so both distance tables agree.
    if v is None: return None
    if isinstance(v, float) and np.isnan(v): return None
    s = str(v).strip()
    if s.lower() in MISSING: return None
    m = re.fullmatch(r"(\d+)\.0+", s)
    return m.group(1) if m else s

vcanon = np.vectorize(canon, otypes=[object])

# --- discovery export ---
x = pd.read_excel(x620_p, dtype=str)
x = x.rename(columns={c: " ".join(str(c).split()) for c in x.columns})
SUF620 = " (S. capitis cgMLST)"
t620 = {c[:-len(SUF620)]: c for c in x.columns
        if c.endswith(SUF620) and not c.startswith(("Complex Type", "Perc."))}
sid620 = [c for c in x.columns if c.startswith("Sample ID")][0]

# --- external export ---
m = pd.read_csv(merged_p, dtype=str, low_memory=False)
m = m[m.SourceBatch == "Scapitis_ExtVal_304thirdBatch"].copy()
SUFM = "_(S._capitis_cgMLST)"
tM = {c[:-len(SUFM)]: c for c in m.columns
      if c.endswith(SUFM) and not c.startswith(("Complex_Type", "Perc."))}
sidM = "Sample_ID_(Sample)"

print("targets 620:", len(t620), " targets merged:", len(tM))
shared = sorted(set(t620) & set(tM))
print("shared target names:", len(shared))
only620 = sorted(set(t620) - set(tM))[:5]
onlyM  = sorted(set(tM) - set(t620))[:5]
print("only in 620 (sample):", only620)
print("only in merged (sample):", onlyM)

# LNZR-1 reference vector over shared targets
r = x.loc[x[sid620] == "GCA_000712995.1"]
assert len(r) == 1
ref = vcanon(np.array([r.iloc[0][t620[t]] for t in shared], dtype=object))

ids  = m[sidM].astype(str).to_numpy()
st   = vcanon(m["ST_(S._capitis_MLST)"].to_numpy(dtype=object))
E    = vcanon(m[[tM[t] for t in shared]].to_numpy(dtype=object))

rows = []
for i in range(len(ids)):
    a = E[i]
    both = np.fromiter(((p is not None and q is not None) for p, q in zip(a, ref)), bool, len(ref))
    diff = np.fromiter(((p != q) for p, q in zip(a, ref)), bool, len(ref))
    rows.append((ids[i], st[i], int((both & diff).sum()), int(both.sum())))
df = pd.DataFrame(rows, columns=["id","ST","dist","n_compared"])
df = df[df.id != "GCF_040739495.1"].reset_index(drop=True)   # drop DSM20326 outgroup
df.to_csv("results/tables/external_dist.csv", index=False)
print("\nexternal genomes after dropping DSM20326:", len(df))
print("shared-target coverage: min %d, max %d, mean %.0f" % (df.n_compared.min(), df.n_compared.max(), df.n_compared.mean()))

st6 = df[df.ST == "6"].sort_values("dist")
print("\nST6 n=%d  distance range: %d - %d" % (len(st6), st6.dist.min(), st6.dist.max()))
print("ST6 distances:", sorted(st6.dist.tolist()))
print("\nall genomes sorted by distance, first 12:")
print(df.sort_values("dist").head(12)[["id","ST","dist"]].to_string(index=False))
print("\nnon-ST6 nearest 6:")
print(df[df.ST != "6"].sort_values("dist").head(6)[["id","ST","dist"]].to_string(index=False))
print("\nunassigned-ST genomes: n=%d, min dist=%d" % ((df.ST.isna()).sum(), df[df.ST.isna()].dist.min()))
print("within 50 alleles:", df[df.dist <= 50][["id","ST","dist"]].to_string(index=False))

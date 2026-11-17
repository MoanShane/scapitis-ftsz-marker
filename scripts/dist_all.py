import sys, re, pandas as pd, numpy as np, warnings
warnings.filterwarnings("ignore")
x620_p, merged_p = sys.argv[1], sys.argv[2]
MISSING = {"", "nan", "none", "-", "?", "na", "n/a"}

def canon(v):
    if v is None: return None
    if isinstance(v, float) and np.isnan(v): return None
    s = str(v).strip()
    if s.lower() in MISSING: return None
    m = re.fullmatch(r"(\d+)\.0+", s)
    return m.group(1) if m else s

vcanon = np.vectorize(canon, otypes=[object])

# ---------- discovery ----------
x = pd.read_excel(x620_p, dtype=str)
x = x.rename(columns={c: " ".join(str(c).split()) for c in x.columns})
S6 = " (S. capitis cgMLST)"
t620 = {c[:-len(S6)]: c for c in x.columns if c.endswith(S6) and not c.startswith(("Complex Type","Perc."))}
sid620 = [c for c in x.columns if c.startswith("Sample ID")][0]
st620  = [c for c in x.columns if c.startswith("ST (")][0]

# ---------- external ----------
m = pd.read_csv(merged_p, dtype=str, low_memory=False)
m = m[m.SourceBatch == "Scapitis_ExtVal_304thirdBatch"].reset_index(drop=True)
SM = "_(S._capitis_cgMLST)"
tM = {c[:-len(SM)]: c for c in m.columns if c.endswith(SM) and not c.startswith(("Complex_Type","Perc."))}
sidM = "Sample_ID_(Sample)"

shared = sorted(set(t620) & set(tM))
assert len(shared) == 1492, len(shared)

D = vcanon(x[[t620[t] for t in shared]].to_numpy(dtype=object))
E = vcanon(m[[tM[t]  for t in shared]].to_numpy(dtype=object))
dids = x[sid620].astype(str).to_numpy(); dst = vcanon(x[st620].to_numpy(dtype=object))
eids = m[sidM].astype(str).to_numpy();   est = vcanon(m["ST_(S._capitis_MLST)"].to_numpy(dtype=object))

ai = int(np.where(dids == "GCA_000712995.1")[0][0])
ref = D[ai]
refok = np.array([v is not None for v in ref])
print("LNZR-1 non-missing targets:", refok.sum(), "of 1492")

def dists(M, ids, sts):
    out = []
    for i in range(len(ids)):
        a = M[i]
        ok = np.fromiter((v is not None for v in a), bool, len(a)) & refok
        df_ = np.fromiter((p != q for p, q in zip(a, ref)), bool, len(a))
        out.append((ids[i], sts[i], int((ok & df_).sum()), int(ok.sum())))
    return pd.DataFrame(out, columns=["id","ST","dist","n_compared"])

dd = dists(D, dids, dst); dd["collection"] = "Discovery"
ee = dists(E, eids, est); ee["collection"] = "External"
ee = ee[ee.id != "GCF_040739495.1"].reset_index(drop=True)

for nm, df in (("DISCOVERY", dd), ("EXTERNAL", ee)):
    print(f"\n===== {nm} (n={len(df)}) =====")
    print("coverage: min %d max %d mean %.0f" % (df.n_compared.min(), df.n_compared.max(), df.n_compared.mean()))
    s6 = df[df.ST == "6"].sort_values("dist")
    print("ST6 n=%d range %d-%d" % (len(s6), s6.dist.min(), s6.dist.max()))
    print("ST6:", sorted(s6.dist.tolist()))
    nn = df[df.ST != "6"].sort_values("dist")
    print("nearest non-ST6:")
    print(nn.head(4)[["id","ST","dist"]].to_string(index=False))
    print("within 50 alleles: n=%d  (ST breakdown: %s)" % (
        (df.dist<=50).sum(), df[df.dist<=50].ST.value_counts(dropna=False).to_dict()))

pd.concat([dd, ee], ignore_index=True).to_csv("results/tables/dist_both.csv", index=False)
print("\nwrote results/tables/dist_both.csv")

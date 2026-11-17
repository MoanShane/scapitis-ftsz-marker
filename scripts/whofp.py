import sys, re, pandas as pd, numpy as np, warnings
warnings.filterwarnings("ignore")
def canon(v):
    if v is None: return None
    if isinstance(v,float) and np.isnan(v): return None
    s=str(v).strip()
    if s.lower() in ("","nan","none","-","?"): return None
    m=re.fullmatch(r"(\d+)\.0+",s); return m.group(1) if m else s
vc=np.vectorize(canon,otypes=[object])
x=pd.read_excel(sys.argv[1],dtype=str)
x=x.rename(columns={c:" ".join(str(c).split()) for c in x.columns})
S6=" (S. capitis cgMLST)"
t620={c[:-len(S6)]:c for c in x.columns if c.endswith(S6) and not c.startswith(("Complex Type","Perc."))}
sid=[c for c in x.columns if c.startswith("Sample ID")][0]
ids=x[sid].astype(str).to_numpy()
cl=pd.read_csv(sys.argv[3],dtype=str)
cmap=dict(zip(cl["Accession"].astype(str).str.strip(),cl["Wang_Cluster"].astype(str).str.strip()))
lab=np.array([cmap.get(i) for i in ids],dtype=object)
keep=np.array([v is not None and v not in ("Not assigned","Outgroup","nan") for v in lab])
isL=np.array([v=="L" for v in lab])
m=pd.read_csv(sys.argv[2],dtype=str,low_memory=False)
m=m[m.SourceBatch=="Scapitis_ExtVal_304thirdBatch"].reset_index(drop=True)
m=m[m["Sample_ID_(Sample)"]!="GCF_040739495.1"].reset_index(drop=True)
SM="_(S._capitis_cgMLST)"
tM={c[:-len(SM)]:c for c in m.columns if c.endswith(SM) and not c.startswith(("Complex_Type","Perc."))}
eids=m["Sample_ID_(Sample)"].astype(str).to_numpy()
est=vc(m["ST_(S._capitis_MLST)"].to_numpy(dtype=object))
eL=np.array([v=="6" for v in est])
full=pd.read_csv("results/tables/cgmlst_target_full.csv")
perf=full[full.disc_FP==0].target.tolist()
shared=sorted(set(t620)&set(tM))
D=vc(x[[t620[t] for t in shared]].to_numpy(dtype=object))
E=vc(m[[tM[t] for t in shared]].to_numpy(dtype=object))
idx={t:j for j,t in enumerate(shared)}
Lrows=np.where(isL&keep)[0]
from collections import Counter
cnt=Counter()
for t in perf:
    j=idx[t]; a={D[i,j] for i in Lrows}.pop()
    fp=[eids[i] for i in range(len(eids)) if (E[i,j]==a) and not eL[i]]
    for f in fp: cnt[f]+=1
print("external genomes flagged by the 26 discovery-perfect targets:")
for g,n in cnt.most_common():
    print(f"   {g}: flagged by {n} of 26 targets")
# ftsZ
j=idx["ftsZ"]; a={D[i,j] for i in Lrows}.pop()
print("\nftsZ (allele %s) external false positives:" % a)
for i in range(len(eids)):
    if E[i,j]==a and not eL[i]:
        print("   ", eids[i], "ST=", est[i])
dist=pd.read_csv("results/tables/dist_both.csv")
dm=dict(zip(dist[dist.collection=="External"].id, dist[dist.collection=="External"].dist))
print("\ncgMLST distance to LNZR-1 of those genomes:")
for g in set(list(cnt)+[eids[i] for i in range(len(eids)) if E[i,j]==a and not eL[i]]):
    print("   ", g, "=", dm.get(g), "alleles | ST =", est[list(eids).index(g)])

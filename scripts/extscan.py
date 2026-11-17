import sys, re, pandas as pd, numpy as np, warnings
warnings.filterwarnings("ignore")
x620_p, merged_p, clus_p = sys.argv[1], sys.argv[2], sys.argv[3]
def canon(v):
    if v is None: return None
    if isinstance(v, float) and np.isnan(v): return None
    s = str(v).strip()
    if s.lower() in ("","nan","none","-","?"): return None
    m = re.fullmatch(r"(\d+)\.0+", s); return m.group(1) if m else s
vc = np.vectorize(canon, otypes=[object])

x = pd.read_excel(x620_p, dtype=str)
x = x.rename(columns={c:" ".join(str(c).split()) for c in x.columns})
S6=" (S. capitis cgMLST)"
t620={c[:-len(S6)]:c for c in x.columns if c.endswith(S6) and not c.startswith(("Complex Type","Perc."))}
sid=[c for c in x.columns if c.startswith("Sample ID")][0]
ids=x[sid].astype(str).to_numpy()
cl=pd.read_csv(clus_p,dtype=str)
cmap=dict(zip(cl["Accession"].astype(str).str.strip(), cl["Wang_Cluster"].astype(str).str.strip()))
lab=np.array([cmap.get(i) for i in ids],dtype=object)
keep=np.array([v is not None and v not in ("Not assigned","Outgroup","nan") for v in lab])
isL=np.array([v=="L" for v in lab])

m=pd.read_csv(merged_p,dtype=str,low_memory=False)
m=m[m.SourceBatch=="Scapitis_ExtVal_304thirdBatch"].reset_index(drop=True)
m=m[m["Sample_ID_(Sample)"]!="GCF_040739495.1"].reset_index(drop=True)
SM="_(S._capitis_cgMLST)"
tM={c[:-len(SM)]:c for c in m.columns if c.endswith(SM) and not c.startswith(("Complex_Type","Perc."))}
est=vc(m["ST_(S._capitis_MLST)"].to_numpy(dtype=object))
eL=np.array([v=="6" for v in est])
print("external n=%d, ST6=%d" % (len(m), int(eL.sum())))

shared=sorted(set(t620)&set(tM))
D=vc(x[[t620[t] for t in shared]].to_numpy(dtype=object))
E=vc(m[[tM[t]  for t in shared]].to_numpy(dtype=object))
Lrows=np.where(isL&keep)[0]
nonL=int((keep&~isL).sum())
rows=[]
for j,t in enumerate(shared):
    col=D[:,j]; lv={col[i] for i in Lrows}
    if len(lv)!=1 or None in lv: continue
    a=lv.pop()
    carr=np.array([v==a for v in col])&keep
    dFP=int((carr&~isL).sum())
    ecol=E[:,j]
    ecarr=np.array([v==a for v in ecol])
    eTP=int((ecarr&eL).sum()); eFP=int((ecarr&~eL).sum()); eFN=int((~ecarr&eL).sum())
    emiss=int(sum(v is None for v in ecol[eL]))
    rows.append((t,a,dFP,eTP,eFN,eFP,emiss))
r=pd.DataFrame(rows,columns=["target","allele","disc_FP","ext_TP","ext_FN","ext_FP","ext_missing_in_ST6"])
r["ext_sens_%"]=(100*r.ext_TP/(r.ext_TP+r.ext_FN)).round(1)
r.to_csv("results/tables/cgmlst_target_external.csv",index=False)

perf=r[r.disc_FP==0]
print("\ntargets perfect in discovery (0 FP):", len(perf))
print("  of those, external sensitivity 100%%: %d" % int((perf.ext_FN==0).sum()))
print("  of those, also 0 external FP:        %d" % int(((perf.ext_FN==0)&(perf.ext_FP==0)).sum()))
print("\ndistribution of external sensitivity among the discovery-perfect targets:")
print(perf["ext_sens_%"].describe().round(1).to_string())
print("\nBEST survivors (0 disc FP, 0 ext FN, fewest ext FP):")
print(perf[(perf.ext_FN==0)].sort_values("ext_FP").head(12).to_string(index=False))
print("\nftsZ:")
print(r[r.target.str.lower()=="ftsz"].to_string(index=False))

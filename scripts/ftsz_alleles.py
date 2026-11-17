import sys, re, pandas as pd, numpy as np, warnings
warnings.filterwarnings("ignore")
def canon(v):
    if v is None: return None
    if isinstance(v,float) and np.isnan(v): return None
    s=str(v).strip()
    if s.lower() in ("","nan","none","-","?"): return None
    m=re.fullmatch(r"(\d+)\.0+",s); return m.group(1) if m else s

x=pd.read_excel(sys.argv[1],dtype=str)
x=x.rename(columns={c:" ".join(str(c).split()) for c in x.columns})
sid=[c for c in x.columns if c.startswith("Sample ID")][0]
fz=[c for c in x.columns if c.strip().startswith("ftsZ")][0]
d=pd.DataFrame({"id":x[sid].astype(str),"ftsZ":x[fz].map(canon),"set":"discovery"})

m=pd.read_csv(sys.argv[2],dtype=str,low_memory=False)
m=m[m.SourceBatch=="Scapitis_ExtVal_304thirdBatch"]
m=m[m["Sample_ID_(Sample)"]!="GCF_040739495.1"]
e=pd.DataFrame({"id":m["Sample_ID_(Sample)"].astype(str),
                "ftsZ":m["ftsZ_(S._capitis_cgMLST)"].map(canon),"set":"external"})
a=pd.concat([d,e],ignore_index=True)
a.to_csv("results/tables/ftsz_allele_assignments.csv",index=False)

print("genomes with a cgMLST ftsZ call:", a.ftsZ.notna().sum(), "of", len(a))
for s in ("discovery","external"):
    sub=a[a.set==s]
    print(f"  {s}: n={len(sub)}, distinct ftsZ alleles={sub.ftsZ.nunique()}")
print("combined distinct cgMLST ftsZ alleles:", a.ftsZ.nunique())
print()
print("top 10 allele frequencies (combined):")
print(a.ftsZ.value_counts(dropna=False).head(10).to_string())
print()
print("genomes carrying the marker allele 4:", int((a.ftsZ=="4").sum()),
      f"(discovery {int(((a.set=='discovery')&(a.ftsZ=='4')).sum())}, "
      f"external {int(((a.set=='external')&(a.ftsZ=='4')).sum())})")

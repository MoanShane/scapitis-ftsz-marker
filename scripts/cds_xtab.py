import sys, re, collections, pandas as pd
fa, assign = sys.argv[1], sys.argv[2]

seqs={}
acc=None; buf=[]
for line in open(fa):
    line=line.rstrip()
    if line.startswith(">"):
        if acc: seqs[acc]="".join(buf)
        acc=line[1:].split()[0]; buf=[]
    else: buf.append(line)
if acc: seqs[acc]="".join(buf)
print("sequences read:", len(seqs))
print("lengths:", collections.Counter(len(s) for s in seqs.values()).most_common())
print("distinct CDS:", len(set(seqs.values())))
print("example accessions:", list(seqs)[:5])

a=pd.read_csv(assign, dtype=str)
a["base"]=a.id.str.replace(r"\.\d+$","",regex=True)
by_exact=dict(zip(a.id, zip(a.ftsZ,a.set)))
by_base =dict(zip(a.base, zip(a.ftsZ,a.set)))

rows=[]
for k,s in seqs.items():
    kb=re.sub(r"\.\d+$","",k)
    hit = by_exact.get(k) or by_base.get(kb)
    rows.append((k, s, hit[0] if hit else None, hit[1] if hit else None))
df=pd.DataFrame(rows, columns=["acc","cds","ftsZ","set"])
print()
print("matched to a cgMLST ftsZ call:", df.ftsZ.notna().sum(), "of", len(df))
print("by set:", df[df.ftsZ.notna()].set.value_counts().to_dict())

ok=df[df.ftsZ.notna()].copy()
cdsid={s:f"C{i+1:02d}" for i,s in enumerate(sorted(ok.cds.unique(), key=lambda z:-(ok.cds==z).sum()))}
ok["cds_type"]=ok.cds.map(cdsid)
print("distinct CDS among matched:", ok.cds_type.nunique(), "| distinct cgMLST alleles among matched:", ok.ftsZ.nunique())
print()
print("=== does each cgMLST allele map to exactly one CDS type? ===")
g=ok.groupby("ftsZ").cds_type.nunique().sort_values(ascending=False)
print("  alleles splitting into >1 CDS type:", int((g>1).sum()), "of", len(g))
print(g.head(8).to_string())
print()
print("=== does each CDS type map to exactly one cgMLST allele? ===")
h=ok.groupby("cds_type").ftsZ.nunique().sort_values(ascending=False)
print("  CDS types covering >1 allele:", int((h>1).sum()), "of", len(h))
print(h.head(8).to_string())
print()
print("=== THE MARKER: cgMLST ftsZ allele 4 ===")
m=ok[ok.ftsZ=="4"]
print("  genomes:", len(m), "| CDS types:", m.cds_type.nunique(), "->", m.cds_type.value_counts().to_dict())
print("  by set:", m.set.value_counts().to_dict())
mt=m.cds_type.value_counts().index[0]
others=ok[ok.cds_type!=mt]
print("  genomes with the marker CDS type but another allele:",
      int((ok[(ok.cds_type==mt)&(ok.ftsZ!='4')]).shape[0]))
ok.drop(columns=["cds"]).to_csv("results/tables/cds_vs_allele.csv", index=False)
import pickle; pickle.dump({"ok":ok, "cdsid":cdsid}, open("results/tables/cds.pkl","wb"))
print("\nwrote results/tables/cds_vs_allele.csv")

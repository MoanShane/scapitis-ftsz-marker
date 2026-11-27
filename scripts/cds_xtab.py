import sys, re, collections, pandas as pd
fa, assign = sys.argv[1], sys.argv[2]

def genome_key(name):
    """One key per genome, whatever the ID style.

    The CDS FASTA stores 134 discovery assemblies twice: once as a bare
    accession (GCA_000174135.1) and once as the download filename
    (GCA_000174135.1_ASM17413v1_genomic). Keying on the bare accession alone
    silently skipped the filename copies; accepting both without collapsing
    them would count those genomes twice. Key on the accession so that each
    genome is counted exactly once.
    """
    m = re.match(r"(GC[AF]_\d+)\.\d+", name)
    return m.group(1) if m else re.sub(r"\.\d+$", "", name)


raw={}
acc=None; buf=[]
for line in open(fa):
    line=line.rstrip()
    if line.startswith(">"):
        if acc: raw[acc]="".join(buf)
        acc=line[1:].split()[0]; buf=[]
    else: buf.append(line)
if acc: raw[acc]="".join(buf)
seqs={}; clash=[]
for name,s in raw.items():
    k=genome_key(name)
    if k in seqs and seqs[k]!=s: clash.append(k)
    seqs[k]=s
if clash:
    raise SystemExit("duplicate FASTA entries disagree for: %s" % sorted(set(clash))[:5])
print("FASTA entries: %d | distinct genomes: %d | duplicate entries collapsed: %d"
      % (len(raw), len(seqs), len(raw)-len(seqs)))
print("sequences read:", len(seqs))
print("lengths:", collections.Counter(len(s) for s in seqs.values()).most_common())
print("distinct CDS:", len(set(seqs.values())))
print("example accessions:", list(seqs)[:5])

a=pd.read_csv(assign, dtype=str)
a["key"]=a.id.map(genome_key)
by_key=dict(zip(a.key, zip(a.ftsZ,a.set)))

rows=[]
for k,s in seqs.items():
    hit = by_key.get(k)
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

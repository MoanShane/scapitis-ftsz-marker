import sys, re, collections, pandas as pd, numpy as np
from math import sqrt


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


def read_cds(path):
    """Read the CDS FASTA, collapsing duplicate entries for the same genome."""
    raw = {}
    acc = None
    buf = []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if acc:
                raw[acc] = "".join(buf)
            acc = line[1:].split()[0]
            buf = []
        else:
            buf.append(line)
    if acc:
        raw[acc] = "".join(buf)
    seqs, clash = {}, []
    for name, s in raw.items():
        k = genome_key(name)
        if k in seqs and seqs[k] != s:
            clash.append(k)
        seqs[k] = s
    if clash:
        raise SystemExit("duplicate FASTA entries disagree for: %s"
                         % sorted(set(clash))[:5])
    return raw, seqs


raw, seqs = read_cds(sys.argv[1])
print("FASTA entries: %d | distinct genomes: %d | duplicate entries collapsed: %d"
      % (len(raw), len(seqs), len(raw) - len(seqs)))
a=pd.read_csv("results/tables/ftsz_allele_assignments.csv",dtype=str)
a["key"]=a.id.map(genome_key)
bykey=dict(zip(a.key,zip(a.ftsZ,a.set)))
rows=[]
for k,s in seqs.items():
    h=bykey.get(k)
    if h: rows.append((k,s,h[0],h[1]))
df=pd.DataFrame(rows,columns=["acc","cds","ftsZ","set"])
cnt=df.cds.value_counts()
cid={s:f"H{i+1:02d}" for i,s in enumerate(cnt.index)}
df["hap"]=df.cds.map(cid)
print("genomes:", len(df), "| CDS haplotypes:", df.hap.nunique(), "| cgMLST alleles:", df.ftsZ.nunique())
print()
print("=== allele -> haplotype ===")
g=df.groupby("ftsZ").hap.nunique()
print("  alleles splitting into >1 haplotype:", int((g>1).sum()), "of", len(g))
print("=== haplotype -> allele ===")
h=df.groupby("hap").ftsZ.nunique()
print("  haplotypes covering >1 allele:", int((h>1).sum()), "of", len(h))
print("  (listing those):")
for k,v in h[h>1].items():
    print("   ", k, "->", df[df.hap==k].ftsZ.value_counts().to_dict())
print()
print("=== THE MARKER, cgMLST allele 4 ===")
m=df[df.ftsZ=="4"]
print("  genomes:", len(m), "| haplotypes:", m.hap.value_counts().to_dict(), "| by set:", m.set.value_counts().to_dict())
MH=m.hap.value_counts().index[0]
print("  marker haplotype:", MH, "| genomes carrying it:", int((df.hap==MH).sum()),
      "| their alleles:", df[df.hap==MH].ftsZ.value_counts().to_dict())

cl=pd.read_csv("data/wang2022_cgmlst_clusters_602.csv",dtype=str)
# keyed the same way as the sequence IDs, so the join cannot silently miss
cmap={genome_key(k.strip()): v.strip() for k, v in zip(cl.Accession, cl.Wang_Cluster)}
df["cluster"]=df.acc.map(cmap.get)
lab=df[(df.set=="discovery") & df.cluster.notna() & ~df.cluster.isin(["Not assigned","Outgroup"])].copy()
lab["isL"]=lab.cluster.eq("L")
def wil(k,n,z=1.959963985):
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n)); h=(z/d)*sqrt(p*(1-p)/n+z*z/(4*n*n)); return (100*max(0,c/d-h),100*min(1,c/d+h))
print(f"\n=== discovery, labeled genomes with a CDS: {len(lab)} (Cluster L {int(lab.isL.sum())}) ===")
for name,carr in (("CDS haplotype "+MH, lab.hap.eq(MH)), ("cgMLST allele 4", lab.ftsZ.eq("4"))):
    TP=int((carr&lab.isL).sum()); FP=int((carr&~lab.isL).sum()); FN=int((~carr&lab.isL).sum()); TN=int((~carr&~lab.isL).sum())
    print(f"  {name:22s} TP {TP:3d} FP {FP:3d} FN {FN:2d} | sens {100*TP/(TP+FN):5.1f}%  spec {100*TN/(TN+FP):6.2f}%  PPV {100*TP/(TP+FP):5.1f}% ({'%.1f-%.1f'%wil(TP,TP+FP)})")

haps={t:df[df.hap==t].cds.iloc[0] for t in df.hap.unique()}
L=len(next(iter(haps.values())))
var=[i for i in range(L) if len({haps[t][i] for t in haps})>1]
uniq=[i for i in var if all(haps[MH][i]!=haps[t][i] for t in haps if t!=MH)]
print(f"\n=== nucleotide structure (CDS {L} bp, {len(haps)} haplotypes) ===")
print("  variable positions:", len(var))
print("  positions where the marker haplotype differs from ALL others:", len(uniq), [p+1 for p in uniq])
near=sorted(((t,sum(1 for i in range(L) if haps[MH][i]!=haps[t][i])) for t in haps if t!=MH), key=lambda z:z[1])
print("  nearest haplotypes (differences from marker):", near[:5])
df.drop(columns=["cds"]).to_csv("results/tables/ftsz_cds_final.csv",index=False)
import pickle; pickle.dump({"df":df,"haps":haps,"MH":MH},open("results/tables/ftsz_final.pkl","wb"))

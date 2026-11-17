import sys, re, collections, pandas as pd
from Bio.Seq import Seq
FRAG=("GAGGTCTTGGTGCAGGTGCAAATCCTGAAATTGGTAAAAAAGCTGCTGAAGAATCTCGTGAACAAATCGAAGATGCTATC"
      "CAAGGTGCAGATATGGTATTCGTTACTGCTGGTATGGGTGGCGGTACTGGTACAGGAGCTGCTCCAGTAGTTGCTAAAATC"
      "GCTAAAGAAATGGGCGCTTTAACAGTAGGTGTTGTAACTCGCCCATTCGGTTTTGAGGGACGTAAACGTCAAACTCAAGCT"
      "GCCGCTGGCGTTGAAGCAATGAAAGCTGCAGTTGATACATTAATTGTTATTCCTAATGATCGTTTATTAGATATCGTAGAC"
      "AAATCTACGCCTATGATGGAAGCATTCAAGGAAGCAGATAATGTGTTACGTCAAGGTGTACAAGGTATCTCAGATTTAATC"
      "GCTGTATCAGGTGAGGTAAACTTAGACTTTGCCGACGTTAAAACTATTATGTCTAACCAAGGTTCAGCATTAATGGGTATC"
      "GGTGTATCTTCTGGTGAAAATAGAGCAGTTGAA")
seqs={}; acc=None; buf=[]
for line in open(sys.argv[1]):
    line=line.rstrip()
    if line.startswith(">"):
        if acc: seqs[acc]="".join(buf)
        acc=line[1:].split()[0]; buf=[]
    else: buf.append(line)
if acc: seqs[acc]="".join(buf)
print("sequences:", len(seqs), "| lengths:", collections.Counter(len(s) for s in seqs.values()).most_common())
print("distinct CDS:", len(set(seqs.values())))

a=pd.read_csv("results/tables/ftsz_allele_assignments.csv", dtype=str)
a["base"]=a.id.str.replace(r"\.\d+$","",regex=True)
ex=dict(zip(a.id, zip(a.ftsZ,a.set))); ba=dict(zip(a.base, zip(a.ftsZ,a.set)))
rows=[]
for k,s in seqs.items():
    h=ex.get(k) or ba.get(re.sub(r"\.\d+$","",k))
    rows.append((k,s,h[0] if h else None,h[1] if h else None))
df=pd.DataFrame(rows,columns=["acc","cds","ftsZ","set"])
print("matched to a cgMLST call:", df.ftsZ.notna().sum(), "| by set:", df[df.ftsZ.notna()].set.value_counts().to_dict())

print("\n=== DECISIVE CHECK: is the 518-bp marker fragment inside the CDS? ===")
m4=df[df.ftsZ=="4"]
print("  genomes with cgMLST allele 4:", len(m4))
hits=sum(FRAG in s for s in m4.cds)
print("  of those, CDS containing the 518-bp marker fragment exactly:", hits, "/", len(m4))
if len(m4):
    s0=m4.cds.iloc[0]; i=s0.find(FRAG)
    print("  offset within the 1,188-bp CDS:", i, "->", (i+1, i+len(FRAG)) if i>=0 else "NOT FOUND")
allhit=sum(FRAG in s for s in df.cds)
print("  across ALL genomes, CDS containing the marker fragment:", allhit, "/", len(df))

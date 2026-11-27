import sys, re, collections, pandas as pd
from Bio.Seq import Seq
FRAG=("GAGGTCTTGGTGCAGGTGCAAATCCTGAAATTGGTAAAAAAGCTGCTGAAGAATCTCGTGAACAAATCGAAGATGCTATC"
      "CAAGGTGCAGATATGGTATTCGTTACTGCTGGTATGGGTGGCGGTACTGGTACAGGAGCTGCTCCAGTAGTTGCTAAAATC"
      "GCTAAAGAAATGGGCGCTTTAACAGTAGGTGTTGTAACTCGCCCATTCGGTTTTGAGGGACGTAAACGTCAAACTCAAGCT"
      "GCCGCTGGCGTTGAAGCAATGAAAGCTGCAGTTGATACATTAATTGTTATTCCTAATGATCGTTTATTAGATATCGTAGAC"
      "AAATCTACGCCTATGATGGAAGCATTCAAGGAAGCAGATAATGTGTTACGTCAAGGTGTACAAGGTATCTCAGATTTAATC"
      "GCTGTATCAGGTGAGGTAAACTTAGACTTTGCCGACGTTAAAACTATTATGTCTAACCAAGGTTCAGCATTAATGGGTATC"
      "GGTGTATCTTCTGGTGAAAATAGAGCAGTTGAA")
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

raw={}; acc=None; buf=[]
for line in open(sys.argv[1]):
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
print("sequences:", len(seqs), "| lengths:", collections.Counter(len(s) for s in seqs.values()).most_common())
print("distinct CDS:", len(set(seqs.values())))

a=pd.read_csv("results/tables/ftsz_allele_assignments.csv", dtype=str)
a["key"]=a.id.map(genome_key)
bykey=dict(zip(a.key, zip(a.ftsZ,a.set)))
rows=[]
for k,s in seqs.items():
    h=bykey.get(k)
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

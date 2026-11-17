import pickle, itertools, pandas as pd
from Bio.Seq import Seq
d=pickle.load(open("results/tables/ftsz_final.pkl","rb")); df=d["df"]; haps=d["haps"]; MH=d["MH"]
L=len(haps[MH]); others=[t for t in haps if t!=MH]
var=[i for i in range(L) if len({haps[t][i] for t in haps})>1]
print(f"CDS {L} bp | {len(haps)} haplotypes | {len(var)} variable positions")

def discriminates(ps):
    sig=tuple(haps[MH][p] for p in ps)
    return all(tuple(haps[t][p] for p in ps)!=sig for t in others)

found=None
for k in (1,2,3,4):
    sets=[c for c in itertools.combinations(var,k) if discriminates(c)]
    print(f"  size {k}: {len(sets)} discriminating set(s)")
    if sets and found is None:
        found=(k,sets); break
k,sets=found
print(f"\nminimum discriminating set size: {k}; {len(sets)} such sets")
from collections import Counter
cc=Counter(p for s in sets for p in s)
print("positions appearing in every minimal set:", [p+1 for p,n in cc.items() if n==len(sets)])
print("most frequent positions:", [(p+1,n) for p,n in cc.most_common(8)])

best=sets[0]
print(f"\nexample minimal set (1-based CDS positions): {[p+1 for p in best]}")
ref=haps[MH]
prot=str(Seq(ref).translate())
for p in sorted({q for s in sets[:40] for q in s} | set(best)):
    codon_i=p//3; pos_in=p%3
    codon=ref[codon_i*3:codon_i*3+3]
    alts={haps[t][p] for t in others if haps[t][p]!=ref[p]}
    effs=set()
    for alt in alts:
        c2=list(codon); c2[pos_in]=alt; c2="".join(c2)
        effs.add(f"{codon}->{c2} {Seq(codon).translate()}->{Seq(c2).translate()}")
    if cc[p] >= len(sets)*0.5:
        print(f"  CDS pos {p+1:4d} (codon {codon_i+1}, base {pos_in+1}): {ref[p]} vs {''.join(sorted(alts))} | {'; '.join(sorted(effs))}")

print("\n=== relation to the 518-bp fragment (CDS 200-717) ===")
print("  fragment position 191 maps to CDS position", 199+191)
print("  is CDS 390 variable?", (390-1) in var)
if (390-1) in var:
    p=390-1
    print("   marker base:", ref[p], "| other haplotypes:", sorted({haps[t][p] for t in others}))
    print("   in how many minimal sets:", cc[p], "of", len(sets))
nearest=sorted(((t,sum(1 for i in range(L) if haps[MH][i]!=haps[t][i])) for t in others), key=lambda z:z[1])[:3]
print("\nnearest haplotypes:", nearest)
for t,nd in nearest:
    dp=[i+1 for i in range(L) if haps[MH][i]!=haps[t][i]]
    n=int((df.hap==t).sum()); al=df[df.hap==t].ftsZ.unique()
    print(f"  {t}: {nd} diff at {dp} | {n} genomes | cgMLST allele {al}")

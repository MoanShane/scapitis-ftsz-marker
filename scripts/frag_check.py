import pickle
FRAG=("GAGGTCTTGGTGCAGGTGCAAATCCTGAAATTGGTAAAAAAGCTGCTGAAGAATCTCGTGAACAAATCGAAGATGCTATC"
      "CAAGGTGCAGATATGGTATTCGTTACTGCTGGTATGGGTGGCGGTACTGGTACAGGAGCTGCTCCAGTAGTTGCTAAAATC"
      "GCTAAAGAAATGGGCGCTTTAACAGTAGGTGTTGTAACTCGCCCATTCGGTTTTGAGGGACGTAAACGTCAAACTCAAGCT"
      "GCCGCTGGCGTTGAAGCAATGAAAGCTGCAGTTGATACATTAATTGTTATTCCTAATGATCGTTTATTAGATATCGTAGAC"
      "AAATCTACGCCTATGATGGAAGCATTCAAGGAAGCAGATAATGTGTTACGTCAAGGTGTACAAGGTATCTCAGATTTAATC"
      "GCTGTATCAGGTGAGGTAAACTTAGACTTTGCCGACGTTAAAACTATTATGTCTAACCAAGGTTCAGCATTAATGGGTATC"
      "GGTGTATCTTCTGGTGAAAATAGAGCAGTTGAA")
print("fragment length:", len(FRAG))
d=pickle.load(open("results/tables/cds.pkl","rb")); ok=d["ok"]; inv={v:k for k,v in d["cdsid"].items()}
c04=inv["C04"]
print("C04 length:", len(c04))
i=c04.find(FRAG)
print("fragment found inside C04 at 0-based offset:", i)
if i<0:
    from Bio.Seq import Seq
    j=str(Seq(c04).reverse_complement()).find(FRAG)
    print("  on reverse complement:", j)
    # partial overlap probe
    for k in (60, 120, 240):
        print(f"  first {k} bp of fragment in C04? ", c04.find(FRAG[:k]))
        print(f"  last  {k} bp of fragment in C04? ", c04.find(FRAG[-k:]))
else:
    print("  => fragment spans CDS positions", i+1, "to", i+len(FRAG))
    # check it in every CDS type
    present={t:(inv[t].find(FRAG)>=0) for t in ok.cds_type.unique()}
    n=sum(present.values())
    print("  CDS types containing this exact fragment:", n, "of", len(present))
    print("  ->", [t for t,v in present.items() if v][:12])
    # how many distinct fragment-sequences across CDS types at that span
    frs={t: inv[t][i:i+len(FRAG)] for t in ok.cds_type.unique()}
    print("  distinct sequences at that span across the 36 CDS types:", len(set(frs.values())))
    ok2=ok.copy(); ok2["frag"]=ok2.cds_type.map(frs)
    print("  distinct fragment sequences across the 603 genomes:", ok2.frag.nunique())

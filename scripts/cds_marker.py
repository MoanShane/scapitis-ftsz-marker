import pickle, pandas as pd, numpy as np
from math import sqrt
d=pickle.load(open("results/tables/cds.pkl","rb")); ok=d["ok"]; inv={v:k for k,v in d["cdsid"].items()}
cl=pd.read_csv(sys.argv[1] if False else "data/wang2022_cgmlst_clusters_602.csv", dtype=str)
cmap=dict(zip(cl.Accession.str.strip(), cl.Wang_Cluster.str.strip()))
ok=ok.copy()
ok["cluster"]=ok.acc.map(lambda a: cmap.get(a) or cmap.get(a.split(".")[0]))
lab=ok[ok.cluster.notna() & ~ok.cluster.isin(["Not assigned","Outgroup"])].copy()
lab["isL"]=lab.cluster.eq("L")
print("labeled genomes with a CDS:", len(lab), "| Cluster L:", int(lab.isL.sum()))

def wilson(k,n,z=1.959963985):
    if n==0: return (float('nan'),)*2
    p=k/n; den=1+z*z/n; c=(p+z*z/(2*n))/den
    h=(z/den)*sqrt(p*(1-p)/n+z*z/(4*n*n)); return (100*max(0,c-h),100*min(1,c+h))

print("\n=== are the Cluster L genomes homogeneous at CDS level? ===")
print(lab[lab.isL].cds_type.value_counts().to_string())
print("\n=== CDS type C04 as a marker, discovery set ===")
carr=lab.cds_type.eq("C04")
TP=int((carr&lab.isL).sum()); FP=int((carr&~lab.isL).sum())
FN=int((~carr&lab.isL).sum()); TN=int((~carr&~lab.isL).sum())
print(f"  TP {TP}  FP {FP}  FN {FN}  TN {TN}")
if TP+FN: print("  sensitivity %.1f%% (%s)" % (100*TP/(TP+FN), "%.1f-%.1f"%wilson(TP,TP+FN)))
if TN+FP: print("  specificity %.2f%%" % (100*TN/(TN+FP)))
if TP+FP: print("  PPV         %.1f%% (%s)" % (100*TP/(TP+FP), "%.1f-%.1f"%wilson(TP,TP+FP)))
print("\n  the false positives and their cgMLST allele / cluster:")
print(lab[carr&~lab.isL][["acc","ftsZ","cluster"]].to_string(index=False))
print("\n=== for comparison, cgMLST ftsZ allele 4 on the SAME genome subset ===")
c2=lab.ftsZ.eq("4")
TP2=int((c2&lab.isL).sum()); FP2=int((c2&~lab.isL).sum()); FN2=int((~c2&lab.isL).sum())
print(f"  TP {TP2}  FP {FP2}  FN {FN2} | sens %.1f%%  PPV %.1f%%" % (
    100*TP2/(TP2+FN2) if TP2+FN2 else 0, 100*TP2/(TP2+FP2) if TP2+FP2 else 0))

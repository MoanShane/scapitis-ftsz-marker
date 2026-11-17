import pandas as pd, numpy as np
r = pd.read_csv("results/tables/cgmlst_target_external.csv")
nonL_d, nL_d = 583, 19
nonL_e, nL_e = 279, 22
r["disc_spec"] = 100*(nonL_d - r.disc_FP)/nonL_d
r["disc_PPV"]  = 100*nL_d/(nL_d + r.disc_FP)
r["ext_spec"]  = 100*(nonL_e - r.ext_FP)/nonL_e
r["ext_PPV"]   = 100*r.ext_TP/(r.ext_TP + r.ext_FP)
r["ext_sens"]  = 100*r.ext_TP/(r.ext_TP + r.ext_FN)

print("candidate targets (100% discovery sensitivity):", len(r))
print("  with 0 discovery FP :", int((r.disc_FP==0).sum()))
print("  with <=11 disc FP   :", int((r.disc_FP<=11).sum()))
print()
print("EXTERNAL behaviour of ALL", len(r), "candidates:")
print("  keep 100% ext sensitivity :", int((r.ext_FN==0).sum()), f"({100*(r.ext_FN==0).mean():.1f}%)")
print("  keep 0 ext false positives:", int((r.ext_FP==0).sum()))
print("  ext FP minimum / median / max:", int(r.ext_FP.min()), int(r.ext_FP.median()), int(r.ext_FP.max()))
print()
p = r[r.disc_FP==0]
print("THE 26 DISCOVERY-PERFECT TARGETS externally:")
print("  ext sensitivity 100%:", int((p.ext_FN==0).sum()), " <100%:", int((p.ext_FN>0).sum()))
print("  ext sensitivity of the 5 that drop:", sorted(p[p.ext_FN>0].ext_sens.round(1).tolist()))
print("  ext FP: min %d, median %d, max %d" % (p.ext_FP.min(), p.ext_FP.median(), p.ext_FP.max()))
print("  ext PPV: min %.1f, median %.1f, max %.1f" % (p.ext_PPV.min(), p.ext_PPV.median(), p.ext_PPV.max()))
print("  how many keep BOTH 100%% sens and 0 FP externally: %d" % int(((p.ext_FN==0)&(p.ext_FP==0)).sum()))
print()
f = r[r.target.str.lower()=="ftsz"].iloc[0]
print("ftsZ: disc FP %d (spec %.1f, PPV %.1f) | ext TP %d FN %d FP %d (sens %.1f, spec %.1f, PPV %.1f)" % (
    f.disc_FP, f.disc_spec, f.disc_PPV, f.ext_TP, f.ext_FN, f.ext_FP, f.ext_sens, f.ext_spec, f.ext_PPV))
print("ftsZ discovery rank by FP:", int((r.disc_FP < f.disc_FP).sum())+1, "of", len(r))
# best overall external
best = r[(r.ext_FN==0)].sort_values(["ext_FP","disc_FP"])
print("\nbest external performers (0 ext FN), top 8:")
print(best.head(8)[["target","allele","disc_FP","ext_FP","ext_PPV"]].round(1).to_string(index=False))
print("\nn candidates with 0 ext FN and exactly 1 ext FP:", int(((r.ext_FN==0)&(r.ext_FP==1)).sum()))
r.to_csv("results/tables/cgmlst_target_full.csv", index=False)

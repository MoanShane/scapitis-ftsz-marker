import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 8, "axes.labelsize": 8.5, "axes.titlesize": 9,
    "xtick.labelsize": 7.5, "ytick.labelsize": 8,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

df = pd.read_csv("results/tables/dist_both.csv", dtype={"ST": str})
LEFT  = (-2, 62)
RIGHT = (385, 1460)
C_L, C_A, C_O = "#0072B2", "#D55E00", "#9A9A9A"

def group(st):
    if st == "6":  return "lineage"
    if st == "46": return "st46"
    return "other"

rng = np.random.default_rng(11)
fig, axes = plt.subplots(2, 2, figsize=(7.0, 3.9), sharey="row",
                         gridspec_kw={"width_ratios": [2.5, 1.7], "wspace": 0.045,
                                      "hspace": 0.40})

panels = [("Discovery collection (n = 620)", "Discovery", 25, 398, 373),
          ("External collection (n = 301)",  "External",  50, 514, 464)]

for r, (title, coll, brk, nextg, gap) in enumerate(panels):
    sub = df[df.collection == coll].copy()
    sub["grp"] = sub.ST.map(group)
    axl, axr = axes[r, 0], axes[r, 1]
    for ax in (axl, axr):
        for g, col, mk, ms, fc in (("other", C_O, "o", 3.4, "none"),
                                   ("st46",  C_A, "^", 5.2, C_A),
                                   ("lineage", C_L, "o", 4.2, C_L)):
            s = sub[sub.grp == g]
            if not len(s): continue
            y = rng.uniform(-0.34, 0.34, len(s))
            ax.plot(s.dist, y, mk, ms=ms, mfc=fc, mec=col, mew=0.8,
                    ls="none", alpha=0.85, zorder=3 if g != "other" else 2)
        ax.axvline(25, color="0.25", lw=0.7, ls=(0, (4, 2)), zorder=1)
        ax.set_ylim(-0.75, 0.75); ax.set_yticks([])
        ax.grid(axis="x", color="0.90", lw=0.4, zorder=0)
        ax.set_axisbelow(True)

    axl.set_xlim(*LEFT); axr.set_xlim(*RIGHT)
    axl.spines[["top", "right"]].set_visible(False)
    axr.spines[["top", "left"]].set_visible(False)
    axl.set_xticks([0, 10, 20, 30, 40, 50, 60])
    axr.set_xticks([400, 700, 1000, 1300])
    # axis break marks
    d = .012
    for ax, xs in ((axl, [1]), (axr, [0])):
        for xx in xs:
            ax.plot([xx, xx], [0, 1], transform=ax.transAxes, color="w", lw=2.2,
                    clip_on=False, zorder=10)
            ax.plot([xx - d*2.6, xx + d*2.6], [-0.015, 0.035], transform=ax.transAxes,
                    color="k", lw=0.7, clip_on=False, zorder=11)
    axl.set_title(title, loc="left", pad=5, fontweight="bold")

    # annotate lineage max and the gap
    lin = sub[sub.grp == "lineage"]
    over = int((lin.dist > 25).sum())
    axl.annotate(f"lineage reaches {lin.dist.max()}", xy=(lin.dist.max(), 0.42),
                 xytext=(lin.dist.max() + 2.0, 0.56), fontsize=7, color=C_L,
                 va="bottom", arrowprops=dict(arrowstyle="-", lw=0.6, color=C_L))
    if over:
        axl.text(0.012, 0.06,
                 f"{over} ST6 genomes exceed the 25-allele cutoff",
                 transform=axl.transAxes, ha="left", va="bottom",
                 fontsize=7, color=C_L)
    axr.annotate(f"nearest non-lineage\ngenome at {nextg}", xy=(nextg, 0.30),
                 xytext=(nextg + 70, 0.58), fontsize=7, color="0.25",
                 arrowprops=dict(arrowstyle="-", lw=0.6, color="0.45"))
    axl.text(0.985, 0.06, f"gap = {gap} alleles", transform=axl.transAxes,
             ha="right", va="bottom", fontsize=7, style="italic", color="0.35")

beyond = (df.dist > RIGHT[1]).sum()
axes[1, 0].set_xlabel("cgMLST allelic distance to LNZR-1", x=0.62)
handles = [Line2D([], [], ls="none", marker="o", ms=4.2, mfc=C_L, mec=C_L, label="Seven-locus ST6 (lineage)"),
           Line2D([], [], ls="none", marker="^", ms=5.2, mfc=C_A, mec=C_A, label="ST46 (boundary genome)"),
           Line2D([], [], ls="none", marker="o", ms=3.4, mfc="none", mec=C_O, mew=0.8, label="Other or unassigned ST"),
           Line2D([], [], color="0.25", lw=0.7, ls=(0, (4, 2)), label="25-allele cutoff from discovery set")]
fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
           fontsize=7, bbox_to_anchor=(0.5, -0.045), handletextpad=0.5, columnspacing=1.4)

fig.savefig("results/figures/Figure2_distance_boundary.png", dpi=400, bbox_inches="tight", facecolor="w")
fig.savefig("results/figures/Figure2_distance_boundary.pdf", bbox_inches="tight", facecolor="w")
print("beyond cap:", beyond)
for coll in ("Discovery", "External"):
    s = df[df.collection == coll]
    print(coll, "n=", len(s), "| in left window:", (s.dist <= LEFT[1]).sum(),
          "| in right window:", ((s.dist >= RIGHT[0]) & (s.dist <= RIGHT[1])).sum())

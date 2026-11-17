import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["DejaVu Sans","Arial"],
    "font.size":8.5,"pdf.fonttype":42,"ps.fonttype":42})

stages=[("All *S. capitis* assemblies in NCBI", 952, None),
        ("Released 2023 onward", 614, "338 released before 2023"),
        ("One record per genome,\nRefSeq copy retained", 307, "307 paired GenBank duplicates"),
        ("Passed independence and\nquality filters", 302, "5 failed quality or independence filters"),
        ("Final external collection", 301, "1 removed on review:\nthe discovery outgroup DSM 20326")]

fig, ax = plt.subplots(figsize=(6.6, 5.6))
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
MAXW, CX = 5.0, 3.3
ys = [8.9, 7.0, 5.1, 3.2, 1.3]
C_BOX, C_EDGE, C_LAST = "#EAF0F6", "#5B7794", "#D6E4F0"

for i,((lab,n,drop), y) in enumerate(zip(stages, ys)):
    w = MAXW * (0.52 + 0.48*(n/952))
    face = C_LAST if i==len(stages)-1 else C_BOX
    lw   = 1.5 if i==len(stages)-1 else 0.9
    ax.add_patch(FancyBboxPatch((CX-w/2, y-0.52), w, 1.04,
                 boxstyle="round,pad=0.02,rounding_size=0.08",
                 facecolor=face, edgecolor=C_EDGE, linewidth=lw, zorder=3))
    ax.text(CX, y+0.17, lab.replace("*",""), ha="center", va="center",
            fontsize=8.2, color="#13293D", zorder=4,
            fontstyle="italic" if i==0 else "normal")
    ax.text(CX, y-0.26, f"n = {n:,}", ha="center", va="center",
            fontsize=10.5, fontweight="bold", color="#13293D", zorder=4)
    if i < len(stages)-1:
        ax.add_patch(FancyArrowPatch((CX, y-0.56), (CX, ys[i+1]+0.56),
                     arrowstyle="-|>", mutation_scale=13, lw=1.0,
                     color="#5B7794", zorder=2))
    if drop:
        ymid = (y + ys[i-1])/2 if i>0 else y
    if i>0 and drop:
        ymid = (ys[i-1] + y)/2
        ax.add_patch(FancyArrowPatch((CX+0.12, ymid), (CX+MAXW/2+0.55, ymid),
                     arrowstyle="-|>", mutation_scale=10, lw=0.8,
                     color="#A33B3B", zorder=2))
        ax.text(CX+MAXW/2+0.70, ymid, drop, ha="left", va="center",
                fontsize=7.2, color="#A33B3B", zorder=4)

ax.text(CX, 9.75, "Selection of the external genome collection",
        ha="center", fontsize=9.5, fontweight="bold", color="#13293D")
ax.text(CX, 0.28, "Enumerated with ncbi-datasets-cli v18.33.1, August 2026",
        ha="center", fontsize=6.8, color="0.45", style="italic")
fig.savefig("results/figures/Figure1_selection_flow.png", dpi=400, bbox_inches="tight", facecolor="w")
fig.savefig("results/figures/Figure1_selection_flow.pdf", bbox_inches="tight", facecolor="w")
print("ok")

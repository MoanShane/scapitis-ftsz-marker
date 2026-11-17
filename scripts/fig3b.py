import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["DejaVu Sans","Arial"],
    "font.size":8,"axes.labelsize":8.5,"xtick.labelsize":7.5,"ytick.labelsize":7.5,
    "axes.linewidth":0.6,"pdf.fonttype":42,"ps.fonttype":42})

r = pd.read_csv("results/tables/cgmlst_target_full.csv")
C_ALL, C_P, C_F, C_NS = "#B9B9B9", "#0072B2", "#D55E00", "#D9D2C5"
rng = np.random.default_rng(7)
def jit(v, s=0.045):
    v = np.asarray(v, float)
    return v * (1 + rng.uniform(-s, s, len(v))) + rng.uniform(-s, s, len(v)) * (v == 0)

fig, ax = plt.subplots(figsize=(6.6, 5.0))

full = r[r.ext_FN == 0]      # kept complete external sensitivity
lost = r[r.ext_FN > 0]       # lost some
ax.scatter(jit(lost.disc_FP), jit(lost.ext_FP), s=7, facecolors="none",
           edgecolors=C_NS, linewidths=0.5, zorder=2, label="_")
ax.scatter(jit(full.disc_FP), jit(full.ext_FP), s=7, color=C_ALL,
           alpha=0.75, linewidths=0, zorder=3, label="_")

perf = r[r.disc_FP == 0]
ax.scatter(jit(perf.disc_FP, .03), jit(perf.ext_FP, .03), s=26, color=C_P,
           edgecolors="white", linewidths=0.5, zorder=6)
f = r[r.target.str.lower() == "ftsz"].iloc[0]
ax.scatter([f.disc_FP], [f.ext_FP], s=70, marker="^", color=C_F,
           edgecolors="white", linewidths=0.7, zorder=7)
g = r[r.ext_FP == 0].iloc[0]
ax.scatter([g.disc_FP], [g.ext_FP], s=42, marker="s", facecolors="none",
           edgecolors="#444444", linewidths=1.0, zorder=7)

ax.set_xscale("symlog", linthresh=1, linscale=0.45)
ax.set_yscale("symlog", linthresh=1, linscale=0.45)
ax.set_xlim(-0.25, 700); ax.set_ylim(-0.25, 400)
ax.set_xticks([0,1,10,100]); ax.set_yticks([0,1,10,100])
ax.set_xticklabels(["0","1","10","100"]); ax.set_yticklabels(["0","1","10","100"])
ax.set_xlabel("Non-lineage genomes flagged in the discovery collection (of 583)")
ax.set_ylabel("Non-lineage genomes flagged externally (of 279)")
ax.grid(color="0.92", lw=0.4, zorder=0); ax.set_axisbelow(True)
ax.spines[["top","right"]].set_visible(False)

# the empty corner
ax.add_patch(Rectangle((-0.25,-0.25), 0.63, 0.63, facecolor="#F6E9E9",
                       edgecolor="#C06C6C", lw=0.7, ls=(0,(3,2)), zorder=1))
ax.annotate("no candidate here: none combined full external\nsensitivity with zero external false positives",
            xy=(0.40, 0.18), xytext=(2.0, 0.44), fontsize=6.6, color="#9A4A4A",
            arrowprops=dict(arrowstyle="->", lw=0.6, color="#9A4A4A"), va="center")
ax.annotate("26 candidates perfectly specific in discovery;\nall flag the same single external genome",
            xy=(0.12, 1.0), xytext=(0.52, 3.2), fontsize=6.8, color=C_P,
            arrowprops=dict(arrowstyle="->", lw=0.7, color=C_P))
ax.annotate("ftsZ", xy=(f.disc_FP, f.ext_FP), xytext=(26, 1.15), fontsize=8,
            color=C_F, style="italic", fontweight="bold",
            arrowprops=dict(arrowstyle="->", lw=0.7, color=C_F))
ax.annotate("$\\it{gyrA}$: no external false positive,\nbut misses a lineage genome",
            xy=(1.25, 0.0), xytext=(26.0, 0.030), fontsize=6.4, color="#444444",
            arrowprops=dict(arrowstyle="->", lw=0.6, color="#666666"), va="center")
ax.text(0.035, 0.80, "Spearman $\\rho$ = 0.94", transform=ax.transAxes, ha="left",
        fontsize=8.5, color="0.2")

handles = [Line2D([],[],ls="none",marker="o",ms=4.6,color=C_P,label="Perfectly specific in discovery (n = 26)"),
           Line2D([],[],ls="none",marker="^",ms=6.0,color=C_F,label="$\\it{ftsZ}$ — the Sanger-readable candidate"),
           Line2D([],[],ls="none",marker="o",ms=3.2,color=C_ALL,label="Other candidates, full external sensitivity"),
           Line2D([],[],ls="none",marker="o",ms=3.2,mfc="none",mec=C_NS,mew=0.7,label="Candidates that lost external sensitivity")]
ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=6.8,
          handletextpad=0.5, borderpad=0.2, labelspacing=0.38)
fig.suptitle("Discovery performance predicts external performance, but no candidate is error-free",
             fontsize=9, fontweight="bold", x=0.02, ha="left", y=0.985)
fig.savefig("results/figures/Figure3_candidate_transfer.png", dpi=400, bbox_inches="tight", facecolor="w")
fig.savefig("results/figures/Figure3_candidate_transfer.pdf", bbox_inches="tight", facecolor="w")
print("candidates:", len(r), "| full ext sens:", len(full), "| lost:", len(lost))
print("ftsZ at", int(f.disc_FP), int(f.ext_FP), "| gyrA at", int(g.disc_FP), int(g.ext_FP))

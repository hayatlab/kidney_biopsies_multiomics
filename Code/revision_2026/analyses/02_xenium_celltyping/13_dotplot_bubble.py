"""(1) scanpy dotplot of the top-3 marker genes per level2_refined, and
(2) a bubble plot of cell counts (x=level2_refined, y=condition; bubble area ~ n cells)
with row (per-condition) and column (per-level2) totals annotated outside the axes."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import numpy as np
import pandas as pd
import scanpy as sc
import config as C

ORDER = ["PT", "aPT", "dPT", "TAL", "DCT", "CNT", "PC", "IC-A", "IC-B",
         "POD", "PEC", "EC", "FIB", "MES", "VSM-P", "Immune"]
COND_ORDER = ["CTRL", "DN", "FSGS", "IgAN", "MN", "TIN"]


def palette(cats):
    present = [g for g in ORDER if g in set(cats)] + [g for g in cats if g not in ORDER]
    base = list(plt.cm.tab20.colors)
    return {g: mcolors.to_hex(base[i % 20]) for i, g in enumerate(present)}


def main():
    a = sc.read_h5ad(C.XENIUM_LABELED)
    a.obs["level2_refined"] = a.obs["level2_refined"].astype("category")
    cats_present = [c for c in ORDER if c in set(a.obs["level2_refined"].cat.categories)]
    a.obs["level2_refined"] = a.obs["level2_refined"].cat.reorder_categories(cats_present)
    cmap = palette(cats_present)

    # ---------- (1) dotplot: top-3 markers per level2_refined ----------
    a.X = a.layers["counts"].copy().astype(np.float32)
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    sc.tl.rank_genes_groups(a, "level2_refined", method="wilcoxon")
    dp = sc.pl.rank_genes_groups_dotplot(
        a, n_genes=3, groupby="level2_refined", standard_scale="var",
        values_to_plot="logfoldchanges", cmap="Reds", vmin=0, return_fig=True)
    dp.savefig(C.FIG / "dotplot_top3_markers_level2_refined.pdf", bbox_inches="tight", dpi=200)
    dp.savefig(C.FIG / "dotplot_top3_markers_level2_refined.png", bbox_inches="tight", dpi=150)
    plt.close("all")

    # ---------- (2) bubble plot: counts by level2_refined x condition ----------
    conds = [c for c in COND_ORDER if c in set(a.obs["condition"].astype(str))]
    ct = pd.crosstab(a.obs["condition"], a.obs["level2_refined"]).reindex(
        index=conds, columns=cats_present, fill_value=0)
    col_tot, row_tot = ct.sum(0), ct.sum(1)
    nX, nY = len(cats_present), len(conds)
    kmax = ct.values.max()
    k = 1500.0 / kmax                               # bubble area ~ n cells
    smin = 12                                       # floor so small counts stay visible

    fig = plt.figure(figsize=(0.72 * nX + 3.5, 0.72 * nY + 3.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[nX, 3.0], height_ratios=[2.8, nY],
                          wspace=0.05, hspace=0.06)
    ax = fig.add_subplot(gs[1, 0])
    axt = fig.add_subplot(gs[0, 0], sharex=ax)      # top: per-type totals (bars)
    axr = fig.add_subplot(gs[1, 1], sharey=ax)      # right: per-condition totals (bars)
    axl = fig.add_subplot(gs[0, 1]); axl.set_axis_off()   # size legend

    for j, t in enumerate(cats_present):
        for i, cond in enumerate(conds):
            n = ct.loc[cond, t]
            if n > 0:
                ax.scatter(j, i, s=max(n * k, smin), color=cmap[t], alpha=0.85,
                           edgecolor="0.25", linewidth=0.4, zorder=3)
    ax.set_xticks(range(nX)); ax.set_xticklabels(cats_present, rotation=90, fontsize=8)
    ax.set_yticks(range(nY)); ax.set_yticklabels(conds, fontsize=9)
    ax.set_xlim(-0.6, nX - 0.4); ax.set_ylim(-0.6, nY - 0.4)
    ax.set_xlabel("level2_refined"); ax.set_ylabel("condition")
    ax.grid(True, color="0.9", linewidth=0.5); ax.set_axisbelow(True)

    axt.bar(range(nX), col_tot.values, width=0.62,
            color=[cmap[t] for t in cats_present], edgecolor="0.3", linewidth=0.3)
    for j, t in enumerate(cats_present):
        axt.text(j, col_tot[t] * 1.03, f"{int(col_tot[t]):,}", ha="center", va="bottom",
                 rotation=90, fontsize=7)
    axt.set_ylim(0, col_tot.max() * 1.55); axt.set_axis_off()
    axt.set_title("Xenium cell counts: level2_refined × condition   (bubble area ∝ n cells)",
                  fontsize=11, pad=8, loc="left")

    axr.barh(range(nY), row_tot.values, height=0.62, color="0.75", edgecolor="0.3", linewidth=0.3)
    for i, cond in enumerate(conds):
        axr.text(row_tot[cond] * 1.03, i, f"{int(row_tot[cond]):,}", va="center", ha="left", fontsize=8)
    axr.set_xlim(0, row_tot.max() * 1.4); axr.set_axis_off()

    ref = [r for r in [500, 3000, 10000] if r <= kmax] or [int(kmax)]
    axl.set_xlim(0, 1); axl.set_ylim(0, 1)
    for yy, n in zip(np.linspace(0.78, 0.18, len(ref)), ref):
        axl.scatter(0.22, yy, s=max(n * k, smin), color="0.6", edgecolor="0.25", linewidth=0.4)
        axl.text(0.45, yy, f"{n:,}", va="center", ha="left", fontsize=8)
    axl.text(0.12, 0.95, "n cells", fontsize=9, ha="left", va="top")

    fig.savefig(C.FIG / "bubble_level2_by_condition.pdf", bbox_inches="tight", dpi=200)
    fig.savefig(C.FIG / "bubble_level2_by_condition.png", bbox_inches="tight", dpi=150)
    plt.close(fig)
    print("wrote dotplot_top3_markers_level2_refined.* and bubble_level2_by_condition.*")
    print("column totals (per level2_refined):\n", col_tot.to_string())
    print("row totals (per condition):\n", row_tot.to_string())


if __name__ == "__main__":
    main()

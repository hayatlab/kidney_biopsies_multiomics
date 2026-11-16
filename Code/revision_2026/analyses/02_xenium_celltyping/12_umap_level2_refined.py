"""UMAP of the Xenium cells coloured by the final level2_refined annotation, using the
same hex palette as the exported label CSV/legend so colours are consistent."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import scanpy as sc
import config as C

plt.rcParams.update({"figure.dpi": 150, "savefig.bbox": "tight", "font.size": 8})
ORDER = ["PT", "aPT", "dPT", "TAL", "DCT", "CNT", "PC", "IC-A", "IC-B",
         "POD", "PEC", "EC", "FIB", "MES", "VSM-P", "Immune"]
# Display names (13 Nov 2026): the co-authors renamed dPT to aPT2 in all figures. The object keeps "dPT"; set
# DISPLAY = {"dPT": "aPT2"} (the Supp. Fig. 4b entry point does) to print the new name. Colours are assigned before
# the rename, so they do not change. Default {} = the 1 July 2026 figure.
DISPLAY = {}


def main():
    a = sc.read_h5ad(C.XENIUM_LABELED)
    assert "X_umap" in a.obsm
    a.obs["level2_refined"] = a.obs["level2_refined"].astype("category")
    cats = list(a.obs["level2_refined"].cat.categories)
    present = [g for g in ORDER if g in set(cats)] + [g for g in cats if g not in ORDER]
    base = list(plt.cm.tab20.colors)                       # same palette as 11_export_labels.py
    cmap = {g: mcolors.to_hex(base[i % 20]) for i, g in enumerate(present)}
    a.uns["level2_refined_colors"] = [cmap[c] for c in cats]   # colours in category order
    if DISPLAY:
        a.obs["level2_refined"] = a.obs["level2_refined"].cat.rename_categories(lambda c: DISPLAY.get(c, c))

    for loc, tag in [("on data", "ondata"), ("right margin", "legend")]:
        sc.pl.umap(a, color="level2_refined", show=False, size=3, frameon=False,
                   legend_loc=loc, legend_fontsize=6, title="Xenium — level2_refined")
        for ax in plt.gcf().axes:
            for coll in ax.collections:
                coll.set_rasterized(True)
        plt.savefig(C.FIG / f"umap_level2_refined_{tag}.pdf", dpi=200)
        plt.savefig(C.FIG / f"umap_level2_refined_{tag}.png", dpi=150)
        plt.close()
    print("wrote umap_level2_refined_{ondata,legend}.{pdf,png} |", len(cats), "types")


if __name__ == "__main__":
    main()

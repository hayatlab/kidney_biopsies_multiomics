"""Heatmap comparing the per-cell transferred level2 (pred_level2) to the final
subcluster-consensus level2_refined. Row-normalised: each pred_level2 row shows the
fraction of those cells that ended up in each level2_refined label. Writes
figures/heatmap_pred_level2_vs_refined.pdf and outputs/metrics/pred_level2_vs_refined.csv."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import metrics as M
import config as C

plt.rcParams.update({"figure.dpi": 150, "savefig.bbox": "tight", "font.size": 8})

# canonical kidney level2 order so the "unchanged" cells fall on the diagonal
ORDER = ["PT", "aPT", "dPT", "TAL", "DCT", "CNT", "PC", "IC-A", "IC-B",
         "POD", "PEC", "EC", "FIB", "MES", "VSM-P", "Immune"]


def main():
    C.ensure_dirs()
    a = sc.read_h5ad(C.XENIUM_LABELED)
    ct = M.crosstab_normalized(a.obs["pred_level2"], a.obs["level2_refined"], normalize="index")
    # order rows/cols by ORDER (keep only present labels)
    rows = [t for t in ORDER if t in ct.index]
    cols = [t for t in ORDER if t in ct.columns]
    ct = ct.loc[rows, cols]
    ct.to_csv(C.MET / "pred_level2_vs_refined.csv")

    plt.figure(figsize=(0.55 * ct.shape[1] + 3, 0.5 * ct.shape[0] + 2))
    sns.heatmap(ct, cmap="magma", vmin=0, vmax=1, annot=True, fmt=".2f",
                annot_kws={"size": 5}, linewidths=0.3, linecolor="0.85",
                cbar_kws={"label": "row fraction"})
    plt.xlabel("level2_refined (subcluster consensus)")
    plt.ylabel("pred_level2 (per-cell transfer)")
    plt.title("Per-cell transferred level2  →  refined level2 (row-normalised)")
    plt.savefig(C.FIG / "heatmap_pred_level2_vs_refined.pdf")
    plt.savefig(C.FIG / "heatmap_pred_level2_vs_refined.png", dpi=140)
    plt.close()
    print("wrote", C.FIG / "heatmap_pred_level2_vs_refined.pdf", "| shape", ct.shape)
    # quick note: where did the dropped dPT cells go?
    if "dPT" in ct.index:
        top = ct.loc["dPT"].sort_values(ascending=False).head(3)
        print("dPT (transferred) -> refined:", {k: round(v, 2) for k, v in top.items()})


if __name__ == "__main__":
    main()

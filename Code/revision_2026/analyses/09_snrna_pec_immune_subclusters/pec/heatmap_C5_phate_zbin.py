#!/usr/bin/env python
"""'Old-style' pseudotime heatmap for the chosen trajectory C_0to9_5cl / phate.

Reproduces the style of `save_heatmap_arrows_only` in
$KB_PROJECT/pathways_traj.py (the "z-scored mean expression per pseudotime bin"
figure), rendered TRANSPOSED to match the reference the colleague showed: pseudotime on X (0..1),
genes on Y ordered by peak, green dashed quintile gridlines, inferno/magma colormap clipped at +-2.

Method faithful to save_heatmap_arrows_only (NOT the MAGIC/plasma one):
  * equal-count binning into n_bins=120 (bin_ids = rank*120//n) -> every bin populated, so PHATE
    pseudotime saturation does NOT collapse bins (the qcut version gave only 91). This is the main
    reason this figure looks smoother than pseudotime_heatmap_magic_peakorder_HVGs.
  * per-bin RAW mean of lognorm X (no MAGIC imputation)
  * z-score per gene across bins (zero-sd guarded)
  * gene order = argmax over gaussian_filter1d(z, sigma=1.2); display gaussian_filter1d(sigma=0.8)
  * clip to +-cap_z (cap_z=2.0, matching the reference colorbar)
Same 1000 HVGs and the SAME delivered pseudotime (obs['pseudotime_via_phate'], job 328231) as the
MAGIC version, so the two are directly comparable. seed 42. No VIA object needed.
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.sparse import issparse
from scipy.ndimage import gaussian_filter1d
import scanpy as sc

DIR = f"{KB_ANALYSIS}/pec_resubcluster_20260717"
CFG = f"{DIR}/trajectories_leiden08/C_0to9_5cl"
OUTD = f"{CFG}/phate"
SRCH5 = f"{CFG}/subset_all_pseudotimes.h5ad"
PT_COL = "pseudotime_via_phate"

N_BINS = 120
CAP_Z = 2.0
SIG_ORDER = 1.2
SIG_DISPLAY = 0.8
N_HVG = 1000
CMAP = "magma"        # black -> purple -> orange -> yellow (matches the reference; 'inferno' near-identical)

# optional marker triangles at the pseudotime edges (up at high-pt end, down at low-pt end).
# set to [] to get the perfectly clean version shown in the reference.
MARK_UP = ["PSD3", "GLS", "NRG3", "CD44", "ITGB1"]
MARK_DN = ["SLC4A11", "CDH1", "PAH"]


def main(out=None):
    """out: output stem; default = next to the trajectory outputs, as originally run."""
    a = sc.read_h5ad(SRCH5)
    assert PT_COL in a.obs, f"{PT_COL} missing"
    # 1000 HVGs on lognorm X (seurat flavor), identical to the MAGIC version
    hvg = sc.pp.highly_variable_genes(a, n_top_genes=N_HVG, inplace=False)
    a = a[:, hvg["highly_variable"].values].copy()
    genes = np.array(a.var_names)

    # sort cells by the delivered phate pseudotime; equal-count bins
    pt = a.obs[PT_COL].to_numpy(dtype=float)
    order_cells = np.argsort(pt, kind="mergesort")
    n = pt.size
    bin_ids = (np.arange(n) * N_BINS // n).astype(int)   # 0..N_BINS-1 over rank
    X = a.X[order_cells, :]
    Xd = X.toarray() if issparse(X) else np.asarray(X)

    expr_bins = np.zeros((N_BINS, a.n_vars), dtype=np.float32)
    for b in range(N_BINS):
        expr_bins[b] = Xd[bin_ids == b].mean(0)

    with np.errstate(invalid="ignore", divide="ignore"):
        mu = expr_bins.mean(0, keepdims=True); sd = expr_bins.std(0, keepdims=True)
        sd[sd == 0] = 1.0
        z = (expr_bins - mu) / sd                         # bins x genes

    z_order = gaussian_filter1d(z, sigma=SIG_ORDER, axis=0, mode="nearest") if SIG_ORDER else z
    order = np.argsort(np.argmax(z_order, axis=0))
    z = z[:, order]; genes = genes[order]
    if SIG_DISPLAY:
        z = gaussian_filter1d(z, sigma=SIG_DISPLAY, axis=0, mode="nearest")
    z = np.clip(z, -CAP_Z, CAP_Z)

    # ---- transposed render: genes (Y) x pseudotime-bin (X) ----
    ng = z.shape[1]
    fig, ax = plt.subplots(figsize=(8.6, 6.2)); fig.patch.set_facecolor("white")
    im = ax.imshow(z.T, origin="lower", aspect="auto", interpolation="nearest",
                   cmap=CMAP, vmin=-CAP_Z, vmax=CAP_Z,
                   extent=[0.0, 1.0, 0, ng])            # X spans pseudotime 0..1 (cell-rank fraction)
    # green dashed quintile gridlines
    for xb in (0.2, 0.4, 0.6, 0.8):
        ax.axvline(xb, color="#2ca02c", lw=1.1, ls="--", alpha=0.85)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel("pseudotime")
    ax.set_yticks([])
    ax.set_ylabel("genes (ordered by peak pseudotime)")
    ax.set_title("z-scored mean expression\nper pseudotime bin", fontsize=13)
    cb = fig.colorbar(im, ax=ax, orientation="vertical", pad=0.02, fraction=0.045)
    cb.set_ticks([-2, -1, 0, 1, 2]); cb.ax.tick_params(labelsize=10)

    # optional marker triangles at the pseudotime edges
    pos = {g: i for i, g in enumerate(genes)}
    xu = [pos[g] for g in MARK_UP if g in pos]           # peak-order index -> y
    xd = [pos[g] for g in MARK_DN if g in pos]
    for yy, col, mk in [(xu, "#D55E00", "<"), (xd, "#0072B2", "<")]:
        pass  # markers drawn below with correct x placement (edge)
    # up markers: small right-pointing triangle at the right (high-pt) edge, at their gene row
    if xu:
        ax.scatter(np.full(len(xu), 0.992), xu, marker="<", s=42, c="#D55E00",
                   edgecolors="white", linewidths=0.5, zorder=6, clip_on=False)
    if xd:
        ax.scatter(np.full(len(xd), 0.008), xd, marker=">", s=42, c="#0072B2",
                   edgecolors="white", linewidths=0.5, zorder=6, clip_on=False)

    fig.tight_layout()
    out = f"{OUTD}/pseudotime_heatmap_zbinmean_HVGs" if out is None else out
    fig.savefig(f"{out}.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(f"{out}.png", dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    pd.Series(genes, name="gene").to_csv(f"{out}_geneorder.csv", index=False)
    print(f"wrote {out}.pdf/.png  ({ng} genes x {N_BINS} equal-count bins)")
    print("marker up present:", [g for g in MARK_UP if g in pos],
          "| down present:", [g for g in MARK_DN if g in pos])


if __name__ == "__main__":
    main()

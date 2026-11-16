#!/usr/bin/env python
"""Gene-trend-along-pseudotime plots for the niche-20 top-5 / bottom-5 markers, on the
C_0to9_5cl PHATE VIA trajectory (root c0 -> terminal c9).

Rebuilds the VIA object with the SAME parameters used by trajectory_leiden08.py for the phate
config (so the trend curves match the existing genetrend_UP/DOWN style), then calls VIA's own
do_impute (MAGIC, 3 steps) + get_gene_expression -- exactly as the original pipeline. The saved
trajectory object already carries the resulting per-cell pseudotime in obs['pseudotime_via_phate'];
we recompute VIA and verify the rebuilt pseudotime matches it (Spearman) before plotting.

Gene sets (user, niche-20 markers):
  UP   (top 5):    CD44, ITGB1, TMSB10, PDLIM7, MMP7
  DOWN (bottom 5): CFH, CDH1, SLC4A11, BMP7, NDRG1

New output files (existing genetrend_UP/DOWN, which used different genes, are left untouched):
  phate/genetrend_UP_top5.{pdf,png}
  phate/genetrend_DOWN_bottom5.{pdf,png}

Env: ViaEnv (pyVIA). Run via SLURM (login node hits a GLIBCXX clash).
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import os
import numpy as np, pandas as pd, scipy.sparse as sp
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
from scipy.stats import spearmanr
try:
    import pyVIA.core as via
except Exception:
    import via  # noqa

BASE = f"{KB_ANALYSIS}/pec_resubcluster_20260717"
H5   = f"{BASE}/trajectories_leiden08/C_0to9_5cl/subset_all_pseudotimes.h5ad"
OUT  = f"{BASE}/trajectories_leiden08/C_0to9_5cl/phate"
ROOT_GROUP, TERM_GROUP = "c0", "c9"

UP   = ["CD44", "ITGB1", "TMSB10", "PDLIM7", "MMP7"]      # niche-20 top 5 (up in niche 20)
DOWN = ["CFH", "CDH1", "SLC4A11", "BMP7", "NDRG1"]        # niche-20 bottom 5 (down in niche 20)


def dense(x):
    return x.toarray() if sp.issparse(x) else np.asarray(x)


def main():
    a = sc.read_h5ad(H5)
    print(f"loaded {a.shape}; X log-norm max={float(a.X.max()):.2f}", flush=True)

    # ---- rebuild VIA on the PHATE embedding (params copied verbatim from trajectory_leiden08.py) ----
    Vdata = a.obsm["X_phate"]
    true_label = list(a.obs["via_label"].astype(str))
    v = via.VIA(data=Vdata, true_label=true_label, edgepruning_clustering_resolution=0.15,
                cluster_graph_pruning=0.5, knn=15, root_user=[ROOT_GROUP], resolution_parameter=1.5,
                dataset='group', random_seed=4, memory=5, do_compute_embedding=False,
                user_defined_terminal_group=[TERM_GROUP])
    v.run_VIA()

    # ---- verify the rebuilt pseudotime matches the stored one (faithful reproduction) ----
    pt_new = np.asarray(v.single_cell_pt_markov, dtype=float)
    pt_new = pt_new / pt_new.max() if pt_new.max() > 0 else pt_new
    pt_old = a.obs["pseudotime_via_phate"].to_numpy(dtype=float)
    rho = spearmanr(pt_new, pt_old, nan_policy="omit").statistic
    n_term = len(np.asarray(getattr(v, "terminal_clusters", [])).ravel())
    print(f"VIA rebuilt: terminal_lineages={n_term}; "
          f"Spearman(rebuilt pt, stored pseudotime_via_phate) = {rho:+.3f}", flush=True)

    # ---- gene trends (VIA do_impute MAGIC + get_gene_expression), same as original pipeline ----
    def gene_trend(genes, stem, title):
        present = [g for g in genes if g in a.var_names]
        missing = [g for g in genes if g not in a.var_names]
        if missing:
            print(f"  [{stem}] not in var_names (skipped): {missing}", flush=True)
        df = pd.DataFrame(dense(a[:, present].X), columns=present, index=a.obs_names)
        dm = v.do_impute(df, magic_steps=3, gene_list=present)
        figg, _ = via.get_gene_expression(via_object=v, gene_exp=dm)
        figg.set_size_inches(14, 7); figg.patch.set_facecolor("white")
        for axx in figg.get_axes():
            axx.set_facecolor("white")
        figg.suptitle(title, y=1.02, fontsize=13)
        figg.tight_layout()
        for ext in ("png", "pdf"):
            figg.savefig(f"{OUT}/{stem}.{ext}", dpi=150 if ext == "png" else None,
                         bbox_inches="tight", facecolor="white")
        plt.close(figg)
        print(f"  wrote {OUT}/{stem}.pdf/.png  genes={present}", flush=True)

    gene_trend(UP, "genetrend_UP_top5",
               "Niche-20 top-5 markers along C_0to9 pseudotime (root c0 -> terminal c9)")
    gene_trend(DOWN, "genetrend_DOWN_bottom5",
               "Niche-20 bottom-5 markers along C_0to9 pseudotime (root c0 -> terminal c9)")
    print("done ->", OUT, flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Fig. 6h — PHATE of the five snRNA-seq PEC subclusters on the aPEC trajectory (Leiden 0.8 clusters 0, 2, 4, 9, 11),
coloured by subcluster and by the crescent (niche-20) ULM score.
Panel        : Main Figures.pdf, Figure 6, panel h
Reads        : $KB_ANALYSIS/pec_resubcluster_20260717/trajectories_leiden08/C_0to9_5cl/subset_all_pseudotimes.h5ad
                 (analyses/09_snrna_pec_immune_subclusters/pec/trajectory_leiden08.py --tags C_0to9_5cl, wrapper
                  run_trajectory_leiden08.sh, job 328231: X_phate, cl, pseudotime_via_phate, pseudotime_dpt,
                  crescent_lesion_score — the layout and pseudotimes are read, never recomputed)
               $KB_ANALYSIS/pec_resubcluster_20260717/new_pec_annotation/res_0.8/cluster_summary_res0.8.csv
                 (legend marker names; analyses/09_snrna_pec_immune_subclusters/pec/pec_leiden_sweep.py)
Writes       : $KB_FIG_OUT/Fig6/Fig6h_pec_phate_clusters_crescent.{pdf,png}
Plotting code: extracted verbatim from analyses/09_snrna_pec_immune_subclusters/pec/trajectory_leiden08.py::run_embedding
               (scatter_panel + figure "(1) the VIA base embedding"; nested in the VIA run, so it cannot be imported);
               constants (CLCOL, cl_label, LABEL, SCORE, KEEP5) imported from the module
Illustrator  : panels 1 (clusters) and 4 (crescent ULM score) of the four-panel figure kept; legend entries
               "c0  SLC8A1/NPHS1" ... renamed "cluster 0" ...; titles edited
Environment  : conda env ViaEnv (LD_LIBRARY_PATH=$CONDA_PREFIX/lib, set by build/verify_run.sh)
Verified     : 2026-11-12, PIXEL-IDENTICAL (PDF, text layer equal; PNG byte-identical) vs
               $KB_ANALYSIS/pec_resubcluster_20260717/trajectories_leiden08/C_0to9_5cl/phate/phate_cluster_pt_crescent.pdf
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import scanpy as sc  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from niches_snrna import ANALYSES_09, load_script  # noqa: E402
from kbpaths import fig_out  # noqa: E402

tl = load_script(ANALYSES_09 / "pec" / "trajectory_leiden08.py")   # module-level: imports, seed, ANN (no outputs)
CFG = next(c for c in tl.CONFIGS if c["tag"] == "C_0to9_5cl")
EMB = "phate"


def main():
    a = sc.read_h5ad(f"{tl.OUT}/{CFG['tag']}/subset_all_pseudotimes.h5ad")
    keep = CFG["keep"]
    root_group, term_group = f"c{CFG['root'][0]}", f"c{CFG['terminal']}"     # single-root config
    pt = a.obs[f"pseudotime_via_{EMB}"].to_numpy()
    dpt = a.obs["pseudotime_dpt"].to_numpy()
    LABEL, SCORE, CLCOL, cl_label = tl.LABEL, tl.SCORE, tl.CLCOL, tl.cl_label

    # ---------------- verbatim from trajectory_leiden08.py::run_embedding ----------------
    Vdata = a.obsm["X_phate"]
    emb = EMB
    E = Vdata[:, :2] if emb != "diffmap" else a.obsm["X_diffmap"][:, 1:3]
    ax_names = {"umap": ("UMAP1", "UMAP2"), "phate": ("PHATE1", "PHATE2"),
                "diffmap": ("DC1", "DC2")}[emb]
    cl = a.obs["cl"].astype(str).to_numpy()

    def scatter_panel(ax, XY, colour_by, cmap=None, title="", lab=("", "")):
        if colour_by is None:                       # categorical clusters
            for c in keep:
                m = cl == f"c{c}"
                ax.scatter(XY[m, 0], XY[m, 1], s=6, lw=0, color=CLCOL[f"c{c}"], label=cl_label(c))
            ax.legend(markerscale=2.5, fontsize=8, frameon=False)
            h = None
        else:
            h = ax.scatter(XY[:, 0], XY[:, 1], s=6, lw=0, c=colour_by, cmap=cmap)
        ax.set_title(title); ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel(lab[0]); ax.set_ylabel(lab[1])
        return h

    # (1) the VIA base embedding: clusters / VIA pt / DPT / crescent
    fig, ax = plt.subplots(1, 4, figsize=(27, 6.2)); fig.patch.set_facecolor("white")
    scatter_panel(ax[0], E, None, title=f"{LABEL} clusters", lab=ax_names)
    h1 = scatter_panel(ax[1], E, pt, "viridis", f"VIA pseudotime ({root_group}->{term_group})", ax_names)
    fig.colorbar(h1, ax=ax[1], shrink=.7)
    h2 = scatter_panel(ax[2], E, dpt, "viridis", "DPT pseudotime (deterministic ref)", ax_names)
    fig.colorbar(h2, ax=ax[2], shrink=.7)
    h3 = scatter_panel(ax[3], E, a.obs[SCORE], "magma", "crescent ULM score", ax_names)
    fig.colorbar(h3, ax=ax[3], shrink=.7)
    fig.suptitle(f"{CFG['tag']} — VIA on {emb.upper()}", fontsize=14)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(fig_out("Fig6", f"Fig6h_pec_phate_clusters_crescent.{ext}"), dpi=150 if ext == "png" else None,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("[Fig6h] wrote", fig_out("Fig6", "Fig6h_pec_phate_clusters_crescent.pdf"))


if __name__ == "__main__":
    main()

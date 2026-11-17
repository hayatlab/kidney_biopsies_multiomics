#!/usr/bin/env python
"""Supp. Fig. 18c-d — snRNA-seq PEC subclusters (Leiden 0.8, 19 clusters, 7,004 nuclei): UMAP (c) and top-5
Wilcoxon marker genes per subcluster (d).
Panel        : All Supp Figures.pdf, Supplementary Figure 18, panels c-d
Reads        : $KB_ANALYSIS/pec_resubcluster_20260717/new_pec_annotation/pec_leiden_sweep.h5ad
                 (analyses/09_snrna_pec_immune_subclusters/pec/pec_leiden_sweep.py, conda base python, log
                  pec_leiden_sweep.log; holds the subset UMAP, leiden_0.8, crescent_lesion_score, provenance, the
                  palettes and the Wilcoxon result uns['rgg_leiden_0.8'] — nothing is recomputed)
Writes       : $KB_FIG_OUT/S18/S18c_pec_subclusters_umap.{pdf,png}, $KB_FIG_OUT/S18/S18d_pec_subclusters_markers_dotplot.{pdf,png}
Plotting code: extracted verbatim from analyses/09_snrna_pec_immune_subclusters/pec/pec_leiden_sweep.py (module-level
               sweep loop, sections "UMAP" and "Wilcoxon markers" at r = 0.8, plus save()); the script runs the whole
               sweep at import, so its code cannot be imported
Illustrator  : c: first of the three UMAP panels (clusters) kept, title replaced ("snRNAseq PEC subclusters");
               d: suptitle replaced ("top 5 marker genes per snRNAseq subcluster"), axis title added
Environment  : conda base python ($CONDA_ROOT/bin/python, as in the original run)
Verified     : 2026-11-12, PIXEL-IDENTICAL (both PDFs, text layers equal; PNGs byte-identical) vs
               $KB_ANALYSIS/pec_resubcluster_20260717/new_pec_annotation/res_0.8/umap_res0.8.pdf and markers_dotplot_res0.8.pdf
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import scanpy as sc  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "config"))
from kbpaths import KB_ANALYSIS, fig_out  # noqa: E402

sc.settings.verbosity = 1
sc.settings.autoshow = False
ROOT = f"{KB_ANALYSIS}/pec_resubcluster_20260717/new_pec_annotation"
r = 0.8


# verbatim from pec_leiden_sweep.py
def save(fig, stem):
    fig.patch.set_facecolor("white")
    for ax in fig.get_axes():
        ax.set_facecolor("white")
    fig.savefig(stem + ".pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(stem + ".png", dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    adata = sc.read_h5ad(f"{ROOT}/pec_leiden_sweep.h5ad")
    key = f"leiden_{r}"
    cats = list(adata.obs[key].cat.categories)
    PROV = "provenance" if "provenance" in adata.obs.columns else None

    # ---- UMAP (verbatim)
    colors = [key, "crescent_lesion_score"] + ([PROV] if PROV else [])
    fig, axs = plt.subplots(1, len(colors), figsize=(6.2 * len(colors), 5.4))
    axs = np.atleast_1d(axs)
    for ax, c in zip(axs, colors):
        sc.pl.umap(adata, color=c, ax=ax, show=False, frameon=False,
                   legend_loc="on data" if c == key else "right margin",
                   legend_fontsize=9, legend_fontoutline=2,
                   cmap="viridis" if c == "crescent_lesion_score" else None,
                   title={key: f"leiden {r} ({len(cats)} clusters)",
                          "crescent_lesion_score": "niche ULM score"}.get(c, c))
    fig.suptitle(f"PEC (level2==PEC, n={adata.n_obs}) - subset HVG/PCA/harmony - leiden {r}",
                 fontsize=13)
    fig.tight_layout()
    save(fig, fig_out("S18", "S18c_pec_subclusters_umap"))

    # ---- Wilcoxon markers (saved result uns['rgg_leiden_0.8'] instead of re-running rank_genes_groups)
    mk = sc.get.rank_genes_groups_df(adata, group=None, key=f"rgg_{key}")
    top5 = {c: mk[mk.group == c].sort_values("scores", ascending=False).names.head(5).tolist()
            for c in cats}
    seen, dot_genes = set(), []
    for c in cats:
        for g in top5[c]:
            if g not in seen:
                seen.add(g); dot_genes.append(g)
    dp = sc.pl.dotplot(adata, var_names=top5, groupby=key, standard_scale="var",
                       return_fig=True,
                       figsize=(max(10, 0.30 * len(dot_genes)), 0.55 * len(cats) + 2.2))
    dp.add_totals().style(cmap="Reds")
    dp.make_figure()                       # dp.fig is None until this is called
    fig = dp.fig
    fig.suptitle(f"top-5 Wilcoxon markers per cluster - leiden {r}", fontsize=12, y=1.02)
    save(fig, fig_out("S18", "S18d_pec_subclusters_markers_dotplot"))
    print("[S18c-d] wrote", fig_out("S18", "S18c_pec_subclusters_umap.pdf"), "and the markers dotplot")


if __name__ == "__main__":
    main()

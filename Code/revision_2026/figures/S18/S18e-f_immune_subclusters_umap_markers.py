#!/usr/bin/env python
"""Supp. Fig. 18e-f — snRNA-seq immune subclusters (Leiden 0.6, 18 clusters): UMAP (e) and top-5 Wilcoxon marker
genes per subcluster (f).
Panel        : All Supp Figures.pdf, Supplementary Figure 18, panels e-f
Reads        : $KB_ANALYSIS/subcluster_level2_20260728/Immune/Immune_leiden_sweep.h5ad
                 (analyses/09_snrna_pec_immune_subclusters/subcluster_sweep/subcluster_sweep.py --celltype Immune,
                  wrapper run_sweep_array.sh task 6, job 328344_6; holds the subset UMAP, leiden_0.6, its palette and the
                  Wilcoxon result uns['rgg_leiden_0.6'] — nothing is recomputed)
Writes       : $KB_FIG_OUT/S18/S18e_immune_subclusters_umap.{pdf,png}, $KB_FIG_OUT/S18/S18f_immune_subclusters_markers_dotplot.{pdf,png}
Plotting code: extracted verbatim from analyses/09_snrna_pec_immune_subclusters/subcluster_sweep/subcluster_sweep.py::main
               (sweep loop at r = 0.6: "UMAP by cluster" and "Wilcoxon markers"); save() and the settings are imported
               from the module (main() subsets, embeds and clusters, so it is not called)
Illustrator  : e: subcluster names added as a legend ("0 T", "1 TOX+ T", ... "17 Immune/epithelial"), title removed;
               f: suptitle removed, axis title "snRNAseq immune subclusters" added. The caption of Supp. Fig. 18e says
               "Leiden resolution 0.8"; the panel is Leiden 0.6.
Environment  : conda env pathfire
Verified     : 2026-11-12, PIXEL-IDENTICAL (both PDFs, text layers equal; PNGs byte-identical) vs
               $KB_ANALYSIS/subcluster_level2_20260728/Immune/res_0.6/umap_res0.6.pdf and markers_dotplot_res0.6.pdf
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

sw = load_script(ANALYSES_09 / "subcluster_sweep" / "subcluster_sweep.py")   # module level: settings, PANELS, save()
CT = "Immune"
r = 0.6


def main():
    ROOT = f"{sw.OUTROOT_DEFAULT}/{CT.replace('/', '_')}"
    adata = sc.read_h5ad(f"{ROOT}/{CT.replace('/', '_')}_leiden_sweep.h5ad")
    save = sw.save
    # ---------------- verbatim from subcluster_sweep.py::main (one resolution) ----------------
    key = f"leiden_{r}"
    cats = list(adata.obs[key].cat.categories)
    # UMAP by cluster
    fig, ax = plt.subplots(figsize=(6.4, 5.6))
    sc.pl.umap(adata, color=key, ax=ax, show=False, frameon=False, legend_loc="on data",
               legend_fontsize=9, legend_fontoutline=2, title=f"{CT} leiden {r} ({len(cats)} clusters)")
    fig.tight_layout(); save(fig, fig_out("S18", "S18e_immune_subclusters_umap"))
    # Wilcoxon markers (saved result uns['rgg_leiden_0.6'] instead of re-running rank_genes_groups)
    mk = sc.get.rank_genes_groups_df(adata, group=None, key=f"rgg_{key}")
    top5 = {c: mk[mk.group == c].sort_values("scores", ascending=False).names.head(5).tolist() for c in cats}
    ndot = len({g for v in top5.values() for g in v})
    dp = sc.pl.dotplot(adata, var_names=top5, groupby=key, standard_scale="var", return_fig=True,
                       figsize=(max(10, 0.30 * ndot), 0.55 * len(cats) + 2.2))
    dp.add_totals().style(cmap="Reds"); dp.make_figure()
    dp.fig.suptitle(f"{CT} top-5 Wilcoxon markers/cluster - leiden {r}", fontsize=12, y=1.02)
    save(dp.fig, fig_out("S18", "S18f_immune_subclusters_markers_dotplot"))
    print("[S18e-f] wrote", fig_out("S18", "S18e_immune_subclusters_umap.pdf"), "and the markers dotplot")


if __name__ == "__main__":
    main()

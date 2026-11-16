#!/usr/bin/env python
"""Fig. 6c — Xenium marker genes per cell type inside glomerular niche 20 (top 5 Wilcoxon markers per cell type,
cell types with >= 5 cells in the niche; dot colour = mean expression scaled 0-1 per gene, dot size = fraction of
cells).
Panel        : Main Figures.pdf, Figure 6, panel c
Reads        : $KB_RESULTS/xenium/leon/new_data/manual_glom_niche_based_on_04_res10_level2/xenium_with_niches.h5ad
               $KB_RESULTS/xenium/leon/new_data/manual_glom_niche_based_on_04_res10_level2/niche_marker_results/
               markers_celltypes_within_cluster_20.csv (both from analyses/01_xenium_preprocessing_niches/
               nichecompass_5regions_downstream.py::per_niche_celltype_markers, nichecompass_replot_level2/
               partH_manualglom10_downstream.sh, job 327992)
Writes       : $KB_FIG_OUT/Fig6/Fig6c_niche20_marker_dotplot.png
Plotting code: extracted from analyses/01_xenium_preprocessing_niches/nichecompass_5regions_downstream.py::
               per_niche_celltype_markers (object preparation and dotplot call verbatim; the Wilcoxon ranking is read
               from the saved marker table instead of being recomputed)
Illustrator  : title "Xenium marker genes in Glom niche 20" added, legends restyled
Environment  : conda env nichecompass
Verified     : 2026-11-12, BYTE-IDENTICAL vs
               $KB_RESULTS/xenium/leon/new_data/manual_glom_niche_based_on_04_res10_level2/
               niche_marker_results/dotplot__cluster_20_dotplot.png
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from niches_maps import use_niche_modules  # noqa: E402

use_niche_modules()
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
from scanpy.plotting._utils import savefig_or_show  # noqa: E402  (scanpy's own save step of sc.pl.dotplot(save=...))
from kbpaths import KB_RESULTS, fig_out  # noqa: E402

RUN = f"{KB_RESULTS}/xenium/leon/new_data/manual_glom_niche_based_on_04_res10_level2"
NICHE_KEY = "latent_leiden_1.0"       # --leiden-res 1.0
CELL_TYPE_KEY = "level2_refined"      # --cell-type-key level2_refined
NICHE = "20"


def main():
    adata = sc.read_h5ad(f"{RUN}/xenium_with_niches.h5ad")
    # --- verbatim from per_niche_celltype_markers: normalise in a copy, X = counts for the dotplot ---
    ad2 = adata.copy()
    if "counts" not in ad2.layers:
        ad2.layers["counts"] = ad2.X.copy()
    sc.pp.normalize_total(ad2, target_sum=1e4)
    sc.pp.log1p(ad2)
    ad2.layers["lognorm"] = ad2.X.copy()
    ad2.X = ad2.layers["counts"]
    sub = ad2[ad2.obs[NICHE_KEY].astype(str) == NICHE].copy()
    # drop tiny cell-type groups
    ct_counts = sub.obs[CELL_TYPE_KEY].astype(str).value_counts()
    keep = ct_counts[ct_counts >= 5].index.tolist()
    sub = sub[sub.obs[CELL_TYPE_KEY].astype(str).isin(keep)].copy()
    # --- saved Wilcoxon ranking (rank_genes_groups output, 200 genes per group, in rank order) ---
    tidy = pd.read_csv(f"{RUN}/niche_marker_results/markers_celltypes_within_cluster_{NICHE}.csv")
    groups = list(dict.fromkeys(tidy["group"].astype(str)))
    assert set(groups) == set(sub.obs[CELL_TYPE_KEY].astype(str)), "marker table and object disagree on cell types"
    # top 5 per group, dedup (verbatim logic)
    top5 = []
    for g in groups:
        names = tidy.loc[tidy["group"].astype(str) == g, "gene"].tolist()
        top5.extend([names[i] for i in range(min(5, len(names)))])
    top5 = list(dict.fromkeys(top5))
    sc.pl.dotplot(sub, var_names=top5, groupby=CELL_TYPE_KEY,
                  standard_scale="var", use_raw=False, show=False)
    out = Path(fig_out("Fig6", "Fig6c_niche20_marker_dotplot.png"))
    sc.settings.figdir = out.parent
    savefig_or_show(out.stem, show=False, save=True, ext="png")


if __name__ == "__main__":
    main()

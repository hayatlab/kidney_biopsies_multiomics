#!/usr/bin/env python
"""Fig. 6a — cell-type composition of the glomerular niches (niche x cell type, row-normalised; manual-glom
NicheCompass model, Leiden 1.0, cell types = level2_refined).
Panel        : Main Figures.pdf, Figure 6, panel a
Reads        : $KB_RESULTS/xenium/leon/new_data/manual_glom_niche_based_on_04_res10_level2/xenium_with_niches.h5ad
               (analyses/01_xenium_preprocessing_niches/nichecompass_5regions_downstream.py --leiden-res 1.0
               --cell-type-key level2_refined, via nichecompass_replot_level2/partH_manualglom10_downstream.sh,
               job 327992; the per-cell niche labels latent_leiden_1.0 in this object are the reference)
Writes       : $KB_FIG_OUT/Fig6/Fig6a_niche_vs_celltype_heatmap.{pdf,png,csv}
Plotting code: analyses/01_xenium_preprocessing_niches/nichecompass_5regions_downstream.py::heatmap_niche_vs_celltype
               (imported unchanged; out_name parameter added by build/patch_rules/fig_niches.py)
Illustrator  : cell values and title removed, axis titles replaced ("Glom. niche"), dPT renamed aPT2
Environment  : conda env nichecompass
Verified     : 2026-11-12, PIXEL-IDENTICAL (text layer equal) vs
               $KB_RESULTS/xenium/leon/new_data/manual_glom_niche_based_on_04_res10_level2/figures/
               niche_vs_celltype_heatmap.pdf (CSV byte-identical)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from niches_maps import use_niche_modules  # noqa: E402

use_niche_modules()
import scanpy as sc  # noqa: E402
from kbpaths import KB_RESULTS, fig_out  # noqa: E402
import nichecompass_5regions_downstream as ds  # noqa: E402

H5AD = f"{KB_RESULTS}/xenium/leon/new_data/manual_glom_niche_based_on_04_res10_level2/xenium_with_niches.h5ad"


def main():
    adata = sc.read_h5ad(H5AD)
    # the downstream run sets these module globals from --leiden-res 1.0 --cell-type-key level2_refined
    ds.LEIDEN_RES = 1.0
    ds.NICHE_KEY = "latent_leiden_1.0"
    ds.CELL_TYPE_KEY = "level2_refined"
    out = Path(fig_out("Fig6", "Fig6a_niche_vs_celltype_heatmap.pdf"))
    ds.heatmap_niche_vs_celltype(adata, out.parent, out_name=out.stem)


if __name__ == "__main__":
    main()

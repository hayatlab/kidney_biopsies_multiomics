#!/usr/bin/env python
"""Extended Data Fig. 4g — cell-type composition of the whole-core NicheCompass niches (overall model with
glomeruli, Leiden 0.4; niche x cell type, row-normalised, cell types = level2_refined).
Panel        : Extended Data Figures.pdf, Extended Data Figure 4, panel g
Reads        : $KB_RESULTS/xenium/leon/new_data/overall_with_glom_04_level2/cell_level2_and_niche.csv
               (per-cell latent_leiden_0.4 + level2_refined, written by analyses/01_xenium_preprocessing_niches/
               nichecompass_replot_level2/replot_partA.py::run_version --version overall in the same run as the
               original heatmap, partE_pecrec_replot.sh; identical to the obs of
               overall_with_glom_04_level2/xenium_with_niches_overall_with_glom.h5ad, checked 2026-11-12)
Writes       : $KB_FIG_OUT/ED4/ED4g_niche_vs_celltype_heatmap_noannot.{pdf,png,csv}
Plotting code: analyses/01_xenium_preprocessing_niches/nichecompass_replot_level2/replot_partA.py::heatmap_noannot
               (imported unchanged)
Illustrator  : transposed (cell types as rows, niches as columns), colour bar horizontal and relabelled
               "row proportion", title removed, dPT renamed aPT2
Environment  : conda env nichecompass
Verified     : 2026-11-12, PIXEL-IDENTICAL (text layer and CSV equal) vs
               $KB_RESULTS/xenium/leon/new_data/overall_with_glom_04_level2/figures/
               niche_vs_celltype_heatmap_noannot.pdf
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from niches_maps import use_niche_modules  # noqa: E402

use_niche_modules()
import anndata as ad  # noqa: E402
import pandas as pd  # noqa: E402
from kbpaths import KB_RESULTS, fig_out  # noqa: E402
import replot_partA as rpa  # noqa: E402

RUN = f"{KB_RESULTS}/xenium/leon/new_data/overall_with_glom_04_level2"
NICHE_KEY = rpa.VERSIONS["overall"]["niche_key"]       # latent_leiden_0.4
CT_KEY = rpa.CT_KEY                                     # level2_refined


def main():
    obs = pd.read_csv(f"{RUN}/cell_level2_and_niche.csv", usecols=[NICHE_KEY, CT_KEY])
    adata = ad.AnnData(obs=obs)          # heatmap_noannot uses obs[niche_key] and obs[ct_key] only
    out = Path(fig_out("ED4", "ED4g_niche_vs_celltype_heatmap_noannot.pdf"))
    rpa.heatmap_noannot(adata, NICHE_KEY, CT_KEY, out.with_suffix(""))


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Fig. 6b — glomerular niches and cell types on the Xenium cell boundaries of two cores: K50 (IgAN with a
crescent; glomerular niche 20) and K18 (TIN, "healthy" glomerulus). Whole core drawn, cells outside the glomerular
NicheCompass object in grey (--full-context); left panel niches (manual-glom model, Leiden 1.0), right panel
level2_refined cell types.
Panel        : Main Figures.pdf, Figure 6, panel b (niche map and cell-type zoom of K50; cell-type zoom of K18)
Reads        : $KB_RESULTS/xenium/leon/new_data/manual_glom_niche_based_on_04_res10_level2/xenium_with_niches.h5ad
               (analyses/01_xenium_preprocessing_niches/nichecompass_5regions_downstream.py via
               nichecompass_replot_level2/partH_manualglom10_downstream.sh)
               $KB_ANALYSIS/celltype_concordance/outputs/xenium_labeled.h5ad (obs only: core of every QC-passed cell,
               for the grey background; analyses/02_xenium_celltyping)
               $KB_XENIUM_RAW/output-<region>/{cell,nucleus}_boundaries.parquet (Xenium Onboard Analysis; external)
Writes       : $KB_FIG_OUT/Fig6/Fig6b_K50_XEN093.{pdf,png}, $KB_FIG_OUT/Fig6/Fig6b_K18_XEN116.{pdf,png}
Plotting code: analyses/01_xenium_preprocessing_niches/nichecompass_replot_level2/replot_boundaries_level2.py::run
               --version manual_glom_10 --full-context (imported unchanged; `cores` argument added by
               build/patch_rules/fig_niches.py: plot only these two cores, title/file named by public core id instead
               of the internal core label) — as partI_manualglom10_boundaries.sh, job 327993
Illustrator  : panels cropped to the glomeruli (K50: niche panel with niche 20 circled + cell-type zoom; K18:
               cell-type zoom), titles and legends replaced, H&E/protein images added alongside (not from code)
Environment  : conda env spatial_env (region boundary files: run through sbatch)
Verified     : 2026-11-12, DIFFERENT (K50 0.27 %, K18 0.24 % of pixels) only in the title row, where the internal
               core label is replaced by the public id; pixel-identical below it and text layers equal after the same
               substitution, vs $KB_RESULTS/xenium/leon/new_data/manual_glom_niche_based_on_04_res10_level2/
               cell_boundary_plots_per_core/<internal label of XEN093>.pdf and <internal label of XEN116>.pdf
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from niches_maps import core_targets, use_niche_modules  # noqa: E402

use_niche_modules()
from kbpaths import fig_out  # noqa: E402
import replot_boundaries_level2 as rbl  # noqa: E402

CORES = {"XEN093": "K50", "XEN116": "K18"}      # public core id -> sample


def main():
    cfg = dict(rbl.VERSIONS["manual_glom_10"])
    cfg["dest"] = Path(fig_out("Fig6", "Fig6b_K50_XEN093.pdf")).parent
    rbl.run(cfg, full_context=True, cores=core_targets("Fig6", "Fig6b", CORES))


if __name__ == "__main__":
    main()

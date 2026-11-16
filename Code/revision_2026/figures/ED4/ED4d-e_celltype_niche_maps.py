#!/usr/bin/env python
"""Extended Data Fig. 4d-e — level2_refined cell types (d) and whole-core NicheCompass niches (e; overall model with
glomeruli, Leiden 0.4) on the Xenium cell boundaries of K82 (FSGS, high eGFR) and K18 (TIN, low eGFR).
Panel        : Extended Data Figures.pdf, Extended Data Figure 4, panels d (cell-type panel of K82 and of K18) and
               e (niche panel of K18)
Reads        : $KB_RESULTS/xenium/leon/archived/overall_with_glom_04/xenium_with_niches_overall_with_glom.h5ad
               (symlink to overall_with_glom_04_level2/xenium_with_niches_overall_with_glom.h5ad; analyses/
               01_xenium_preprocessing_niches/nichecompass_overall_with_glom.py -> nichecompass_5regions_downstream.py,
               level2_refined baked in by nichecompass_replot_level2/bake_level2_into_niche_objects.py; the original
               run merged the same labels from celltype_concordance/outputs/xenium_labeled.h5ad — identical for all
               258,375 cells, checked 2026-11-12)
               $KB_XENIUM_RAW/output-<region>/{cell,nucleus}_boundaries.parquet (Xenium Onboard Analysis; external)
Writes       : $KB_FIG_OUT/ED4/ED4d-e_K82_XEN119.{pdf,png}, $KB_FIG_OUT/ED4/ED4d-e_K18_XEN116.{pdf,png}
Plotting code: analyses/01_xenium_preprocessing_niches/nichecompass_replot_level2/replot_boundaries_level2.py::run
               --version overall (imported unchanged; `cores` argument added by build/patch_rules/fig_niches.py:
               plot only these two cores, title/file named by public core id instead of the internal core label) —
               as partE_pecrec_boundaries.sh, job 327951
Illustrator  : single panels taken from the two-panel PDFs, titles "FSGS with high eGFR" / "TIN with low eGFR" /
               "niches", legends restyled, dPT renamed aPT2, box "enriched for niche 1&4" drawn on the niche panel
Environment  : conda env spatial_env (1 GB object + region boundary files: run through sbatch)
Verified     : 2026-11-12, DIFFERENT (K82 0.16 %, K18 0.16 % of pixels) only in the title row, where the internal
               core label is replaced by the public id; pixel-identical below it and text layers equal after the same
               substitution, vs $KB_RESULTS/xenium/leon/new_data/overall_with_glom_04_level2/
               cell_boundary_plots_per_core/<internal label of XEN119>.pdf and <internal label of XEN116>.pdf
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from niches_maps import core_targets, use_niche_modules  # noqa: E402

use_niche_modules()
from kbpaths import fig_out  # noqa: E402
import replot_boundaries_level2 as rbl  # noqa: E402

CORES = {"XEN119": "K82", "XEN116": "K18"}      # public core id -> sample


def main():
    cfg = dict(rbl.VERSIONS["overall"])
    cfg["dest"] = Path(fig_out("ED4", "ED4d-e_K82_XEN119.pdf")).parent
    rbl.run(cfg, full_context=False, cores=core_targets("ED4", "ED4d-e", CORES))


if __name__ == "__main__":
    main()

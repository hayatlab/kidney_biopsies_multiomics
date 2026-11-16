#!/usr/bin/env python
"""Extended Data Fig. 4h — enrichment of eGFR groups (<30, 30-60, 60-90, >90 mL/min/1.73 m2, eGFR at biopsy of the
core's donor) in the whole-core NicheCompass niches (overall model with glomeruli, Leiden 0.4): one-sided Fisher test
per niche x group on cells, BH-adjusted; colour = sign(log2 enrichment) x -log10(FDR); niches with >= 20 cells.
Panel        : Extended Data Figures.pdf, Extended Data Figure 4, panel h
Reads        : $KB_RESULTS/xenium/leon/new_data/overall_with_glom_04_level2/xenium_with_niches_overall_with_glom.h5ad
               (analyses/01_xenium_preprocessing_niches/nichecompass_overall_with_glom.py ->
               nichecompass_5regions_downstream.py; level2_refined baked in by bake_level2_into_niche_objects.py)
               $KB_ANALYSIS/celltype_concordance/outputs/per_sample_metadata_updated.csv (per-core clinical table,
               identifier join; analyses/02_xenium_celltyping/15_update_metadata_slope.py)
Writes       : $KB_FIG_OUT/ED4/ED4h_egfr_group_enrichment_signedlog10fdr_clustermap.{pdf,png} (panel) plus the other
               outputs of run_categorical_enrichment with the same prefix (counts, proportions, Fisher/binomial tables)
Plotting code: analyses/01_xenium_preprocessing_niches/nichecompass_5regions_downstream.py::run_categorical_enrichment
               (imported unchanged; prefix argument = output file prefix), on the object prepared by
               nichecompass_replot_level2/regenerate_metadata_figures.py::load_config("overall04") (load_config
               factored out of run_config by build/patch_rules/fig_niches.py) — as in partK_joinfix_regen.sh, job 328177
Illustrator  : dendrograms removed, columns put in clinical order (30-60 and 60-90 swapped), column labels moved to
               the top with the title "eGFR (mL/min per 1.73m2)", axis titles "Niche", colour bar relabelled
               "signed log10(FDR)"
Environment  : conda env nichecompass (1 GB object: run through sbatch)
Verified     : 2026-11-12, PIXEL-IDENTICAL (text layer equal) vs
               $KB_RESULTS/xenium/leon/new_data/overall_with_glom_04_level2/figures/
               egfr_group_enrichment_signedlog10fdr_clustermap.pdf
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from niches_maps import use_niche_modules  # noqa: E402

use_niche_modules()
from kbpaths import fig_out  # noqa: E402
import regenerate_metadata_figures as rmf  # noqa: E402

CONFIG = "overall04"
CAT = "egfr_group"
PREFIX = "ED4h_egfr_group"


def main():
    adata = rmf.load_config(rmf.CONFIGS[CONFIG])
    fig_dir = Path(fig_out("ED4", f"{PREFIX}_enrichment_signedlog10fdr_clustermap.pdf")).parent
    # as regenerate_metadata_figures.regen_dir for cat = egfr_group (min_niche_cells=20)
    rmf.ds.run_categorical_enrichment(adata, CAT, fig_dir, PREFIX, min_niche_cells=20)


if __name__ == "__main__":
    main()

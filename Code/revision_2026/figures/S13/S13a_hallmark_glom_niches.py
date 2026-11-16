#!/usr/bin/env python
"""Supp. Fig. 13a — MSigDB Hallmark pathway activity across the glomerular niches (manual-glom NicheCompass model,
Leiden 1.0): decoupler ULM on log-normalised counts, top-3 pathways per niche (rankby_group, t-test), per-niche mean
score min-max scaled per pathway, clustered heatmap without the top-marker stars.
Panel        : All Supp Figures.pdf, Supplementary Figure 13, panel a
Reads        : $KB_RESULTS/xenium/leon/new_data/manual_glom_niche_based_on_04_res10_level2/xenium_with_niches.h5ad
               (analyses/01_xenium_preprocessing_niches/nichecompass_5regions_downstream.py via
               nichecompass_replot_level2/partH_manualglom10_downstream.sh); MSigDB Hallmark gene sets fetched at run
               time from OmniPath by decoupler (dc.op.hallmark) — external, not versioned
Writes       : $KB_FIG_OUT/S13/S13a_hallmark_marker_pathways_dedup_heatmap_nostars.{pdf,png} (panel) plus the
               starred heatmap, dotplot, S13a_hallmark_score_means.csv and S13a_hallmark_top3_per_niche.csv
Plotting code: analyses/01_xenium_preprocessing_niches/glom_postprocess.py::hallmark_pathway_dotplot (imported
               unchanged; prefix parameter added by build/patch_rules/fig_niches.py), called as in
               nichecompass_replot_level2/regenerate_hallmark_dedup.py::run_config --config manual_glom10
Illustrator  : colour bar moved and relabelled "min-max scaled scores", dendrograms removed, x-axis title removed,
               column labels re-typed (the 16th column is labelled "10" a second time; it is niche 16)
Environment  : conda env nichecompass
Verified     : 2026-11-12, PIXEL-IDENTICAL (text layer equal) vs
               $KB_RESULTS/xenium/leon/new_data/manual_glom_niche_based_on_04_res10_level2/
               figures/hallmark_marker_pathways_dedup_heatmap_nostars.pdf
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from niches_maps import use_niche_modules  # noqa: E402

use_niche_modules()
import scanpy as sc  # noqa: E402
from kbpaths import fig_out  # noqa: E402
import regenerate_hallmark_dedup as rhd  # noqa: E402

CONFIG = "manual_glom10"
PREFIX = "S13a_"


def main():
    cfg = rhd.CONFIGS[CONFIG]
    adata = sc.read_h5ad(cfg["adata"])                    # as regenerate_hallmark_dedup.run_config
    out_dir = Path(fig_out("S13", f"{PREFIX}hallmark_marker_pathways_dedup_heatmap_nostars.pdf")).parent
    rhd.gp.hallmark_pathway_dotplot(adata, niche_key=cfg["niche_key"], out_dir=out_dir,
                                    n_markers=rhd.N_MARKERS, prefix=PREFIX)


if __name__ == "__main__":
    main()

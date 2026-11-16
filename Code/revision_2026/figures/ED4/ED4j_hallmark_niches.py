#!/usr/bin/env python
"""Extended Data Fig. 4j — MSigDB Hallmark pathway activity across the whole-core NicheCompass niches (overall model
with glomeruli, Leiden 0.4): decoupler ULM on log-normalised counts, top-3 pathways per niche (rankby_group,
t-test), per-niche mean score min-max scaled per pathway, clustered heatmap without the top-marker stars.
Panel        : Extended Data Figures.pdf, Extended Data Figure 4, panel j
Reads        : $KB_RESULTS/xenium/leon/new_data/overall_with_glom_04_level2/xenium_with_niches_overall_with_glom.h5ad
               (analyses/01_xenium_preprocessing_niches/nichecompass_overall_with_glom.py ->
               nichecompass_5regions_downstream.py; level2_refined baked in by
               nichecompass_replot_level2/bake_level2_into_niche_objects.py); MSigDB Hallmark gene sets fetched at
               run time from OmniPath by decoupler (dc.op.hallmark) — external, not versioned
Writes       : $KB_FIG_OUT/ED4/ED4j_hallmark_marker_pathways_dedup_heatmap_nostars.{pdf,png} (panel) plus the
               starred heatmap, dotplot, ED4j_hallmark_score_means.csv and ED4j_hallmark_top3_per_niche.csv
Plotting code: analyses/01_xenium_preprocessing_niches/glom_postprocess.py::hallmark_pathway_dotplot (imported
               unchanged; prefix parameter added by build/patch_rules/fig_niches.py), called as in
               nichecompass_replot_level2/regenerate_hallmark_dedup.py::run_config --config overall04
Illustrator  : transposed (niches as rows, pathways as columns), dendrograms removed,
               EPITHELIAL_MESENCHYMAL_TRANSITION shortened to EMT, colour bar relabelled "Z-scaled scores"
               (the values are min-max scaled, not z-scores)
Environment  : conda env nichecompass (1 GB object: run through sbatch)
Verified     : 2026-11-12, PIXEL-IDENTICAL (text layer equal) vs
               $KB_RESULTS/xenium/leon/new_data/overall_with_glom_04_level2/figures/
               hallmark_marker_pathways_dedup_heatmap_nostars.pdf
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from niches_maps import use_niche_modules  # noqa: E402

use_niche_modules()
import scanpy as sc  # noqa: E402
from kbpaths import fig_out  # noqa: E402
import regenerate_hallmark_dedup as rhd  # noqa: E402

CONFIG = "overall04"
PREFIX = "ED4j_"


def main():
    cfg = rhd.CONFIGS[CONFIG]
    adata = sc.read_h5ad(cfg["adata"])                    # as regenerate_hallmark_dedup.run_config
    out_dir = Path(fig_out("ED4", f"{PREFIX}hallmark_marker_pathways_dedup_heatmap_nostars.pdf")).parent
    rhd.gp.hallmark_pathway_dotplot(adata, niche_key=cfg["niche_key"], out_dir=out_dir,
                                    n_markers=rhd.N_MARKERS, prefix=PREFIX)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Extended Data Fig. 4i — condition composition of the whole-core NicheCompass niches (overall model with glomeruli,
Leiden 0.4): cells per niche stacked by condition (left panel of the original figure; niches with >= 50 cells,
ordered by size).
Panel        : Extended Data Figures.pdf, Extended Data Figure 4, panel i
Reads        : $KB_RESULTS/xenium/leon/new_data/overall_with_glom_04_level2/xenium_with_niches_overall_with_glom.h5ad
               (analyses/01_xenium_preprocessing_niches/nichecompass_overall_with_glom.py ->
               nichecompass_5regions_downstream.py; level2_refined baked in by bake_level2_into_niche_objects.py)
               $KB_ANALYSIS/celltype_concordance/outputs/per_sample_metadata_updated.csv (applied by load_config;
               `condition` itself comes from the object and does not depend on it)
Writes       : $KB_FIG_OUT/ED4/ED4i_condition_stacked_bars.{pdf,png}
Plotting code: analyses/01_xenium_preprocessing_niches/nichecompass_5regions_downstream.py::plot_stacked_proportions
               (imported unchanged), on the object prepared by nichecompass_replot_level2/
               regenerate_metadata_figures.py::load_config("overall04") (load_config factored out of run_config by
               build/patch_rules/fig_niches.py) — as in partK_joinfix_regen.sh, job 328177
Illustrator  : left panel (absolute counts) only; "n=" labels removed, title "Condition composition of niches",
               legend moved into the panel, y axis "cell counts"
Environment  : conda env nichecompass (1 GB object: run through sbatch)
Verified     : 2026-11-12, PIXEL-IDENTICAL (text layer equal) vs
               $KB_RESULTS/xenium/leon/new_data/overall_with_glom_04_level2/figures/
               condition_stacked_bars.pdf
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from niches_maps import use_niche_modules  # noqa: E402

use_niche_modules()
from kbpaths import fig_out  # noqa: E402
import regenerate_metadata_figures as rmf  # noqa: E402

CONFIG = "overall04"
CAT = "condition"


def main():
    adata = rmf.load_config(rmf.CONFIGS[CONFIG])
    out = Path(fig_out("ED4", "ED4i_condition_stacked_bars.pdf"))
    # as regenerate_metadata_figures.regen_dir for cat = condition (min_total=50)
    rmf.ds.plot_stacked_proportions(adata, CAT, out.with_suffix(""), min_total=50)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Fig. 6j (heatmap) — z-scored mean expression of 1,000 highly variable genes per pseudotime bin along the aPEC
trajectory (VIA pseudotime on PHATE, PEC subclusters 0, 2, 4, 9, 11); genes ordered by peak pseudotime.
Panel        : Main Figures.pdf, Figure 6, panel j (right, "z-scored mean expression per pseudotime bin")
Reads        : $KB_ANALYSIS/pec_resubcluster_20260717/trajectories_leiden08/C_0to9_5cl/subset_all_pseudotimes.h5ad
                 (analyses/09_snrna_pec_immune_subclusters/pec/trajectory_leiden08.py --tags C_0to9_5cl, job 328231;
                  obs pseudotime_via_phate and log-normalised X are read, nothing is recomputed upstream)
Writes       : $KB_FIG_OUT/Fig6/Fig6j_pseudotime_heatmap.{pdf,png}, Fig6j_pseudotime_heatmap_geneorder.csv
Plotting code: analyses/09_snrna_pec_immune_subclusters/pec/heatmap_C5_phate_zbin.py::main (imported; patch rule in
               build/patch_rules/fig_niches.py adds the out= keyword, default = original location). main() itself
               selects the 1,000 HVGs (scanpy seurat flavour), bins cells into 120 equal-count pseudotime bins and
               z-scores per gene — light, deterministic steps that are part of the plotting script.
Illustrator  : y-axis label reworded ("genes (ordered by peak pseudotime)" -> "genes ordered by peak pseudotime")
Environment  : conda env ViaEnv
Verified     : 2026-11-12, PIXEL-IDENTICAL (PDF, text layer equal; PNG and gene-order CSV byte-identical) vs
               $KB_ANALYSIS/pec_resubcluster_20260717/trajectories_leiden08/C_0to9_5cl/phate/pseudotime_heatmap_zbinmean_HVGs.pdf
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from niches_snrna import ANALYSES_09, load_script  # noqa: E402
from kbpaths import fig_out  # noqa: E402

hm = load_script(ANALYSES_09 / "pec" / "heatmap_C5_phate_zbin.py")

if __name__ == "__main__":
    out = fig_out("Fig6", "Fig6j_pseudotime_heatmap")      # stem: .pdf/.png/_geneorder.csv are appended
    hm.main(out=out)

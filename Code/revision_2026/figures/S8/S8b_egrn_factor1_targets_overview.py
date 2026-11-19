"""Supplementary Fig. 8b — the Factor-1 top-positive PT genes that are targets in the r > 0.7 scMEGA eGRN (43 of 80),
ranked by the number of TFs regulating each gene; bars coloured by network module.
Panel        : All Supp Figures.pdf, Supplementary Figure 8, panel b
Reads        : $KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/ATAC/scmega/egrn_tutorial/archr_dm/{grn_edges.csv,
                 grn_edges_all.csv, tf_centrality.csv, meta.json} (analyses/10_atac_pt_apt_egrn/ATAC/scripts/
                 43_scmega_egrn_tutorial.R, TRAJ=archr_dm)
               $KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/cis_regulatory_pt_apt/outputs/task3/
                 C4_factor1_top200_positive_genes.csv (columns variable, Loading; rebuilt by
                 analyses/10_atac_pt_apt_egrn/make_factor1_pt_list.py)
Writes       : $KB_FIG_OUT/S8/S8b_egrn_factor1_targets_overview/egrn_factor1_targets_overview.{pdf,png} (the panel) and
               the other outputs of the script (network view = ED Fig. 3m, per-gene views, tables)
Plotting code: analyses/10_atac_pt_apt_egrn/ATAC/scripts/45_egrn_paper_views.py::main (imported unchanged; run with
               --run archr_dm as in job 339554 of 26 Oct 2026; module-level OUT redirected to $KB_FIG_OUT)
Illustrator  : title re-worded ("Factor 1 top positive genes in the eGRN (43 of 80 are targets)"); module legend
               removed (all 43 genes are in module I); panel letter
Environment  : conda env agent_env
Verified     : 2026-11-12, PIXEL-IDENTICAL (PNG byte-identical) vs
               $KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/ATAC/scmega/figures_paper/egrn_factor1_targets_overview.pdf
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from atac_egrn import run_main  # noqa: E402  (also puts config/ on sys.path)
from kbpaths import fig_out  # noqa: E402

OUT = fig_out("S8", "S8b_egrn_factor1_targets_overview")
run_main("45_egrn_paper_views.py", ["--run", "archr_dm"], OUT)
print(f"[S8b] panel: {OUT}/egrn_factor1_targets_overview.pdf", flush=True)

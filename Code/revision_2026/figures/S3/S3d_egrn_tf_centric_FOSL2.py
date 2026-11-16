"""Supplementary Fig. 3d — FOSL2 regulon in the r > 0.7 scMEGA eGRN: FOSL2 at the centre, its 114 target genes on a
ring ordered by TF-gene correlation; the 19 Factor-1 top-positive targets in orange and labelled (edge width = r).
Panel        : All Supp Figures.pdf, Supplementary Figure 3, panel d
Reads        : $KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/ATAC/scmega/figures_paper/egrn_paper_views_edges.csv
                 (analyses/10_atac_pt_apt_egrn/ATAC/scripts/45_egrn_paper_views.py)
               $KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/ATAC/scmega/figures_paper/egrn_factor1_tf_side.csv
                 (analyses/10_atac_pt_apt_egrn/ATAC/scripts/47_egrn_factor1_tf_side.py)
               $KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/cis_regulatory_pt_apt/outputs/task3/
                 C4_factor1_top200_positive_genes.csv (column variable; rebuilt by
                 analyses/10_atac_pt_apt_egrn/make_factor1_pt_list.py)
Writes       : $KB_FIG_OUT/S3/S3d_egrn_tf_centric_FOSL2/egrn_tf_centric/FOSL2.{pdf,png} (the panel) and the script's
               other TF-centred views (12 TFs; REL is not in the network), the multipage PDF, the composite and
               48_meta.json
Plotting code: analyses/10_atac_pt_apt_egrn/ATAC/scripts/48_egrn_tf_centric_panels.py::main (imported unchanged; default
               TF list as in job 339554 of 26 Oct 2026; module-level OUT redirected to $KB_FIG_OUT)
Illustrator  : title shortened to "FOSL2: 114 targets in the r>0.7 eGRN; 19 of them factor-1 top-positive genes";
               bottom legend removed; panel letter
Environment  : conda env agent_env
Verified     : 2026-11-12, PIXEL-IDENTICAL (PNG byte-identical) vs
               $KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/ATAC/scmega/figures_paper/egrn_tf_centric/FOSL2.pdf
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from atac_egrn import run_main  # noqa: E402  (also puts config/ on sys.path)
from kbpaths import fig_out  # noqa: E402

OUT = fig_out("S3", "S3d_egrn_tf_centric_FOSL2")
run_main("48_egrn_tf_centric_panels.py", [], OUT)
print(f"[S3d] panel: {OUT}/egrn_tf_centric/FOSL2.pdf", flush=True)

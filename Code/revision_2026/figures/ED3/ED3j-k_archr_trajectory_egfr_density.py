"""Extended Data Fig. 3j-k — PT->aPT scMEGA trajectory: joint diffusion map coloured by the ArchR pseudotime (j) and
eGFR groups along that pseudotime (k; per-donor density, mean +/- 95% t-interval across donors).
Panel        : Extended Data Figures.pdf, Extended Data Figure 3, panels j and k (j = panel "a", k = panel "b′" of the
               output file; its panel "b", disease groups, is not in the figure)
Reads        : $KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/ATAC/tables/HANDOFF_joint_pseudotime.csv
                 (analyses/10_atac_pt_apt_egrn/ATAC/scripts/12_joint_pseudotime.py)
               $KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/ATAC/tables/archr_trajectories.csv
                 (analyses/10_atac_pt_apt_egrn/ATAC/scripts/31_archr_trajectory.R)
               $KB_ANALYSIS/vivien_takeover/figures/milo_donor_composition/sample_metadata_milo.csv (only the columns
                 sample, disease, egfr_group are read; producer vivien_takeover/scripts/22_sample_metadata_table.py, not
                 part of this repository)
Writes       : $KB_FIG_OUT/ED3/ED3j-k_archr_trajectory_egfr_density/panel_a_b_archr_dm.{pdf,png} and
               panel_b_donor_densities_archr_dm.csv (the per-donor curves behind panels b / b′)
Plotting code: analyses/10_atac_pt_apt_egrn/ATAC/scripts/34_panels_a_b_archr.py, run unchanged with the options of the
               published version (TRAJ=archr_dm ONLY_USED=1, job 339466 of 26 Oct 2026); FIGDIR and TABLEDIR (release
               patch, default = original folder) point to $KB_FIG_OUT
Illustrator  : j: title replaced by "scMEGA on GLUE - ArchR trajectory from PT to aPT / joint diffusion map & ArchR
               pseudotime", PT / aPT labels re-set in plain type, colour-bar label removed, axis titles DC1 / DC2, panel
               letter; k: only panel b′ is used, title "eGFR group along the pseudotime", donor counts removed from the
               legend; panel b (disease) dropped.
               Sample K68 has no disease / eGFR group in sample_metadata_milo.csv and is therefore not in k (printed by
               the script; see analyses/10_atac_pt_apt_egrn/README.md)
Environment  : conda env agent_env
Verified     : 2026-11-12, PIXEL-IDENTICAL (PNG byte-identical; density table byte-identical) vs
               $KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/ATAC/scmega/figures_paper/panel_a_b_archr_dm.pdf
"""
import os
import runpy
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from atac_egrn import SCRIPTS  # noqa: E402  (also puts config/ on sys.path)
from kbpaths import fig_out  # noqa: E402

OUT = fig_out("ED3", "ED3j-k_archr_trajectory_egfr_density")
os.environ.update(TRAJ="archr_dm", ONLY_USED="1", FIGDIR=OUT, TABLEDIR=OUT)
os.makedirs(OUT, exist_ok=True)
runpy.run_path(str(SCRIPTS / "34_panels_a_b_archr.py"), run_name="__main__")
print(f"[ED3j-k] panel j = left, panel k = right ('b′') of {OUT}/panel_a_b_archr_dm.pdf", flush=True)

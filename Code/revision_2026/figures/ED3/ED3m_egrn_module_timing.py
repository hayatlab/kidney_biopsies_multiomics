"""Extended Data Fig. 3m — the two eGRN modules along the PT->aPT pseudotime: z-scored TF motif activity, TF expression
and target expression per TF (thin lines) and module means with 95% CI (thick lines); grey band = both modules
mid-transition; dashed line = crossing of the module means.
Panel        : Extended Data Figures.pdf, Extended Data Figure 3, panel m (top: TF activity = panel A of the output
               file; middle: TF expression = panel B; bottom: Target expression = panel C)
Reads        : $KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/ATAC/scmega/egrn_tutorial/archr_dm/{grn_edges.csv,
                 pseudotime_curves.csv.gz, selected_tfs.csv} (analyses/10_atac_pt_apt_egrn/ATAC/scripts/
                 43_scmega_egrn_tutorial.R, TRAJ=archr_dm)
Writes       : $KB_FIG_OUT/ED3/ED3m_egrn_module_timing/egrn_module_timing.{pdf,png} (panels A-C are ED 3m) and the
               script's tables (per-TF timing, module curves, permutation tests, cross-correlation, meta json)
Plotting code: analyses/10_atac_pt_apt_egrn/ATAC/scripts/46_egrn_module_timing.py::main (imported unchanged; default
               options as in job 339554 of 26 Oct 2026; module-level OUT redirected to $KB_FIG_OUT)
Illustrator  : panels A-C of the output transposed into one column; titles shortened to "TF activity", "TF expression",
               "Target expression"; per-panel subtitles, panels D-F, suptitle and caveat text removed; module labels
               placed on the curves; x-axis titles re-typed ("pseudotime bin (0=PT like to 100 = aPT like)").
               NOTE: the printed bottom panel (Target expression) is panel C of the SEPTEMBER version of this figure
               (glue_v2 embedding, before the AKI / pt_C5 exclusion:
               $KB_ANALYSIS/ATAC/scmega/figures_paper/egrn_module_timing.pdf, 25 Sep 2026; module crossing at bin 60).
               This script draws the rerun version (glue_v3), whose panel C differs (crossing at bin 59, transition
               window 56-73 instead of 58-72)
Environment  : conda env agent_env
Verified     : 2026-11-12, PIXEL-IDENTICAL (PNG and all tables identical) vs
               $KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/ATAC/scmega/figures_paper/egrn_module_timing.pdf (rerun;
               source of the printed top and middle panels); DIFFERENT (15% of pixels) vs the September file that
               supplied the printed bottom panel, as expected
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from atac_egrn import run_main  # noqa: E402  (also puts config/ on sys.path)
from kbpaths import fig_out  # noqa: E402

OUT = fig_out("ED3", "ED3m_egrn_module_timing")
run_main("46_egrn_module_timing.py", [], OUT)
print(f"[ED3m] panels A-C of {OUT}/egrn_module_timing.pdf", flush=True)

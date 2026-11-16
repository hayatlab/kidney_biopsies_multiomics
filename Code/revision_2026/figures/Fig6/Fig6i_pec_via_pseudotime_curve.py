#!/usr/bin/env python
"""Fig. 6i — VIA pseudotime of the aPEC trajectory (root PEC subcluster 0 -> terminal subcluster 9) on the PHATE layout,
with the VIA lineage curve.
Panel        : Main Figures.pdf, Figure 6, panel i
Reads        : nothing (pointer script, see below)
Writes       : nothing
Plotting code: analyses/09_snrna_pec_immune_subclusters/pec/trajectory_leiden08.py::run_embedding, step "(3) VIA lineage
               curves" (pyVIA 0.2.4 plot_trajectory_curves on the PHATE coordinates; right panel of
               trajectories_leiden08/C_0to9_5cl/phate/via_lineage_curves.pdf, job 328231)
Illustrator  : right panel ("pseudotime") of via_lineage_curves.pdf kept; title "pseudotime snRNA-seq"
Environment  : conda env ViaEnv
Verified     : 2026-11-12, NOT REGENERATED — this script only points at the analysis step (reason below)

Why there is no re-draw: plot_trajectory_curves() needs the fitted VIA object (VIA clusters, cluster graph, root and
terminal states, raw Markov pseudotime). The object was never saved, and VIA is not seed-deterministic (job logs of
three identical runs differ; the 2026-07-30 VIA rebuild in genetrend_niche_markers.py correlates rho = 0.985 with the
stored pseudotime), so the curve cannot be reproduced exactly. What IS saved — the PHATE layout and the per-cell VIA
pseudotime of job 328231 — is drawn by figures/Fig6/Fig6h_pec_phate_clusters_crescent.py (panel 2 of that figure is
this pseudotime without the curve).

Regenerating the panel means re-running the analysis step, which recomputes the subset embeddings and VIA and writes
into the analysis workspace (not into $KB_FIG_OUT):
    source config/env.sh
    sbatch analyses/09_snrna_pec_immune_subclusters/pec/run_trajectory_leiden08.sh --tags C_0to9_5cl --embeddings phate
    -> $KB_ANALYSIS/pec_resubcluster_20260717/trajectories_leiden08/C_0to9_5cl/phate/via_lineage_curves.pdf (right panel)
"""
import sys

if __name__ == "__main__":
    sys.exit(__doc__.split("Why there is no re-draw:")[0].splitlines()[0] + "\n\nNOT REGENERATED from saved outputs."
             "\nWhy there is no re-draw:" + __doc__.split("Why there is no re-draw:")[1])

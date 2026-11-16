#!/usr/bin/env python
"""Fig. 6j (trends) — expression of niche-20 marker genes along the aPEC VIA pseudotime (MAGIC-imputed, VIA gene
trends): CD44, ITGB1 (up in niche 20) and SLC4A11, NDRG1 (down in niche 20) of the top-5 / bottom-5 panels.
Panel        : Main Figures.pdf, Figure 6, panel j (left, four gene-trend plots)
Reads        : nothing (pointer script, see below)
Writes       : nothing
Plotting code: analyses/09_snrna_pec_immune_subclusters/pec/genetrend_niche_markers.py::main (VIA do_impute, 3 MAGIC
               steps + pyVIA get_gene_expression) -> trajectories_leiden08/C_0to9_5cl/phate/genetrend_UP_top5.pdf and
               genetrend_DOWN_bottom5.pdf (job 328398, 2026-07-30)
Illustrator  : CD44 and ITGB1 taken from genetrend_UP_top5.pdf, SLC4A11 and NDRG1 from genetrend_DOWN_bottom5.pdf;
               rearranged as a 2 x 2 grid; axis labels renamed ("Intensity" -> "expression", "time" -> "pseudotime");
               lineage legend and suptitle removed
Environment  : conda env ViaEnv
Verified     : 2026-11-12, NOT REGENERATED — this script only points at the analysis step (reason below)

Why there is no re-draw: the trend curves are computed from a fitted VIA object (lineage probabilities, terminal state,
Markov pseudotime) that was never saved. genetrend_niche_markers.py therefore rebuilds VIA on the saved PHATE layout
with the parameters of trajectory_leiden08.py; VIA is not seed-deterministic, and the rebuild used for the published
panels correlated rho = +0.985 (Spearman) with the stored pseudotime of Fig. 6h/i (log genetrend_niche.328398.out), i.e.
the trends come from a second, slightly different VIA realisation. A new run gives another realisation.

Regenerating the panels means re-running the analysis step (it reads the saved PHATE layout, rebuilds VIA and writes
into the analysis workspace, not into $KB_FIG_OUT):
    source config/env.sh
    sbatch analyses/09_snrna_pec_immune_subclusters/pec/run_genetrend_niche.sh
    -> $KB_ANALYSIS/pec_resubcluster_20260717/trajectories_leiden08/C_0to9_5cl/phate/genetrend_UP_top5.pdf, genetrend_DOWN_bottom5.pdf
"""
import sys

if __name__ == "__main__":
    sys.exit(__doc__.split("Why there is no re-draw:")[0].splitlines()[0] + "\n\nNOT REGENERATED from saved outputs."
             "\nWhy there is no re-draw:" + __doc__.split("Why there is no re-draw:")[1])

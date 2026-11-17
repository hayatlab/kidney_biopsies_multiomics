"""Extended Data Fig. 3l (network) — scMEGA eGRN (r > 0.7) on the ArchR PT->aPT diffusion-map trajectory: two disjoint
modules (I: AP-1 / NF-kB / STAT injury programme; II: HNF4A / RXRA / VDR differentiated-PT programme), TFs coloured by
the pseudotime bin of peak motif activity, Factor-1 top-positive target genes in orange.
Panel        : Extended Data Figures.pdf, Extended Data Figure 3, panel l (network; the JUN / HNF4A insets are drawn by
               ED3l_insets_scMEGA_TF_JUN_HNF4A.R)
Reads        : $KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/ATAC/scmega/egrn_tutorial/archr_dm/{grn_edges.csv,
                 grn_edges_all.csv, tf_centrality.csv, meta.json} (analyses/10_atac_pt_apt_egrn/ATAC/scripts/
                 43_scmega_egrn_tutorial.R, TRAJ=archr_dm: peak-to-gene FDR < 1e-50)
               $KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/cis_regulatory_pt_apt/outputs/task3/
                 C4_factor1_top200_positive_genes.csv (columns variable, Loading; rebuilt by
                 analyses/10_atac_pt_apt_egrn/make_factor1_pt_list.py)
Writes       : $KB_FIG_OUT/ED3/ED3l_egrn_two_modules/egrn_archr_dm_two_modules.{pdf,png} (the panel) and the other
               outputs of the script (per-gene views, node / edge tables, egrn_factor1_targets_overview = Supp. Fig. 8b)
Plotting code: analyses/10_atac_pt_apt_egrn/ATAC/scripts/45_egrn_paper_views.py::main (imported unchanged; run with
               --run archr_dm as in job 339554 of 26 Oct 2026; module-level OUT redirected to $KB_FIG_OUT)
Revision     : 2026-11-13, TF and target labels 1.4x larger (--label-scale, default here 1.4; writes
               egrn_archr_dm_two_modules_labels1.4x.{pdf,png}); --label-scale 1 draws the published version
Illustrator  : suptitle replaced by "scMEGA eGRN on the ArchR PT to aPT diffusion map trajectory"; module titles
               shortened to "module I - AP-1/NF-kB/STAT injury programme" and "module II - HNF4A/RXRA/VDR differentiated
               PT programme" (node / edge / TF counts removed); legend re-worded and colour bar re-labelled ("0=PT like
               to 100 = aPT like"; the scale is bins 1-100) below the network; JUN / HNF4A insets placed with leader
               lines; panel letter
Environment  : conda env agent_env
Verified     : 2026-11-12, PIXEL-IDENTICAL (PNG byte-identical; node / edge tables identical) vs
               $KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/ATAC/scmega/figures_paper/egrn_archr_dm_two_modules.pdf
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from atac_egrn import run_main  # noqa: E402  (also puts config/ on sys.path)
from kbpaths import fig_out  # noqa: E402

ap = argparse.ArgumentParser(description="Extended Data Fig. 3l network")
ap.add_argument("--label-scale", type=float, default=1.4, help="label font scale (1 = published version)")
args = ap.parse_args()
OUT = fig_out("ED3", "ED3l_egrn_two_modules")
run_main("45_egrn_paper_views.py", ["--run", "archr_dm", "--label-scale", f"{args.label_scale:g}"], OUT)
stem = "egrn_archr_dm_two_modules" + ("" if args.label_scale == 1 else f"_labels{args.label_scale:g}x")
print(f"[ED3l] network panel: {OUT}/{stem}.pdf", flush=True)

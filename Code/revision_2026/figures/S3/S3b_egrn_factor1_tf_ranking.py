"""Supplementary Fig. 3b — TFs upstream of the 43 Factor-1 top-positive genes that are targets in the r > 0.7 scMEGA
eGRN: number of those genes each TF regulates (bars coloured by the pseudotime bin of peak motif activity), and all
targets of the TF with the Factor-1 share and the hypergeometric over-representation (BH-FDR < 0.05, asterisks).
Panel        : All Supp Figures.pdf, Supplementary Figure 3, panel b
Reads        : $KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/ATAC/scmega/figures_paper/egrn_paper_views_edges.csv
                 (analyses/10_atac_pt_apt_egrn/ATAC/scripts/45_egrn_paper_views.py)
               $KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/ATAC/scmega/egrn_tutorial/archr_dm/{tf_centrality.csv,
                 grn_edges_all.csv, peak_to_gene_links.csv} (analyses/10_atac_pt_apt_egrn/ATAC/scripts/
                 43_scmega_egrn_tutorial.R, TRAJ=archr_dm)
               $KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/cis_regulatory_pt_apt/outputs/task3/
                 C4_factor1_top200_positive_genes.csv (columns variable, Loading; rebuilt by
                 analyses/10_atac_pt_apt_egrn/make_factor1_pt_list.py)
Writes       : $KB_FIG_OUT/S3/S3b_egrn_factor1_tf_ranking/egrn_factor1_tf_ranking.{pdf,png} (the panel) and the other
               outputs of the script (TF / gene tables, sub-network, specificity and adjacency views)
Plotting code: analyses/10_atac_pt_apt_egrn/ATAC/scripts/47_egrn_factor1_tf_side.py::main (imported unchanged; run with
               the published options --highlight RELA,REL,NFKB1,NFKB2,RELB,FOS,JUN,JUNB,JDP2,FOSL2,HNF4A,RXRA,RXRB,
               RXRG,PPARA,PPARG as in job 339554 of 26 Oct 2026; module-level OUT redirected to $KB_FIG_OUT)
Illustrator  : titles and x-axis titles re-typed ("TFs upstream of the 43 factor 1-genes", "number of 43 F1 target genes
               the TF regulates", "all targets vs factor 1 share"); figure-level legend at the bottom removed (dotted
               line and asterisks are explained in the figure legend); colour bar moved into the bar panel and
               re-labelled; panel letter. NOTE: the removed legend (and 47's meta json) say "BH-FDR < 0.05 over 86 TFs"
               and "population = 388 targets, 46 Factor-1 targets": stale text from the September run. The test itself
               uses the rerun network (80 TFs, 337 targets, 43 Factor-1 targets), so the asterisks are right
Environment  : conda env agent_env
Verified     : 2026-11-12, PIXEL-IDENTICAL (PNG byte-identical; TF / gene tables identical) vs
               $KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/ATAC/scmega/figures_paper/egrn_factor1_tf_ranking.pdf
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from atac_egrn import run_main  # noqa: E402  (also puts config/ on sys.path)
from kbpaths import fig_out  # noqa: E402

HIGHLIGHT = "RELA,REL,NFKB1,NFKB2,RELB,FOS,JUN,JUNB,JDP2,FOSL2,HNF4A,RXRA,RXRB,RXRG,PPARA,PPARG"
OUT = fig_out("S3", "S3b_egrn_factor1_tf_ranking")
run_main("47_egrn_factor1_tf_side.py", ["--highlight", HIGHLIGHT], OUT)
print(f"[S3b] panel: {OUT}/egrn_factor1_tf_ranking.pdf", flush=True)

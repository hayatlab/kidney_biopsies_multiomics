#!/usr/bin/env python
"""Fig. 6k — volcano plot of the snRNA-seq MAST comparison aPEC (PEC subcluster 9) vs PEC subclusters 0, 2, 4 and 11,
after removal of canonical non-PEC lineage markers; top 10 up- and top 10 down-regulated genes by log2FC labelled.
Panel        : Main Figures.pdf, Figure 6, panel k
Reads        : $KB_ANALYSIS/pec_resubcluster_20260717/de_9vs0241/de_results/9_vs_0_2_4_11/DE/markers_MAST_9_vs_0_2_4_11.csv
                 (analyses/09_snrna_pec_immune_subclusters/pec/de_9vs0241/de_analysis_9vs0241.R, wrapper run_de_9vs0241.sh;
                  input built by pec/prep_de_9vs0241.py)
               analyses/09_snrna_pec_immune_subclusters/subcluster_sweep/subcluster_sweep.py (PANELS literal = the
                 canonical lineage marker panels whose non-PEC genes are removed)
Writes       : $KB_FIG_OUT/Fig6/Fig6k_volcano_aPEC_vs_PEC.{pdf,png}, Fig6k_volcano_aPEC_vs_PEC_labeledgenes.csv
Plotting code: analyses/09_snrna_pec_immune_subclusters/pec/de_9vs0241/filter_and_volcano_pecfiltered.py::main (imported
               unchanged; its module-level output folder DE_DIR is pointed at a scratch folder under $KB_FIG_OUT and
               N_TOP_LIST at [10]; main() re-applies the deterministic non-PEC marker filter, draws its volcanos there and
               the top10/bottom10 volcano is kept as the panel; the scratch folder is removed)
Illustrator  : title replaced ("aPEC (subcluster 9) vs PEC subclusters 0, 2, 4 and 11 / top 10 up + 10 down genes by
               log2FC (padj<0.05) (snRNA-seq)"); legend reduced to the UP/DOWN entries and reworded ("UP in cl9
               (2134)" -> "UP in cluster 9 (2134 genes)"); x label "(cl9 vs cl0/2/4/11)" suffix removed
Environment  : conda env nichecompass
Verified     : 2026-11-12, PIXEL-IDENTICAL (PDF, text layer equal; PNG and labelled-gene CSV byte-identical) vs
               $KB_ANALYSIS/pec_resubcluster_20260717/de_9vs0241/de_results/9_vs_0_2_4_11/DE/
               volcano_top10bottom10_9_vs_0_2_4_11_pecfiltered.pdf
"""
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from niches_snrna import ANALYSES_09, load_script  # noqa: E402
from kbpaths import fig_out  # noqa: E402

vol = load_script(ANALYSES_09 / "pec" / "de_9vs0241" / "filter_and_volcano_pecfiltered.py")
N = 10
STEM = f"volcano_top{N}bottom{N}_9_vs_0_2_4_11_pecfiltered"


def main():
    fig_dir = Path(fig_out("Fig6", "Fig6k_volcano_aPEC_vs_PEC.pdf")).parent
    work = Path(tempfile.mkdtemp(prefix="_work_Fig6k_", dir=fig_dir))
    try:
        vol.DE_DIR = work            # outputs only; the input table path (vol.CSV) was bound at import
        vol.N_TOP_LIST = [N]
        vol.main()
        for ext in ("pdf", "png"):
            shutil.move(str(work / f"{STEM}.{ext}"), fig_out("Fig6", f"Fig6k_volcano_aPEC_vs_PEC.{ext}"))
        shutil.move(str(work / f"volcano_top{N}bottom{N}_labeledgenes.csv"),
                    fig_out("Fig6", "Fig6k_volcano_aPEC_vs_PEC_labeledgenes.csv"))
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print("[Fig6k] wrote", fig_out("Fig6", "Fig6k_volcano_aPEC_vs_PEC.pdf"))


if __name__ == "__main__":
    main()

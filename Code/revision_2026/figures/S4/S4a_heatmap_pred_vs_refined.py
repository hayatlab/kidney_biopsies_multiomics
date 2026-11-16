#!/usr/bin/env python
"""Supp. Fig. 4a — concordance of the per-cell scANVI label transfer from snRNA-seq (pred_level2, rows) with the final
Xenium annotation (level2_refined, columns), row-normalised: fraction of the cells of each transferred label that end
up in each final label.
Panel        : All Supp Figures.pdf, Supplementary Fig. 4, panel a
Reads        : $KB_ANALYSIS/celltype_concordance/outputs/xenium_labeled.h5ad, obs (analyses/02_xenium_celltyping/:
                 pred_level2 from 02_train_transfer.py, level2_refined from 09_level2_subcluster.py and
                 17_recover_pec_from_apt.py)
Writes       : $KB_FIG_OUT/S4/S4a_heatmap_pred_level2_vs_refined.pdf (+ .png, S4a_pred_level2_vs_refined.csv)
Plotting code: analyses/02_xenium_celltyping/10_heatmap_pred_vs_refined.py::main (imported unchanged; its output
               folders config.FIG / config.MET are redirected and the files renamed)
               (original run: python 10_heatmap_pred_vs_refined.py, 2026-07-01, no log kept)
Illustrator  : axis labels replaced ("cell types snRNAseq" for the rows, "cell types Xenium" for the columns); dPT
               renamed aPT2
Environment  : conda env kpmp_repurposing
Verified     : 2026-11-12, DIFFERENT (1.09% of pixels) vs $KB_ANALYSIS/celltype_concordance/figures/
               heatmap_pred_level2_vs_refined.pdf on the final object (post PEC recovery: 6 matrix cells change, max 0.018;
               4 printed values); PIXEL-IDENTICAL with --xenium $KB_ANALYSIS/celltype_concordance/outputs/
               xenium_labeled.PRE_METADATA_UPDATE_20260702.h5ad (crosstab CSV byte-identical to outputs/metrics/pred_level2_vs_refined.csv)

Data version: the published panel was drawn on 2026-07-01, before 17_recover_pec_from_apt.py (2026-07-09) moved 182
cells from aPT to PEC (aPT 15,072 -> 14,890; PEC 2,395 -> 2,577). On the final object this script draws the
post-recovery labels; `--xenium <archived copy of 2026-07-01>` redraws the published version.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from xenium_figs import load_script, rename  # noqa: E402
from kbpaths import fig_out  # noqa: E402

m10 = load_script("02_xenium_celltyping/10_heatmap_pred_vs_refined.py", on_path=True)


def main():
    ap = argparse.ArgumentParser(description="Supp. Fig. 4a")
    ap.add_argument("--xenium", default=None, help="read this object instead of outputs/xenium_labeled.h5ad")
    ap.add_argument("--tag", default="", help="suffix for the output file names")
    args = ap.parse_args()
    C = m10.C
    outdir = Path(fig_out("S4", "S4a_heatmap_pred_level2_vs_refined.pdf")).parent
    C.FIG = C.MET = outdir
    C.ensure_dirs = lambda: None
    if args.xenium:
        C.XENIUM_LABELED = Path(args.xenium)
    print(f"[S4a] object: {C.XENIUM_LABELED}", flush=True)
    m10.main()
    t = args.tag
    for f in rename(outdir, {"heatmap_pred_level2_vs_refined.pdf": f"S4a_heatmap_pred_level2_vs_refined{t}.pdf",
                             "heatmap_pred_level2_vs_refined.png": f"S4a_heatmap_pred_level2_vs_refined{t}.png",
                             "pred_level2_vs_refined.csv": f"S4a_pred_level2_vs_refined{t}.csv"}):
        print(f"[S4a] wrote {f}", flush=True)


if __name__ == "__main__":
    main()

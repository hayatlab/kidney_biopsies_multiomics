#!/usr/bin/env python
"""Supp. Fig. 9c-d — c, dot plot of the top-3 marker genes per Xenium cell type (level2_refined; scanpy Wilcoxon
ranking on log1p(CP10k) of the raw counts; dot size = fraction of cells expressing; colour = log fold change, scaled
per gene); d, number of Xenium cells per cell type and condition (bubble area ~ n cells) with per-type and
per-condition totals.
Panel        : All Supp Figures.pdf, Supplementary Fig. 9, panels c and d
Reads        : $KB_ANALYSIS/celltype_concordance/outputs/xenium_labeled.h5ad, layers["counts"] + obs (level2_refined,
                 condition) (analyses/02_xenium_celltyping/)
Writes       : $KB_FIG_OUT/S9/S9c_dotplot_top3_markers_level2_refined.pdf, $KB_FIG_OUT/S9/S9d_bubble_level2_by_condition.pdf
               (+ .png)
Plotting code: analyses/02_xenium_celltyping/13_dotplot_bubble.py::main (imported unchanged; its output folder
               config.FIG is redirected and the files renamed). The marker ranking (sc.tl.rank_genes_groups,
               Wilcoxon) is computed on the fly, as in the original (a few minutes on 258,375 cells).
               (original run: python 13_dotplot_bubble.py, 2026-07-01, no log kept)
Illustrator  : c: dPT renamed aPT2; d: title and axis titles removed, "CTRL" written "CTRL_Xenium", dPT renamed aPT2
Environment  : conda env kpmp_repurposing
Verified     : 2026-11-12, on the final object DIFFERENT vs $KB_ANALYSIS/celltype_concordance/figures/
               dotplot_top3_markers_level2_refined.pdf (0.47% of pixels, aPT/PEC rows; same genes) and bubble_level2_by_condition.pdf
               (0.13%; aPT 14,890 / PEC 2,577 instead of 15,072 / 2,395); both PIXEL-IDENTICAL with --xenium $KB_ANALYSIS/
               celltype_concordance/outputs/xenium_labeled.PRE_METADATA_UPDATE_20260702.h5ad

Data version: drawn on 2026-07-01, before 17_recover_pec_from_apt.py moved 182 cells from aPT to PEC (the published
bubble plot shows aPT 15,072 and PEC 2,395; the final object gives 14,890 and 2,577). `--xenium <archived copy of
2026-07-01>` redraws the published version. Colour in c is the log fold change, not the scaled mean expression.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from xenium_figs import load_script, rename  # noqa: E402
from kbpaths import fig_out  # noqa: E402

m13 = load_script("02_xenium_celltyping/13_dotplot_bubble.py", on_path=True)


def main():
    ap = argparse.ArgumentParser(description="Supp. Fig. 9c-d")
    ap.add_argument("--xenium", default=None, help="read this object instead of outputs/xenium_labeled.h5ad")
    ap.add_argument("--tag", default="", help="suffix for the output file names")
    args = ap.parse_args()
    C = m13.C
    outdir = Path(fig_out("S9", "S9c_dotplot_top3_markers_level2_refined.pdf")).parent
    C.FIG = outdir
    if args.xenium:
        C.XENIUM_LABELED = Path(args.xenium)
    print(f"[S9c-d] object: {C.XENIUM_LABELED}", flush=True)
    m13.main()
    t = args.tag
    names = {"dotplot_top3_markers_level2_refined": "S9c", "bubble_level2_by_condition": "S9d"}
    for f in rename(outdir, {f"{s}.{e}": f"{p}_{s}{t}.{e}" for s, p in names.items() for e in ("pdf", "png")}):
        print(f"[S9c-d] wrote {f}", flush=True)


if __name__ == "__main__":
    main()

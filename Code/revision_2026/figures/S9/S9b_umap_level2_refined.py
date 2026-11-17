#!/usr/bin/env python
"""Supp. Fig. 9b — UMAP of all 258,375 Xenium cells (X_umap stored in the object: NicheCompass latent space)
coloured by the final annotation level2_refined, labels on the data, tab20 palette as in the exported label legend.
Panel        : All Supp Figures.pdf, Supplementary Fig. 9, panel b
Reads        : $KB_ANALYSIS/celltype_concordance/outputs/xenium_labeled.h5ad, obs + obsm["X_umap"]
                 (analyses/02_xenium_celltyping/; X_umap carried over from the NicheCompass object,
                  analyses/01_xenium_preprocessing_niches/nichecompass_overall_with_glom.py)
Writes       : $KB_FIG_OUT/S9/S9b_umap_level2_refined_ondata.pdf (+ .png; and the right-margin-legend version
               S9b_umap_level2_refined_legend.pdf/.png, not used)
Plotting code: analyses/02_xenium_celltyping/12_umap_level2_refined.py::main (imported unchanged; its output folder
               config.FIG is redirected and the files renamed)
               (original run: python 12_umap_level2_refined.py, 2026-07-01, no log kept)
Illustrator  : title "Xenium — level2_refined" replaced by "N=258,375 cells" (dPT is NOT renamed in this panel)
Revision     : 2026-11-13, dPT is printed as aPT2 (co-authors' rename; same colour). --original-labels draws "dPT"
Environment  : conda env kpmp_repurposing
Verified     : 2026-11-12, DIFFERENT (0.29% of pixels; 182 cells recoloured aPT -> PEC) vs $KB_ANALYSIS/celltype_concordance/
               figures/umap_level2_refined_ondata.pdf on the final object; PIXEL-IDENTICAL with --xenium $KB_ANALYSIS/
               celltype_concordance/outputs/xenium_labeled.PRE_METADATA_UPDATE_20260702.h5ad

Data version: drawn on 2026-07-01, before 17_recover_pec_from_apt.py moved 182 cells from aPT to PEC; see
S9a_heatmap_pred_vs_refined.py. `--xenium <archived copy of 2026-07-01>` redraws the published version.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from xenium_figs import load_script, rename  # noqa: E402
from kbpaths import fig_out  # noqa: E402

m12 = load_script("02_xenium_celltyping/12_umap_level2_refined.py", on_path=True)


def main():
    ap = argparse.ArgumentParser(description="Supp. Fig. 9b")
    ap.add_argument("--xenium", default=None, help="read this object instead of outputs/xenium_labeled.h5ad")
    ap.add_argument("--tag", default="", help="suffix for the output file names")
    ap.add_argument("--original-labels", action="store_true", help='print "dPT" as in the object (published version)')
    args = ap.parse_args()
    C = m12.C
    outdir = Path(fig_out("S9", "S9b_umap_level2_refined_ondata.pdf")).parent
    C.FIG = outdir
    if args.xenium:
        C.XENIUM_LABELED = Path(args.xenium)
    m12.DISPLAY = {} if args.original_labels else {"dPT": "aPT2"}
    print(f"[S9b] object: {C.XENIUM_LABELED}; display names {m12.DISPLAY}", flush=True)
    m12.main()
    t = args.tag
    for f in rename(outdir, {f"umap_level2_refined_{v}.{e}": f"S9b_umap_level2_refined_{v}{t}.{e}"
                             for v in ("ondata", "legend") for e in ("pdf", "png")}):
        print(f"[S9b] wrote {f}", flush=True)


if __name__ == "__main__":
    main()

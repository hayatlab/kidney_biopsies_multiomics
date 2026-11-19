#!/usr/bin/env python
"""Supp. Fig. 22a-f — Venn diagrams of glomerular niche-20 PEC DE genes vs the top-200 MCFA factor 1 loadings of the PT and PEC views (restricted universe).
Panel        : All Supp Figures.pdf, Supplementary Figure 22, panels a-f (a PT UP, b PEC UP, c 3-way UP, d PT DOWN, e PEC DOWN, f 3-way DOWN)
Reads        : $KB_ANALYSIS/niche20_pec/DE/niche20_vs_others_pecOnly/markers_MAST_niche20_vs_others_pecOnly.csv (analyses/08_glomerular_niche20_crescent/niche20_volcano_gsea.R: copy of markers_manualglom10.csv from de_crescent_pec_mast.R)
               $KB_ANALYSIS/loadings~hvg~5.csv (MCFA factor loadings; producer outside this folder, see analyses/08_glomerular_niche20_crescent/README.md)
Writes       : $KB_FIG_OUT/S22/S22{a..f}_venn_{PT_UP,PEC_UP,3way_UP,PT_DOWN,PEC_DOWN,3way_DOWN}.pdf
Plotting code: analyses/08_glomerular_niche20_crescent/venn_niche20_mofa_factor1.py::load_de, load_factor1, draw_venn2, draw_venn3
               (imported unchanged); the set logic of main(), branch main_restricted_universe/top200, single-panel files,
               is extracted verbatim (the original loop also writes the other three variants and the tables).
Illustrator  : "PT" expanded to "PT (PT, aPT, aPT2)" in titles and set labels; typography restyled (bold titles,
               regular-weight counts)
Environment  : conda env decoupler
Verified     : 2026-11-12, PIXEL-IDENTICAL (all six) vs $KB_ANALYSIS/niche20_pec/venn_mofa_factor1/main_restricted_universe/top200/venn_*.pdf
"""
import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "config"))
from kbpaths import fig_out  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# importing the module sets the Agg backend, its rcParams (pdf.fonttype 42, savefig.bbox tight) and the
# equal-size Venn layouts
v = _load(REPO / "analyses" / "08_glomerular_niche20_crescent" / "venn_niche20_mofa_factor1.py", "venn_niche20_mofa_factor1")

PANEL = {("PT", "UP"): "S22a", ("PEC", "UP"): "S22b", ("3way", "UP"): "S22c",
         ("PT", "DOWN"): "S22d", ("PEC", "DOWN"): "S22e", ("3way", "DOWN"): "S22f"}


def main() -> None:
    up, dn, panel = v.load_de()
    F = v.load_factor1()
    key_pos, key_neg = "pos_top", "neg_top"                          # loading cut "top200"
    U_PT = panel & F["PT"]["all_genes"]                              # restricted universes
    U_PEC = panel & F["PEC"]["all_genes"]
    U_3 = panel & F["PT"]["all_genes"] & F["PEC"]["all_genes"]

    def R(s, U):
        return s & U

    for direction, de_set, key in [("UP", up, key_pos), ("DOWN", dn, key_neg)]:
        sign_word = "positive" if direction == "UP" else "negative"
        for view, U, colour in [("PT", U_PT, v.COL_PT), ("PEC", U_PEC, v.COL_PEC)]:
            A = R(F[view][key], U)
            B = R(de_set, U)
            lab_a = f"{view} F1 {sign_word}\n(n={len(A)})"
            lab_b = f"niche20 {direction}\n(n={len(B)})"
            f1_, a1 = plt.subplots(figsize=(4.6, 3.9))
            v.draw_venn2(a1, A, B, lab_a, lab_b, colour, v.COL_N20,
                         f"{view} Factor1 {sign_word} vs niche20 {direction}")
            out = fig_out("S22", f"{PANEL[(view, direction)]}_venn_{view}_{direction}.pdf")
            f1_.savefig(out)
            plt.close(f1_)
            print(f"[S22] {out}: |A|={len(A)} |B|={len(B)} overlap={len(A & B)}")

        A = R(F["PT"][key], U_3)
        B = R(F["PEC"][key], U_3)
        C = R(de_set, U_3)
        labs = (f"PT F1 {sign_word}\n(n={len(A)})",
                f"PEC F1 {sign_word}\n(n={len(B)})",
                f"niche20 {direction}\n(n={len(C)})")
        f3_, a3 = plt.subplots(figsize=(5.0, 4.7))
        v.draw_venn3(a3, A, B, C, labs, f"PT vs PEC Factor1 {sign_word} vs niche20 {direction}")
        out = fig_out("S22", f"{PANEL[('3way', direction)]}_venn_3way_{direction}.pdf")
        f3_.savefig(out)
        plt.close(f3_)
        print(f"[S22] {out}: |PT|={len(A)} |PEC|={len(B)} |niche20|={len(C)} all three={len(A & B & C)}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""
Venn diagrams: niche-20 (crescent PEC) DE genes  vs  MOFA Factor 1 loadings (PT and PEC views).

Species human (GRCh38), HGNC symbols throughout. Deterministic (no sampling / no seed needed).

Sets
----
niche20 UP / DOWN : MAST crescent-PEC vs other-PEC, |avg_log2FC| > 0.25 & p_val_adj < 0.05
                    (positive avg_log2FC = up in the crescent niche)
PT / PEC Factor1  : loadings~hvg~5.csv, Factor == "Factor1", view in {PT, PEC}
                    sign-matched to the DE direction (UP <-> positive loading, DOWN <-> negative)

Two loading cuts (both produced, per user request):
  top200     : the 200 largest-|loading| genes of that view/factor, then split by sign
  all_signed : every gene of that view/factor with loading > 0 (resp. < 0)

Two gene universes (both produced):
  main_restricted_universe : intersect every set with (468-gene Xenium panel) n (genes of that MOFA
                             view). A gene outside a circle then means "tested in both, not a member",
                             so the hypergeometric enrichment test is valid.
  supp_unrestricted        : raw sets, no intersection. Overlap counts only - the hypergeometric
                             p-value is reported for reference but is NOT valid, because the two sets
                             are drawn from different measurable gene spaces.

Six Venns per variant: PT-UP, PT-DOWN, PEC-UP, PEC-DOWN, 3-way-UP, 3-way-DOWN.
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)

import os
import itertools
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from matplotlib_venn.layout.venn2 import DefaultLayoutAlgorithm as Venn2Layout
from matplotlib_venn.layout.venn3 import DefaultLayoutAlgorithm as Venn3Layout
from scipy.stats import hypergeom

# circles of equal size: the sets differ by up to ~25x (123 vs 2904) and an area-proportional
# layout collapses the small circle to a dot. Counts are printed inside every region.
UNWEIGHTED2 = Venn2Layout(fixed_subset_sizes=(1, 1, 1))
UNWEIGHTED3 = Venn3Layout(fixed_subset_sizes=(1, 1, 1, 1, 1, 1, 1))

BASE = f"{KB_ANALYSIS}"
DE_CSV = os.path.join(
    BASE, "niche20_pec/DE/niche20_vs_others_pecOnly/markers_MAST_niche20_vs_others_pecOnly.csv"
)
LOADINGS_CSV = os.path.join(BASE, "loadings~hvg~5.csv")
OUTDIR = os.path.join(BASE, "niche20_pec/venn_mofa_factor1")

LFC_CUT = 0.25
PADJ_CUT = 0.05
FACTOR = "Factor1"
TOPN = 200

# colour-blind safe (Okabe-Ito)
COL_PT = "#0072B2"      # blue   - PT Factor1
COL_PEC = "#E69F00"     # orange - PEC Factor1
COL_N20 = "#009E73"     # green  - niche 20

plt.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 9, "axes.titlesize": 10,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


# ----------------------------------------------------------------------------- input sets
def load_de():
    de = pd.read_csv(DE_CSV)
    panel = set(de["gene"])
    sig = de[(de["p_val_adj"] < PADJ_CUT) & (de["avg_log2FC"].abs() > LFC_CUT)]
    up = set(sig.loc[sig["avg_log2FC"] > 0, "gene"])
    dn = set(sig.loc[sig["avg_log2FC"] < 0, "gene"])
    return up, dn, panel


def load_factor1():
    L = pd.read_csv(LOADINGS_CSV)
    f1 = L[L["Factor"] == FACTOR]
    out = {}
    for view in ("PT", "PEC"):
        s = f1[f1["view"] == view]
        assert s["variable"].is_unique, f"duplicate genes in {view}/{FACTOR}"
        out[view] = {
            "all_genes": set(s["variable"]),
            "pos_all": set(s.loc[s["Loading"] > 0, "variable"]),
            "neg_all": set(s.loc[s["Loading"] < 0, "variable"]),
        }
        top = s.reindex(s["Loading"].abs().sort_values(ascending=False).index).head(TOPN)
        out[view]["pos_top"] = set(top.loc[top["Loading"] > 0, "variable"])
        out[view]["neg_top"] = set(top.loc[top["Loading"] < 0, "variable"])
        out[view]["top_min_abs"] = float(top["Loading"].abs().min())
    return out


# ----------------------------------------------------------------------------- stats
def hypergeom_pair(a, b, M):
    """One-sided (over-representation) hypergeometric test for |a n b| given universe size M."""
    k = len(a & b)
    n, K = len(a), len(b)
    if M == 0 or n == 0 or K == 0:
        return k, np.nan, np.nan
    expected = n * K / M
    # P(X >= k)
    p = hypergeom.sf(k - 1, M, K, n)
    return k, expected, p


def triple_expected(a, b, c, M):
    if M == 0:
        return np.nan
    return M * (len(a) / M) * (len(b) / M) * (len(c) / M)


# ----------------------------------------------------------------------------- plotting
def draw_venn2(ax, A, B, labA, labB, colA, colB, title):
    v = venn2([A, B], set_labels=(labA, labB), set_colors=(colA, colB), alpha=0.55, ax=ax,
              layout_algorithm=UNWEIGHTED2)
    for t in v.set_labels:
        if t:
            t.set_fontsize(9)
    for t in v.subset_labels:
        if t:
            t.set_fontsize(9)
            t.set_fontweight("bold")
    ax.set_title(title, fontsize=10, pad=6)


def draw_venn3(ax, A, B, C, labs, title):
    v = venn3([A, B, C], set_labels=labs, set_colors=(COL_PT, COL_PEC, COL_N20),
              alpha=0.55, ax=ax, layout_algorithm=UNWEIGHTED3)
    for t in v.set_labels:
        if t:
            t.set_fontsize(9)
    for t in v.subset_labels:
        if t:
            t.set_fontsize(8.5)
            t.set_fontweight("bold")
    ax.set_title(title, fontsize=10, pad=13)


# ----------------------------------------------------------------------------- membership export
def membership_rows(sets_named, variant, panel_label):
    """Long table: one row per gene per Venn, with which circles it falls in."""
    rows = []
    for venn_name, d in sets_named.items():
        names = list(d.keys())
        allg = set().union(*d.values())
        for g in sorted(allg):
            memb = [n for n in names if g in d[n]]
            rows.append({
                "variant": variant, "venn": venn_name, "gene": g,
                "n_sets": len(memb), "sets": "|".join(memb),
                "universe": panel_label,
            })
    return rows


def main():
    up, dn, panel = load_de()
    F = load_factor1()

    print(f"[input] niche20 panel={len(panel)} UP={len(up)} DOWN={len(dn)} "
          f"(|lfc|>{LFC_CUT}, padj<{PADJ_CUT})")
    for v in ("PT", "PEC"):
        print(f"[input] {FACTOR}/{v}: {len(F[v]['all_genes'])} genes | "
              f"all pos {len(F[v]['pos_all'])} neg {len(F[v]['neg_all'])} | "
              f"top{TOPN} pos {len(F[v]['pos_top'])} neg {len(F[v]['neg_top'])} "
              f"(min |loading| {F[v]['top_min_abs']:.3f})")

    stats_rows, memb_rows = [], []

    for universe_tag, restrict in [("main_restricted_universe", True),
                                   ("supp_unrestricted", False)]:
        for cut_tag, key_pos, key_neg in [("top200", "pos_top", "neg_top"),
                                          ("all_signed", "pos_all", "neg_all")]:
            variant = f"{universe_tag}/{cut_tag}"
            odir = os.path.join(OUTDIR, universe_tag, cut_tag)
            os.makedirs(odir, exist_ok=True)

            # ---- universes
            U_PT = (panel & F["PT"]["all_genes"]) if restrict else None
            U_PEC = (panel & F["PEC"]["all_genes"]) if restrict else None
            U_3 = (panel & F["PT"]["all_genes"] & F["PEC"]["all_genes"]) if restrict else None

            def R(s, U):
                return (s & U) if U is not None else s

            # M for the hypergeometric test
            M_PT = len(U_PT) if restrict else len(panel | F["PT"]["all_genes"])
            M_PEC = len(U_PEC) if restrict else len(panel | F["PEC"]["all_genes"])
            M_3 = len(U_3) if restrict else len(panel | F["PT"]["all_genes"] | F["PEC"]["all_genes"])

            fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.2))
            sets_named = {}

            for row, (direction, de_set, key) in enumerate(
                    [("UP", up, key_pos), ("DOWN", dn, key_neg)]):
                sign_word = "positive" if direction == "UP" else "negative"

                # ---------- panels 1 & 2 : PT and PEC, two-way
                for col, (view, U, M, colour) in enumerate(
                        [("PT", U_PT, M_PT, COL_PT), ("PEC", U_PEC, M_PEC, COL_PEC)]):
                    A = R(F[view][key], U)
                    B = R(de_set, U)
                    k, exp, p = hypergeom_pair(A, B, M)
                    jac = k / len(A | B) if len(A | B) else np.nan
                    lab_a = f"{view} F1 {sign_word}\n(n={len(A)})"
                    lab_b = f"niche20 {direction}\n(n={len(B)})"
                    draw_venn2(axes[row, col], A, B, lab_a, lab_b, colour, COL_N20,
                               f"{view} Factor1 {sign_word} vs niche20 {direction}")

                    name = f"{view}_{direction}"
                    sets_named[name] = {f"{view}_F1_{sign_word}": A,
                                        f"niche20_{direction}": B}
                    stats_rows.append(dict(
                        variant=variant, universe=universe_tag, loading_cut=cut_tag,
                        venn=name, type="2-way", direction=direction, view=view,
                        n_loading_set=len(A), n_niche20=len(B), overlap=k,
                        expected=exp, fold=k / exp if exp else np.nan,
                        jaccard=jac, universe_M=M, hypergeom_p=p,
                        p_valid=bool(restrict)))

                    # single-panel file
                    f1_, a1 = plt.subplots(figsize=(4.6, 3.9))
                    draw_venn2(a1, A, B, lab_a, lab_b, colour, COL_N20,
                               f"{view} Factor1 {sign_word} vs niche20 {direction}")
                    f1_.savefig(os.path.join(odir, f"venn_{view}_{direction}.pdf"))
                    f1_.savefig(os.path.join(odir, f"venn_{view}_{direction}.png"), dpi=300)
                    plt.close(f1_)

                # ---------- panel 3 : three-way
                A = R(F["PT"][key], U_3)
                B = R(F["PEC"][key], U_3)
                C = R(de_set, U_3)
                labs = (f"PT F1 {sign_word}\n(n={len(A)})",
                        f"PEC F1 {sign_word}\n(n={len(B)})",
                        f"niche20 {direction}\n(n={len(C)})")
                triple = len(A & B & C)
                exp3 = triple_expected(A, B, C, M_3)
                # PT and PEC are two views of the same factor and are NOT independent, so exp3
                # (mutual independence) understates chance. Condition on PT n PEC instead:
                # "given the genes both views load, is niche20 enriched among them?"
                n_ptpec = len(A & B)
                exp3_cond = n_ptpec * len(C) / M_3 if M_3 else np.nan
                p3_cond = (hypergeom.sf(triple - 1, M_3, n_ptpec, len(C))
                           if M_3 and n_ptpec and len(C) else np.nan)
                draw_venn3(axes[row, 2], A, B, C, labs,
                           f"PT vs PEC Factor1 {sign_word} vs niche20 {direction}")

                sets_named[f"3way_{direction}"] = {
                    f"PT_F1_{sign_word}": A, f"PEC_F1_{sign_word}": B,
                    f"niche20_{direction}": C}
                stats_rows.append(dict(
                    variant=variant, universe=universe_tag, loading_cut=cut_tag,
                    venn=f"3way_{direction}", type="3-way", direction=direction, view="PT+PEC",
                    n_loading_set=f"PT={len(A)};PEC={len(B)}", n_niche20=len(C),
                    overlap=triple, expected=exp3_cond,
                    fold=triple / exp3_cond if exp3_cond else np.nan,
                    jaccard=triple / len(A | B | C) if len(A | B | C) else np.nan,
                    universe_M=M_3, hypergeom_p=p3_cond, p_valid=bool(restrict),
                    n_PT_and_PEC=n_ptpec,
                    expected_indep_null=exp3,
                    fold_indep_null=triple / exp3 if exp3 else np.nan))

                f3_, a3 = plt.subplots(figsize=(5.0, 4.7))
                draw_venn3(a3, A, B, C, labs,
                           f"PT vs PEC Factor1 {sign_word} vs niche20 {direction}")
                f3_.savefig(os.path.join(odir, f"venn_3way_{direction}.pdf"))
                f3_.savefig(os.path.join(odir, f"venn_3way_{direction}.png"), dpi=300)
                plt.close(f3_)

            head = ("restricted universe (Xenium panel ∩ MOFA view)" if restrict
                    else "unrestricted (raw sets; panel / HVG spaces differ)")
            fig.suptitle(f"niche 20 (crescent PEC) DE vs MOFA {FACTOR} loadings — "
                         f"{cut_tag} · {head}\ncircles not area-proportional",
                         fontsize=11.5, y=0.995)
            fig.tight_layout(rect=[0, 0, 1, 0.945], h_pad=1.5, w_pad=1.0)
            fig.savefig(os.path.join(odir, "figure_all_venns.pdf"))
            fig.savefig(os.path.join(odir, "figure_all_venns.png"), dpi=300)
            plt.close(fig)

            memb_rows += membership_rows(sets_named, variant, universe_tag)
            pd.DataFrame(membership_rows(sets_named, variant, universe_tag)).to_csv(
                os.path.join(odir, "gene_set_membership.csv"), index=False)
            print(f"[done] {variant} -> {odir}")

    os.makedirs(OUTDIR, exist_ok=True)
    st = pd.DataFrame(stats_rows)
    st.to_csv(os.path.join(OUTDIR, "overlap_stats_all_variants.tsv"), sep="\t", index=False)
    pd.DataFrame(memb_rows).to_csv(
        os.path.join(OUTDIR, "gene_set_membership_all_variants.csv"), index=False)

    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print("\n=== overlap summary ===")
        print(st[["variant", "venn", "n_loading_set", "n_niche20", "overlap",
                  "expected", "fold", "hypergeom_p", "p_valid"]].to_string(index=False))


if __name__ == "__main__":
    main()

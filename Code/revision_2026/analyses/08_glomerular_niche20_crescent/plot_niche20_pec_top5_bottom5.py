"""Top-5-up / top-5-down marker dotplot for niche-20 (manual_glom) PEC vs other-niche PEC.

Mirrors the niche-9 figure (`niche9_marker_plots_extra.py`, plot A) but points at the
manual_glom NicheCompass run where the crescent = niche 20.

Data flow (kept identical to the DE so cells + labels match exactly):
  - Expression  : celltype_concordance/outputs/crescent_check/de/pec_manualglom10.h5ad
                  (PEC-only, X = raw integer counts, obs['grp'] in {'crescent','other'};
                   'crescent' = niche-20 PEC, 'other' = all other-niche PEC).
  - DE ranking  : celltype_concordance/outputs/crescent_check/de/markers_manualglom10.csv
                  (MAST p-value, Seurat linear-space avg_log2FC; ident.1='crescent' vs 'other').

Two figures are written:
  * v1 (audit)   : padj-ranked top-5 up / top-5 down, standard_scale='var'  -- the faithful
                   niche-9-recipe reproduction. Includes AFDN in the down set.
                   -> dotplot_niche20_vs_other_pec_top5_bottom5_v1_padjrank_varscale.*
  * final (canon): DIRECTION-CONSISTENT selection (keep only genes whose mean log-norm
                   direction agrees with the Seurat avg_log2FC sign) + standard_scale=None
                   (colour = raw mean log-norm expression on a shared scale).
                   -> dotplot_niche20_vs_other_pec_top5_bottom5.*

Why the two differ: Seurat's avg_log2FC is a LINEAR-space statistic (log2 of the ratio of
mean expm1(lognorm)), computed independently of the MAST p-value, whereas the dotplot colour
and MAST's own hurdle model live in LOG space. For weak-effect genes the two can disagree in
sign. AFDN is the only panel gene where they do (Seurat log2FC=-0.72 "down", but mean log-norm
is +0.21 HIGHER in niche-20). The direction-consistency filter drops AFDN and pulls in NDRG1
(= niche-9's 5th down gene). standard_scale=None additionally avoids the per-gene 0..1 rescale
that stretches a trivial log-mean difference into full colour contrast.
"""
from __future__ import annotations
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from matplotlib import font_manager as _fm

for _p in (
    "/usr/share/fonts/liberation-sans/LiberationSans-Regular.ttf",
    "/usr/share/fonts/liberation-sans/LiberationSans-Bold.ttf",
    "/usr/share/fonts/liberation-sans/LiberationSans-Italic.ttf",
    "/usr/share/fonts/liberation-sans/LiberationSans-BoldItalic.ttf",
):
    try:
        _fm.fontManager.addfont(_p)
    except Exception:
        pass

mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "DejaVu Sans"]
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"

BASE = Path(f"{KB_ANALYSIS}/celltype_concordance/outputs/crescent_check/de")
ADATA_PATH = BASE / "pec_manualglom10.h5ad"
DE_CSV = BASE / "markers_manualglom10.csv"
OUT_DIR = BASE

TARGET = "niche-20 PEC"
OTHER = "other-niche PEC"
N = 5
PADJ = 0.05
EXTS = ("pdf", "png", "svg")


def group_logmean_diff(adata: sc.AnnData) -> pd.Series:
    """Mean log-norm(niche-20) - mean log-norm(other) per gene (the dotplot colour basis)."""
    is_c = (adata.obs["grp"].astype(str) == "crescent").values
    X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
    diff = X[is_c].mean(0) - X[~is_c].mean(0)
    return pd.Series(np.asarray(diff).ravel(), index=adata.var_names)


def select_top_bottom(de: pd.DataFrame, panel: set[str], n: int, padj: float,
                      require_consistent: bool, logmean_diff: pd.Series | None):
    df = de[de["gene"].isin(panel)].copy()
    df["abs_lfc"] = df["avg_log2FC"].abs()
    df = df[df["p_val_adj"] < padj]
    if require_consistent:
        assert logmean_diff is not None
        lmd = df["gene"].map(logmean_diff)
        df = df[np.sign(df["avg_log2FC"]) == np.sign(lmd)]
    def pick(sign):
        s = df[np.sign(df["avg_log2FC"]) == sign]
        return s.sort_values(["p_val_adj", "abs_lfc"], ascending=[True, False]).head(n)
    up, down = pick(1.0), pick(-1.0)
    assert not up.empty and not down.empty, (len(up), len(down))
    return up["gene"].tolist(), down["gene"].tolist(), pd.concat(
        [up.assign(direction="up"), down.assign(direction="down")], ignore_index=True
    )[["gene", "direction", "avg_log2FC", "p_val_adj", "pct.1", "pct.2"]]


def make_dotplot(adata, up, down, scale, title, tag):
    vg = {f"Up in {TARGET}": up, f"Down in {TARGET}": down}
    dp = sc.pl.dotplot(adata, var_names=vg, groupby="grp_label",
                       standard_scale=scale, return_fig=True, show=False)
    fig = dp.get_axes()["mainplot_ax"].get_figure() if hasattr(dp, "get_axes") else dp.fig
    for ax in fig.axes:
        if ax.get_title():
            ax.set_title("")
    fig.suptitle(title, y=1.06, fontsize=12)
    for ext in EXTS:
        fig.savefig(OUT_DIR / f"{tag}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {tag}.{{{','.join(EXTS)}}}")


def main() -> None:
    print(f"Loading {ADATA_PATH}")
    adata = sc.read_h5ad(ADATA_PATH)
    assert "grp" in adata.obs.columns
    print("  grp:", dict(adata.obs["grp"].value_counts()))

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    lab = np.where(adata.obs["grp"].astype(str) == "crescent", TARGET, OTHER)
    adata.obs["grp_label"] = pd.Categorical(lab, categories=[TARGET, OTHER], ordered=True)

    de = pd.read_csv(DE_CSV)
    panel = set(adata.var_names)
    lmd = group_logmean_diff(adata)

    # ---- v1: faithful padj-ranked reproduction (incl. AFDN), standard_scale='var' ----
    up1, down1, sm1 = select_top_bottom(de, panel, N, PADJ, require_consistent=False,
                                        logmean_diff=None)
    print(f"[v1  ] up={up1}\n       down={down1}")
    sm1.to_csv(OUT_DIR / "niche20_vs_other_pec_top5_bottom5_v1.csv", index=False)
    make_dotplot(adata, up1, down1, "var",
                 f"Top {N} up / Top {N} down — {TARGET} vs {OTHER}  (padj-ranked, scaled)",
                 "dotplot_niche20_vs_other_pec_top5_bottom5_v1_padjrank_varscale")

    # ---- final: direction-consistent selection + standard_scale=None ----
    up2, down2, sm2 = select_top_bottom(de, panel, N, PADJ, require_consistent=True,
                                        logmean_diff=lmd)
    print(f"[final] up={up2}\n       down={down2}")
    sm2 = sm2.assign(logmean_diff_niche20_minus_other=sm2["gene"].map(lmd).values)
    sm2.to_csv(OUT_DIR / "niche20_vs_other_pec_top5_bottom5.csv", index=False)
    make_dotplot(adata, up2, down2, None,
                 f"Top {N} up / Top {N} down — {TARGET} vs {OTHER}  (direction-consistent)",
                 "dotplot_niche20_vs_other_pec_top5_bottom5")

    print("\nDone ->", OUT_DIR)


if __name__ == "__main__":
    main()

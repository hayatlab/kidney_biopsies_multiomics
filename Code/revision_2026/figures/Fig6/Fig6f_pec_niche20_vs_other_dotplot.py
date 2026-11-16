#!/usr/bin/env python
"""Fig. 6f — Xenium marker genes in PECs: top-5 up / top-5 down genes of glomerular niche-20 PEC vs PEC of all other niches.
Panel        : Main Figures.pdf, Figure 6, panel f
Reads        : $KB_ANALYSIS/celltype_concordance/outputs/crescent_check/de/pec_manualglom10.h5ad (analyses/08_glomerular_niche20_crescent/00_build_pec_manualglom10.py)
               $KB_ANALYSIS/celltype_concordance/outputs/crescent_check/de/markers_manualglom10.csv (analyses/08_glomerular_niche20_crescent/de_crescent_pec_mast.R)
Writes       : $KB_FIG_OUT/Fig6/Fig6f_pec_niche20_vs_other_dotplot.pdf
Plotting code: analyses/08_glomerular_niche20_crescent/plot_niche20_pec_top5_bottom5.py::select_top_bottom + make_dotplot (imported
               unchanged). The published panel is an ad hoc variant of that script (28 Jul 2026, code not saved): the gene
               selection of its 'final' figure (padj < 0.05, direction-consistent, ranked by padj -> NDRG1 replaces AFDN)
               drawn with the colour scaling and title of its 'v1' figure (standard_scale='var': mean log-normalised
               expression min-max scaled 0-1 per gene). This script reproduces exactly that combination.
Illustrator  : title replaced ("Xenium marker genes in PECs"); rows renamed "niche 20" / "all other niches"; up/down
               brackets removed; legend titles re-added and legends moved
Environment  : conda env nichecompass
Verified     : 2026-11-12, PIXEL-IDENTICAL vs $KB_ANALYSIS/celltype_concordance/outputs/crescent_check/de/dotplot_niche20_vs_other_pec_top5_bottom5_ndrg1_varscale.pdf
"""
import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "config"))
from kbpaths import fig_out  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# importing the module applies its font registration and rcParams (Arial/Liberation Sans, pdf.fonttype 42)
m = _load(REPO / "analyses" / "08_glomerular_niche20_crescent" / "plot_niche20_pec_top5_bottom5.py", "plot_niche20_pec_top5_bottom5")


def main() -> None:
    # data preparation exactly as in plot_niche20_pec_top5_bottom5.py::main
    adata = sc.read_h5ad(m.ADATA_PATH)
    assert "grp" in adata.obs.columns
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    lab = np.where(adata.obs["grp"].astype(str) == "crescent", m.TARGET, m.OTHER)
    adata.obs["grp_label"] = pd.Categorical(lab, categories=[m.TARGET, m.OTHER], ordered=True)
    de = pd.read_csv(m.DE_CSV)
    panel = set(adata.var_names)
    lmd = m.group_logmean_diff(adata)

    # 'final' selection (direction-consistent) ...
    up, down, table = m.select_top_bottom(de, panel, m.N, m.PADJ, require_consistent=True, logmean_diff=lmd)
    print(f"up={up}\ndown={down}")
    # ... drawn with the 'v1' scaling and title
    out = Path(fig_out("Fig6", "Fig6f_pec_niche20_vs_other_dotplot.pdf"))
    m.OUT_DIR, m.EXTS = out.parent, ("pdf",)
    m.make_dotplot(adata, up, down, "var",
                   f"Top {m.N} up / Top {m.N} down — {m.TARGET} vs {m.OTHER}  (padj-ranked, scaled)", out.stem)


if __name__ == "__main__":
    main()

"""Shared helpers for re-annotating NicheCompass objects with the refined
`level2_refined` cell typing (from the celltype-concordance work) and for a
project-consistent colour palette.

The palette matches celltype_concordance/11_export_labels.py and 12_umap_* :
tab20 indexed by the fixed ORDER, so a given cell type has the SAME colour in
every figure (UMAP, dotplot, exported legend, and now the niche plots).
"""
from __future__ import annotations
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from matplotlib import cm, colors as mcolors

# Canonical order used across the concordance figures
ORDER = ["PT", "aPT", "dPT", "TAL", "DCT", "CNT", "PC", "IC-A", "IC-B",
         "POD", "PEC", "EC", "FIB", "MES", "VSM-P", "Immune"]

# Fixed hex colour per cell type (global ORDER index -> tab20), so colours are
# stable across the overall / glom-subset objects even if some types are absent.
_BASE = list(cm.tab20.colors)
L2_COLORS = {g: mcolors.to_hex(_BASE[i % 20]) for i, g in enumerate(ORDER)}

LABELED_DEFAULT = Path(
    f"{KB_ANALYSIS}/celltype_concordance/"
    "outputs/xenium_labeled.h5ad"
)


def load_level2_map(labeled_path: Path = LABELED_DEFAULT,
                    key: str = "level2_refined") -> dict:
    """(sample, cell_id) -> level2_refined, read from the labeled object's obs only."""
    lab = ad.read_h5ad(str(labeled_path), backed="r")
    lo = lab.obs
    return dict(zip(zip(lo["sample"].astype(str), lo["cell_id"].astype(str)),
                    lo[key].astype(str)))


def merge_level2(adata, labeled_path: Path = LABELED_DEFAULT,
                 key: str = "level2_refined", verbose: bool = True):
    """Add `level2_refined` to adata.obs by joining on (sample, cell_id).

    Sets an ordered categorical (ORDER, filtered to those present) and writes
    adata.uns[f"{key}_colors"] with the fixed palette in category order.
    Reports the number/percent of unmatched (NaN) cells explicitly and returns
    (adata, n_nan). Does NOT drop unmatched cells.
    """
    m = load_level2_map(labeled_path, key)
    keys = list(zip(adata.obs["sample"].astype(str),
                    adata.obs["cell_id"].astype(str)))
    vals = np.array([m.get(k, None) for k in keys], dtype=object)
    n = len(keys)
    n_nan = int(sum(v is None for v in vals))
    if verbose:
        print(f"[merge_level2] {n - n_nan}/{n} cells matched on (sample,cell_id); "
              f"{n_nan} unmatched -> NaN ({100 * n_nan / max(n, 1):.3f}%)", flush=True)
    present = [c for c in ORDER if c in set(v for v in vals if v is not None)]
    adata.obs[key] = pd.Categorical(vals, categories=present, ordered=True)
    adata.uns[f"{key}_colors"] = [L2_COLORS[c] for c in present]
    if verbose:
        print(f"[merge_level2] categories ({len(present)}): {present}", flush=True)
    return adata, n_nan


def level2_colour_dict(adata, key: str = "level2_refined") -> dict:
    """cat -> hex colour for the categories present in this object (fixed palette)."""
    if pd.api.types.is_categorical_dtype(adata.obs[key]):
        cats = list(adata.obs[key].cat.categories)
    else:
        cats = [c for c in ORDER if c in set(adata.obs[key].astype(str))]
    return {c: L2_COLORS.get(c, "#999999") for c in cats}

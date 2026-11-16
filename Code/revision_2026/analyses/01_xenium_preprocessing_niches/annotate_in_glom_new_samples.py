"""Add `in_glom` annotation for the two Region_1 Xenium samples that were
added during review (Combi 1734 Region 1, FSGS 42348 Region 1).

The parent enriched adata (`xenium_with_niches.h5ad`) already carries
`in_glom` for the three original ROI samples (`1734__Region_3`,
`42348__Region_2`, `42348__Region_3`) but `in_glom == False` for every
cell in the two newly-added Region_1 samples because no glom polygons
existed for them yet.

This script reads the two polygon CSVs in
`review_analyses/in_glom_for_new_samples/`, builds glom polygons (any
`Selection` containing "glom", case insensitive), tests cells in the
matching Xenium region against the union of those polygons using
`matplotlib.path.Path.contains_points` (same convention used in
`$KB_XENIUM_RAW/only_glom_roi/new_way_of_tying.ipynb`,
which produced the original `in_glom` labels), and writes a small
parquet with `cell_id` -> `in_glom` for the two Region_1 samples.

Outputs
-------
`review_analyses/in_glom_for_new_samples/in_glom_new_samples.parquet`
columns: `cell_id` (str), `sample` (str), `in_glom` (bool)

Notes
-----
- The polygon CSV X/Y are in the same coordinate system as
  `adata.obsm["spatial"]` for the corresponding sample (Xenium pixel
  coords); confirmed by matching bounding boxes.
- `cell_id` alone is not globally unique across regions (Xenium re-uses
  per-region cell_ids), so polygons are restricted to the cells of
  their own sample before the point-in-polygon test.
"""
from __future__ import annotations
from kbpaths import KB_ANALYSIS, KB_PROJECT  # repository path config (config/kbpaths.py)

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.path import Path as MplPath


PARENT_H5AD = Path(
    f"{KB_PROJECT}/niche_compass_results/artifacts/"
    "single_sample/22042026_221229/xenium_with_niches.h5ad"
)
ANNOT_DIR = Path(
    f"{KB_ANALYSIS}/in_glom_for_new_samples"
)
OUT_PARQUET = ANNOT_DIR / "in_glom_new_samples.parquet"

# Mapping: parent-adata sample name -> polygon CSV path
SAMPLE_CSV = {
    "XETG00229__0041734__Region_1__20250416__104148":
        ANNOT_DIR / "coordinates_combi_region1_480genepanel_incl_glom.csv",
    "XETG00229__0042348__Region_1__20250416__104148":
        ANNOT_DIR / "FSGS 42348 region1 with Glom coordinates.csv",
}


def load_polygons(csv_path: Path) -> pd.DataFrame:
    """Read QuPath-style polygon CSV (Selection,X,Y[,Class,Color]).

    The exported CSV starts with `#`-comment metadata lines (selection
    names + areas) that pandas treats as comments when `comment="#"`.
    """
    df = pd.read_csv(csv_path, comment="#")
    needed = {"Selection", "X", "Y"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"{csv_path}: missing columns {missing}")
    return df[["Selection", "X", "Y"]]


def annotate_in_glom(sample: str, csv_path: Path,
                     points: np.ndarray) -> np.ndarray:
    """Return a bool mask of length `len(points)` marking cells inside
    any glom polygon for this sample."""
    polygons = load_polygons(csv_path)
    is_glom_sel = polygons["Selection"].astype(str).str.lower().str.contains("glom")
    glom_polys = polygons.loc[is_glom_sel]
    n_polys = glom_polys["Selection"].nunique()
    print(f"  [{sample}] {len(polygons)} polygon vertices total, "
          f"{n_polys} glom polygons, {len(glom_polys)} glom vertices",
          flush=True)

    in_glom = np.zeros(len(points), dtype=bool)
    for sel, grp in glom_polys.groupby("Selection", sort=False):
        coords = grp[["X", "Y"]].to_numpy()
        if coords.shape[0] < 3:
            print(f"    [skip] {sel}: only {coords.shape[0]} vertices",
                  flush=True)
            continue
        if not np.array_equal(coords[0], coords[-1]):
            coords = np.vstack([coords, coords[0]])
        path = MplPath(coords)
        inside = path.contains_points(points)
        n_new = int((inside & ~in_glom).sum())
        in_glom |= inside
        if n_new:
            print(f"    [glom] {sel!r}: {int(inside.sum())} cells "
                  f"({n_new} newly tagged)", flush=True)
    return in_glom


def main() -> int:
    print(f"[load] {PARENT_H5AD}", flush=True)
    a = sc.read_h5ad(PARENT_H5AD)
    print(f"[load] adata {a.shape}", flush=True)

    rows = []
    for sample, csv_path in SAMPLE_CSV.items():
        if not csv_path.exists():
            print(f"[error] CSV missing for {sample}: {csv_path}", flush=True)
            return 1
        mask = (a.obs["sample"].astype(str) == sample).to_numpy()
        n = int(mask.sum())
        if n == 0:
            print(f"[warn] no cells found for sample {sample}", flush=True)
            continue
        print(f"[run] {sample}: {n} cells, csv={csv_path.name}", flush=True)

        points = a.obsm["spatial"][mask]
        cell_ids = a.obs.loc[mask, "cell_id"].astype(str).to_numpy()

        in_glom = annotate_in_glom(sample, csv_path, points)
        pct = 100.0 * in_glom.sum() / max(n, 1)
        print(f"  -> {int(in_glom.sum())}/{n} cells in_glom ({pct:.1f}%)",
              flush=True)

        rows.append(pd.DataFrame({
            "cell_id": cell_ids,
            "sample": sample,
            "in_glom": in_glom,
        }))

    if not rows:
        print("[error] no per-sample annotations produced", flush=True)
        return 1

    out = pd.concat(rows, ignore_index=True)
    print(f"[save] {OUT_PARQUET}  ({len(out)} rows, "
          f"{int(out['in_glom'].sum())} in_glom==True)", flush=True)
    out.to_parquet(OUT_PARQUET, index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())

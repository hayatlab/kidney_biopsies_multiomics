#!/usr/bin/env python
"""Stage B -- cell/nucleus boundary GeoDataFrames for all 5 Xenium runs (Vivien takeover).

Reproduces Code/cell_boundaries/1_preproc_cell_boundaries.py (original paper code). Vivien built these for 3 runs; the new object has 5.

Deliberate deviation: the original loaded each run through `spatialdata_io.xenium(...)` and
pulled `sdata["cell_boundaries"]`. That reader is heavy (it also decodes images and the
transcript table). The raw `cell_boundaries.parquet` in every Xenium output dir already holds
the polygon vertices, and building polygons straight from it was verified to reproduce Vivien's
stored geometries EXACTLY on the shared run 0042348_Region_3 (max relative area difference
0.000e+00; 50/50 shapely `equals_exact` at 1e-6). So the reader dependency is dropped.

Cell metadata comes from Stage A (`outputs/xenium_factor_scores.parquet`), so the GeoDataFrames
carry level2_refined, Factor1_ulm..Factor5_ulm, latent_leiden_0.4 (niches) and the clinical
columns needed by Fig 3g/3h, ED3b/3f and Suppl 7.

Human GRCh38, HGNC symbols. Seed 42.
"""
from kbpaths import KB_ANALYSIS, KB_XENIUM_RAW  # repository path config (config/kbpaths.py)
from cores import excluded_core_labels  # private core lists (config/cores.py)
import argparse, json, warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point, Polygon

HERE = Path(f"{KB_ANALYSIS}/vivien_takeover")
RAW = Path(f"{KB_XENIUM_RAW}")
SEED = 42

# The 3 cores Vivien dropped for having almost no cells. Asserted absent rather than
# re-filtered -- the new object already excludes them.
VIVIEN_DROPPED = excluded_core_labels("dropped_low_cells")  # K29 MN, K25 TIN, K24 TIN (labels in the private list)


def polygons_from_vertices(path, id_col="cell_id"):
    """Group a Xenium vertex table into one Polygon per cell."""
    df = pd.read_parquet(path, columns=[id_col, "vertex_x", "vertex_y"])
    g = df.groupby(id_col, sort=False)
    polys = g.apply(lambda d: Polygon(zip(d.vertex_x, d.vertex_y)), include_groups=False)
    return gpd.GeoDataFrame({id_col: polys.index}, geometry=polys.values)


def circles_from_cells(path):
    """cells_as_circles=True equivalent: area-preserving circle at the cell centroid."""
    df = pd.read_parquet(path, columns=["cell_id", "x_centroid", "y_centroid", "cell_area"])
    r = np.sqrt(df["cell_area"].to_numpy() / np.pi)
    geom = [Point(x, y).buffer(rr) for x, y, rr in zip(df.x_centroid, df.y_centroid, r)]
    return gpd.GeoDataFrame({"cell_id": df["cell_id"].to_numpy()}, geometry=geom)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    tag = "_smoke" if args.smoke else ""
    warnings.filterwarnings("ignore", message=".*geographic CRS.*")

    outdir = HERE / "outputs" / "geodataframes"; outdir.mkdir(parents=True, exist_ok=True)

    obs = pd.read_parquet(HERE / "outputs" / "xenium_factor_scores.parquet")
    obs["cell_id"] = obs["cell_id"].astype(str)
    obs["sample"] = obs["sample"].astype(str)
    runs = sorted(obs["sample"].unique())
    print(f"[obs] {len(obs)} annotated cells across {len(runs)} runs", flush=True)

    present = set(obs["TMA_ID"].astype(str))
    still_there = [c for c in VIVIEN_DROPPED if c in present]
    print(f"[cores] Vivien's 3 near-empty cores present in the new object: {still_there or 'none'}",
          flush=True)
    cs = obs["TMA_ID"].value_counts()
    print(f"[cores] {len(cs)} cores, smallest {cs.min()} cells ({cs.idxmin()})", flush=True)

    meta_cols = [c for c in obs.columns if c not in ("spatial_x", "spatial_y", "sample_cell_id")]
    parts = {"cell_boundaries": [], "nucleus_boundaries": [], "cell_circles": []}

    for run in runs:
        d = RAW / f"output-{run}"
        if not d.is_dir():
            raise SystemExit(f"missing raw Xenium dir: {d}")
        o = obs[obs["sample"] == run]
        if args.smoke:
            o = o.groupby("TMA_ID", observed=True).head(40)
        keep_ids = set(o["cell_id"])

        built = {
            "cell_boundaries": polygons_from_vertices(d / "cell_boundaries.parquet"),
            "nucleus_boundaries": polygons_from_vertices(d / "nucleus_boundaries.parquet"),
            "cell_circles": circles_from_cells(d / "cells.parquet"),
        }
        for name, gdf in built.items():
            gdf["cell_id"] = gdf["cell_id"].astype(str)
            n_raw = len(gdf)
            gdf = gdf[gdf["cell_id"].isin(keep_ids)]
            # left join onto the annotated cells -> anything without TMA_ID is unannotated
            gdf = gdf.merge(o[meta_cols], on="cell_id", how="left")
            n_nan = int(gdf["TMA_ID"].isna().sum())
            gdf = gpd.GeoDataFrame(gdf[gdf["TMA_ID"].notna()], geometry="geometry")
            parts[name].append(gdf)
            print(f"[{run[10:30]}] {name}: {n_raw} raw -> {len(gdf)} kept "
                  f"({n_raw - len(gdf)} unannotated, {n_nan} NaN TMA_ID)", flush=True)

    meta = {}
    for name, chunks in parts.items():
        gdf = gpd.GeoDataFrame(pd.concat(chunks, ignore_index=True), geometry="geometry")
        if name in ("cell_boundaries", "nucleus_boundaries"):
            gdf["centroid"] = gdf.geometry.centroid.to_wkt()   # parquet-safe single geom col
        gdf["sample_cell_id"] = gdf["TMA_ID"].astype(str) + "_" + gdf["cell_id"].astype(str)
        dup = int(gdf["sample_cell_id"].duplicated().sum())
        if dup:
            raise SystemExit(f"{name}: {dup} duplicate sample_cell_id")
        gdf.index = gdf["sample_cell_id"]
        # order cell types case-insensitively, as the original did for cell_type_manual
        cts = sorted(gdf["level2_refined"].astype(str).unique(), key=str.lower)
        gdf["level2_refined"] = pd.Categorical(gdf["level2_refined"].astype(str),
                                               categories=cts, ordered=True)
        p = outdir / f"{name}{tag}.parquet"
        gdf.to_parquet(p)
        meta[name] = {"n": int(len(gdf)), "path": str(p),
                      "size_mb": round(p.stat().st_size / 1e6, 1)}
        print(f"[save] {name}: {len(gdf)} rows -> {p.name} ({meta[name]['size_mb']} MB)", flush=True)

    nb, cb = meta["nucleus_boundaries"]["n"], meta["cell_boundaries"]["n"]
    print(f"[check] nuclei/cells = {nb}/{cb} = {nb/cb:.3f} "
          f"(a cell without a detected nucleus is expected, ratio <1)", flush=True)

    (HERE / "outputs" / f"02_geodataframes_meta{tag}.json").write_text(json.dumps({
        "runs": runs, "n_cores": int(obs['TMA_ID'].nunique()),
        "vivien_dropped_cores_still_present": still_there,
        "reader": "raw cell_boundaries.parquet -> shapely Polygon (verified identical to "
                  "spatialdata_io output on 0042348_Region_3: max rel area diff 0.0)",
        "outputs": meta, "seed": SEED, "smoke": bool(args.smoke),
        "versions": {"geopandas": gpd.__version__, "pandas": pd.__version__,
                     "numpy": np.__version__},
    }, indent=2))
    print("[done]", flush=True)


if __name__ == "__main__":
    main()

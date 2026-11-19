#!/usr/bin/env python
"""Stage C -- distance-bounded Voronoi graphs + 2-hop neighbourhood composition.

Reproduces the first-submission script nh_analysis/1_compute_nh.py (original paper code, no longer in this
repository) on the new 5-run Xenium object.

Algorithm, unchanged from the original:
  1. per TMA core, Voronoi diagram on cell CENTROIDS
  2. graph edges from `vor.ridge_points` (Voronoi adjacency)
  3. prune every edge whose POLYGON-to-POLYGON distance exceeds `--threshold` um
     (default 10 -- this is the "distance-bounded" in "distance-bounded Voronoi graph")
  4. per cell, the 2-hop neighbourhood = {root} u {nodes within graph distance 2},
     and count its composition (root INCLUDED, as in the original)

Two deliberate implementation changes, both verified not to alter the result:
  * The original precomputed a full N x N polygon distance matrix per core and pickled it
    (`distances.pkl`, 4.1 GB). Those distances are only ever READ for Voronoi ridge pairs
    (~6 per cell), so only those pairs are computed here. Same edges, no 4 GB artifact.
  * The original mutated `voronoi_graph` while iterating `.edges()`; edges to drop are
    collected first and removed afterwards.

Compositions are computed over BOTH `level2_refined` (16 cell types -> Fig 3h, ED3b) and
`latent_leiden_0.4` (10 niches -> ED3f), replacing the old cell_type_manual (21) / niche (12).

Human GRCh38, HGNC symbols. Seed 42.
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import argparse, json, pickle, time
from pathlib import Path

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
from scipy.spatial import Voronoi

HERE = Path(f"{KB_ANALYSIS}/vivien_takeover")
SEED = 42
CT_COL, NICHE_COL = "level2_refined", "latent_leiden_0.4"


def core_graph(gdf_core, threshold):
    """Distance-bounded Voronoi adjacency graph for one core, nodes = sample_cell_id."""
    cent = gdf_core.geometry.centroid
    coords = np.column_stack([cent.x.to_numpy(), cent.y.to_numpy()])
    if len(coords) < 5:
        return nx.Graph(), 0, 0
    vor = Voronoi(coords)
    ridge = np.asarray(vor.ridge_points)          # (n_ridges, 2) positional indices

    geoms = gdf_core.geometry.to_numpy()
    # polygon-to-polygon distance, ridge pairs only (the original built the full N x N matrix)
    d = np.fromiter((geoms[i].distance(geoms[j]) for i, j in ridge),
                    dtype=float, count=len(ridge))
    keep = ridge[d <= threshold]

    labels = gdf_core.index.to_numpy()
    g = nx.Graph()
    g.add_nodes_from(labels)
    g.add_edges_from(zip(labels[keep[:, 0]], labels[keep[:, 1]]))
    return g, len(ridge), len(keep)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--threshold", type=float, default=10.0)
    args = ap.parse_args()
    tag = "_smoke" if args.smoke else ""
    thr = args.threshold
    np.random.seed(SEED)

    outdir = HERE / "outputs" / "neighborhoods"; outdir.mkdir(parents=True, exist_ok=True)

    # smoke uses the FULL-density boundaries and limits the number of CORES -- subsetting
    # cells within a core would thin the Voronoi graph and give meaningless neighbourhoods
    src = HERE / "outputs" / "geodataframes" / "cell_boundaries.parquet"
    gdf = gpd.read_parquet(src)
    gdf.index = gdf["sample_cell_id"]
    print(f"[load] {src.name}: {len(gdf)} cells, {gdf['TMA_ID'].nunique()} cores", flush=True)

    cell_types = list(pd.Categorical(gdf[CT_COL]).categories)
    niches = sorted(gdf[NICHE_COL].astype(str).unique(), key=lambda s: int(s))
    print(f"[levels] {len(cell_types)} cell types {cell_types}", flush=True)
    print(f"[levels] {len(niches)} niches {niches}", flush=True)

    cores = list(pd.unique(gdf["TMA_ID"]))
    if args.smoke:
        cores = cores[:4]
        gdf = gdf[gdf["TMA_ID"].isin(cores)]

    graphs, ct_rows, niche_rows = {}, [], []
    n_ridge_tot = n_kept_tot = 0
    t0 = time.time()

    for k, core in enumerate(cores, 1):
        gc = gdf[gdf["TMA_ID"] == core]
        g, n_ridge, n_kept = core_graph(gc, thr)
        graphs[core] = g
        n_ridge_tot += n_ridge; n_kept_tot += n_kept

        ct = gc[CT_COL].astype(str)
        nz = gc[NICHE_COL].astype(str)
        for idx in gc.index:
            nb = nx.single_source_shortest_path_length(g, idx, cutoff=2).keys()
            members = list(nb)                       # includes the root, as in the original
            ct_rows.append({"TMA_ID": core, "sample_cell_id": idx,
                            "root_cell_type": ct.loc[idx],
                            **ct.loc[members].value_counts().reindex(cell_types, fill_value=0).to_dict()})
            niche_rows.append({"TMA_ID": core, "sample_cell_id": idx,
                               "root_niche": nz.loc[idx],
                               **nz.loc[members].value_counts().reindex(niches, fill_value=0).to_dict()})
        if k % 20 == 0 or k == len(cores):
            print(f"[core {k}/{len(cores)}] {core[:40]:40s} n={len(gc):5d} "
                  f"edges {n_kept}/{n_ridge}  ({time.time()-t0:.0f}s)", flush=True)

    print(f"[graph] {n_kept_tot}/{n_ridge_tot} Voronoi edges kept "
          f"({n_kept_tot/max(n_ridge_tot,1)*100:.1f}%) at threshold {thr} um", flush=True)

    with open(outdir / f"voronoi_graphs_{int(thr)}{tag}.pkl", "wb") as f:
        pickle.dump(graphs, f)

    ct_df = pd.DataFrame(ct_rows).set_index("sample_cell_id")
    nz_df = pd.DataFrame(niche_rows).set_index("sample_cell_id")
    dis = gdf.groupby("TMA_ID", observed=True)["condition"].first()
    ct_df["disease"] = ct_df["TMA_ID"].map(dis); nz_df["disease"] = nz_df["TMA_ID"].map(dis)
    ct_df.to_csv(outdir / f"nh_2hop_compo_{int(thr)}{tag}.csv")
    nz_df.to_csv(outdir / f"nh_2hop_niche_compo_{int(thr)}{tag}.csv")

    size = ct_df[cell_types].sum(axis=1)
    iso = int((size == 1).sum())
    print(f"[nh] 2-hop size: median {size.median():.0f}, mean {size.mean():.1f}, "
          f"p5 {size.quantile(.05):.0f}, p95 {size.quantile(.95):.0f}, max {size.max()}", flush=True)
    print(f"[nh] isolated cells (neighbourhood = itself only): {iso} ({iso/len(size)*100:.2f}%)", flush=True)
    print(f"[save] {outdir}", flush=True)

    (HERE / "outputs" / f"03_neighborhoods_meta{tag}.json").write_text(json.dumps({
        "source": str(src), "threshold_um": thr, "n_cells": int(len(gdf)),
        "n_cores": len(cores), "cell_types": cell_types, "niches": niches,
        "voronoi_edges_total": int(n_ridge_tot), "voronoi_edges_kept": int(n_kept_tot),
        "nh_size_median": float(size.median()), "nh_size_mean": float(size.mean()),
        "isolated_cells": iso, "seed": SEED, "smoke": bool(args.smoke),
        "versions": {"geopandas": gpd.__version__, "networkx": nx.__version__,
                     "numpy": np.__version__, "pandas": pd.__version__},
    }, indent=2))
    print("[done]", flush=True)


if __name__ == "__main__":
    main()

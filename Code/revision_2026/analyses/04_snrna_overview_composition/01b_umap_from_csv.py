#!/usr/bin/env python
"""Build umap_and_obs.parquet from Leon's supplied UMAP embeddings.

Supersedes 01_compute_umap.py, which computed a UMAP from X_harmony because the object
shipped without one. Leon provided the real embedding on 2026-09-02, so the computed one
is discarded -- the figures must show the same coordinates as the rest of the project.

  embeddings : review_analyses/umap_embeddings.csv  (cellid, umap_1, umap_2)
               363,763 rows; covers all 363,656 cells of the object plus 107 cells that
               were filtered out of it. Joined on cellid, reindexed to the object's order.

Output schema is identical to 01_compute_umap.py's, so 02_overview_figures.py is unchanged.
Human GRCh38, HGNC symbols.
"""
from kbpaths import KB_ANALYSIS, KB_DATA  # repository path config (config/kbpaths.py)
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

HERE = Path(f"{KB_ANALYSIS}/snrna_overview_20260831")
SRC = f"{KB_DATA}/charlotte_new_data/adata_20260717_with_clinical.h5ad"
UMAP_CSV = f"{KB_ANALYSIS}/umap_embeddings.csv"
OBS_KEYS = ["sample", "patient_ID", "sex_clinical", "disease", "Study",
            "level1", "level2", "nCount_RNA", "nFeature_RNA", "pct_counts_mt"]


def read_obs():
    with h5py.File(SRC, "r") as f:
        g = f["obs"]
        index = [x.decode() if isinstance(x, bytes) else x for x in g[g.attrs["_index"]][:]]
        out = {}
        for k in OBS_KEYS:
            n = g[k]
            if isinstance(n, h5py.Group) and "categories" in n:
                cats = [x.decode() if isinstance(x, bytes) else x for x in n["categories"][:]]
                out[k] = pd.Categorical.from_codes(n["codes"][:], cats)
            else:
                v = n[:]
                out[k] = np.array([x.decode() for x in v]) if v.dtype.kind == "S" else v
    return pd.DataFrame(out, index=pd.Index(index, name="cellid"))


def main():
    obs = read_obs()
    print(f"[obs] {obs.shape} from {Path(SRC).name}", flush=True)

    um = pd.read_csv(UMAP_CSV, index_col=0)
    print(f"[umap] {um.shape} from {Path(UMAP_CSV).name}", flush=True)
    missing = obs.index.difference(um.index)
    extra = um.index.difference(obs.index)
    print(f"[join] object cells without coordinates: {len(missing)}; "
          f"coordinates for cells not in the object (ignored): {len(extra)}", flush=True)
    if len(missing):
        raise SystemExit(f"{len(missing)} cells have no UMAP coordinates, e.g. {list(missing[:5])}")

    um = um.reindex(obs.index)
    obs.insert(0, "UMAP1", um["umap_1"].to_numpy())
    obs.insert(1, "UMAP2", um["umap_2"].to_numpy())
    print(f"[umap] UMAP1 [{obs.UMAP1.min():.2f}, {obs.UMAP1.max():.2f}]  "
          f"UMAP2 [{obs.UMAP2.min():.2f}, {obs.UMAP2.max():.2f}]  "
          f"NaN {int(obs[['UMAP1','UMAP2']].isna().sum().sum())}", flush=True)

    out = HERE / "umap_and_obs.parquet"
    obs.to_parquet(out)
    print(f"[save] {out} {obs.shape}", flush=True)

    (HERE / "umap_run_meta.json").write_text(json.dumps({
        "source_h5ad": SRC, "umap_source": UMAP_CSV,
        "umap_origin": "supplied by Leon 2026-09-02; replaces the X_harmony UMAP "
                       "computed by 01_compute_umap.py",
        "n_cells": int(len(obs)), "n_coords_not_in_object": int(len(extra)),
        "species": "human GRCh38", "gene_id_format": "HGNC symbols",
        "versions": {"numpy": np.__version__, "pandas": pd.__version__, "h5py": h5py.__version__},
    }, indent=2))
    print("[done]", flush=True)


if __name__ == "__main__":
    main()

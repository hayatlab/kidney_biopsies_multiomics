"""Bake the newest cell-type annotation into the three NicheCompass niche objects so they
are self-contained (show level2_refined when opened directly, not just at plot time).

Adds, by joining MAIN xenium_labeled.h5ad on (sample, cell_id):
  * level2_refined   (ordered categorical + fixed palette in uns, via merge_level2)
  * subcluster_level2
  * in_glom          (the COMPLETED 5-sample flag, overwriting the stale/incomplete copy)

Backs up each object first. Writes in place (paths the replot pipeline already uses).
"""
from __future__ import annotations
from kbpaths import KB_RESULTS  # repository path config (config/kbpaths.py)
import shutil
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from level2_common import merge_level2, LABELED_DEFAULT   # noqa: E402

NEW = Path(f"{KB_RESULTS}/xenium/leon/new_data")
ARCH = Path(f"{KB_RESULTS}/xenium/leon/archived")
TS = "20260709"

OBJS = {
    "overall": ARCH / "overall_with_glom_04/xenium_with_niches_overall_with_glom.h5ad",
    "in_glom_niche": NEW / "in_glom_niche_based_on_04/xenium_with_niches_in_glom_niche.h5ad",
    "manual_glom": ARCH / "manual_glom_niche_based_on_04/xenium_with_niches_manual_glom_niche.h5ad",
}


def main() -> int:
    mo = ad.read_h5ad(str(LABELED_DEFAULT), backed="r").obs
    mkey = list(zip(mo["sample"].astype(str), mo["cell_id"].astype(str)))
    scl_map = dict(zip(mkey, mo["subcluster_level2"].astype(str)))
    ig_map = dict(zip(mkey, mo["in_glom"].astype(bool)))
    print(f"[main] {len(mkey)} cells; subcluster_level2 + in_glom maps built", flush=True)

    for name, path in OBJS.items():
        backup = Path(str(path)[:-5] + f".PRE_LEVEL2_{TS}.h5ad")
        if not backup.exists():
            shutil.copy2(path, backup)
            print(f"\n[{name}] backup -> {backup.name}", flush=True)
        else:
            print(f"\n[{name}] backup already exists: {backup.name}", flush=True)

        a = sc.read_h5ad(str(path))
        a, n_nan = merge_level2(a, key="level2_refined")     # level2_refined + colors
        if n_nan:
            raise SystemExit(f"[{name}] {n_nan} unmatched cells for level2_refined -- aborting")
        keys = list(zip(a.obs["sample"].astype(str), a.obs["cell_id"].astype(str)))
        a.obs["subcluster_level2"] = pd.Categorical([scl_map.get(k) for k in keys])
        a.obs["in_glom"] = np.array([ig_map.get(k, False) for k in keys], dtype=bool)

        vc = a.obs["level2_refined"].astype(str).value_counts()
        print(f"[{name}] baked: PEC={vc.get('PEC',0)} aPT={vc.get('aPT',0)} "
              f"| in_glom True={int(a.obs['in_glom'].sum())} "
              f"| subcluster_level2 recovered tags={int((a.obs['subcluster_level2'].astype(str)=='aPTrec_PEC').sum())}",
              flush=True)
        a.write_h5ad(str(path))
        print(f"[{name}] wrote self-contained object -> {path}", flush=True)

    print("\n[ALL DONE] level2_refined + subcluster_level2 + completed in_glom baked into all 3 objects", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

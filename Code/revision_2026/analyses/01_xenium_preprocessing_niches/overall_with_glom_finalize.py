"""Resume the post-training steps that were cut off when job 321201 hit
the 8h SLURM time limit on 2026-05-15 05:10.

Training itself finished (early-stopped at epoch 36) and
`<run_dir>/model/adata.h5ad` was written successfully — what timed out
was the `sc.pp.neighbors` + `sc.tl.umap` + `sc.tl.leiden` + standalone
adata write block.

This script:
1. Loads `model/adata.h5ad` from the latest overall_with_glom run.
2. Re-runs neighbours / UMAP / leiden(LEIDEN_RES) on `nichecompass_latent`.
3. Writes the niche-9 recovery diagnostic CSV under
   `<run_dir>/figures/niche9_match_res_<LEIDEN_RES>.csv`.
4. Writes the enriched `<run_dir>/adata.h5ad` so downstream tools
   (`glom_postprocess.py`, `copy_results_to_kramann.py`) pick it up.

Usage: `python -u review_analyses/overall_with_glom_finalize.py [--run-dir <path>]`
"""
from __future__ import annotations
from kbpaths import KB_PROJECT  # repository path config (config/kbpaths.py)

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


OVERALL_BASE = Path(
    f"{KB_PROJECT}/niche_compass_results/artifacts/overall_with_glom"
)
LATENT_KEY = "nichecompass_latent"
LEIDEN_RES = 0.4
IN_GLOM_KEY = "in_glom"


def latest_run_dir(base: Path) -> Path:
    runs = [p for p in base.iterdir() if p.is_dir()
            and (p / "model" / "adata.h5ad").exists()]
    if not runs:
        raise FileNotFoundError(f"no run dirs with model/adata.h5ad under {base}")
    runs.sort(key=lambda p: p.stat().st_mtime)
    return runs[-1]


def report_niche9_match(adata, niche_key: str, out_csv: Path) -> None:
    ct = (
        adata.obs.groupby(niche_key, observed=True)[IN_GLOM_KEY]
        .agg(n_total="count", n_in_glom="sum")
    )
    ct["pct_in_glom"] = 100.0 * ct["n_in_glom"] / ct["n_total"]
    total_in_glom = int(ct["n_in_glom"].sum())
    ct["share_of_in_glom"] = 100.0 * ct["n_in_glom"] / max(total_in_glom, 1)
    ct = ct.sort_values("pct_in_glom", ascending=False)
    ct.to_csv(out_csv)
    print(f"\n[niche9_match] saved {out_csv}", flush=True)
    print(ct.to_string(), flush=True)
    if not ct.empty:
        top = ct.index[0]
        print(f"[niche9_match] top niche by pct_in_glom: {top!r} "
              f"({ct.loc[top, 'pct_in_glom']:.1f}% of its cells, "
              f"{ct.loc[top, 'share_of_in_glom']:.1f}% of all in_glom)",
              flush=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", type=Path, default=None)
    ap.add_argument("--leiden-res", type=float, default=LEIDEN_RES)
    args = ap.parse_args()

    run_dir = args.run_dir if args.run_dir is not None else latest_run_dir(OVERALL_BASE)
    print(f"[run] run_dir={run_dir}", flush=True)
    src = run_dir / "model" / "adata.h5ad"
    if not src.exists():
        print(f"[error] missing {src}", flush=True)
        return 1

    print(f"[load] {src}", flush=True)
    adata = sc.read_h5ad(src)
    print(f"[load] shape={adata.shape}  obsm={list(adata.obsm.keys())}", flush=True)

    if LATENT_KEY not in adata.obsm:
        print(f"[error] {LATENT_KEY} missing from obsm", flush=True)
        return 1

    print(f"[neighbors] use_rep={LATENT_KEY}", flush=True)
    sc.pp.neighbors(adata, use_rep=LATENT_KEY, key_added=LATENT_KEY)

    print(f"[umap] running", flush=True)
    sc.tl.umap(adata, neighbors_key=LATENT_KEY)

    leiden_key = f"latent_leiden_{args.leiden_res}"
    print(f"[leiden] res={args.leiden_res} -> {leiden_key}", flush=True)
    sc.tl.leiden(adata=adata, resolution=args.leiden_res,
                 key_added=leiden_key, neighbors_key=LATENT_KEY)
    print(f"[leiden] niches: {adata.obs[leiden_key].value_counts().to_dict()}",
          flush=True)

    fig_dir = run_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    report_niche9_match(
        adata, leiden_key,
        out_csv=fig_dir / f"niche9_match_res_{args.leiden_res}.csv",
    )

    standalone = run_dir / "adata.h5ad"
    print(f"[save] writing standalone adata -> {standalone}", flush=True)
    for col in adata.obs.columns:
        if adata.obs[col].dtype == object:
            adata.obs[col] = adata.obs[col].astype(str)
    adata.write_h5ad(standalone, compression="gzip")
    size_mb = standalone.stat().st_size / 1e6
    print(f"[save] wrote {standalone} ({size_mb:.1f} MB)", flush=True)

    print(f"\n=== DONE run_dir={run_dir} ===", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

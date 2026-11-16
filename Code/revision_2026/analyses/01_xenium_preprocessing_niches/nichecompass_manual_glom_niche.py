"""NicheCompass on the manual glomerulus subset (in_glom == True).

Sibling of `nichecompass_glom_subniche.py`, which filtered to the
NicheCompass parent niche 5. This variant uses the human-annotated
glomerulus polygons (`adata.obs['in_glom']`) instead.

Pipeline
--------
1. Load the parent enriched adata (`xenium_with_niches.h5ad`).
2. Merge in the `in_glom` values for the two Region_1 samples that
   were added during review. The merge is driven by
   `in_glom_new_samples.parquet`, produced by
   `annotate_in_glom_new_samples.py` (which point-in-polygon-tests the
   QuPath CSVs in `review_analyses/in_glom_for_new_samples/`).
3. Subset to `in_glom == True`.
4. Apply the existing TMA blacklist (plus the Combi D2 IgAN QC-fail core).
5. Rebuild the spatial graph + GP masks; retrain NicheCompass with the
   exact same hyperparameters used for `nichecompass_glom_subniche.py`.
6. Save model + adata.

Downstream analyses (composition, GP marker plots, etc.) can be done
with `nichecompass_5regions_downstream.py --leiden-res 0.6`, pointed at
the new run directory.
"""
from __future__ import annotations
from kbpaths import KB_PROJECT  # repository path config (config/kbpaths.py)
from cores import excluded_core_labels  # private core lists (config/cores.py)

import argparse
import sys
import warnings
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

import nichecompass as nc
from nichecompass.models import NicheCompass
from nichecompass.utils import (
    add_gps_from_gp_dict_to_adata,
    extract_gp_dict_from_mebocost_ms_interactions,
    extract_gp_dict_from_nichenet_lrt_interactions,
    extract_gp_dict_from_omnipath_lr_interactions,
    filter_and_combine_gp_dict_gps_v2,
)

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

# ---------------------------------------------------------------------------
PROJECT_DIR = Path(f"{KB_PROJECT}")
PARENT_RUN_DIR = Path(
    f"{KB_PROJECT}/niche_compass_results/artifacts/single_sample/22042026_221229"
)
PARENT_H5AD = PARENT_RUN_DIR / "xenium_with_niches.h5ad"
GP_DATA_DIR = PROJECT_DIR / "niche_compass_results/data/gene_programs"
ARTIFACTS_DIR = PROJECT_DIR / "niche_compass_results/artifacts/manual_glom_niche"

NEW_IN_GLOM_PARQUET = PROJECT_DIR / (
    "review_analyses/in_glom_for_new_samples/in_glom_new_samples.parquet"
)

PARENT_NICHE_KEY = "latent_leiden_0.4"
IN_GLOM_KEY = "in_glom"

CELL_TYPE_KEY = "predic"
SAMPLE_KEY = "sample"
CONDITION_KEY = "condition"
TMA_KEY = "TMA_ID"
SPATIAL_KEY = "spatial"
ADJ_KEY = "spatial_connectivities"
COUNTS_KEY = "counts"

GP_TARGETS_MASK_KEY = "nichecompass_gp_targets"
GP_TARGETS_CATEGORIES_MASK_KEY = "nichecompass_gp_targets_categories"
GP_SOURCES_MASK_KEY = "nichecompass_gp_sources"
GP_SOURCES_CATEGORIES_MASK_KEY = "nichecompass_gp_sources_categories"
GP_NAMES_KEY = "nichecompass_gp_names"
ACTIVE_GP_NAMES_KEY = "nichecompass_active_gp_names"
LATENT_KEY = "nichecompass_latent"

# graph + model + training hyperparameters (identical to parent / sibling)
N_NEIGHBORS = 4
N_EPOCHS = 400
N_EPOCHS_ALL_GPS = 25
LR = 1e-3
LAMBDA_EDGE_RECON = 5e5
LAMBDA_GENE_EXPR_RECON = 300.0
LAMBDA_L1_MASKED = 0.0
LAMBDA_L1_ADDON = 30.0
EDGE_BATCH_SIZE = 1024
N_SAMPLED_NEIGHBORS = 4
LEIDEN_RES = 0.6  # finer sub-clustering inside the glom subset

# Cores excluded before NicheCompass: unassigned cells, the Apollo d0 A2 control core and ten disease cores that
# failed QC (Utrecht cores of K4 DN, K24 TIN, K29 MN, K23 TIN, K25 TIN, K71 FSGS, K46 MN, K78 FSGS and K60 IgAN,
# and the Combi D2 IgAN core). Their internal labels contain pathology numbers and are read from a private list
# (config/cores.py); none of these cores is part of the de-identified Zenodo object.
TO_EXCLUDE = set(excluded_core_labels("nichecompass_blacklist"))


def merge_new_in_glom(adata: ad.AnnData) -> ad.AnnData:
    """Overlay the in_glom values for the two new Region_1 samples on
    top of the parent adata. The parent already has in_glom for the 3
    original ROI samples; for the two new ones it's all False until we
    overlay the polygon-based annotation produced by
    annotate_in_glom_new_samples.py."""
    if not NEW_IN_GLOM_PARQUET.exists():
        raise FileNotFoundError(
            f"Missing per-sample in_glom annotations parquet: "
            f"{NEW_IN_GLOM_PARQUET}\n"
            f"Run review_analyses/annotate_in_glom_new_samples.py first."
        )
    annot = pd.read_parquet(NEW_IN_GLOM_PARQUET)
    print(f"[merge] in_glom_new_samples.parquet: {len(annot)} rows, "
          f"{int(annot['in_glom'].sum())} in_glom==True", flush=True)

    # Build (sample, cell_id) key — Xenium cell_ids re-use per-region
    annot_key = (
        annot["sample"].astype(str) + "::" + annot["cell_id"].astype(str)
    )
    annot_map = pd.Series(annot["in_glom"].to_numpy(), index=annot_key)

    obs_key = (
        adata.obs[SAMPLE_KEY].astype(str) + "::"
        + adata.obs["cell_id"].astype(str)
    ).to_numpy()
    in_glom_new = pd.Series(obs_key).map(annot_map)

    # only overlay where we actually have an annotation (i.e. the two new samples)
    have = in_glom_new.notna().to_numpy()
    print(f"[merge] cells in parent that match parquet rows: {int(have.sum())}",
          flush=True)
    before = int(adata.obs[IN_GLOM_KEY].sum())
    new_vals = adata.obs[IN_GLOM_KEY].to_numpy().copy()
    new_vals[have] = in_glom_new[have].astype(bool).to_numpy()
    adata.obs[IN_GLOM_KEY] = new_vals.astype(bool)
    after = int(adata.obs[IN_GLOM_KEY].sum())
    print(f"[merge] in_glom True: {before} -> {after} (+{after - before})",
          flush=True)

    # per-sample sanity
    summary = (
        adata.obs.groupby(SAMPLE_KEY, observed=True)[IN_GLOM_KEY]
        .agg(["sum", "count"])
        .rename(columns={"sum": "n_in_glom", "count": "n_total"})
    )
    print(f"[merge] per-sample summary:\n{summary}", flush=True)
    return adata


def load_glom_subset() -> ad.AnnData:
    print(f"[load] reading {PARENT_H5AD}", flush=True)
    a = sc.read_h5ad(PARENT_H5AD)
    print(f"[load] full: {a.shape}", flush=True)

    a = merge_new_in_glom(a)

    # subset to in_glom == True
    mask = a.obs[IN_GLOM_KEY].astype(bool).to_numpy()
    n_in_glom = int(mask.sum())
    print(f"[subset] cells with in_glom==True: {n_in_glom}", flush=True)
    a = a[mask].copy()

    # apply TMA blacklist
    before = a.n_obs
    keep = (
        (~a.obs[TMA_KEY].astype(str).isin(TO_EXCLUDE))
        & a.obs[TMA_KEY].notna()
        & (a.obs[TMA_KEY].astype(str) != "nan")
    )
    a = a[keep].copy()
    if pd.api.types.is_categorical_dtype(a.obs[TMA_KEY]):
        a.obs[TMA_KEY] = a.obs[TMA_KEY].cat.remove_unused_categories()
    print(f"[subset] after TMA blacklist: {a.shape} (dropped {before - a.n_obs})",
          flush=True)

    # ensure counts layer exists
    if COUNTS_KEY not in a.layers:
        print(f"[load] adding layers['{COUNTS_KEY}'] from .X", flush=True)
        a.layers[COUNTS_KEY] = a.X.copy()

    # drop parent leiden / parent latent / etc. so the new model learns fresh
    drop_obsm = [k for k in list(a.obsm.keys())
                 if k in (LATENT_KEY, "X_pca_harmony", "X_umap")]
    for k in drop_obsm:
        del a.obsm[k]
    drop_obs = [k for k in list(a.obs.columns)
                if k.startswith("latent_leiden_")]
    if drop_obs:
        a.obs = a.obs.drop(columns=drop_obs)
    drop_uns = [k for k in list(a.uns.keys())
                if k.startswith("nichecompass_") or k == "neighbors" or k == LATENT_KEY]
    for k in drop_uns:
        a.uns.pop(k, None)
    drop_obsp = [k for k in list(a.obsp.keys())
                 if k.startswith("nichecompass_") or k.startswith(LATENT_KEY)]
    for k in drop_obsp:
        del a.obsp[k]

    # sanitize obs keys for h5ad writes
    bad_obs = {c: c.replace("/", "_").replace("?", "_").replace("#", "_").replace("%", "_pct")
               for c in a.obs.columns if any(ch in c for ch in "/?#%")}
    if bad_obs:
        print(f"[load] sanitizing obs columns: {bad_obs}", flush=True)
        a.obs = a.obs.rename(columns=bad_obs)

    print(f"[load] samples: {sorted(a.obs[SAMPLE_KEY].astype(str).unique())}",
          flush=True)
    print(f"[load] conditions: {a.obs[CONDITION_KEY].astype(str).value_counts().to_dict()}",
          flush=True)
    print(f"[load] TMAs: {a.obs[TMA_KEY].nunique()} unique", flush=True)
    return a


def build_gp_dict() -> dict:
    GP_DATA_DIR.mkdir(parents=True, exist_ok=True)
    print("[gp] loading omnipath / nichenet / mebocost gene programs", flush=True)

    omnipath_csv = GP_DATA_DIR / "omnipath_lr_network.csv"
    omnipath_gp_dict = extract_gp_dict_from_omnipath_lr_interactions(
        species="human",
        load_from_disk=omnipath_csv.exists(),
        save_to_disk=not omnipath_csv.exists(),
        lr_network_file_path=str(omnipath_csv),
        plot_gp_gene_count_distributions=False,
    )
    print(f"[gp] omnipath gps: {len(omnipath_gp_dict)}", flush=True)

    nichenet_files = {
        "lr_network_file_path": str(GP_DATA_DIR / "nichenet_lr_network_v2_human.csv"),
        "ligand_target_matrix_file_path": str(GP_DATA_DIR / "nichenet_ligand_target_matrix_v2_human.csv"),
    }
    nichenet_gp_dict = extract_gp_dict_from_nichenet_lrt_interactions(
        species="human",
        version="v2",
        keep_target_genes_ratio=1.0,
        max_n_target_genes_per_gp=250,
        load_from_disk=all(Path(p).exists() for p in nichenet_files.values()),
        save_to_disk=not all(Path(p).exists() for p in nichenet_files.values()),
        plot_gp_gene_count_distributions=False,
        **nichenet_files,
    )
    print(f"[gp] nichenet gps: {len(nichenet_gp_dict)}", flush=True)

    mebocost_dir = GP_DATA_DIR / "metabolite_enzyme_sensor_gps"
    mebocost_dir.mkdir(parents=True, exist_ok=True)
    mebocost_gp_dict = extract_gp_dict_from_mebocost_ms_interactions(
        species="human",
        dir_path=str(mebocost_dir),
        plot_gp_gene_count_distributions=False,
    )
    print(f"[gp] mebocost gps: {len(mebocost_gp_dict)}", flush=True)

    combined = filter_and_combine_gp_dict_gps_v2(
        gp_dicts=[omnipath_gp_dict, nichenet_gp_dict, mebocost_gp_dict],
        verbose=False,
    )
    print(f"[gp] combined filtered gps: {len(combined)}", flush=True)
    return combined


def build_spatial_graph(adata: ad.AnnData) -> None:
    print(f"[graph] spatial_neighbors(n_neighs={N_NEIGHBORS})", flush=True)
    sq.gr.spatial_neighbors(
        adata, coord_type="generic", spatial_key=SPATIAL_KEY,
        n_neighs=N_NEIGHBORS,
    )
    adata.obsp[ADJ_KEY] = adata.obsp[ADJ_KEY].maximum(adata.obsp[ADJ_KEY].T)


def fit_nichecompass(adata: ad.AnnData, gp_dict: dict) -> NicheCompass:
    print("[nc] add_gps_from_gp_dict_to_adata", flush=True)
    add_gps_from_gp_dict_to_adata(
        gp_dict=gp_dict, adata=adata,
        gp_targets_mask_key=GP_TARGETS_MASK_KEY,
        gp_targets_categories_mask_key=GP_TARGETS_CATEGORIES_MASK_KEY,
        gp_sources_mask_key=GP_SOURCES_MASK_KEY,
        gp_sources_categories_mask_key=GP_SOURCES_CATEGORIES_MASK_KEY,
        gp_names_key=GP_NAMES_KEY,
        min_genes_per_gp=2, min_source_genes_per_gp=1, min_target_genes_per_gp=1,
        max_genes_per_gp=None, max_source_genes_per_gp=None, max_target_genes_per_gp=None,
    )
    print("[nc] init model", flush=True)
    model = NicheCompass(
        adata, counts_key=COUNTS_KEY, adj_key=ADJ_KEY,
        gp_names_key=GP_NAMES_KEY, active_gp_names_key=ACTIVE_GP_NAMES_KEY,
        gp_targets_mask_key=GP_TARGETS_MASK_KEY,
        gp_targets_categories_mask_key=GP_TARGETS_CATEGORIES_MASK_KEY,
        gp_sources_mask_key=GP_SOURCES_MASK_KEY,
        gp_sources_categories_mask_key=GP_SOURCES_CATEGORIES_MASK_KEY,
        latent_key=LATENT_KEY, conv_layer_encoder="gatv2conv",
        active_gp_thresh_ratio=0.01,
    )
    print(f"[nc] train  n_epochs={N_EPOCHS}  n_epochs_all_gps={N_EPOCHS_ALL_GPS}",
          flush=True)
    model.train(
        n_epochs=N_EPOCHS, n_epochs_all_gps=N_EPOCHS_ALL_GPS, lr=LR,
        lambda_edge_recon=LAMBDA_EDGE_RECON,
        lambda_gene_expr_recon=LAMBDA_GENE_EXPR_RECON,
        lambda_l1_masked=LAMBDA_L1_MASKED, lambda_l1_addon=LAMBDA_L1_ADDON,
        edge_batch_size=EDGE_BATCH_SIZE, n_sampled_neighbors=N_SAMPLED_NEIGHBORS,
        use_cuda_if_available=True, verbose=False,
    )
    return model


def main() -> int:
    global N_EPOCHS
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-epochs", type=int, default=N_EPOCHS)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    N_EPOCHS = args.n_epochs

    ts = datetime.now().strftime("%d%m%Y_%H%M%S")
    run_dir = ARTIFACTS_DIR / ts
    model_dir = run_dir / "model"
    fig_dir = run_dir / "figures"
    model_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    print(f"[run] timestamp={ts}", flush=True)
    print(f"[run] model_dir={model_dir}", flush=True)
    print(f"[run] fig_dir={fig_dir}", flush=True)

    adata = load_glom_subset()

    if args.smoke:
        rng = np.random.default_rng(0)
        keep = rng.choice(adata.n_obs, min(5_000, adata.n_obs), replace=False)
        adata = adata[keep].copy()
        N_EPOCHS = 3
        print(f"[smoke] adata -> {adata.shape}, n_epochs={N_EPOCHS}", flush=True)

    gp_dict = build_gp_dict()
    build_spatial_graph(adata)
    model = fit_nichecompass(adata, gp_dict)

    # neighbours + UMAP + leiden 0.6 BEFORE saving so the saved adata already has them
    sc.pp.neighbors(model.adata, use_rep=LATENT_KEY, key_added=LATENT_KEY)
    sc.tl.umap(model.adata, neighbors_key=LATENT_KEY)
    leiden_key = f"latent_leiden_{LEIDEN_RES}"
    sc.tl.leiden(adata=model.adata, resolution=LEIDEN_RES,
                 key_added=leiden_key, neighbors_key=LATENT_KEY)
    print(f"[leiden res={LEIDEN_RES}] {model.adata.obs[leiden_key].value_counts().to_dict()}",
          flush=True)

    # save adata first as a standalone file so we don't lose ~60 min training
    # if model.save trips downstream
    standalone_adata = run_dir / "adata.h5ad"
    print(f"[save] writing standalone adata -> {standalone_adata}", flush=True)
    try:
        for col in model.adata.obs.columns:
            if model.adata.obs[col].dtype == object:
                model.adata.obs[col] = model.adata.obs[col].astype(str)
        model.adata.write_h5ad(standalone_adata, compression="gzip")
        print(f"[save] standalone adata written "
              f"({standalone_adata.stat().st_size / 1e6:.1f} MB)", flush=True)
    except Exception as exc:
        print(f"[save] standalone adata write FAILED: {exc!r}", flush=True)
        import traceback as tb; tb.print_exc()

    print(f"[save] model.save -> {model_dir}", flush=True)
    try:
        model.save(dir_path=str(model_dir), overwrite=True,
                   save_adata=True, adata_file_name="adata.h5ad")
        print(f"[save] model + adata saved to {model_dir}", flush=True)
    except Exception as exc:
        print(f"[save] model.save FAILED: {exc!r}", flush=True)
        import traceback as tb
        tb.print_exc()
        try:
            import torch
            torch.save(model.model.state_dict(), model_dir / "model_state_dict.pt")
            print(f"[save] state_dict fallback written -> "
                  f"{model_dir}/model_state_dict.pt", flush=True)
        except Exception as exc2:
            print(f"[save] state_dict fallback FAILED: {exc2!r}", flush=True)

    print(f"\n=== DONE ts={ts} run_dir={run_dir} ===", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

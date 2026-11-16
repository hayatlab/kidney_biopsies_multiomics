"""Re-run of the overall (5-region) NicheCompass training, but with the
`in_glom` annotation carried into the output adata.

Mirrors `nichecompass_5regions.py` end-to-end (same input h5ad, same TMA
blacklist, same GPs + hyperparameters, same leiden resolution = 0.4) so
the resulting niches are comparable to the original 22042026_221229
run that produced niche 5 (a.k.a. "niche 9" in earlier numbering). The
only difference: an `in_glom` column is merged into `adata.obs` before
the model is trained, with values pulled from

* the existing parent adata `xenium_with_niches.h5ad`
  (the three original ROI samples already carry `in_glom`)
* `review_analyses/in_glom_for_new_samples/in_glom_new_samples.parquet`
  (the two newly-added Region_1 samples)

Writes a fresh timestamped directory under
`niche_compass_results/artifacts/overall_with_glom/<ts>/` — does NOT
touch the existing single_sample/22042026_221229 run.
"""
from __future__ import annotations
from kbpaths import KB_PROJECT, KB_XENIUM_RAW  # repository path config (config/kbpaths.py)
from cores import excluded_core_labels  # private core lists (config/cores.py)

import argparse
import sys
import traceback
import warnings
from datetime import datetime
from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from anndata import AnnData

import nichecompass as nc
from nichecompass.models import NicheCompass
from nichecompass.utils import (
    add_gps_from_gp_dict_to_adata,
    compute_communication_gp_network,
    create_new_color_dict,
    extract_gp_dict_from_mebocost_ms_interactions,
    extract_gp_dict_from_nichenet_lrt_interactions,
    extract_gp_dict_from_omnipath_lr_interactions,
    filter_and_combine_gp_dict_gps_v2,
    generate_enriched_gp_info_plots,
    visualize_communication_gp_network,
)

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)


# ----- defaults (match nichecompass_5regions.py) -------------------------
INPUT_H5AD = Path(
    f"{KB_XENIUM_RAW}/h5ad_files_only_transcript/"
    "marker_genes/clustered_annotated_final_file_5regions.h5ad"
)
PROJECT_DIR = Path(f"{KB_PROJECT}")
GP_DATA_DIR = PROJECT_DIR / "niche_compass_results/data/gene_programs"
ARTIFACTS_DIR = PROJECT_DIR / "niche_compass_results/artifacts/overall_with_glom"

# Sources for the in_glom annotation
PARENT_H5AD = Path(
    f"{KB_PROJECT}/niche_compass_results/artifacts/"
    "single_sample/22042026_221229/xenium_with_niches.h5ad"
)
NEW_IN_GLOM_PARQUET = PROJECT_DIR / (
    "review_analyses/in_glom_for_new_samples/in_glom_new_samples.parquet"
)

# obs / obsm / obsp / var / layer keys
CELL_TYPE_KEY = "predic"
SAMPLE_KEY = "sample"
CONDITION_KEY = "condition"
TMA_KEY = "TMA_ID"
SPATIAL_KEY = "spatial"
ADJ_KEY = "spatial_connectivities"
COUNTS_KEY = "counts"
IN_GLOM_KEY = "in_glom"

# NicheCompass mask keys
GP_TARGETS_MASK_KEY = "nichecompass_gp_targets"
GP_TARGETS_CATEGORIES_MASK_KEY = "nichecompass_gp_targets_categories"
GP_SOURCES_MASK_KEY = "nichecompass_gp_sources"
GP_SOURCES_CATEGORIES_MASK_KEY = "nichecompass_gp_sources_categories"
GP_NAMES_KEY = "nichecompass_gp_names"
ACTIVE_GP_NAMES_KEY = "nichecompass_active_gp_names"
LATENT_KEY = "nichecompass_latent"

# training hyperparameters (identical to nichecompass_5regions.py)
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
LEIDEN_RES = 0.4

# Cores excluded before NicheCompass: unassigned cells, the Apollo d0 A2 control core and ten disease cores that
# failed QC (Utrecht cores of K4 DN, K24 TIN, K29 MN, K23 TIN, K25 TIN, K71 FSGS, K46 MN, K78 FSGS and K60 IgAN,
# and the Combi D2 IgAN core). Their internal labels contain pathology numbers and are read from a private list
# (config/cores.py); none of these cores is part of the de-identified Zenodo object.
TO_EXCLUDE = set(excluded_core_labels("nichecompass_blacklist"))


def build_in_glom_map() -> pd.Series:
    """Build a Series indexed by `<sample>::<cell_id>` -> bool in_glom,
    combining the existing parent adata (3 original ROI samples) and
    the new Region_1 annotations parquet (2 new samples)."""
    print(f"[in_glom] reading parent {PARENT_H5AD}", flush=True)
    parent = sc.read_h5ad(PARENT_H5AD)
    pkey = (parent.obs[SAMPLE_KEY].astype(str) + "::"
            + parent.obs["cell_id"].astype(str))
    pser = pd.Series(parent.obs[IN_GLOM_KEY].astype(bool).to_numpy(), index=pkey)
    print(f"[in_glom]   parent rows: {len(pser)}  ({int(pser.sum())} True)",
          flush=True)
    del parent

    print(f"[in_glom] reading new-samples parquet {NEW_IN_GLOM_PARQUET}", flush=True)
    annot = pd.read_parquet(NEW_IN_GLOM_PARQUET)
    nkey = annot["sample"].astype(str) + "::" + annot["cell_id"].astype(str)
    nser = pd.Series(annot[IN_GLOM_KEY].astype(bool).to_numpy(), index=nkey)
    print(f"[in_glom]   parquet rows: {len(nser)}  ({int(nser.sum())} True)",
          flush=True)

    # The 2 Region_1 samples appear in `parent` with in_glom == False (no
    # polygons existed yet). The parquet overrides those with the
    # polygon-derived values. Concat then dedup keeping the parquet entry
    # for keys that appear in both.
    combined = pd.concat([pser, nser])
    combined = combined[~combined.index.duplicated(keep="last")]
    print(f"[in_glom] combined unique keys: {len(combined)}  "
          f"({int(combined.sum())} True)", flush=True)
    return combined


def load_adata() -> AnnData:
    print(f"[load] reading {INPUT_H5AD}", flush=True)
    a = sc.read_h5ad(INPUT_H5AD)
    print(f"[load] full: {a.shape}", flush=True)

    if COUNTS_KEY not in a.layers:
        print(f"[load] adding layers['{COUNTS_KEY}'] from .X", flush=True)
        a.layers[COUNTS_KEY] = a.X.copy()

    before = a.n_obs
    mask = (~a.obs[TMA_KEY].astype(str).isin(TO_EXCLUDE)) & a.obs[TMA_KEY].notna() \
        & (a.obs[TMA_KEY].astype(str) != "nan")
    a = a[mask].copy()
    if pd.api.types.is_categorical_dtype(a.obs[TMA_KEY]):
        a.obs[TMA_KEY] = a.obs[TMA_KEY].cat.remove_unused_categories()
    print(f"[load] after TMA filter: {a.shape} (dropped {before - a.n_obs})",
          flush=True)

    nnull = a.obs[CELL_TYPE_KEY].isna().sum()
    if nnull:
        print(f"[load] WARNING: {nnull} cells with null {CELL_TYPE_KEY}, dropping",
              flush=True)
        a = a[a.obs[CELL_TYPE_KEY].notna()].copy()

    print(f"[load] samples: {sorted(a.obs[SAMPLE_KEY].astype(str).unique())}",
          flush=True)

    # merge in in_glom
    in_glom_map = build_in_glom_map()
    obs_key = (a.obs[SAMPLE_KEY].astype(str) + "::"
               + a.obs["cell_id"].astype(str)).to_numpy()
    in_glom_vals = pd.Series(obs_key).map(in_glom_map)
    n_missing = int(in_glom_vals.isna().sum())
    if n_missing:
        print(f"[load] WARN: {n_missing} cells without in_glom annotation -> "
              f"defaulted to False", flush=True)
    a.obs[IN_GLOM_KEY] = in_glom_vals.fillna(False).astype(bool).to_numpy()
    print(f"[load] in_glom True: {int(a.obs[IN_GLOM_KEY].sum())} / {a.n_obs} "
          f"({100*a.obs[IN_GLOM_KEY].mean():.2f}%)", flush=True)
    per_sample = (
        a.obs.groupby(SAMPLE_KEY, observed=True)[IN_GLOM_KEY]
        .agg(["sum", "count"])
        .rename(columns={"sum": "n_in_glom", "count": "n_total"})
    )
    print(f"[load] per-sample in_glom:\n{per_sample}", flush=True)

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


def build_spatial_graph(adata: AnnData) -> None:
    print(f"[graph] spatial_neighbors(n_neighs={N_NEIGHBORS}, coord_type=generic)",
          flush=True)
    sq.gr.spatial_neighbors(
        adata, coord_type="generic", spatial_key=SPATIAL_KEY,
        n_neighs=N_NEIGHBORS,
    )
    adata.obsp[ADJ_KEY] = adata.obsp[ADJ_KEY].maximum(adata.obsp[ADJ_KEY].T)


def fit_nichecompass(adata: AnnData, gp_dict: dict) -> NicheCompass:
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


def report_niche9_match(adata: AnnData, niche_key: str, out_csv: Path) -> None:
    """For each niche, report % of its cells with in_glom == True and
    its share of all in_glom cells. The niche that maximises both metrics
    is the analog of the original niche 9 (a.k.a. niche 5)."""
    ct = (
        adata.obs.groupby(niche_key, observed=True)[IN_GLOM_KEY]
        .agg(n_total="count", n_in_glom="sum")
    )
    ct["pct_in_glom"] = 100.0 * ct["n_in_glom"] / ct["n_total"]
    total_in_glom = int(ct["n_in_glom"].sum())
    ct["share_of_in_glom"] = 100.0 * ct["n_in_glom"] / max(total_in_glom, 1)
    ct = ct.sort_values("pct_in_glom", ascending=False)
    ct.to_csv(out_csv)
    print(f"\n[niche9_match] in_glom enrichment per niche (saved to {out_csv}):",
          flush=True)
    print(ct.to_string(), flush=True)
    if not ct.empty:
        top = ct.index[0]
        print(f"[niche9_match] top niche by pct_in_glom: {top!r} "
              f"({ct.loc[top, 'pct_in_glom']:.1f}% of its cells, "
              f"{ct.loc[top, 'share_of_in_glom']:.1f}% of all in_glom)",
              flush=True)


def main() -> int:
    global LEIDEN_RES, N_EPOCHS
    ap = argparse.ArgumentParser()
    ap.add_argument("--leiden-res", type=float, default=LEIDEN_RES)
    ap.add_argument("--n-epochs", type=int, default=N_EPOCHS)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    LEIDEN_RES = args.leiden_res
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

    adata = load_adata()
    if args.smoke:
        rng = np.random.default_rng(0)
        keep = rng.choice(adata.n_obs, min(10_000, adata.n_obs), replace=False)
        adata = adata[keep].copy()
        N_EPOCHS = 3
        print(f"[smoke] adata -> {adata.shape}, n_epochs={N_EPOCHS}", flush=True)

    gp_dict = build_gp_dict()
    build_spatial_graph(adata)
    model = fit_nichecompass(adata, gp_dict)
    model.save(dir_path=str(model_dir), overwrite=True,
               save_adata=True, adata_file_name="adata.h5ad")
    print(f"[nc] saved model + adata to {model_dir}", flush=True)

    latent_cluster_key = f"latent_leiden_{LEIDEN_RES}"
    sc.pp.neighbors(model.adata, use_rep=LATENT_KEY, key_added=LATENT_KEY)
    sc.tl.umap(model.adata, neighbors_key=LATENT_KEY)
    sc.tl.leiden(adata=model.adata, resolution=LEIDEN_RES,
                 key_added=latent_cluster_key, neighbors_key=LATENT_KEY)
    print(f"[leiden] niches: {model.adata.obs[latent_cluster_key].value_counts().to_dict()}",
          flush=True)

    # niche-9 recovery diagnostic
    report_niche9_match(
        model.adata, latent_cluster_key,
        out_csv=fig_dir / f"niche9_match_res_{LEIDEN_RES}.csv",
    )

    # save adata explicitly with the new niche labels + in_glom so downstream
    # tools can pick it up directly (mirror in_glom_niche layout)
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
        traceback.print_exc()

    print(f"\n=== DONE ts={ts} run_dir={run_dir} ===", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

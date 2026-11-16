#!/usr/bin/env python
"""Map snRNA immune subclusters (leiden_0.6, 18 clusters) onto Xenium immune cells via scANVI
reference mapping / query surgery -- the same scVI->scANVI "scArches" workflow used for PEC
(mapping_pec_xenium.ipynb / mapping_immune_xenium.ipynb).

The `scarches` conda env was deleted, but scArches's `sca.models.SCVI/SCANVI` are just re-exports of
scvi-tools models and query surgery (`load_query_data`) is native to scvi-tools. This script uses
`scvi.model` directly (scvi 1.1.6, env `kpmp_repurposing` -- the env the project already uses for its
scANVI label transfer, review_analyses/celltype_concordance/02_train_transfer.py). Equivalent result.

Changed vs the PEC/immune notebooks (reference + label only):
  reference : subcluster_level2_20260728/Immune/Immune_leiden_sweep.h5ad  (19486 cells)
  label     : leiden_0.6  (18 clusters "0".."17")   [was pec_subcluster_0.80 / immune_subcluster]
  batch     : patient_ID   [sweep object has no donor_id]
  query     : xenium_with_niches.h5ad, cell_type_manual in {Immune, PL/B, B}  (32158 cells)
  output    : review_analyses/immune_mapping_leiden06/prediction/   (fresh dir; old prediction_immune kept)

Representation note (faithful to PEC): reference `.X` is log1p and Xenium query `.X` is raw counts --
fed to scVI/scANVI as-is, exactly as the PEC pipeline did (the immune reference has no raw-counts
layer, so raw-counts-for-both isn't available). Env: kpmp_repurposing. CPU only (no GPUs).

Run:  python map_immune_leiden06_xenium.py            # full (SLURM)
      python map_immune_leiden06_xenium.py --smoke    # quick API/data smoke test
"""
from kbpaths import KB_ANALYSIS, KB_RESULTS  # repository path config (config/kbpaths.py)
import os, argparse
import numpy as np, pandas as pd, scanpy as sc, scvi, torch
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, seaborn as sns

scvi.settings.seed = 42
np.random.seed(42); torch.manual_seed(42)
sc.settings.verbosity = 1

BASE = f"{KB_ANALYSIS}/immune_mapping_leiden06"
condition_key = "patient_ID"        # batch (sweep object has no donor_id)
cell_type_key = "leiden_0.6"        # transferred label -> becomes 'predictions'
ARCH = dict(n_layers=2, encode_covariates=True, deeply_inject_covariates=False,
            use_layer_norm="both", use_batch_norm="none")
SRC = f"{KB_ANALYSIS}/subcluster_level2_20260728/Immune/Immune_leiden_sweep.h5ad"
XEN = f"{KB_RESULTS}/xenium/leon/xenium_with_niches.h5ad"
XEN_KEEP = ["Immune", "PL/B", "B"]


def main(smoke=False):
    EP = 2 if smoke else 100
    QEP = 2 if smoke else 50
    OUT = f"{BASE}/prediction_smoke" if smoke else f"{BASE}/prediction"
    os.makedirs(OUT, exist_ok=True); sc.settings.figdir = OUT

    # ---------- reference (snRNA immune) ----------
    source_adata = sc.read_h5ad(SRC)
    source_adata.obs[cell_type_key] = source_adata.obs[cell_type_key].astype(str).astype("category")
    if smoke:
        sc.pp.subsample(source_adata, n_obs=3000, random_state=42)
        source_adata.obs[cell_type_key] = source_adata.obs[cell_type_key].cat.remove_unused_categories()
    print(f"[ref] {source_adata.shape} | {cell_type_key}: {source_adata.obs[cell_type_key].nunique()} clusters "
          f"| {condition_key}: {source_adata.obs[condition_key].nunique()} batches", flush=True)

    # ---------- query (Xenium immune) ----------
    target_adata = sc.read_h5ad(XEN)
    target_adata.obs["cell_type_manual"] = target_adata.obs["cell_type_manual"].replace(
        {"EC/Pericyte/VSMC": "MC/Pericyte/VSMC"})
    target_adata.obs["cell_type_manual"] = target_adata.obs["cell_type_manual"].astype("category")
    keep = np.isin(target_adata.obs["cell_type_manual"], XEN_KEEP)
    xen_immune = target_adata[keep].copy()          # 32158 x 468, keeps cell_type_manual + cell_id + spatial
    if smoke:
        sc.pp.subsample(xen_immune, n_obs=3000, random_state=42)
    xen_immune.obs["cell_type_manual"] = xen_immune.obs["cell_type_manual"].cat.remove_unused_categories()
    print(f"[query] {xen_immune.shape} | cell_type_manual: {list(xen_immune.obs['cell_type_manual'].cat.categories)}", flush=True)

    # ---------- common genes (Xenium panel) ----------
    common_genes = np.intersect1d(xen_immune.var_names, source_adata.var_names)
    source_adata = source_adata[:, common_genes].copy()
    print(f"[genes] common = {len(common_genes)}", flush=True)

    # ---------- scVI base ----------
    scvi.model.SCVI.setup_anndata(source_adata, batch_key=condition_key)
    vae = scvi.model.SCVI(source_adata, **ARCH)
    vae.train(max_epochs=EP)

    # ---------- scANVI reference model ----------
    scanvae = scvi.model.SCANVI.from_scvi_model(
        vae, unlabeled_category="Unknown", labels_key=cell_type_key, adata=source_adata)
    scanvae.train(max_epochs=EP)
    ref_path = f"{BASE}/scanvi_ref_smoke" if smoke else f"{BASE}/scanvi_ref"
    scanvae.save(ref_path, overwrite=True)

    reference_latent = sc.AnnData(scanvae.get_latent_representation())
    reference_latent.obs["cell_type"] = source_adata.obs[cell_type_key].to_numpy()
    reference_latent.obs["predictions"] = scanvae.predict()
    ref_acc = float(np.mean(reference_latent.obs.predictions.to_numpy() == reference_latent.obs.cell_type.to_numpy()))
    print(f"[ref] scANVI training accuracy = {ref_acc:.4f}", flush=True)

    # ---------- query surgery onto Xenium immune ----------
    humanized_adata = xen_immune[:, source_adata.var_names].copy()
    humanized_adata.obs.rename(columns={"TMA_ID": condition_key}, inplace=True)
    humanized_adata.obs["orig_cell_types"] = humanized_adata.obs["cell_type_manual"].to_numpy()
    humanized_adata.obs[cell_type_key] = scanvae.unlabeled_category_   # all query cells unlabeled

    scvi.model.SCANVI.prepare_query_anndata(humanized_adata, ref_path)
    model = scvi.model.SCANVI.load_query_data(humanized_adata, ref_path)
    model._unlabeled_indices = np.arange(humanized_adata.n_obs)
    model._labeled_indices = []
    print(f"[query] labelled={len(model._labeled_indices)} unlabelled={len(model._unlabeled_indices)}", flush=True)
    model.train(max_epochs=QEP, plan_kwargs=dict(weight_decay=0.0), check_val_every_n_epoch=10)
    if not smoke:
        model.save(f"{BASE}/scanvi_query", overwrite=True)

    query_latent = sc.AnnData(model.get_latent_representation())
    query_latent.obs["cell_type"] = humanized_adata.obs["orig_cell_types"].to_numpy()
    query_latent.obs["batch"] = humanized_adata.obs[condition_key].to_numpy()
    query_latent.obs["predictions"] = model.predict()
    query_latent.obs["pred_confidence"] = model.predict(soft=True).max(axis=1).to_numpy()

    # ---------- outputs ----------
    df = query_latent.obs.groupby(["cell_type", "predictions"], observed=True).size().unstack(fill_value=0)
    (df / df.sum(axis=0)).to_csv(f"{OUT}/resolution_prediction_normalised_df.csv")
    df.to_csv(f"{OUT}/resolution_prediction_df.csv")
    pd.DataFrame(common_genes).to_csv(f"{OUT}/resolution_common_genes.csv")
    query_latent.write_h5ad(f"{OUT}/file_after_predictions.h5ad")
    query_latent.obs.to_csv(f"{OUT}/predictions_incsv.csv")
    print("[predictions] value_counts:\n" +
          query_latent.obs["predictions"].value_counts().sort_index().to_string(), flush=True)
    print(f"[predictions] median confidence = {float(query_latent.obs['pred_confidence'].median()):.3f} | "
          f"frac <0.5 = {float((query_latent.obs['pred_confidence'] < 0.5).mean()):.3f}", flush=True)

    # positional assign onto full-gene Xenium immune object + confusion heatmap
    assert query_latent.n_obs == xen_immune.n_obs, (query_latent.n_obs, xen_immune.n_obs)
    xen_immune.obs["predictions"] = query_latent.obs["predictions"].to_numpy()
    xen_immune.obs["pred_confidence"] = query_latent.obs["pred_confidence"].to_numpy()

    heatmap_data = pd.crosstab(xen_immune.obs["cell_type_manual"], xen_immune.obs["predictions"])
    plt.figure(figsize=(18, 5)); sns.heatmap(heatmap_data, annot=True, fmt="d", cmap="Blues")
    plt.title("Xenium cell_type_manual vs predicted immune leiden_0.6"); plt.tight_layout()
    plt.savefig(f"{OUT}/heatmap_pred_vs_manual_annotation.pdf"); plt.close()

    xen_immune.write_h5ad(f"{OUT}/xenium_immune_with_prediction.h5ad")
    print(f"[done] wrote {OUT}/xenium_immune_with_prediction.h5ad  ({xen_immune.n_obs} cells; "
          f"predictions = {xen_immune.obs['predictions'].nunique()} leiden_0.6 labels)", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    main(**vars(ap.parse_args()))

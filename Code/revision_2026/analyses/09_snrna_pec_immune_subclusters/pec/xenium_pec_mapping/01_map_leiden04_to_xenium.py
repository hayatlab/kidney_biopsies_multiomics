"""Transfer new PEC leiden subcluster labels (snRNA) onto the Xenium PEC compartment
via scArches surgery (SCVI -> SCANVI -> load_query_data).

Parameterised by --label (default leiden_0.4). All outputs are tagged by resolution
so different resolutions never clobber each other (leiden_0.4 -> tag "leiden04",
leiden_0.1 -> "leiden01", ...).

Mirrors review_analyses/celltype_concordance/02_train_transfer.py:
  * reference = the 7,004-cell new PEC subset, labelled by <label>
  * query     = xenium_labeled.h5ad restricted to level2_refined == "PEC" (2,577 cells)
Reference expression = ambient-corrected counts from pec_new_20260717_raw.h5ad, subset
to the 468 Xenium panel genes (Xenium order) and rounded to integers. <label> is joined
from pec_leiden_sweep.h5ad by barcode. Query integer counts come from layers["counts"].

seed 42.  Outputs (in ./, TAG = e.g. leiden04):
  xenium_pec_<TAG>.h5ad                 query PEC cells + pred_<TAG> + conf_<TAG>
  humanized_pec_ref_<TAG>.h5ad          the 468-gene integer reference actually trained on
  pec_<TAG>_predictions.csv             cell_id,sample,TMA_ID,in_glom,pred_<TAG>,conf_<TAG>
  transfer_confidence_summary_<TAG>.csv per predicted cluster: n, mean/median confidence
  models_<TAG>/{scvi_base,scanvi,query}
"""
from __future__ import annotations
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scvi

from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
HERE = Path(f"{KB_ANALYSIS}/pec_resubcluster_20260717/xenium_pec_mapping")  # analysis workspace (outputs live here)
BASE = HERE.parent  # pec_resubcluster_20260717/

PEC_RAW   = BASE / "pec_new_20260717_raw.h5ad"
PEC_SWEEP = BASE / "new_pec_annotation" / "pec_leiden_sweep.h5ad"
XENIUM_LABELED = Path(f"{KB_ANALYSIS}/celltype_concordance/"
                      "outputs/xenium_labeled.h5ad")
SHARED_GENES = Path(f"{KB_ANALYSIS}/celltype_concordance/"
                    "outputs/shared_genes.csv")

SEED = 42
SC_BATCH = "sample"
XEN_BATCH = "sample"
PEC_QUERY_COL = "level2_refined"
PEC_QUERY_VAL = "PEC"
ARCH = dict(n_layers=2, encode_covariates=True, deeply_inject_covariates=False,
            use_layer_norm="both", use_batch_norm="none")


def tag_for(label: str) -> str:
    return "leiden" + label.split("_", 1)[1].replace(".", "")


def load_shared_genes():
    g = pd.read_csv(SHARED_GENES).iloc[:, 0].tolist()
    assert len(g) == 468, len(g)
    return g


def build_reference(genes, label, smoke=False):
    ref = sc.read_h5ad(PEC_RAW)
    lab = sc.read_h5ad(PEC_SWEEP, backed="r").obs[label]
    assert set(ref.obs_names) == set(lab.index), "barcode mismatch raw<->sweep"
    ref.obs[label] = lab.reindex(ref.obs_names).astype(str).astype("category")

    missing = [g for g in genes if g not in ref.var_names]
    assert not missing, f"genes missing from PEC ref: {missing[:10]}"
    ref = ref[:, genes].copy()
    assert list(ref.var_names) == genes

    X = ref.X.tocsr() if hasattr(ref.X, "tocsr") else ref.X
    data = X.data if hasattr(X, "data") else X
    frac_noninteger = float(np.mean(np.abs(data - np.rint(data)) > 1e-6))
    if hasattr(X, "data"):
        r = np.rint(X.data); r[r < 0] = 0
        X.data = r.astype(np.float32); X.eliminate_zeros(); ref.X = X
    else:
        r = np.rint(X); r[r < 0] = 0; ref.X = r.astype(np.float32)
    print(f"[ref] {ref.shape}  frac_nonzero_noninteger_before_round={frac_noninteger:.4f}")
    print(f"[ref] {label} counts:\n{ref.obs[label].value_counts().sort_index().to_string()}")
    print(f"[ref] {SC_BATCH} nunique={ref.obs[SC_BATCH].nunique()}")

    if smoke:
        keep = ref.obs[SC_BATCH].isin(pd.Index(ref.obs[SC_BATCH].unique())[:6])
        ref = ref[keep].copy()
        print(f"[ref/smoke] -> {ref.shape}")
    ref.obs[SC_BATCH] = ref.obs[SC_BATCH].astype(str).astype("category")
    return ref


def prep_query(genes, smoke=False):
    xen = sc.read_h5ad(XENIUM_LABELED)
    keep = (xen.obs[PEC_QUERY_COL].astype(str) == PEC_QUERY_VAL).to_numpy()
    xen = xen[keep].copy()
    print(f"[query] PEC cells: {xen.n_obs}")
    xen = xen[:, genes].copy()
    assert list(xen.var_names) == genes, "Xenium var order must equal reference"
    if "counts" not in xen.layers:
        raise ValueError("xenium_labeled missing layers['counts']")
    xen.X = xen.layers["counts"].copy().astype(np.float32)
    assert XEN_BATCH in xen.obs, f"query missing batch column '{XEN_BATCH}'"
    xen.obs[XEN_BATCH] = xen.obs[XEN_BATCH].astype(str).astype("category")
    if smoke:
        sc.pp.subsample(xen, n_obs=400, random_state=SEED)
        print(f"[query/smoke] -> {xen.n_obs}")
    return xen


def main(label="leiden_0.4", smoke=False):
    tag = tag_for(label)
    pcol, ccol = f"pred_{tag}", f"conf_{tag}"
    print(f"[cfg] label={label}  tag={tag}")
    scvi.settings.seed = SEED
    mdir = HERE / (f"models_{tag}_smoke" if smoke else f"models_{tag}")
    mdir.mkdir(parents=True, exist_ok=True)
    max_epochs = 2 if smoke else None
    q_epochs = 2 if smoke else 40

    genes = load_shared_genes()
    ref = build_reference(genes, label, smoke=smoke)
    if not smoke:
        ref.write_h5ad(HERE / f"humanized_pec_ref_{tag}.h5ad")

    scvi.model.SCVI.setup_anndata(ref, batch_key=SC_BATCH)
    vae = scvi.model.SCVI(ref, **ARCH)
    vae.train(max_epochs=max_epochs, early_stopping=not smoke)
    vae.save(str(mdir / "scvi_base"), overwrite=True)

    scan = scvi.model.SCANVI.from_scvi_model(
        vae, unlabeled_category="Unknown", labels_key=label, adata=ref)
    scan.train(max_epochs=max_epochs, early_stopping=not smoke)
    scan.save(str(mdir / "scanvi"), overwrite=True)

    xen = prep_query(genes, smoke=smoke)
    scvi.model.SCANVI.prepare_query_anndata(xen, str(mdir / "scanvi"))
    q = scvi.model.SCANVI.load_query_data(xen, str(mdir / "scanvi"))
    q.train(max_epochs=q_epochs, plan_kwargs=dict(weight_decay=0.0),
            check_val_every_n_epoch=10)
    q.save(str(mdir / "query"), overwrite=True)

    pred = np.asarray(q.predict())
    conf = q.predict(soft=True).max(axis=1).to_numpy()
    xen.obs[pcol] = pd.Categorical(
        pred, categories=[str(c) for c in sorted({int(x) for x in np.unique(pred)})])
    xen.obs[ccol] = conf

    suffix = "_smoke" if smoke else ""
    out_h5 = HERE / f"xenium_pec_{tag}{suffix}.h5ad"
    xen.write_h5ad(out_h5)
    cols = [c for c in ["cell_id", "sample", "TMA_ID", "in_glom", pcol, ccol] if c in xen.obs.columns]
    xen.obs[cols].to_csv(HERE / f"pec_{tag}_predictions{suffix}.csv", index_label="obs_name")
    g = xen.obs.groupby(pcol, observed=True)[ccol]
    summ = pd.DataFrame({"n": g.size(), "mean_conf": g.mean().round(4),
                         "median_conf": g.median().round(4)})
    summ.to_csv(HERE / f"transfer_confidence_summary_{tag}{suffix}.csv")

    print(f"\nwrote {out_h5}  {xen.shape}")
    print(f"{pcol} distribution:")
    print(xen.obs[pcol].value_counts().sort_index().to_string())
    print(f"\nlow-confidence (<0.5) fraction: {float((conf < 0.5).mean()):.3f}")
    print("per-cluster confidence:\n" + summ.to_string())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", default="leiden_0.4")
    ap.add_argument("--smoke", action="store_true")
    main(**vars(ap.parse_args()))

# review_analyses/celltype_concordance/02_train_transfer.py
"""Train one SCVI base on the humanized reference, derive two SCANVI models
(level1, level2), and transfer onto Xenium via scArches surgery (load_query_data)."""
import argparse, csv
import numpy as np
import scanpy as sc
import scvi
import config as C

ARCH = dict(n_layers=2, encode_covariates=True, deeply_inject_covariates=False,
            use_layer_norm="both", use_batch_norm="none")


def _shared():
    with open(C.SHARED_GENES_CSV) as fh:
        r = csv.reader(fh); next(r)
        return [row[0] for row in r]


def prep_xenium(smoke=False):
    xen = sc.read_h5ad(C.XENIUM)
    genes = _shared()
    xen = xen[:, genes].copy()
    assert list(xen.var_names) == genes
    # use integer counts in .X for the query
    if "counts" not in xen.layers:
        raise ValueError("Xenium object missing layers['counts'] — integer counts required for scANVI")
    xen.X = xen.layers["counts"].copy()
    xen.X = xen.X.astype(np.float32)
    assert C.XEN_BATCH in xen.obs, f"Xenium query missing batch column '{C.XEN_BATCH}'"
    if smoke:
        sc.pp.subsample(xen, n_obs=4000, random_state=C.SEED)
    return xen


def main(smoke=False):
    C.ensure_dirs()
    mdir = (C.MODELS / "smoke") if smoke else C.MODELS
    mdir.mkdir(parents=True, exist_ok=True)
    scvi.settings.seed = C.SEED
    max_epochs = 2 if smoke else None          # None => scvi auto-heuristic (low for big data)
    q_epochs = 2 if smoke else 40

    ref = sc.read_h5ad(C.OUT / ("humanized_ref_smoke.h5ad" if smoke else "humanized_ref.h5ad"))

    # --- SCVI base ---
    scvi.model.SCVI.setup_anndata(ref, batch_key=C.SC_BATCH)
    vae = scvi.model.SCVI(ref, **ARCH)
    vae.train(max_epochs=max_epochs, early_stopping=not smoke)
    vae.save(str(mdir / "scvi_base"), overwrite=True)

    genes = list(ref.var_names)
    xen = prep_xenium(smoke=smoke)
    assert list(xen.var_names) == genes, "Xenium var order must equal reference"

    summary_rows = []
    # NOTE: `xen` is reused across both SCANVI surgeries; safe only because both heads share an
    # identical gene set + batch registry. If a 3rd label level with a different registry is added, copy xen per level.
    for level in C.LABELS:  # ["level1","level2"]
        scan = scvi.model.SCANVI.from_scvi_model(
            vae, unlabeled_category="Unknown", labels_key=level, adata=ref)
        scan.train(max_epochs=max_epochs, early_stopping=not smoke)
        scan.save(str(mdir / f"scanvi_{level}"), overwrite=True)

        # surgery onto Xenium
        scvi.model.SCANVI.prepare_query_anndata(xen, str(mdir / f"scanvi_{level}"))
        q = scvi.model.SCANVI.load_query_data(xen, str(mdir / f"scanvi_{level}"))
        q.train(max_epochs=q_epochs, plan_kwargs=dict(weight_decay=0.0),
                check_val_every_n_epoch=10)
        q.save(str(mdir / f"query_{level}"), overwrite=True)

        pred = q.predict()                       # hard labels
        proba = q.predict(soft=True)             # DataFrame of class probs
        conf = proba.max(axis=1).to_numpy()
        xen.obs[f"pred_{level}"] = np.asarray(pred)
        xen.obs[f"conf_{level}"] = conf

        df = xen.obs.groupby(f"pred_{level}", observed=True)[f"conf_{level}"]
        for cls, g in df:
            summary_rows.append((level, cls, int(g.size),
                                 round(float(g.mean()), 4), round(float(g.median()), 4)))

    conf_path = C.MET / ("transfer_confidence_summary_smoke.csv" if smoke else "transfer_confidence_summary.csv")
    with open(conf_path, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["label_level","pred_class","n","mean_conf","median_conf"])
        w.writerows(summary_rows)

    out = C.XENIUM_LABELED if not smoke else C.OUT / "xenium_labeled_smoke.h5ad"
    xen.write_h5ad(out)
    print("wrote", out, xen.shape)
    print("low-confidence (<0.5) fractions:",
          {lv: float((xen.obs[f'conf_{lv}'] < 0.5).mean()) for lv in C.LABELS})


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    main(**vars(ap.parse_args()))

# review_analyses/celltype_concordance/01_prep_reference.py
"""Humanize the snRNA reference: subset to the 468 Xenium panel genes (Xenium order),
round ambient-corrected counts to integers, report the rounding footprint."""
import argparse, csv
import numpy as np
import scanpy as sc
import config as C


def load_shared_genes():
    with open(C.SHARED_GENES_CSV) as fh:
        r = csv.reader(fh); next(r)
        return [row[0] for row in r]


def main(smoke=False):
    C.ensure_dirs()
    genes = load_shared_genes()
    assert len(genes) == 468

    adata = sc.read_h5ad(C.SC_REF)
    if smoke:
        # preserve composition: stratify-ish by taking first 3 samples
        keep = adata.obs[C.SC_BATCH].isin(adata.obs[C.SC_BATCH].cat.categories[:3])
        adata = adata[keep].copy()

    # subset + reorder to Xenium panel order (master)
    missing = [g for g in genes if g not in adata.var_names]
    assert not missing, f"genes missing from ref: {missing[:10]}"
    adata = adata[:, genes].copy()
    assert list(adata.var_names) == genes

    X = adata.X
    X = X.tocsr() if hasattr(X, "tocsr") else X
    data = X.data
    n_nonzero = data.size
    frac_noninteger = float(np.mean(np.abs(data - np.rint(data)) > 1e-6))
    min_val = float(data.min()) if n_nonzero else 0.0
    libsize_before = np.asarray(X.sum(axis=1)).ravel()

    # round to integers, clip at 0
    rounded = np.rint(data)
    rounded[rounded < 0] = 0
    X.data = rounded.astype(np.float32)
    X.eliminate_zeros()
    adata.X = X
    libsize_after = np.asarray(adata.X.sum(axis=1)).ravel()

    rel_libsize_shift = float(np.median(np.abs(libsize_after - libsize_before)
                                        / np.maximum(libsize_before, 1)))

    rows = [
        ("n_cells", adata.n_obs), ("n_genes", adata.n_vars),
        ("n_nonzero_entries", n_nonzero),
        ("frac_nonzero_noninteger", round(frac_noninteger, 4)),
        ("min_value_before_round", round(min_val, 4)),
        ("median_rel_libsize_shift", round(rel_libsize_shift, 5)),
        ("any_negative_after", int((adata.X.data < 0).any())),
    ]
    met_path = C.MET / ("rounding_report_smoke.csv" if smoke else "rounding_report.csv")
    with open(met_path, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["metric", "value"])
        for k, v in rows:
            w.writerow([k, v])
    print("ROUNDING REPORT:")
    for k, v in rows:
        print(f"  {k}: {v}")

    out = C.HUMANIZED_REF if not smoke else C.OUT / "humanized_ref_smoke.h5ad"
    adata.write_h5ad(out)
    print("wrote", out, adata.shape)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    main(**vars(ap.parse_args()))

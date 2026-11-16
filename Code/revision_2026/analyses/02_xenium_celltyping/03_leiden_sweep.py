"""Recompute neighbors on X_pca_harmony and run a Leiden resolution sweep."""
import argparse
import scanpy as sc
import config as C


def main(smoke=False):
    path = C.OUT / "xenium_labeled_smoke.h5ad" if smoke else C.XENIUM_LABELED
    adata = sc.read_h5ad(path)
    assert "X_pca_harmony" in adata.obsm, "missing X_pca_harmony"
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", random_state=C.SEED)
    res_list = [0.4, 1.0] if smoke else C.LEIDEN_RES
    for res in res_list:
        key = f"leiden_harmony_{res}"
        sc.tl.leiden(adata, resolution=res, key_added=key,
                     flavor="igraph", n_iterations=2, directed=False,
                     random_state=C.SEED)
        print(key, "->", adata.obs[key].nunique(), "clusters")
    adata.write_h5ad(path)
    print("updated", path)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    main(**vars(ap.parse_args()))

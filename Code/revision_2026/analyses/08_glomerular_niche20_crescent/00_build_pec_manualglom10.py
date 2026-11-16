#!/usr/bin/env python
"""Build pec_manualglom10.h5ad: the PEC cells of the manual-glomerulus NicheCompass object (Leiden 1.0), labelled
glomerular niche 20 ("crescent") versus all other glomerular niches. It is the input of de_crescent_pec_mast.R (MAST,
crescent PEC vs other PEC) and of the Fig. 6f dotplot.

Hand-written reconstruction (12 Nov 2026): the script that produced the file used in the paper was not saved. With
--verify-against <file> the result is compared with an existing file cell by cell (only counts and booleans are
printed). Against the paper's file this gives: same 1,711 cells in the same order, identical raw counts, gene order,
obs names, grp and donor columns (values and category order).

Input  : $KB_RESULTS/xenium/leon/new_data/manual_glom_niche_based_on_04_res10_level2/xenium_with_niches.h5ad
         (nichecompass_manual_glom_niche.py, then nichecompass_5regions_downstream.py at Leiden 1.0 with level2_refined,
         wrapper nichecompass_replot_level2/partH_manualglom10_downstream.sh; folder analyses/01_xenium_preprocessing_niches)
Output : --out FILE. The paper's copy is $KB_ANALYSIS/celltype_concordance/outputs/crescent_check/de/pec_manualglom10.h5ad
         (that path is where de_crescent_pec_mast.R and plot_niche20_pec_top5_bottom5.py read it from).
         AnnData 1,711 PEC cells x 468 panel genes, nothing else stored:
           X            raw counts (layer 'counts' of the input, CSR float32)
           obs['grp']   'crescent' if latent_leiden_1.0 == '20' else 'other' (categories ['crescent', 'other'])
           obs['donor'] the input's H.number column; missing values are the string 'nan' (category appended last)
           obs_names    '0' ... '1710' (the input's cell names are not kept)
Caveats: H.number is filled only for the two regions added during review (_origin 'query', 386 PEC). All 1,325 PEC of the
         three original regions, including all 134 niche-20 PEC, have donor 'nan', so de_crescent_pec_mast.R sees one donor
         level in the crescent group and drops the donor term: the MAST contrast is not adjusted for donor or sample.
         The donor values are pathology numbers, so the output file is private: never share it or print its donor column.
         Niche labels are those stored in the input object; a Leiden re-run can renumber niches.
Environment: conda env nichecompass (anndata 0.10.9, pandas 2.3.3, numpy 1.26.4)
Usage  : python 00_build_pec_manualglom10.py --out pec_manualglom10.h5ad [--verify-against <existing pec_manualglom10.h5ad>]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "config"))
from kbpaths import KB_RESULTS  # repository path config (config/kbpaths.py)

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

IN_H5AD = f"{KB_RESULTS}/xenium/leon/new_data/manual_glom_niche_based_on_04_res10_level2/xenium_with_niches.h5ad"
CELLTYPE_KEY, CELLTYPE = "level2_refined", "PEC"
NICHE_KEY, CRESCENT_NICHE = "latent_leiden_1.0", "20"
DONOR_KEY = "H.number"


def build(in_h5ad: str) -> ad.AnnData:
    a = ad.read_h5ad(in_h5ad)
    pec = a[(a.obs[CELLTYPE_KEY].astype(str) == CELLTYPE).values].copy()  # .copy() drops unused categories
    counts = pec.layers["counts"]
    assert sp.issparse(counts) and np.all(np.mod(counts.data, 1) == 0), "layer 'counts' must hold raw integer counts"
    grp = np.where(pec.obs[NICHE_KEY].astype(str).values == CRESCENT_NICHE, "crescent", "other")
    donor = pec.obs[DONOR_KEY].cat.add_categories("nan").fillna("nan")
    obs = pd.DataFrame({"grp": pd.Categorical(grp, categories=["crescent", "other"]),
                        "donor": donor.values},
                       index=pd.Index([str(i) for i in range(pec.n_obs)]))
    out = ad.AnnData(X=sp.csr_matrix(counts, dtype=np.float32), obs=obs, var=pd.DataFrame(index=pec.var_names.copy()))
    return out


def summary(x: ad.AnnData) -> str:
    d = x.obs["donor"].astype(str)
    return (f"{x.n_obs} cells x {x.n_vars} genes; grp {x.obs['grp'].value_counts().to_dict()}; "
            f"{d.nunique()} donor levels ({int((d == 'nan').sum())} cells 'nan'; "
            f"{d[x.obs['grp'] == 'crescent'].nunique()} level(s) in the crescent group)")


def verify(new: ad.AnnData, ref_path: str) -> bool:
    """Compare with an existing file; prints only counts and booleans (never donor values)."""
    ref = ad.read_h5ad(ref_path)
    checks = {
        "shape": new.shape == ref.shape,
        "obs_names": list(new.obs_names) == list(ref.obs_names),
        "var_names": list(new.var_names) == list(ref.var_names),
        "obs columns": list(new.obs.columns) == list(ref.obs.columns),
        "X (raw counts)": new.X.shape == ref.X.shape and (new.X != ref.X).nnz == 0 and new.X.dtype == ref.X.dtype,
        "grp values": bool((new.obs["grp"].astype(str).values == ref.obs["grp"].astype(str).values).all()),
        "grp categories": list(new.obs["grp"].cat.categories) == list(ref.obs["grp"].cat.categories),
        "donor values": bool((new.obs["donor"].astype(str).values == ref.obs["donor"].astype(str).values).all()),
        "donor categories (order)": list(new.obs["donor"].cat.categories) == list(ref.obs["donor"].cat.categories),
        "no layers/obsm/uns": not (ref.layers.keys() or ref.obsm.keys() or ref.uns.keys() or len(ref.var.columns)),
    }
    for k, v in checks.items():
        print(f"[verify] {k:26s} {'OK' if v else 'DIFFERS'}")
    print(f"[verify] reference: {summary(ref)}")
    return all(checks.values())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--in", dest="in_h5ad", default=IN_H5AD, help="manual_glom Leiden-1.0 object (default: %(default)s)")
    ap.add_argument("--out", required=True, help="output h5ad (private: contains pathology numbers in obs['donor'])")
    ap.add_argument("--verify-against", help="existing pec_manualglom10.h5ad to compare with (read only)")
    args = ap.parse_args()
    new = build(args.in_h5ad)
    print(f"[build] {summary(new)}")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    new.write_h5ad(args.out)
    print(f"[build] wrote {args.out}")
    if args.verify_against:  # compare the file as written (round trip included)
        ok = verify(ad.read_h5ad(args.out), args.verify_against)
        print(f"[verify] {'IDENTICAL CONTENT' if ok else 'DIFFERENCES FOUND'}")
        sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

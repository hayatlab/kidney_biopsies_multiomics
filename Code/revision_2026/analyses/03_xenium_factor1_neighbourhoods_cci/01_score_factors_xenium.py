#!/usr/bin/env python
"""Stage A -- score MCFA factors onto the NEW Xenium object (Vivien takeover).

Reproduces Code/plot_egfr_f2_on_spatial/1_preprocess_xenium.py (Vivien Goepp) on:
  Xenium   : celltype_concordance/outputs/xenium_labeled.h5ad  (258,375 x 468, 5 runs, 119 cores)
  loadings : review_analyses/loadings~hvg~5.csv                (5 factors x 10 views, 2026-08-04)

Differences from the original, all deliberate (see PLAN.md):
  * Factor 1 is the factor of interest now, not Factor 2. All five are scored anyway.
  * Views are the 10 new snRNA level1 classes; the old model used 12 different views and
    needed an IC-A/IC-B -> IC special case. The new view set matches level1 exactly, so the
    special case is gone. The level2 -> level1 map is DERIVED from the snRNA object at
    runtime and asserted to be 1:1 -- never hard-coded.
  * The original scored on log-normalised X. The new object's .X is raw integer counts, so
    normalize_total() + log1p is applied here (the original used scanpy defaults:
    target_sum=None -> median library size; that was 193 in the old object).
  * `view` is derived from level2_refined (Leon's chosen annotation, 2026-09-01). Scores are
    ALSO computed with view=pred_level1 as a sensitivity check and the correlation reported.

Human GRCh38, HGNC symbols. Seed 42.
"""
from kbpaths import KB_ANALYSIS, KB_DATA  # repository path config (config/kbpaths.py)
import argparse, json, sys
from pathlib import Path

import anndata as ad
import decoupler as dc
import h5py
import numpy as np
import pandas as pd
import scanpy as sc

HERE = Path(f"{KB_ANALYSIS}/vivien_takeover")
XENIUM = f"{KB_ANALYSIS}/celltype_concordance/outputs/xenium_labeled.h5ad"
SNRNA = f"{KB_DATA}/charlotte_new_data/adata_20260717_with_clinical.h5ad"
LOADINGS = f"{KB_ANALYSIS}/loadings~hvg~5.csv"
FACTORS = [f"Factor{i}" for i in range(1, 6)]
SEED = 42


def level2_to_level1_map():
    """Derive the level2 -> level1 hierarchy from the snRNA object; assert it is 1:1."""
    with h5py.File(SNRNA, "r") as f:
        def cat(c):
            g = f["obs"][c]
            cc = [x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]]
            return pd.Categorical.from_codes(g["codes"][:], cc)
        l1, l2 = pd.Series(cat("level1")), pd.Series(cat("level2"))
    ct = pd.crosstab(l2, l1)
    ambiguous = {lv2: list(ct.columns[ct.loc[lv2] > 0]) for lv2 in ct.index
                 if (ct.loc[lv2] > 0).sum() != 1}
    if ambiguous:
        raise SystemExit(f"level2 -> level1 is not 1:1: {ambiguous}")
    return {lv2: ct.columns[ct.loc[lv2].argmax()] for lv2 in ct.index}


def score_view(adata, net, view, view_col):
    """ULM + manual weighted sum for one view. Mirrors Vivien's score_factor/_wsum."""
    sub = adata[adata.obs[view_col] == view].copy()
    net_v = net[net["view"] == view]
    net_v = net_v[net_v["target"].isin(sub.var_names)]
    if sub.n_obs == 0 or net_v.empty:
        return None, None, 0

    # drop all-zero cells -- ULM is undefined for them (original did the same)
    nz = np.asarray((sub.X != 0).sum(axis=1)).ravel() > 0
    n_empty = int((~nz).sum())
    sub = sub[nz].copy()

    dc.mt.ulm(data=sub, net=net_v, tmin=1, raw=False, verbose=False)
    key = next(k for k in sub.obsm if k.startswith("score_"))
    ulm = pd.DataFrame(sub.obsm[key], index=sub.obs_names)

    # weighted sum: X . W, W aligned on var_names
    wsum = {}
    X = sub.X.toarray() if hasattr(sub.X, "toarray") else np.asarray(sub.X)
    for factor in FACTORS:
        w = net_v[net_v["source"] == factor].set_index("target")["weight"]
        W = np.zeros(sub.n_vars)
        mask = sub.var_names.isin(w.index)
        W[mask] = w.reindex(sub.var_names[mask]).to_numpy()
        wsum[factor] = X.dot(W)
    return ulm, pd.DataFrame(wsum, index=sub.obs_names), n_empty


def score_all(adata, net, view_col, tag):
    views = sorted(set(net["view"]) & set(adata.obs[view_col].astype(str).unique()))
    missing = sorted(set(adata.obs[view_col].astype(str).unique()) - set(net["view"]))
    if missing:
        print(f"[{tag}] views in Xenium but NOT in loadings (cells left NaN): {missing}", flush=True)
    ulms, wsums, empties = [], [], 0
    for v in views:
        u, w, ne = score_view(adata, net, v, view_col)
        empties += ne
        if u is not None:
            ulms.append(u); wsums.append(w)
        print(f"[{tag}] {v}: {0 if u is None else len(u)} cells scored ({ne} all-zero dropped)", flush=True)
    ulm = pd.concat(ulms); wsum = pd.concat(wsums)
    print(f"[{tag}] total {len(ulm)}/{adata.n_obs} cells scored, {empties} all-zero dropped", flush=True)
    return ulm, wsum


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    tag = "_smoke" if args.smoke else ""
    np.random.seed(SEED)
    sc.settings.verbosity = 1

    outdir = HERE / "outputs"; outdir.mkdir(parents=True, exist_ok=True)

    l2l1 = level2_to_level1_map()
    print("[map] level2 -> level1:", json.dumps(l2l1), flush=True)

    print("[load] reading Xenium ...", flush=True)
    adata = sc.read_h5ad(XENIUM)
    print(f"[load] {adata.shape}", flush=True)

    if args.smoke:  # stratified by TMA core, never a flat random subsample
        idx = (adata.obs.groupby("TMA_ID", observed=True)
                        .apply(lambda d: d.sample(min(len(d), 60), random_state=SEED),
                               include_groups=False)
                        .index.get_level_values(-1))
        adata = adata[idx].copy()
        print(f"[smoke] stratified to {adata.shape} ({adata.obs.TMA_ID.nunique()} cores)", flush=True)

    # --- normalise exactly as the original did (scanpy defaults) ---------------
    adata.X = adata.layers["counts"].copy()
    libsize = np.asarray(adata.X.sum(axis=1)).ravel()
    med_lib = float(np.median(libsize))
    sc.pp.normalize_total(adata)          # target_sum=None -> median library size
    sc.pp.log1p(adata)
    print(f"[norm] normalize_total(target_sum=None): median library size = {med_lib:.1f}, "
          f"then log1p (original object was 193.0)", flush=True)

    # --- views ----------------------------------------------------------------
    adata.obs["view"] = adata.obs["level2_refined"].astype(str).map(l2l1).astype("category")
    n_unmapped = int(adata.obs["view"].isna().sum())
    print(f"[view] from level2_refined -> {sorted(adata.obs['view'].cat.categories)} "
          f"({n_unmapped} unmapped)", flush=True)

    net = pd.read_csv(LOADINGS, index_col=0).rename(
        columns={"Factor": "source", "variable": "target", "Loading": "weight"})
    on_panel = net[net["target"].isin(adata.var_names)]
    print(f"[net] {len(net)} loadings, {len(on_panel)} on the 468-gene panel "
          f"({on_panel.groupby('view')['target'].nunique().to_dict()} genes/view)", flush=True)

    ulm, wsum = score_all(adata, net, "view", "level2_refined")
    for factor in FACTORS:
        adata.obs[f"{factor}_ulm"] = adata.obs_names.map(ulm[factor])
        adata.obs[f"{factor}_wsum"] = adata.obs_names.map(wsum[factor])

    # sensitivity: same scoring but views taken straight from pred_level1
    ulm_p1, _ = score_all(adata, net, "pred_level1", "pred_level1")
    for factor in FACTORS:
        adata.obs[f"{factor}_ulm_predlevel1"] = adata.obs_names.map(ulm_p1[factor])

    # per-core mean of the factor of interest (Vivien's Factor2_ulm_mean equivalent)
    for factor in ["Factor1", "Factor2"]:
        m = adata.obs.groupby("TMA_ID", observed=True)[f"{factor}_ulm"].mean()
        adata.obs[f"{factor}_ulm_mean"] = adata.obs["TMA_ID"].map(m).astype(float)

    # --- report ---------------------------------------------------------------
    miss = {f"{f}_ulm": round(float(adata.obs[f"{f}_ulm"].isna().mean() * 100), 3) for f in FACTORS}
    print("[missing %] " + json.dumps(miss), flush=True)
    both = adata.obs[["Factor1_ulm", "Factor1_ulm_predlevel1"]].dropna()
    r = float(both.corr(method="pearson").iloc[0, 1]) if len(both) > 2 else float("nan")
    rs = float(both.corr(method="spearman").iloc[0, 1]) if len(both) > 2 else float("nan")
    print(f"[sensitivity] Factor1_ulm  level2_refined-view vs pred_level1-view: "
          f"pearson r={r:.4f} spearman={rs:.4f} on n={len(both)}", flush=True)
    print("[Factor1_ulm by view]\n" + str(
        adata.obs.groupby("view", observed=True)["Factor1_ulm"].describe()[["count","mean","50%","std"]]), flush=True)

    # --- save -----------------------------------------------------------------
    # NB: use EGFR (float, 2.7% NaN), NOT eGFR_biopsy -- the latter is a decimal-comma
    # STRING column that is 88% empty. supp_egfr_at_time_of_biopsy is the same as EGFR.
    keep = ["cell_id", "sample", "TMA_ID", "region", "condition", "level2_refined",
            "pred_level1", "pred_level2", "predic", "view", "latent_leiden_0.4",
            "in_glom", "EGFR", "supp_egfr_at_time_of_biopsy", "egfr_group",
            "ifta_group", "IFTA_pct", "egfr_slope", "supp_jm_slope",
            "supp_classification", "H.number", "Disease", "total_counts", "cell_area"]
    keep = [c for c in keep if c in adata.obs.columns]
    obs = adata.obs[keep + [c for c in adata.obs.columns if c.startswith("Factor")]].copy()
    if "supp_classification" in obs:
        norm = {"fast progressor": "Fast Progressor", "stable": "stable",
                "improver": "regenerate", "regenerate": "regenerate",
                "not calculable": np.nan, "nan": np.nan, "": np.nan}
        raw = obs["supp_classification"].astype(str).str.strip().str.lower()
        obs["classification"] = raw.map(norm).astype("category")
        print("[classification] raw levels %d -> normalised %s"
              % (raw.nunique(), obs["classification"].value_counts(dropna=False).to_dict()), flush=True)
    obs["spatial_x"] = adata.obsm["spatial"][:, 0]
    obs["spatial_y"] = adata.obsm["spatial"][:, 1]
    obs["sample_cell_id"] = obs["TMA_ID"].astype(str) + "_" + obs["cell_id"].astype(str)
    obs.to_parquet(outdir / f"xenium_factor_scores{tag}.parquet")
    print(f"[save] {outdir / f'xenium_factor_scores{tag}.parquet'}  {obs.shape}", flush=True)

    (outdir / f"01_factor_scoring_meta{tag}.json").write_text(json.dumps({
        "xenium": XENIUM, "loadings": LOADINGS, "snrna_for_hierarchy": SNRNA,
        "factor_of_interest": "Factor1", "view_column": "level2_refined -> level1",
        "level2_to_level1": l2l1, "median_library_size": med_lib,
        "normalisation": "counts -> normalize_total(target_sum=None) -> log1p",
        "n_cells": int(adata.n_obs), "n_cores": int(adata.obs.TMA_ID.nunique()),
        "missing_pct": miss, "sensitivity_pearson_r": r, "sensitivity_spearman": rs,
        "seed": SEED, "smoke": bool(args.smoke),
        "species": "human GRCh38", "gene_id_format": "HGNC symbols",
        "versions": {"scanpy": sc.__version__, "anndata": ad.__version__,
                     "decoupler": dc.__version__, "numpy": np.__version__,
                     "pandas": pd.__version__, "python": sys.version.split()[0]},
    }, indent=2, default=str))
    print("[done]", flush=True)


if __name__ == "__main__":
    main()

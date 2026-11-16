#!/usr/bin/env python
"""NF-kB pathway activity per Xenium cell, for an ED 3f-style niche ribbon (Leon, 2026-09-10).

Leon: "make an additional figure like ed3f_niche_vs_factor1.pdf but with a pathway score for
nf kappa b signaling on the x axis". Score chosen by Leon (2026-09-10): PROGENy NFkB.

  primary      PROGENy NFkB footprint (human, top 500 genes per pathway), decoupler ULM
  sensitivity  CollecTRI NF-kB transcription-factor activity (source "NFKB"), decoupler ULM
  robustness   leave-one-gene-out: PROGENy NFkB rescored without each on-panel gene in turn

Only a slice of either gene set is on the 468-gene panel. The slice is counted and written
out, never assumed, and the score is a panel-restricted estimate. Normalisation mirrors
01_score_factors_xenium.py exactly (raw counts -> normalize_total(target_sum=None) -> log1p),
so the NF-kB axis and the Factor-1 axis come from the same expression values.

Unlike Factor 1 (per-view MCFA loadings), PROGENy is one footprint for all cells, so cell
identity feeds the score (immune cells express TNF, CD40, IRF4, JAK3). The per-cell-type
distribution and the within-cell-type correlation with Factor 1 are written out so that can
be judged.

Networks: compute nodes may have no internet, so `--fetch-nets` (run once on the login node)
caches PROGENy + CollecTRI from OmniPath to outputs/nets/ with the retrieval date; scoring
reads only the cache. All-zero cells stay NaN (ULM is undefined for them), as in 01.
Human GRCh38, HGNC symbols. Seed 42.
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import argparse, json, sys
from datetime import datetime
from pathlib import Path

import anndata as ad
import decoupler as dc
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

HERE = Path(f"{KB_ANALYSIS}/vivien_takeover")
XENIUM = f"{KB_ANALYSIS}/celltype_concordance/outputs/xenium_labeled.h5ad"
OUT = HERE / "outputs"
NETS = OUT / "nets"
SEED, TOP, TMIN = 42, 500, 5
CT_COL = "level2_refined"
PROGENY_SRC, COLLECTRI_SRC = "NFkB", "NFKB"


def fetch_nets():
    NETS.mkdir(parents=True, exist_ok=True)
    pg = dc.op.progeny(organism="human", top=TOP)
    ct = dc.op.collectri(organism="human")
    pg.to_csv(NETS / f"progeny_human_top{TOP}.csv", index=False)
    ct.to_csv(NETS / "collectri_human.csv", index=False)
    (NETS / "nets_meta.json").write_text(json.dumps({
        "retrieved": datetime.now().isoformat(timespec="seconds"), "source": "OmniPath via decoupler",
        "decoupler": dc.__version__, "progeny_top": TOP,
        "progeny_rows": int(len(pg)), "collectri_rows": int(len(ct)),
        "progeny_nfkb_targets": int((pg["source"] == PROGENY_SRC).sum()),
        "collectri_nfkb_targets": int((ct["source"] == COLLECTRI_SRC).sum())}, indent=2))
    print(f"[nets] cached to {NETS}: progeny {len(pg)} rows, collectri {len(ct)} rows", flush=True)


def ulm(adata, net, source):
    n = net.loc[net["source"] == source, ["source", "target", "weight"]]
    dc.mt.ulm(data=adata, net=n, tmin=TMIN, raw=False, verbose=False)
    key = next(k for k in adata.obsm if k.startswith("score_"))
    s = pd.DataFrame(adata.obsm[key], index=adata.obs_names)[source].astype(float)
    for k in [k for k in list(adata.obsm) if k.startswith(("score_", "padj_", "pvals_"))]:
        del adata.obsm[k]
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fetch-nets", action="store_true", help="cache the networks and exit")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    if args.fetch_nets:
        fetch_nets(); return
    tag = "_smoke" if args.smoke else ""
    np.random.seed(SEED)

    pg = pd.read_csv(NETS / f"progeny_human_top{TOP}.csv")
    ct = pd.read_csv(NETS / "collectri_human.csv")
    nets_meta = json.loads((NETS / "nets_meta.json").read_text())

    adata = sc.read_h5ad(XENIUM)
    print(f"[load] {adata.shape}", flush=True)
    if args.smoke:      # stratified by TMA core, as in 01
        idx = (adata.obs.groupby("TMA_ID", observed=True)
                        .apply(lambda d: d.sample(min(len(d), 60), random_state=SEED),
                               include_groups=False).index.get_level_values(-1))
        adata = adata[idx].copy()
        print(f"[smoke] stratified to {adata.shape}", flush=True)

    adata.X = adata.layers["counts"].copy()
    med_lib = float(np.median(np.asarray(adata.X.sum(axis=1)).ravel()))
    sc.pp.normalize_total(adata)                 # target_sum=None -> median library size
    sc.pp.log1p(adata)
    adata.obs_names = (adata.obs["TMA_ID"].astype(str) + "_"
                       + adata.obs["cell_id"].astype(str)).to_numpy()
    nz = np.asarray((adata.X != 0).sum(axis=1)).ravel() > 0
    print(f"[norm] median library size {med_lib:.1f}; {int((~nz).sum())} all-zero cells "
          f"left NaN", flush=True)
    sub = adata[nz].copy()

    panel = set(adata.var_names)
    on = {}
    for name, net, src in [("progeny", pg, PROGENY_SRC), ("collectri", ct, COLLECTRI_SRC)]:
        n = net[net["source"] == src]
        on[name] = n[n["target"].isin(panel)].sort_values("weight", key=np.abs, ascending=False)
        print(f"[net] {name} {src}: {len(n)} targets, {len(on[name])} on the panel "
              f"(+{int((on[name].weight > 0).sum())} / -{int((on[name].weight < 0).sum())})",
              flush=True)
    on["progeny"].to_csv(OUT / f"nfkb_progeny_panel_genes{tag}.csv", index=False)
    print("[net] PROGENy NFkB on panel: " + ", ".join(
        f"{r.target} {r.weight:+.2f}" for r in on["progeny"].itertuples()), flush=True)

    res = pd.DataFrame(index=adata.obs_names)
    res["NFkB_progeny"] = ulm(sub, pg, PROGENY_SRC).reindex(res.index)
    res["NFkB_collectri"] = ulm(sub, ct, COLLECTRI_SRC).reindex(res.index)
    genes = on["progeny"]["target"].tolist()
    res["n_progeny_nfkb_genes_detected"] = np.asarray(
        (adata[:, genes].layers["counts"] > 0).sum(axis=1)).ravel()
    for c in ["TMA_ID", CT_COL]:
        res[c] = adata.obs[c].astype(str).to_numpy()

    # --- robustness: leave each on-panel PROGENy gene out ----------------------
    full = res["NFkB_progeny"].dropna()
    loo = []
    for g in genes:
        s = ulm(sub, pg[pg["target"] != g], PROGENY_SRC).reindex(full.index)
        loo.append({"gene_left_out": g, "spearman_vs_full": float(spearmanr(full, s).statistic)})
    loo = pd.DataFrame(loo).sort_values("spearman_vs_full")
    loo.to_csv(OUT / f"nfkb_progeny_leave_one_gene_out{tag}.csv", index=False)
    print(f"[loo] Spearman(full, leave-one-out) min {loo.spearman_vs_full.min():.3f} "
          f"({loo.gene_left_out.iloc[0]}), median {loo.spearman_vs_full.median():.3f}", flush=True)

    # --- how it relates to Factor 1 and to cell identity -----------------------
    f1 = pd.read_parquet(OUT / "xenium_factor_scores.parquet",
                         columns=["sample_cell_id", "Factor1_ulm"]).set_index("sample_cell_id")
    res["Factor1_ulm"] = f1["Factor1_ulm"].reindex(res.index)
    miss = {c: round(float(res[c].isna().mean() * 100), 3)
            for c in ["NFkB_progeny", "NFkB_collectri", "Factor1_ulm"]}
    print(f"[missing %] {miss}", flush=True)
    both = res.dropna(subset=["NFkB_progeny", "NFkB_collectri", "Factor1_ulm"])
    rho = {"progeny_vs_collectri": float(spearmanr(both.NFkB_progeny, both.NFkB_collectri).statistic),
           "progeny_vs_factor1": float(spearmanr(both.NFkB_progeny, both.Factor1_ulm).statistic)}
    per_ct = both.groupby(CT_COL).apply(lambda d: pd.Series({
        "n_cells": len(d), "median_NFkB_progeny": d.NFkB_progeny.median(),
        "iqr_NFkB_progeny": d.NFkB_progeny.quantile(.75) - d.NFkB_progeny.quantile(.25),
        "pct_no_nfkb_gene_detected": (d.n_progeny_nfkb_genes_detected == 0).mean() * 100,
        "spearman_progeny_vs_factor1": spearmanr(d.NFkB_progeny, d.Factor1_ulm).statistic,
        "spearman_progeny_vs_collectri": spearmanr(d.NFkB_progeny, d.NFkB_collectri).statistic,
    }), include_groups=False).sort_values("median_NFkB_progeny", ascending=False)
    per_ct.to_csv(OUT / f"nfkb_by_celltype{tag}.csv")
    print(f"[rho] {rho}\n[per cell type]\n{per_ct.round(3).to_string()}", flush=True)

    res = res.drop(columns=["Factor1_ulm"]).rename_axis("sample_cell_id").reset_index()
    res.to_parquet(OUT / f"xenium_nfkb_scores{tag}.parquet", index=False)
    print(f"[save] {OUT / f'xenium_nfkb_scores{tag}.parquet'} {res.shape}", flush=True)

    (OUT / f"26_nfkb_meta{tag}.json").write_text(json.dumps({
        "xenium": XENIUM, "nets": nets_meta, "score_primary": f"PROGENy {PROGENY_SRC} top{TOP}, ULM",
        "score_sensitivity": f"CollecTRI {COLLECTRI_SRC}, ULM", "tmin": TMIN,
        "normalisation": "counts -> normalize_total(target_sum=None) -> log1p (as 01)",
        "median_library_size": med_lib, "n_cells": int(adata.n_obs),
        "n_all_zero_cells": int((~nz).sum()), "missing_pct": miss,
        "progeny_targets_on_panel": int(len(on["progeny"])),
        "collectri_targets_on_panel": int(len(on["collectri"])),
        "spearman": rho, "loo_min": float(loo.spearman_vs_full.min()),
        "loo_min_gene": str(loo.gene_left_out.iloc[0]), "seed": SEED, "smoke": bool(args.smoke),
        "species": "human GRCh38", "gene_id_format": "HGNC symbols",
        "versions": {"scanpy": sc.__version__, "anndata": ad.__version__,
                     "decoupler": dc.__version__, "numpy": np.__version__,
                     "pandas": pd.__version__, "python": sys.version.split()[0]},
    }, indent=2, default=str))
    print("[done]", flush=True)


if __name__ == "__main__":
    main()

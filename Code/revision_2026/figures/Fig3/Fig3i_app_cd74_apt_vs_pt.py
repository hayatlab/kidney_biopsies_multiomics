#!/usr/bin/env python
"""Fig. 3i — APP-CD74 communication from aPT vs PT to immune cells: APP in the sender (the more abundant of PT / aPT
in the core) and CD74 in immune cells on an FSGS core with higher eGFR (K82, eGFR 63) and a TIN core with low eGFR
(K18, eGFR 23); right, the local CellPhoneDB v5 APP-CD74 score with immune cells as target averaged per core, aPT vs PT
sender (box plot + one point per core; paired Wilcoxon over the 116 cores with both senders).
Panel        : Main Figures.pdf, Figure 3, panel i
Reads        : $KB_ANALYSIS/vivien_takeover/outputs/local_cci_results/interaction_scores_*.parquet
                 (analyses/03_xenium_factor1_neighbourhoods_cci/08_run_local_cci.py via 08_run_local_cci_array.sh)
               $KB_ANALYSIS/vivien_takeover/outputs/geodataframes/cell_boundaries.parquet (02_build_geodataframes.py)
               $KB_ANALYSIS/vivien_takeover/outputs/authoritative_core_egfr.csv (10_build_authoritative_egfr.py)
               $KB_ANALYSIS/celltype_concordance/outputs/xenium_labeled.h5ad, layers["counts"] (final object,
                 analyses/02_xenium_celltyping/, steps 02 -> 17b)
Writes       : $KB_FIG_OUT/Fig3/Fig3i_app_cd74_apt_vs_pt.pdf (+ .png, .svg, _meta.json)
Plotting code: analyses/03_xenium_factor1_neighbourhoods_cci/34_cci_panel.py::draw_pair (imported unchanged); the data
               assembly is extracted verbatim from 34_cci_panel.py::main, for the pair APP_CD74 only (original run:
               34b_cci_top_pairs.sh = --pairs <18 pairs, APP_CD74 first> --outsub cci_lr_panels, SLURM job 330299,
               2026-09-28 -> fig3j_app_cd74_apt_pt; it supersedes 29_fig3j_app_cd74.py). The two map cores were given
               as TMA_ID prefixes ("Utrecht2 F5", "Utrecht2 F2"); here they are selected by public id (XEN119 = K82,
               XEN116 = K18). The per-core score table of main() (rows = internal core labels) is not written.
Illustrator  : map titles "FSGS, eGFR 63" / "TIN, eGFR 23" rewritten "FSGS with higher eGFR" / "TIN with low eGFR" and
               the gene names coloured; box-plot title "(average per core, 116 cores)" relabelled "(average per
               sample)"; "(aPT higher)" removed from the significance label; colour bars and legend moved
Environment  : conda env spatialdata
Verified     : 2026-11-12, PIXEL-IDENTICAL vs $KB_ANALYSIS/vivien_takeover/figures/xenium/cci_lr_panels/fig3j_app_cd74_apt_pt.pdf
               (same test: 116 cores, median paired difference 0.655, bootstrap 95% CI 0.131-1.016, Wilcoxon p = 0.0022)
"""
import json
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from xenium_figs import load_script  # noqa: E402
from kbpaths import fig_out  # noqa: E402
from cores import core_label  # noqa: E402

PAIR = "APP_CD74"
SENDERS = ["aPT", "PT"]
MAP_CORES = ["XEN119", "XEN116"]   # left: K82 FSGS (eGFR 63); right: K18 TIN (eGFR 23)

cc = load_script("03_xenium_factor1_neighbourhoods_cci/34_cci_panel.py")
xf = load_script("03_xenium_factor1_neighbourhoods_cci/05_xenium_figures.py")


def main():
    tag = ""
    rng = np.random.default_rng(cc.SEED)
    out = Path(fig_out("Fig3", "Fig3i_app_cd74_apt_vs_pt.pdf"))
    fig_dir = out_dir = out.parent
    stem = out.stem

    # ---- extracted from 34_cci_panel.py::main (one pair; same order of random draws) ----
    gdf = gpd.read_parquet(cc.HERE / f"outputs/geodataframes/cell_boundaries{tag}.parquet")
    egfr, _, _, _ = xf.corrected_core_egfr(gdf[["TMA_ID", "EGFR"]])
    egfr.index = egfr.index.astype(str)
    cond = gdf.groupby("TMA_ID", observed=True)["condition"].first().astype(str)
    cells = gdf.groupby("TMA_ID", observed=True).size()
    n_type = gdf.groupby(["TMA_ID", "level2_refined"], observed=True).size().unstack(fill_value=0)

    pairs = [PAIR]
    d = cc.load_scores(pairs, SENDERS, False)
    lr = d.groupby("interacting_pair")[["gene_a", "gene_b"]].first()

    picks = [core_label(x) for x in MAP_CORES]
    missing = [x for x, c in zip(MAP_CORES, picks) if c not in set(gdf["TMA_ID"].astype(str))]
    if missing:
        raise SystemExit(f"map cores not in the object: {missing}")
    for c in picks:
        print(f"[core] {cond.get(c)}  eGFR {float(egfr.get(c, np.nan)):.0f}  {cells[c]} cells  "
              f"aPT {n_type.loc[c, 'aPT']}  PT {n_type.loc[c, 'PT']}  "
              f"{cc.TARGET} {n_type.loc[c, cc.TARGET]}", flush=True)

    genes = sorted({g for pr in pairs for g in (lr.loc[pr, "gene_a"], lr.loc[pr, "gene_b"])})
    import anndata as ad
    a = ad.read_h5ad(cc.XENIUM)
    X = a.layers["counts"] if "counts" in a.layers else a.X
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
    lib = np.asarray(X.sum(axis=1)).ravel()
    key = a.obs["TMA_ID"].astype(str).str.strip() + "_" + a.obs["cell_id"].astype(str)
    expr = pd.DataFrame({g: np.log1p(np.asarray(X[:, a.var_names.get_loc(g)].todense()).ravel()
                                     / np.where(lib > 0, lib, 1) * 1e4) for g in genes},
                        index=key.to_numpy())
    del a, X
    joined = gdf.join(expr, how="left")
    frac = joined[genes[0]].notna().mean() * 100
    print(f"[expr] {len(genes)} genes joined onto {frac:.1f} % of cells", flush=True)
    assert frac > 95, "cell ids of the boundaries and the h5ad do not match"

    lig, rec = lr.loc[PAIR, "gene_a"], lr.loc[PAIR, "gene_b"]
    dp = d[d.interacting_pair == PAIR]
    per = dp.groupby(["core", "source"]).score.mean().reset_index()
    per["disease"] = per["core"].map(cond)
    per["egfr"] = per["core"].map(egfr)
    wide = per.pivot(index="core", columns="source", values="score")
    tests = [cc.paired_stats(wide, SENDERS[0], b, rng) for b in SENDERS[1:] if b in wide]
    print(f"\n=== {PAIR} ({lig} -> {rec}) ===\n{json.dumps(tests, indent=1)}", flush=True)
    cc.draw_pair(PAIR, lig, rec, SENDERS, per, tests, wide, joined, picks, gdf, egfr, cond, cells, n_type,
                 rng, fig_dir, out_dir, stem, tag)


if __name__ == "__main__":
    main()

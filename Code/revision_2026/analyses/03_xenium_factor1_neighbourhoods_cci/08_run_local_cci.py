#!/usr/bin/env python
"""Stage E.1 -- local CellPhoneDB CCI within 2-hop Voronoi neighbourhoods (Vivien takeover).

Reproduces the first-submission script spatial_cci/1_run_local_cci.py (no longer in this repository) on the
new 5-run Xenium object.
For every root cell, its 2-hop neighbourhood is treated as a mini tissue and
`cpdb_statistical_analysis_complex_method` is run on it (200 permutations, as in the
original's `work()`), giving a ligand-receptor score and p-value per cell-type pair.

Inputs (all produced by earlier stages):
  outputs/neighborhoods/voronoi_graphs_10.pkl   Stage C -- distance-bounded Voronoi graphs
  outputs/geodataframes/cell_boundaries.parquet Stage B -- TMA_ID + level2_refined per cell
  outputs/cpdb_db/cellphonedb.zip               CellPhoneDB v5.0.0 (re-downloaded; Vivien's
                                                copy sat under his locked /data home)
  the Xenium h5ad                               raw integer counts

Deviations from the original, all deliberate:
  * cell types are `level2_refined` (16), not `cell_type_manual` (21)  [Leon, 2026-09-01]
  * 258,375 roots instead of 198,361
  * `sample_cell_id` is built once in Stage B and reused, so the root ids match the graphs
  * chunk bounds are passed explicitly per array task and stamped into the output filename,
    so a partial run is resumable and chunks can never silently overlap

Human GRCh38, HGNC symbols. Seed 42 (CellPhoneDB permutation seed is its own; set via
--seed so a rerun is reproducible).
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import argparse, contextlib, io, json, logging, os, pickle, sys, time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(f"{KB_ANALYSIS}/vivien_takeover")
XENIUM = f"{KB_ANALYSIS}/celltype_concordance/outputs/xenium_labeled.h5ad"
CPDB_ZIP = HERE / "outputs" / "cpdb_db" / "cellphonedb.zip"
CT_COL = "level2_refined"
SEED = 42

ID_VARS = ['id_cp_interaction', 'interacting_pair', 'partner_a', 'partner_b',
           'gene_a', 'gene_b', 'secreted', 'receptor_a', 'receptor_b',
           'annotation_strategy', 'is_integrin', 'directionality', 'classification', 'root']


def melt_result(df, root, value_name):
    df = df.dropna(subset=["gene_a", "gene_b"]).copy()
    df["root"] = root
    idv = [c for c in ID_VARS if c in df.columns]
    out = df.melt(id_vars=idv, var_name="interaction", value_name=value_name)
    parts = out["interaction"].str.split("|", n=1, expand=True)
    if parts.shape[1] != 2:
        return None
    out["source"], out["target"] = parts[0], parts[1]
    return out.drop(columns=["interaction"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--end", type=int, default=-1, help="-1 = all roots")
    ap.add_argument("--iterations", type=int, default=200)   # the original's work()
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    t0 = time.time()

    import scanpy as sc
    from joblib import Parallel, delayed
    # CellPhoneDB prints a per-permutation percentage to stdout for EVERY root; with 258k
    # roots that is millions of lines. Silence its logger and swallow its stdout per call.
    logging.getLogger("cellphonedb").setLevel(logging.ERROR)
    for name in list(logging.root.manager.loggerDict):
        if "cellphonedb" in name.lower() or name.upper() == "CORE":
            logging.getLogger(name).setLevel(logging.ERROR)
    from cellphonedb.src.core.methods import cpdb_statistical_analysis_complex_method as cpdb
    from cellphonedb.utils import db_utils

    outdir = HERE / "outputs" / "local_cci_results"; outdir.mkdir(parents=True, exist_ok=True)

    # --- neighbourhoods -------------------------------------------------------
    with open(HERE / "outputs" / "neighborhoods" / "voronoi_graphs_10.pkl", "rb") as f:
        graphs = pickle.load(f)
    import networkx as nx
    adj = {}
    for core, g in graphs.items():
        for n in g.nodes():
            nb = list(nx.single_source_shortest_path_length(g, n, cutoff=2).keys())
            adj[n] = [x for x in nb if x != n]
    roots = sorted(adj)
    print(f"[nh] {len(roots)} roots across {len(graphs)} cores "
          f"(median 2-hop size {int(np.median([len(v) for v in adj.values()]))+1}) "
          f"({time.time()-t0:.0f}s)", flush=True)

    # --- counts + metadata ----------------------------------------------------
    gdf = pd.read_parquet(HERE / "outputs" / "geodataframes" / "cell_boundaries.parquet",
                          columns=["sample_cell_id", "cell_id", CT_COL, "TMA_ID"])
    adata = sc.read_h5ad(XENIUM)
    adata.obs["sample_cell_id"] = (adata.obs["TMA_ID"].astype(str) + "_"
                                   + adata.obs["cell_id"].astype(str))
    adata.obs[CT_COL] = adata.obs[CT_COL].astype(str)
    adata = adata[adata.obs["sample_cell_id"].isin(set(gdf["sample_cell_id"]))].copy()
    adata.X = adata.layers["counts"]          # CellPhoneDB wants raw integer counts
    print(f"[counts] {adata.shape}", flush=True)

    interactions, genes, complex_composition, complex_expanded, _, receptor2tfs = \
        db_utils.get_interactions_genes_complex(str(CPDB_ZIP))
    idmap = genes.set_index("hgnc_symbol")["id_multidata"].to_dict()
    counts_relations = genes[["id_multidata", "ensembl", "gene_name", "hgnc_symbol"]].drop_duplicates()

    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    counts_full = pd.DataFrame(X.T.astype(np.float64),
                               index=adata.var_names.astype(str),
                               columns=adata.obs["sample_cell_id"].astype(str).values)
    db_genes = set(genes["hgnc_symbol"].dropna().astype(str))
    counts_full = counts_full[counts_full.index.isin(db_genes)]
    counts_full.index = counts_full.index.to_series().map(idmap)
    print(f"[db] {counts_full.shape[0]} of the {adata.n_vars} panel genes are in CellPhoneDB "
          f"v5.0.0", flush=True)

    meta = adata.obs[["sample_cell_id", CT_COL]].copy()
    meta.columns = ["Cell", "cell_type"]
    meta = meta.set_index("Cell")

    def run_one(root):
        members = [root] + adj[root]
        m = meta.loc[meta.index.isin(members)].copy()
        if m.shape[0] <= 1:
            return None
        m["cell_type"] = m["cell_type"].astype("category").cat.remove_unused_categories()
        if m["cell_type"].isnull().any():
            return None
        cs = counts_full.loc[:, m.index.values]
        cs = pd.DataFrame(np.ascontiguousarray(cs.values, dtype=np.float64),
                          index=cs.index, columns=cs.columns)
        micro = pd.DataFrame({"cell_type": m["cell_type"], "microenvironment": "all"})
        try:
            res = cpdb.call(
                meta=m.copy(), counts=cs, counts_relations=counts_relations,
                counts_data="hgnc_symbol", active_tf2cell_types=None,
                interactions=interactions, genes=genes, complexes=complex_expanded,
                complex_compositions=complex_composition, microenvs=micro,
                receptor2tfs=receptor2tfs, pvalue=0.05, separator="|",
                iterations=args.iterations, threshold=0.0, threads=1,
                debug_seed=args.seed, result_precision=3, debug=False, output_path=None)
        except Exception as e:
            print(f"[warn] root {root} failed: {type(e).__name__}: {e}", flush=True)
            return None
        if not res or res.get("means") is None:
            return None
        s = melt_result(res["means"], root, "score")
        p = melt_result(res["pvalues"], root, "pvalue")
        return (s, p) if s is not None and p is not None else None

    end = len(roots) if args.end < 0 else min(args.end, len(roots))
    sel = roots[args.start:end]
    if args.smoke:
        sel = sel[:12]
    n_jobs = args.threads or int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count()))
    print(f"[run] roots {args.start}..{end} ({len(sel)} to process), "
          f"{args.iterations} permutations, n_jobs={n_jobs}", flush=True)

    t1 = time.time()
    # one redirect around the whole block: CellPhoneDB writes a per-permutation percentage
    # to stdout for every root, which is millions of lines over 258k roots
    _sink = io.StringIO()
    with contextlib.redirect_stdout(_sink):
        out = Parallel(n_jobs=n_jobs, prefer="threads", verbose=0)(
            delayed(run_one)(r) for r in sel)
    for line in _sink.getvalue().splitlines():
        if "[warn]" in line:          # keep genuine per-root failures
            print(line, flush=True)
    ok = [r for r in out if r is not None]
    dt = time.time() - t1
    print(f"[run] {len(ok)}/{len(sel)} roots produced results in {dt:.0f}s "
          f"({dt/max(len(sel),1):.2f}s per root, wall)", flush=True)
    if not ok:
        raise SystemExit("no roots produced results")

    tag = "_smoke" if args.smoke else f"_{args.start}_{end}"
    sdf = pd.concat([r[0] for r in ok], ignore_index=True)
    pdf = pd.concat([r[1] for r in ok], ignore_index=True)
    sdf.to_parquet(outdir / f"interaction_scores{tag}.parquet")
    pdf.to_parquet(outdir / f"interaction_pvalues{tag}.parquet")
    print(f"[save] scores {sdf.shape} / pvalues {pdf.shape} -> {outdir}", flush=True)

    (outdir / f"08_cci_meta{tag}.json").write_text(json.dumps({
        "start": args.start, "end": end, "n_roots_selected": len(sel),
        "n_roots_with_results": len(ok), "iterations": args.iterations,
        "cell_type_col": CT_COL, "cpdb_version": "v5.0.0", "seed": args.seed,
        "seconds_per_root_wall": dt / max(len(sel), 1), "n_jobs": n_jobs,
        "genes_in_db": int(counts_full.shape[0]),
    }, indent=2))
    print(f"[done] {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

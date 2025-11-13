#
# Created on Sun Sep 21 2025
#
# Copyright (c) 2025 Vivien Goepp
#

# --------------------------------------------------------------------------------------------
# This script performs **local CellPhoneDB-based cell–cell interaction (CCI) analysis**
# within 2-hop neighborhoods derived from Xenium kidney biopsy data.
#
# Datasets used/loaded:
# - Cell, nucleus, and circular boundaries:
#       geodataframes/{cell_boundaries, nucleus_boundaries, cell_circles}.parquet
# - Processed Xenium transcriptomic data:
#       xenium_processed.h5ad
# - Precomputed Voronoi adjacency graphs:
#       neighborhoods/voronoi_graphs_10.pkl
#
# Datasets generated/saved:
# - Interaction means and p-values (per 2-hop root neighborhood):
#       interaction_scores_<start>_<end>.parquet
#       interaction_pvalues_<start>_<end>.parquet
#
# External dependencies:
# - CellPhoneDB database v5.0.0 (downloaded locally via `db_utils.download_database`)
#
# The analysis pipeline:
# 1. Load Xenium spatial and transcriptomic data, filtering to annotated tissue cores.
# 2. Build adjacency and 2-hop neighborhood graphs from Voronoi connectivity.
# 3. Prepare gene–protein mapping and count matrices compatible with CellPhoneDB.
# 4. For each root cell (and its 2-hop neighborhood):
#      - Extract the corresponding submatrix of counts and metadata.
#      - Map genes to `id_multidata` and compute ligand–receptor interactions using
#        `cpdb_statistical_analysis_complex_method`.
#      - Extract both mean interaction strengths and p-values, expanding to long format.
# 5. Run per-root analysis in parallel (with optional SLURM array integration).
# 6. Save consolidated parquet tables of interaction scores and p-values for downstream use.
#
# This pipeline produces spatially resolved ligand–receptor networks reflecting
# **local intercellular communication** patterns in kidney biopsies, at single-cell
# resolution, across thousands of 2-hop microenvironments.
# --------------------------------------------------------------------------------------------


import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import sys
import seaborn as sns
import os
import sys
import networkx as nx
import numpy as np
import geopandas as gpd
from matplotlib.colors import ListedColormap
import pickle
import warnings
from joblib import Parallel, delayed
from cellphonedb.src.core.methods import cpdb_statistical_analysis_complex_method
from cellphonedb.utils import db_utils
import argparse
warnings.filterwarnings(
    "ignore", message="nopython is set for njit and is ignored")


# Set directories ---------------------------
data_dir = "/data/zd775156/kidney_biopsy/data/"
xenium_path = "/data/Kramann_lab/From_Sidrah/new_data/output-XETG00229__0042348__Region_3__20250416__104148/"
res_dir = "/data/Kramann_lab/kidney_biopsy/data/processed/xenium/"
cci_res_dir = "/data/Kramann_lab/From_Sidrah/new_data/only_glom_roi/cci_nongrp/csv_output_fiels_all_cells_plotting/"
fig_dir = "/data/Kramann_lab/From_Vivien/kidney_biopsy/spatial_cci/local_cci_results/"

# Load helper functions ----------------------------
sys.path.append(os.path.abspath("../utils/"))
from running_on_slurm import running_on_slurm

# -------------------------------------------------------------------------------
# Load and prepare data ---------------------------------------------------------
# -------------------------------------------------------------------------------

# Load xenium with boundaries ---------------------------------------------------
gdf_cell_boundaries = gpd.read_parquet(
    res_dir+"geodataframes/cell_boundaries.parquet")
gdf_nucleus_boundaries = gpd.read_parquet(
    res_dir+"geodataframes/nucleus_boundaries.parquet")
gdf_cell_circles = gpd.read_parquet(
    res_dir+"geodataframes/cell_circles.parquet")

cores = gdf_cell_boundaries["TMA_ID"].unique().tolist()

cell_types = gdf_cell_boundaries["cell_type_manual"].cat.categories.tolist()

# Create custom palette for cell_types ---------------------------------
colors_tab20 = plt.get_cmap("tab20")(range(20))
colors_husl = sns.color_palette('husl', len(cell_types) - 20)
my_pal = list(colors_tab20) + colors_husl
my_cmap = ListedColormap(my_pal, name="my_custom")


# Load adata ---------------------------------------------------
adata = sc.read_h5ad(res_dir + "xenium_processed.h5ad")

# Set .X to raw counts ---------------------------------------------------
adata.X = adata.layers["counts"]


adata.obs["sample_cell_id"] = adata.obs["TMA_ID"].astype(
    str) + "_" + adata.obs["cell_id"].astype(str)

# Print cores in adata.obs["TMA_ID"] and not in cores -----
set(adata.obs["TMA_ID"].unique().tolist()) - set(cores)
# {'Utrecht 2 A5 T03-04110 K29 MN',
#  'Utrecht1 B5 T18-01374 K25 TIN',
#  'Utrecht2 C3 T17-11176 K24 TIN'}

# Remove them from adata ------
adata = adata[adata.obs["TMA_ID"].isin(cores), :].copy()


# Load Voronoi graph ---------------------------------------------------
max_threshold = 10
with open(f"{res_dir}neighborhoods/voronoi_graphs_{max_threshold}.pkl", "rb") as f:
    voronoi_graphs = pickle.load(f)

cores = list(voronoi_graphs.keys())

# Create adj list ---------------------------------------------------
adj_list = {
    core: {n: list(voronoi_graphs[core].neighbors(n))
           for n in voronoi_graphs[core].nodes()}
    for core in cores
}

# Compute adj list of two hop neighbors (excluding self) ------------------------------
adj_2hop_list = {
    core: {
        n: list(
            set(nx.single_source_shortest_path_length(
                voronoi_graphs[core], n, cutoff=2).keys()) - {n}
        )
        for n in voronoi_graphs[core].nodes()
    }
    for core in cores
}

# Flatten adj_2hop_list to a single dict {root: neighbors} ------------------------------
adj_2hop_dict = {
    root: neighbors
    for core in adj_2hop_list.keys()
    for root, neighbors in adj_2hop_list[core].items()
}


# Test that adj_list is subset of adj_2hop_list ------------------------------
for core in cores:
    assert all([set(adj_list[core][n]).issubset(set(adj_2hop_list[core][n]))
                for n in adj_list[core].keys()])


# Download cellphoneDB database -----------------------------------------------------------
cpdb_version = "v5.0.0"
cpdb_target_dir = "cci_tmp_files/"
db_utils.download_database(cpdb_target_dir, cpdb_version)


# -------------------------------------------------------------------------------
# Setup variables for cellphone run ---------------------------------------------
# -------------------------------------------------------------------------------


# Set root list ---
root = "Utrecht2 F2 T05-15961 K18 TIN_oepkbcng-1"
core = "Utrecht2 F2 T05-15961 K18 TIN"

roots = list(adj_2hop_dict.keys())

# cellphonedb imports (same modules used in shipped code)
cpdb_file_path = "/data/zd775156/kidney_biopsy/code/cci_utils/cci_tmp_files/cellphonedb.zip"

# 1) preload DB once
interactions, genes, complex_composition, complex_expanded, gene_synonym2gene_name, receptor2tfs = \
    db_utils.get_interactions_genes_complex(cpdb_file_path)
microenvs = pd.DataFrame(columns=['cell_type'])  # Include all cell types in the analysis

# 2) build a mapping counts_data -> id_multidata (do once)
# counts_data is e.g. 'hgnc_symbol' or 'ensembl'
idmap = genes.set_index("hgnc_symbol")['id_multidata'].to_dict()

# 3) build a single in-memory counts matrix: rows = gene IDs (counts_data), cols = sample_cell_id
#    (avoid creating many small AnnData objects)
X = adata.X.astype(np.int64)
if hasattr(X, "toarray"):
    X = X.toarray()
gene_ids = adata.var.index.astype(str).values
cell_ids = adata.obs['sample_cell_id'].astype(str).values
counts_full = pd.DataFrame(X.T.astype(np.int64), index=gene_ids, columns=cell_ids)

# optional: filter to genes present in DB to reduce size
db_genes = set(genes["hgnc_symbol"].dropna().astype(str).unique())
counts_full = counts_full[counts_full.index.isin(db_genes)]

# helper to create counts_relations (simple, used by complex_method)
counts_relations = genes[['id_multidata', 'ensembl', 'gene_name', 'hgnc_symbol']].drop_duplicates()

meta_df = adata.obs.loc[:, ['sample_cell_id', 'cell_type_manual']].copy()
meta_df.columns = ['Cell', 'cell_type']
meta_df = meta_df.set_index("Cell")

counts_full = counts_full.astype(np.float64)


# 4) fast per-neighborhood processing using in-memory DataFrames (no disk)


def run_one_root(root, neighbors, meta_df, method="liana_cpdb", **cpdb_kwargs):

    print(f"--------------- {root} ---------------")

    # make meta DataFrame expected by CPDB
    members = [root] + neighbors
    meta_df_sub = meta_df.loc[meta_df.index.isin(members), :].copy()

    if meta_df_sub.shape[0] <= 1:
        return None

    # subset counts quickly (no copy of AnnData)
    counts_sub = counts_full.loc[:, meta_df_sub.index.values]

    counts_sub = pd.DataFrame(
        np.ascontiguousarray(counts_sub.values, dtype=np.float64),
        index=counts_sub.index, columns=counts_sub.columns
    )

    # ensure categorical cell_type has no NA and only integer codes
    meta_df_sub['cell_type'] = meta_df_sub['cell_type'].astype('category')
    if meta_df_sub['cell_type'].isnull().any():
        raise ValueError("NA in cell_type after subsetting")
    meta_df_sub['cell_type'] = meta_df_sub['cell_type'].cat.remove_unused_categories()

    # map gene ids -> id_multidata
    mapped = counts_sub.index.to_series().map(idmap)
    if mapped.isna().sum() > 0:
        raise ValueError(f"{mapped.isna().sum()} genes were not mapped")

    counts_sub.index = mapped
    counts_sub = counts_sub[counts_sub.index.notna()]

    # Select only cells in counts
    meta_df_sub['cell_type'] = meta_df_sub['cell_type'].cat.remove_unused_categories()

    microenvs = pd.DataFrame({
        'cell_type': meta_df_sub['cell_type'],
        'microenvironment': 'all'  # all cells in one environment
    })

    analysis_result = cpdb_statistical_analysis_complex_method.call(
        meta=meta_df_sub.copy(),
        counts=counts_sub,
        counts_relations=counts_relations,
        counts_data="hgnc_symbol",
        active_tf2cell_types=None,
        interactions=interactions,
        genes=genes,
        complexes=complex_expanded,
        complex_compositions=complex_composition,
        microenvs=microenvs,
        receptor2tfs=receptor2tfs,
        pvalue=cpdb_kwargs.get("pvalue", 0.05),
        separator="|",
        iterations=cpdb_kwargs.get("iterations", 1000),
        threshold=cpdb_kwargs.get("threshold", 0.0),
        threads=cpdb_kwargs.get("threads", 1),
        debug_seed=-1,
        result_precision=cpdb_kwargs.get("result_precision", 3),
        debug=False,
        output_path=None
    )

    # Postprocessing of results
    if 'means' in analysis_result and analysis_result['means'] is not None:

        # Extract and pivot longer the means DataFrame ---
        df_score = analysis_result["means"].dropna(
            subset=["gene_a", "gene_b"]).copy()
        df_score['root'] = root

        df_score = df_score.melt(
            id_vars=['id_cp_interaction', 'interacting_pair', 'partner_a', 'partner_b',
                     'gene_a', 'gene_b', 'secreted', 'receptor_a', 'receptor_b',
                     'annotation_strategy', 'is_integrin', 'directionality', 'classification', 'root'],
            var_name='interaction',
            value_name='score'
        )

        split_cols = df_score['interaction'].str.split('|', n=1, expand=True)
        if split_cols.shape[1] == 0:
            print(
                f"Warning: interaction column could not be split for root {root}. Result: {split_cols}")
            return None
        elif split_cols.shape[1] > 2:
            print(
                f"Warning: interaction column split into more than 2 parts for root {root}. Result: {split_cols}")
        df_score['source'] = split_cols[0]
        df_score['target'] = split_cols[1]
        df_score = df_score.drop(columns=['interaction'])

        # df_score[['source', 'target']] = pd.DataFrame(
        #     df_score['interaction'].str.split('|', n=1).to_list(),
        #     index=df_score.index
        # )

        # Extract and pivot longer the pvalues DataFrame -----
        df_pval = analysis_result["pvalues"].dropna(
            subset=["gene_a", "gene_b"]).copy()
        df_pval['root'] = root

        df_pval = df_pval.melt(
            id_vars=['id_cp_interaction', 'interacting_pair', 'partner_a', 'partner_b',
                     'gene_a', 'gene_b', 'secreted', 'receptor_a', 'receptor_b',
                     'annotation_strategy', 'is_integrin', 'directionality', 'classification', 'root'],
            var_name='interaction',
            value_name='pvalue'
        )

        if split_cols.shape[1] == 0:
            print(
                f"Warning: interaction column could not be split for root {root}. Result: {split_cols}")
            return None
        elif split_cols.shape[1] > 2:
            print(
                f"Warning: interaction column split into more than 2 parts for root {root}. Result: {split_cols}")

        df_pval['source'] = split_cols[0]
        df_pval['target'] = split_cols[1]
        df_pval = df_pval.drop(columns=['interaction'])

        return df_score, df_pval

    return None


# -------------------------------------------------------------------------------
# Parallel cellphone run --------------------------------------------------------
# -------------------------------------------------------------------------------


def work(root):
    return run_one_root(root, adj_2hop_dict[root], meta_df,
                        iterations=200, threads=1)


if running_on_slurm():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nroots", type=int, help="Number of roots to process",
                        default=5000)
    parser.add_argument("--start_root", type=int, help="Index of first root to process",
                        default=0)
    args = parser.parse_args()
    nroots = args.nroots
    start_root = args.start_root
    if nroots == -1:
        nroots = len(roots) - start_root

    # Split across array jobs AFTER reading args
    slurm_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))
    slurm_n = int(os.environ.get("SLURM_ARRAY_TASK_COUNT", 1))

    # Divide total range into equal chunks
    total_roots = roots[start_root:start_root + nroots]
    chunk_size = int(np.ceil(len(total_roots) / slurm_n))
    chunk_start = start_root + slurm_id * chunk_size
    chunk_end = min(start_root + nroots, chunk_start + chunk_size)

    start_root, nroots = chunk_start, chunk_end - chunk_start
    print(f"SLURM_ARRAY_TASK_ID={slurm_id}, processing roots {start_root} to {start_root+nroots} (total {nroots})")

else:
    nroots = 1000
    start_root = 0

end_root = start_root + nroots
if nroots == -1 or end_root > len(roots):
    end_root = len(roots)
    nroots = end_root - start_root

print(f"Running for {nroots} roots...")
print(f"Starting at root index {start_root} and ending at {end_root}...")

n_jobs = int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count()))
results = Parallel(n_jobs=n_jobs, prefer="threads")(
    delayed(work)(r) for r in roots[start_root:(end_root)])

print("Done computing interaction scores and p-values.")
print(f"Results for {len([r for r in results if r is not None])} roots.")


# Postprocess and save results --------------------------------------------------------

res_score = [r[0] if r is not None else None for r in results]
res_pval = [r[1] if r is not None else None for r in results]

# Combine and write *once* (fast formats)
if res_score:
    cci_df = pd.concat(res_score, ignore_index=True)
    cci_df.to_parquet(fig_dir+f"interaction_scores_{start_root}_{end_root}.parquet", index=False)

if res_pval:
    cci_pval_df = pd.concat(res_pval, ignore_index=True)
    cci_pval_df.to_parquet(fig_dir+f"interaction_pvalues_{start_root}_{end_root}.parquet", index=False)


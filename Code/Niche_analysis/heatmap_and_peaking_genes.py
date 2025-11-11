import os
import numpy as np
import pandas as pd
import scanpy as sc
import decoupler as dc
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.stats import combine_pvalues
import scipy.sparse as sp

# --------------------------
# Config: paths and output
# --------------------------
adata_path = "/data/ep924610/kidney_biopsy/trajectory_pec/pec_with_pseudotime.h5ad"
results_root = "/data/ep924610/kidney_biopsy/trajectory_pec/reactome_results"
os.makedirs(results_root, exist_ok=True)

# --------------------------
# Load data
# --------------------------
adata = sc.read_h5ad(adata_path)

# --------------------------
# Load MSigDB Reactome and filter
# --------------------------
reactome = dc.op.resource(name="MSigDB")
reactome = reactome[reactome["collection"].astype(str) == "reactome_pathways"]
reactome = reactome[~reactome.duplicated(["geneset", "genesymbol"])]
reactome = reactome.rename(columns={"geneset": "source", "genesymbol": "target"})
reactome = reactome[["source", "target"]]

# Intersect with genes present in data
genes_present = set(adata.var_names.astype(str))
reactome = reactome[reactome["target"].astype(str).isin(genes_present)]

# --------------------------
# Prepare pseudotime quintiles
# --------------------------
if "pseudotime_via" not in adata.obs:
    raise KeyError("`adata.obs['pseudotime_via']` is required for binning.")

adata = adata[~adata.obs["pseudotime_via"].isna()].copy()
adata = adata[adata.obs.sort_values("pseudotime_via").index].copy()

adata.obs["pseudotime_quintile"] = pd.qcut(
    adata.obs["pseudotime_via"],
    q=5,
    duplicates="drop",
    labels=[f"Q{i}" for i in range(1, 6)]
)

# All pathways for later column naming if needed
all_pathways = np.array(sorted(reactome["source"].astype(str).unique()))

# --------------------------
# Helpers
# --------------------------
def _ensure_df(matrix, obs_names, all_cols=None):
    # Already a DF
    if isinstance(matrix, pd.DataFrame):
        df = matrix.copy()
    # Sparse
    elif sp.issparse(matrix):
        df = pd.DataFrame(matrix.toarray(), index=obs_names)
    # Numpy array
    elif isinstance(matrix, np.ndarray):
        if matrix.ndim == 1:
            df = pd.DataFrame(matrix.reshape(-1, 1), index=obs_names)
        elif matrix.ndim == 2:
            df = pd.DataFrame(matrix, index=obs_names)
        else:
            raise ValueError(f"Expected 1D/2D array, got {matrix.ndim}D with shape {matrix.shape}.")
    # Dict-like
    elif isinstance(matrix, dict):
        for k in ("estimate", "z", "nes", "scores", "acts", "activity"):
            if k in matrix:
                return _ensure_df(matrix[k], obs_names, all_cols=all_cols)
        for k in ("pvals", "pval", "pvalues", "fdr", "padj"):
            if k in matrix:
                return _ensure_df(matrix[k], obs_names, all_cols=all_cols)
        raise TypeError(
            f"Got a dict from obsm with keys {list(matrix.keys())}. "
            "None look like a 2D matrix (e.g. 'estimate', 'z', 'nes', 'pvals'). "
        )
    else:
        raise TypeError(f"Unsupported type for matrix: {type(matrix)}")

    n_obs = len(obs_names)
    if df.shape[0] != n_obs and df.shape[1] == n_obs:
        df = df.T
        df.index = obs_names
    elif df.shape[0] != n_obs:
        raise ValueError(
            f"Matrix shape {df.shape} is incompatible with number of observations {n_obs}."
        )

    if all_cols is not None:
        try:
            if df.shape[1] == len(all_cols) and (
                isinstance(df.columns, pd.RangeIndex) or
                all(isinstance(c, (int, np.integer)) for c in df.columns)
            ):
                df.columns = all_cols
        except Exception:
            pass
    return df

def _combine_pvals(series):
    vals = series.dropna().astype(float)
    vals = vals[vals > 0]
    if len(vals) == 0:
        return np.nan
    _, combined = combine_pvalues(vals, method="fisher")
    return combined

def _detect_gsea_keys(ad):
    avail = list(ad.obsm_keys())
    lower_map = {k: k.lower() for k in avail}
    gsea_keys = [k for k in avail if "gsea" in lower_map[k]]
    if not gsea_keys:
        raise KeyError(f"No GSEA outputs found in .obsm. Available keys: {avail}")

    def is_pval_key(k):
        lk = lower_map[k]
        return any(s in lk for s in ("pval", "pvals", "pvalue", "padj", "fdr"))

    score_candidates = [k for k in gsea_keys if not is_pval_key(k)]
    pval_candidates = [k for k in gsea_keys if is_pval_key(k)]

    if not score_candidates:
        raise KeyError(f"Could not find a GSEA score matrix in .obsm. GSEA-related keys: {gsea_keys}")
    if not pval_candidates:
        raise KeyError(f"Could not find a GSEA p-value matrix in .obsm. GSEA-related keys: {gsea_keys}")

    prefer_order = ["gsea_z", "gsea_estimate", "gsea_nes", "gsea_scores", "gsea"]
    for pref in prefer_order:
        for k in score_candidates:
            if k.lower() == pref:
                score_key = k
                break
        else:
            continue
        break
    else:
        score_key = score_candidates[0]

    prefer_p = ["gsea_pvals", "gsea_pval", "gsea_pvalues", "gsea_fdr", "gsea_padj"]
    for pref in prefer_p:
        for k in pval_candidates:
            if k.lower() == pref:
                pval_key = k
                break
        else:
            continue
        break
    else:
        pval_key = pval_candidates[0]

    return score_key, pval_key


import os
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.sparse as sp

# Where to store ORA-specific combined outputs
ora_root = os.path.join(results_root, "ora_peak_genes")
os.makedirs(ora_root, exist_ok=True)
group_key = "pseudotime_quintile"
# Ensure quintiles exist and are ordered
if group_key not in adata.obs:
    adata.obs[group_key] = pd.qcut(
        adata.obs[pseudotime_key], q=5,
        labels=[f"Q{i}" for i in range(1, 6)],
        duplicates="drop"
    )

if not pd.api.types.is_categorical_dtype(adata.obs[group_key]):
    adata.obs[group_key] = pd.Categorical(
        adata.obs[group_key],
        categories=[f"Q{i}" for i in range(1, 6) if f"Q{i}" in adata.obs[group_key].astype(str).unique()],
        ordered=True
    )

ordered_groups = list(adata.obs[group_key].cat.categories)

def _detect_ora_keys(ad):
    """
    Detect score and p-value keys for ORA in ad.obsm.
    Returns (score_key, pval_key).
    """
    avail = list(ad.obsm_keys())
    lower_map = {k: k.lower() for k in avail}
    ora_keys = [k for k in avail if "ora" in lower_map[k]]
    if not ora_keys:
        raise KeyError(f"No ORA outputs found in .obsm. Available keys: {avail}")

    def is_pval_key(k):
        lk = lower_map[k]
        return any(s in lk for s in ("pval", "pvals", "pvalue", "padj", "fdr"))

    score_candidates = [k for k in ora_keys if not is_pval_key(k)]
    pval_candidates = [k for k in ora_keys if is_pval_key(k)]

    if not score_candidates:
        raise KeyError(f"Could not find an ORA score matrix in .obsm. ORA-related keys: {ora_keys}")
    if not pval_candidates:
        raise KeyError(f"Could not find an ORA p-value matrix in .obsm. ORA-related keys: {ora_keys}")

    # Heuristic preferences
    prefer_score = ["ora_z", "ora_estimate", "ora_odds", "ora_scores", "ora"]
    for pref in prefer_score:
        for k in score_candidates:
            if k.lower() == pref:
                score_key = k
                break
        else:
            continue
        break
    else:
        score_key = score_candidates[0]

    prefer_p = ["ora_pvals", "ora_pval", "ora_pvalues", "ora_fdr", "ora_padj"]
    for pref in prefer_p:
        for k in pval_candidates:
            if k.lower() == pref:
                pval_key = k
                break
        else:
            continue
        break
    else:
        pval_key = pval_candidates[0]

    return score_key, pval_key

# --------------------------
# Identify peaking genes per quintile
# Definition: a gene peaks in Qk if Qk attains both:
# - the highest mean expression across quintiles (ties allowed), and
# - the highest single-cell expression across quintiles (ties allowed).
# Genes with zero expression across all cells are excluded.
# --------------------------
group_key = "pseudotime_quintile"
if group_key not in adata.obs:
    raise KeyError(f"Missing {group_key} in adata.obs.")

groups = list(adata.obs[group_key].cat.categories) if hasattr(adata.obs[group_key], "cat") else sorted(adata.obs[group_key].unique())
var_names = adata.var_names.astype(str).to_numpy()
n_genes = adata.shape[1]
G = len(groups)

import scipy.sparse as sp
import numpy as np

def axis_stat_1d(X, stat="max", axis=0):
    # Ensure consistent sparse type
    if sp.issparse(X):
        X = X.tocsr()

    # Compute reduction
    if stat == "max":
        out = X.max(axis=axis)
    elif stat == "min":
        out = X.min(axis=axis)
    elif stat == "mean":
        out = X.mean(axis=axis)
    elif stat == "sum":
        out = X.sum(axis=axis)
    else:
        raise ValueError(f"Unsupported stat: {stat}")

    # Convert to 1D numeric ndarray
    if sp.issparse(out):
        vec = out.toarray().ravel()
    elif hasattr(out, "A1"):  # numpy.matrix path
        vec = out.A1
    else:
        vec = np.asarray(out).ravel()
    return vec

mean_mat = np.zeros((G, n_genes), dtype=np.float64)
max_mat  = np.zeros((G, n_genes), dtype=np.float64)

for gi, g in enumerate(groups):
    idx = np.where(adata.obs[group_key].to_numpy() == g)[0]
    if idx.size == 0:
        continue
    Xg = adata.X[idx, :]
    if sp.issparse(Xg):
        mean_mat[gi, :] = axis_stat_1d(Xg, stat="mean", axis=0)
        max_mat[gi, :] = axis_stat_1d(Xg, stat="max", axis=0)
    else:
        mean_mat[gi, :] = Xg.mean(axis=0)
        max_mat[gi, :]  = Xg.max(axis=0)

mean_max = mean_mat.max(axis=0)
max_max  = max_mat.max(axis=0)

# Exclude genes that are zero across all cells (never expressed)
nonzero_mask = (mean_max > 0) | (max_max > 0)
eps = 1e-12

peak_genes_by_group = {}
peak_counts = []

for gi, g in enumerate(groups):
    peak_mask = (
        (mean_mat[gi, :] >= (mean_max - eps)) &
        (max_mat[gi, :]  >= (max_max  - eps)) &
        nonzero_mask
    )
    peak_genes = var_names[peak_mask]
    peak_genes_by_group[g] = peak_genes.tolist()
    peak_counts.append({"group": g, "n_peak_genes": int(peak_mask.sum())})

# Save peak gene lists per quintile and a summary table
pd.DataFrame(peak_counts).to_csv(os.path.join(ora_root, "peak_genes_counts_per_quintile.csv"), index=False)
for g in groups:
    gdir = os.path.join(results_root, g)
    os.makedirs(gdir, exist_ok=True)
    with open(os.path.join(gdir, f"{g}_peak_genes.txt"), "w") as fh:
        fh.write("\n".join(peak_genes_by_group[g]))

# --------------------------
# Run ORA per quintile on peak gene lists
# Using dc.mt.ora with a one-sample AnnData (binary indicator vector).
# --------------------------
ora_long_list = []
n_top_ora = 10  # how many top pathways to display per group

for g in groups:
    gdir = os.path.join(results_root, g)
    os.makedirs(gdir, exist_ok=True)

    peak_genes = set(peak_genes_by_group.get(g, []))
    if len(peak_genes) == 0:
        # Save empty artifacts and continue
        empty_df = pd.DataFrame(columns=["pathway", "score", "pval", "group"])
        empty_df.to_csv(os.path.join(gdir, f"{g}_ora_summary.csv"), index=False)
        continue

    indicator = np.isin(var_names, list(peak_genes)).astype(float).reshape(1, -1)
    ad_ora = sc.AnnData(
        X=indicator,
        var=pd.DataFrame(index=var_names),
        obs=pd.DataFrame(index=[f"{g}_peaks"])
    )

    # Run ORA; same API pattern as GSEA
    dc.mt.ora(ad_ora, net=reactome, verbose=False)

    # Detect returned matrices
    available_obsm = list(ad_ora.obsm_keys())
    try:
        ora_score_key, ora_pval_key = _detect_ora_keys(ad_ora)
    except KeyError as e:
        raise KeyError(
            f"{e}\nAvailable `.obsm` keys after ORA for group {g}: {available_obsm}"
        )

    # Normalize to DataFrames
    ora_scores_wide = _ensure_df(ad_ora.obsm[ora_score_key], ad_ora.obs_names, all_cols=None)
    ora_pvals_wide  = _ensure_df(ad_ora.obsm[ora_pval_key],  ad_ora.obs_names, all_cols=None)

    # Attach pathway names if necessary
    if (ora_scores_wide.shape[1] == len(all_pathways) and
        (isinstance(ora_scores_wide.columns, pd.RangeIndex) or
         all(isinstance(c, (int, np.integer)) for c in ora_scores_wide.columns))):
        ora_scores_wide.columns = all_pathways
    if (ora_pvals_wide.shape[1] == len(all_pathways) and
        (isinstance(ora_pvals_wide.columns, pd.RangeIndex) or
         all(isinstance(c, (int, np.integer)) for c in ora_pvals_wide.columns))):
        ora_pvals_wide.columns = all_pathways

    # Save raw ORA matrices
    ora_scores_wide.to_parquet(os.path.join(gdir, f"{g}_ora_scores.parquet"))
    ora_pvals_wide.to_parquet(os.path.join(gdir, f"{g}_ora_pvals.parquet"))
    ora_scores_wide.to_csv(os.path.join(gdir, f"{g}_ora_scores.csv"))
    ora_pvals_wide.to_csv(os.path.join(gdir, f"{g}_ora_pvals.csv"))

    # Long format (one sample row -> straightforward melt)
    scores_long = (
        ora_scores_wide
        .reset_index(names="sample")
        .melt(id_vars="sample", var_name="pathway", value_name="score")
    )
    pvals_long = (
        ora_pvals_wide
        .reset_index(names="sample")
        .melt(id_vars="sample", var_name="pathway", value_name="pval")
    )

    merged = scores_long.merge(pvals_long, on=["sample", "pathway"], how="left")
    merged["group"] = g
    merged["neg_log10_p"] = -np.log10(pd.to_numeric(merged["pval"], errors="coerce").clip(lower=1e-300))

    # Save per-group ORA summary
    merged.to_parquet(os.path.join(gdir, f"{g}_ora_summary.parquet"), index=False)
    merged.to_csv(os.path.join(gdir, f"{g}_ora_summary.csv"), index=False)

    # Per-group plot: top enrichment scores (score column)
    top_g = merged.sort_values("score", ascending=False).head(n_top_ora).copy()
    if not top_g.empty:
        fig = plt.figure(figsize=(6.5, max(3, 0.4 * len(top_g))))
        sns.barplot(
            data=top_g,
            y="pathway", x="score",
            color="#4477AA", edgecolor="black"
        )
        plt.title(f"ORA on peaking genes in {g}: top enrichment scores")
        plt.xlabel("ORA enrichment score")
        plt.ylabel("")
        plt.tight_layout()
        plt.savefig(os.path.join(gdir, f"{g}_ora_top{n_top_ora}_scores_bar.png"), dpi=300)
        plt.close(fig)

    ora_long_list.append(merged[["pathway", "score", "pval", "neg_log10_p", "group"]])

# --------------------------
# Combine ORA across quintiles and plot comparison heatmap
# --------------------------
if len(ora_long_list) > 0:
    ora_df = pd.concat(ora_long_list, ignore_index=True)
    ora_df.to_parquet(os.path.join(ora_root, "ora_all_quintiles.parquet"), index=False)
    ora_df.to_csv(os.path.join(ora_root, "ora_all_quintiles.csv"), index=False)

    # Choose pathways to show: union of top n per group by score
    top_union = (
        ora_df.sort_values(["group", "score"], ascending=[True, False])
              .groupby("group", as_index=False)
              .head(n_top_ora)["pathway"]
              .unique()
    )
    # Preserve the order of pseudotime quintiles for plotting
    group_key = "pseudotime_quintile"

    # Ensure it's categorical and ordered
    if not pd.api.types.is_categorical_dtype(adata.obs[group_key]):
        # If you created Q1..Q5 labels above, make them ordered here
        adata.obs[group_key] = pd.Categorical(adata.obs[group_key], 
                                            categories=[f"Q{i}" for i in range(1, 6) if f"Q{i}" in adata.obs[group_key].astype(str).unique()],
                                            ordered=True)

    ordered_groups = list(adata.obs[group_key].cat.categories)

    # Use the same order for any iteration per group
    groups = ordered_groups
    heat_scores = (
        ora_df[ora_df["pathway"].isin(top_union)]
        .pivot(index="pathway", columns="group", values="score")
        .reindex(columns=ordered_groups)
    )
    # Order rows by maximum score
    heat_scores = heat_scores.loc[heat_scores.max(axis=1).sort_values(ascending=False).index]
    import textwrap

    def wrap_pathway_label(s, width=22):
        s = str(s)
        if s.startswith("REACTOME_"):
            s = s[len("REACTOME_"):]
        s = s.replace("_", " ")
        return textwrap.fill(s, width=width, break_long_words=False)

    # Example:
    # heat_scores is your DataFrame used for the heatmap (pathways as index)
    yticklabels_wrapped = [wrap_pathway_label(p, width=22) for p in heat_scores.index]
    sns.heatmap(
        heat_scores,
        yticklabels=yticklabels_wrapped,
        cmap="mako",
        linewidths=0.5,
        linecolor="white",
        annot=True,
        fmt=".2f",
        cbar_kws={"label": "ORA enrichment score"}
    )
    plt.title("ORA of peaking genes across pseudotime quintiles (score)")
    plt.xlabel("Pseudotime quintile")
    plt.ylabel("")
    plt.tight_layout()
    plt.savefig(os.path.join(ora_root, f"ora_peak_genes_heatmap_top{n_top_ora}_by_group_scores.png"), dpi=300)
    plt.close(fig)
else:
    # No ORA results (e.g., empty peak lists everywhere)
    with open(os.path.join(ora_root, "NO_ORA_RESULTS.txt"), "w") as fh:
        fh.write("No peaking genes detected for any quintile; ORA was skipped.")
import os
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.sparse as sp
import textwrap

# --------------------------
# Inputs
# --------------------------
import os
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.ndimage import gaussian_filter1d
import textwrap
from typing import List, Optional

# --------------------------
# User inputs
# --------------------------
adata_path = "/data/ep924610/kidney_biopsy/trajectory_pec/pec_with_pseudotime.h5ad"
pseudotime_key = "pseudotime_via"
group_key = "pseudotime_quintile"

ora_csv = "/data/ep924610/kidney_biopsy/trajectory_pec/reactome_results/ora_peak_genes/ora_all_quintiles.csv"
output_path = "/data/ep924610/kidney_biopsy/trajectory_pec/reactome_results/ora_peak_genes/pseudotime_heatmap_plus_ORA.png"

# Provide peaking genes either as a CSV/TSV list or from adata.var
# Option A: a CSV/TSV with one column (e.g., "gene") of gene names
peaking_genes_file: Optional[str] = None  # e.g., "/path/to/peaking_genes.csv"
peaking_genes_column_in_file: Optional[str] = None  # column name that holds gene names; if None, use first column

# Option B: a boolean column in adata.var that flags peaking genes
peaking_gene_flag_in_var: Optional[str] = None  # e.g., "is_peaking"

# Option C: a rank column in adata.var; take top N by rank (ascending = more peaky)
peaking_gene_rank_in_var: Optional[str] = None  # e.g., "peak_rank"
top_n_from_rank: int = 300

# Plot and analysis parameters
n_bins = 120                    # pseudotime bins (y-axis)
smooth_sigma_for_order = 1.2    # smoothing along pseudotime to determine peak position
smooth_sigma_for_display = 0.8  # smoothing for display (heatmap), set 0 to disable
cap_z = 2.5                     # clamp z-scores in heatmap to [-cap_z, cap_z]
max_genes_to_plot = 100000        # optional cap to avoid huge plots
alpha = 0.01                    # FDR cutoff (BH within each quintile)
top_k_pathways = 5              # per quintile
cmap_name = "vlag"              # heatmap colormap

# --------------------------
# Helpers
# --------------------------
def bh_adjust(pvals: pd.Series) -> pd.Series:
    p = np.asarray(pvals, dtype=float)
    m = np.isfinite(p).sum()
    out = np.full_like(p, np.nan, dtype=float)
    if m == 0:
        return pd.Series(out, index=pvals.index)
    idx = np.where(np.isfinite(p))[0]
    p_nonan = p[idx]
    order = np.argsort(p_nonan)
    ranks = np.arange(1, order.size + 1)
    q = p_nonan[order] * m / ranks
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.minimum(q, 1.0)
    tmp = np.empty_like(q)
    tmp[order] = q
    out[idx] = tmp
    return pd.Series(out, index=pvals.index)

def wrap_pathway(s: str, width: int = 22) -> str:
    s = str(s)
    if s.upper().startswith("REACTOME_"):
        s = s[len("REACTOME_"):]
    s = s.replace("_", " ")
    return textwrap.fill(s, width=width, break_long_words=False)

def choose_score_column(df: pd.DataFrame) -> str:
    candidates = ["odds_ratio", "OR", "or", "score"]
    for c in candidates:
        if c in df.columns:
            return c
    # Fallback: if neither present, synthesize from neg_log10_p (not an OR, but preserves ranking)
    if "neg_log10_p" in df.columns:
        df["_synthetic_score"] = df["neg_log10_p"].astype(float)
        return "_synthetic_score"
    raise KeyError("No odds ratio/score column found. Expected one of: odds_ratio, OR, score, or neg_log10_p.")

def get_peaking_genes(
    adata,
    peaking_genes_file: Optional[str],
    peaking_genes_column_in_file: Optional[str],
    peaking_gene_flag_in_var: Optional[str],
    peaking_gene_rank_in_var: Optional[str],
    top_n_from_rank: int,
    max_genes_to_plot: int
) -> List[str]:
    var_names = pd.Index(adata.var_names.astype(str))

    # A) File
    if peaking_genes_file is not None:
        ext = os.path.splitext(peaking_genes_file)[1].lower()
        if ext in [".csv"]:
            df = pd.read_csv(peaking_genes_file)
        elif ext in [".tsv", ".txt"]:
            df = pd.read_csv(peaking_genes_file, sep="\t")
        else:
            df = pd.read_csv(peaking_genes_file)
        if peaking_genes_column_in_file is None:
            gene_col = df.columns[0]
        else:
            gene_col = peaking_genes_column_in_file
        genes = df[gene_col].astype(str).tolist()
        genes = [g for g in genes if g in var_names]
        if not genes:
            raise ValueError("No peaking genes from file match adata.var_names.")
        return genes[:max_genes_to_plot]

    # B) Boolean flag in var
    if peaking_gene_flag_in_var is not None and peaking_gene_flag_in_var in adata.var.columns:
        mask = adata.var[peaking_gene_flag_in_var].astype(bool).to_numpy()
        genes = var_names[mask].tolist()
        if not genes:
            raise ValueError(f"No genes with {peaking_gene_flag_in_var}=True.")
        return genes[:max_genes_to_plot]

    # C) Rank column in var
    if peaking_gene_rank_in_var is not None and peaking_gene_rank_in_var in adata.var.columns:
        ser = pd.to_numeric(adata.var[peaking_gene_rank_in_var], errors="coerce")
        genes = ser.sort_values(ascending=True).index.tolist()[:top_n_from_rank]
        genes = [g for g in genes if g in var_names]
        if not genes:
            raise ValueError(f"No genes selected from rank column {peaking_gene_rank_in_var}.")
        return genes[:max_genes_to_plot]

    raise ValueError(
        "Please provide peaking genes via one of:\n"
        "- peaking_genes_file (CSV/TSV) + column name, or\n"
        "- peaking_gene_flag_in_var (boolean column in adata.var), or\n"
        "- peaking_gene_rank_in_var (numeric rank in adata.var)."
    )

# --------------------------
# Load AnnData and prepare
# --------------------------
adata = sc.read_h5ad(adata_path)
if pseudotime_key not in adata.obs:
    raise KeyError(f"Missing {pseudotime_key} in adata.obs.")
adata = adata[~adata.obs[pseudotime_key].isna()].copy()
adata = adata[adata.obs.sort_values(pseudotime_key).index].copy()

# Ensure quintiles exist and are ordered
if group_key not in adata.obs:
    adata.obs[group_key] = pd.qcut(
        adata.obs[pseudotime_key],
        q=5,
        labels=[f"Q{i}" for i in range(1, 6)],
        duplicates="drop",
    )
if not pd.api.types.is_categorical_dtype(adata.obs[group_key]):
    adata.obs[group_key] = pd.Categorical(
        adata.obs[group_key],
        categories=[f"Q{i}" for i in range(1, 6) if f"Q{i}" in adata.obs[group_key].astype(str).unique()],
        ordered=True,
    )
ordered_groups = list(adata.obs[group_key].cat.categories)

from pathlib import Path
from pandas.api.types import CategoricalDtype  # optional: for deprecation fix

# --------------------------
# Peaking genes: use the ones computed for ORA (union of Q1..Q5 lists)
# --------------------------
results_root_guess = Path(ora_csv).resolve().parent.parent  # .../reactome_results
peak_files = [(g, results_root_guess / g / f"{g}_peak_genes.txt") for g in ordered_groups]
use_union_from_files = any(p.is_file() for _, p in peak_files)

if use_union_from_files:
    peak_union = []
    for g, p in peak_files:
        if p.is_file():
            with open(p) as fh:
                peak_union.extend([ln.strip() for ln in fh if ln.strip()])
    # keep only genes present in adata
    peak_union = pd.Index(peak_union).intersection(adata.var_names.astype(str)).unique().tolist()
    if len(peak_union) == 0:
        raise RuntimeError(
            f"Found peak-gene files under {results_root_guess}, but none matched adata.var_names."
        )
    peaking_genes = peak_union[:max_genes_to_plot]
else:
    # Fall back to user-provided options (file/flag/rank)
    peaking_genes = get_peaking_genes(
        adata=adata,
        peaking_genes_file=peaking_genes_file,
        peaking_genes_column_in_file=peaking_genes_column_in_file,
        peaking_gene_flag_in_var=peaking_gene_flag_in_var,
        peaking_gene_rank_in_var=peaking_gene_rank_in_var,
        top_n_from_rank=top_n_from_rank,
        max_genes_to_plot=max_genes_to_plot,
    )

if len(peaking_genes) == 0:
    raise ValueError("Peaking gene list is empty after matching to adata.var_names.")
gene_idx = pd.Index(adata.var_names.astype(str)).get_indexer(peaking_genes)
if np.any(gene_idx < 0):
    raise ValueError("Some peaking genes were not found in adata.var_names.")

# --------------------------
# Bin cells along pseudotime with equal-count bins (no empty bins)
# --------------------------
pt = adata.obs[pseudotime_key].to_numpy()
n_cells = pt.size
# adata is already sorted by pseudotime above
cell_order = np.arange(n_cells)  # already sorted
bin_ids = (cell_order * n_bins // n_cells).astype(int)  # 0..n_bins-1, all bins populated

# Compute per-bin mean expression for selected genes (shape: n_bins x n_genes)
X = adata.X
n_plot_genes = len(peaking_genes)
expr_bins = np.zeros((n_bins, n_plot_genes), dtype=np.float32)

for b in range(n_bins):
    idx = np.where(bin_ids == b)[0]
    Xb = X[idx, :][:, gene_idx]
    if sp.issparse(Xb):
        expr_bins[b, :] = np.asarray(Xb.mean(axis=0)).ravel()
    else:
        expr_bins[b, :] = Xb.mean(axis=0)

# Z-score per gene across bins (along y-axis)
with np.errstate(invalid="ignore", divide="ignore"):
    mu = expr_bins.mean(axis=0, keepdims=True)
    sd = expr_bins.std(axis=0, keepdims=True)
    sd[sd == 0] = 1.0
    z_expr = (expr_bins - mu) / sd

# Smooth for ordering (optional)
if smooth_sigma_for_order and smooth_sigma_for_order > 0:
    z_for_order = gaussian_filter1d(z_expr, sigma=smooth_sigma_for_order, axis=0, mode="nearest")
else:
    z_for_order = z_expr

# Determine peak bin for each gene and order columns
peak_bins = np.argmax(z_for_order, axis=0)
order = np.argsort(peak_bins)
z_ordered = z_expr[:, order]
genes_ordered = [peaking_genes[i] for i in order]

# Optional smoothing for display
if smooth_sigma_for_display and smooth_sigma_for_display > 0:
    z_disp = gaussian_filter1d(z_ordered, sigma=smooth_sigma_for_display, axis=0, mode="nearest")
else:
    z_disp = z_ordered

# Clamp extremes
z_disp = np.clip(z_disp, -cap_z, cap_z)

# --------------------------
# Quintile bands aligned to bins (by majority label per bin)
# --------------------------
bin_groups = []
obs_groups = adata.obs[group_key].to_numpy()
for b in range(n_bins):
    idx = np.where(bin_ids == b)[0]
    if idx.size == 0:
        bin_groups.append(None)
    else:
        # majority vote of the quintile within the bin
        g = pd.Series(obs_groups[idx]).mode().iat[0]
        bin_groups.append(g)
bin_groups = np.array(bin_groups, dtype=object)

# Count bins per quintile and boundary positions
bins_per_group = (
    pd.Series(bin_groups, dtype="category")
      .astype(pd.CategoricalDtype(categories=ordered_groups, ordered=True))
      .value_counts(sort=False)
      .reindex(ordered_groups, fill_value=0)
      .to_numpy()
)

boundary_bins = np.where(pd.Series(bin_groups).shift(-1) != pd.Series(bin_groups))[0][:-1]  # edges between groups
q_centers = []
for g in ordered_groups:
    pos = np.where(bin_groups == g)[0]
    q_centers.append(pos.mean() if pos.size else np.nan)
# --------------------------
# Load ORA results, FDR filter, select top-K per group
# --------------------------
ora = pd.read_csv(ora_csv)
if "group" not in ora.columns or "pval" not in ora.columns:
    raise KeyError("ORA CSV must contain at least 'group' and 'pval' columns.")
ora = ora[ora["group"].isin(ordered_groups)].copy()
# Compute BH within each group
ora["padj"] = ora.groupby("group", group_keys=False)["pval"].apply(bh_adjust)

score_col = choose_score_column(ora)

top_paths_by_group = {}
for g in ordered_groups:
    gdf = ora[(ora["group"] == g) & (ora["padj"] < alpha)].copy()
    if gdf.empty:
        top_paths_by_group[g] = pd.DataFrame(columns=["pathway", score_col, "padj"])
        continue
    gdf = gdf.sort_values(score_col, ascending=False).head(top_k_pathways)
    top_paths_by_group[g] = gdf[["pathway", score_col, "padj"]].copy()

# Global x-limit for barplots
all_disp = pd.concat(top_paths_by_group.values(), ignore_index=True) if len(top_paths_by_group) else pd.DataFrame()
global_max_score = float(all_disp[score_col].max()) if not all_disp.empty else 1.0
global_max_score = max(global_max_score, 1.0)

import numpy as np
from matplotlib.lines import Line2D

# ==========================
# Separate heatmap (arrows only) + individual barplots
# ==========================
import os
import re
import textwrap
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D

sns.set_theme(context="talk", style="white")
expr_cmap = sns.color_palette("inferno", as_cmap=True)

# Colors
marker_up_color = "#D55E00"   # orange
marker_dn_color = "#0072B2"   # blue

# Color per group (adjust palette if you want unique colors)
okabe_ito = ["#033D8F", "#E69F00", "#009E73", "#CC79A7", "#56B4E9"]
colors_by_group = {g: okabe_ito[i % len(okabe_ito)] for i, g in enumerate(ordered_groups)}
marker_genes_up = ['CD44','TMSB10','SPECC1','ITGB1','MMP7','PDLIM7','EPHB2','ITGA11','RUNX1','SOX4']
marker_genes_down = ["SLC4A11", "CDH1", "FGFR2", "PAH", "PIGR"]
# Reverse order so top barplot would correspond to top quintile if you later recombine
groups_plot = list(ordered_groups)[::-1]

# Global score max for consistent x-limits across per-group barplots
if hasattr(all_disp, "empty") and not all_disp.empty:
    global_max_score = float(all_disp[score_col].max())
else:
    # fall back to scanning the dict if needed
    vals = []
    for g in groups_plot:
        gdf = top_paths_by_group.get(g, None)
        if gdf is not None and hasattr(gdf, "empty") and not gdf.empty and score_col in gdf:
            vals.append(float(gdf[score_col].max()))
    global_max_score = max(vals) if len(vals) else 1.0
global_max_score = max(global_max_score, 1.0)
label_for_x = "Odds ratio" if score_col.lower() in ["odds_ratio", "or"] else "Enrichment score"

# Output naming
base_no_ext, _ = os.path.splitext(output_path)
outdir = os.path.dirname(base_no_ext)
if outdir == "":
    outdir = "."
os.makedirs(outdir, exist_ok=True)

def slugify(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", str(s)).strip("_")

def clean_pathway_label(s):
    s = str(s)
    s = re.sub(r'^(REACTOME[_\s]+)', '', s, flags=re.IGNORECASE)
    return s.replace("_", " ").strip()

def wrap_pathway(s, width=30):
    return "\n".join(textwrap.wrap(str(s), width=width))

# --------------------------
# HEATMAP (arrows only)
# --------------------------
def save_heatmap_arrows_only():
    # Figure size: wide; height scales mildly with number of bins
    fig_w = 10
    fig_h = float(np.clip(5 + 0.01 * n_bins, 4, 10))
    fig, ax_heat = plt.subplots(figsize=(fig_w, fig_h), constrained_layout=True)

    im = ax_heat.imshow(
        z_disp,
        origin="lower",
        aspect="auto",
        interpolation="nearest",
        cmap=expr_cmap,
        vmin=-cap_z, vmax=cap_z,
    )

    # Y ticks as % pseudotime
    n_yticks = 6
    yticks = np.linspace(0, n_bins - 1, n_yticks)
    ax_heat.set_yticks(yticks)
    ax_heat.set_yticklabels([f"{int(v)}%" for v in np.linspace(0, 100, n_yticks)], fontsize=11)
    ax_heat.set_ylabel("Pseudotime (binned)")

    # Remove any x tick labels entirely (we keep only arrows)
    ax_heat.set_xticks([])
    ax_heat.set_xlabel("")

    # Quintile boundaries and labels
    for yb in boundary_bins:
        ax_heat.axhline(y=yb + 0.5, color="k", lw=1.2, ls="--", alpha=0.7)
    for y, g in zip(q_centers, ordered_groups):
        if np.isfinite(y):
            ax_heat.text(-0.8, y, g, va="center", ha="right", fontsize=12, clip_on=False)

    # Colorbar
    cb = fig.colorbar(im, ax=ax_heat, orientation="vertical", pad=0.01, fraction=0.03)
    cb.set_label("Z-scored mean expression per pseudotime bin", fontsize=12)
    cb.ax.tick_params(labelsize=11)

    # Markers (triangles) at top/bottom edges
    pos_map = {g: i for i, g in enumerate(genes_ordered)}
    marker_genes_up_local = list(marker_genes_up)
    marker_genes_down_local = list(marker_genes_down)

    x_up = [pos_map[g] for g in marker_genes_up_local if g in pos_map]
    x_dn = [pos_map[g] for g in marker_genes_down_local if g in pos_map]

    if len(x_up):
        ax_heat.scatter(
            x_up, np.full(len(x_up), n_bins - 1 - 0.2),
            marker="v", s=70, c=marker_up_color,
            edgecolors="white", linewidths=0.6, zorder=6
        )
    if len(x_dn):
        ax_heat.scatter(
            x_dn, np.full(len(x_dn), 0 + 0.2),
            marker="^", s=70, c=marker_dn_color,
            edgecolors="white", linewidths=0.6, zorder=6
        )

    legend_handles = []
    if len(x_up):
        legend_handles.append(Line2D([], [], marker="v", linestyle="None",
                                     markerfacecolor=marker_up_color, markeredgecolor="white",
                                     markersize=8, label="Markers (up)"))
    if len(x_dn):
        legend_handles.append(Line2D([], [], marker="^", linestyle="None",
                                     markerfacecolor=marker_dn_color, markeredgecolor="white",
                                     markersize=8, label="Markers (down)"))
    if legend_handles:
        ax_heat.legend(handles=legend_handles, loc="upper left", frameon=True,
                       framealpha=0.9, fontsize=10, title="Genes of interest", title_fontsize=11)

    # Report marker order (left->right along x), plus full gene order
    marker_set = set(marker_genes_up_local) | set(marker_genes_down_local)
    markers_in_order = [g for g in genes_ordered if g in marker_set]

    # Print to console
    print("\nArrowed marker genes (left → right):")
    for i, g in enumerate(markers_in_order, 1):
        direction = "up" if g in marker_genes_up_local else ("down" if g in marker_genes_down_local else "?")
        print(f"{i:>3}. {g}\t[{direction}]")

    # Write TSVs
    arrow_order_path = f"{base_no_ext}_marker_arrow_order.tsv"
    with open(arrow_order_path, "w", encoding="utf-8") as fh:
        fh.write("rank\tgene\tdirection\tx_index\n")
        for i, g in enumerate(markers_in_order, 1):
            direction = "up" if g in marker_genes_up_local else ("down" if g in marker_genes_down_local else "")
            fh.write(f"{i}\t{g}\t{direction}\t{pos_map.get(g, '')}\n")

    full_gene_order_path = f"{base_no_ext}_full_gene_order.tsv"
    with open(full_gene_order_path, "w", encoding="utf-8") as fh:
        fh.write("rank\tgene\n")
        for i, g in enumerate(genes_ordered, 1):
            fh.write(f"{i}\t{g}\n")

    missing_up = [g for g in marker_genes_up_local if g not in pos_map]
    missing_dn = [g for g in marker_genes_down_local if g not in pos_map]
    if missing_up or missing_dn:
        print("\nMarkers not found among genes_ordered:")
        print({"up": missing_up, "down": missing_dn})

    # Save
    heat_png = f"{base_no_ext}_heatmap_arrows_only.png"
    heat_pdf = f"{base_no_ext}_heatmap_arrows_only.pdf"
    fig.savefig(heat_png, dpi=600, bbox_inches="tight")
    fig.savefig(heat_pdf, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved heatmap: {heat_png}\nSaved heatmap: {heat_pdf}")
    print(f"Saved marker order TSV: {arrow_order_path}")
    print(f"Saved full gene order TSV: {full_gene_order_path}")

# --------------------------
# Individual BARPLOTS per group
# --------------------------
def save_barplot_for_group(g):
    gdf = top_paths_by_group.get(g, None)

    # Build labels and sizing
    if gdf is None or (hasattr(gdf, "empty") and gdf.empty):
        labels = []
        n = 0
    else:
        gdf = gdf.sort_values(score_col, ascending=True)
        labels_clean = [clean_pathway_label(p) for p in gdf["pathway"]]
        labels = [wrap_pathway(p, width=30) for p in labels_clean]
        n = len(labels)

    # Dynamic figure height based on number of bars
    # Tweak these two numbers to get more/less vertical real estate per bar
    fig_w = 12
    fig_h = float(np.clip(1.6 + n * 0.55, 2.8, 18))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), constrained_layout=True)

    # If no data for the group, just write a message
    if n == 0:
        ax.text(0.5, 0.5, "No enriched\npathways\n(FDR < 0.05)", ha="center", va="center", fontsize=14)
        ax.set_axis_off()
    else:
        # bar geometry
        bar_thickness = 1.40
        bar_step = 1.58
        y = np.arange(n)
        y_pos = y * bar_step + bar_thickness / 2

        ax.barh(
            y_pos,
            gdf[score_col],
            color=colors_by_group.get(g, "#033D8F"),
            edgecolor="black",
            height=bar_thickness,
            alpha=0.95,
        )

        ax.set_xlim(0, global_max_score * 1.20)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=13)
        top = y_pos[-1] + bar_thickness * 0.8
        ax.set_ylim(-bar_thickness * 0.6, top)

        ax.grid(axis="x", color="0.9", linewidth=0.8)
        ax.tick_params(axis="x", labelsize=12)
        ax.set_xlabel(label_for_x, fontsize=13)

    ax.set_title(str(g), loc="left", fontsize=16, weight="bold")
    ax.margins(y=0.01)

    slug = slugify(g)
    png = f"{base_no_ext}_bar_{slug}.png"
    pdf = f"{base_no_ext}_bar_{slug}.pdf"
    fig.savefig(png, dpi=600, bbox_inches="tight")
    fig.savefig(pdf, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved barplot for {g}: {png}")
    print(f"Saved barplot for {g}: {pdf}")

# --------------------------
# Run the separate saves
# --------------------------
save_heatmap_arrows_only()

for g in groups_plot:
    save_barplot_for_group(g)

sc.pp.normalize_total(adata)
# Logarithmize the data
sc.pp.log1p(adata)
sc.tl.rank_genes_groups(adata, groupby=group_key, method="wilcoxon")
dp = sc.pl.rank_genes_groups_dotplot(adata,groupby=group_key, n_genes=5, dendrogram=False, return_fig = True)
dp.savefig("/data/ep924610/kidney_biopsy/trajectory_pec/marker_genes.pdf")

marking_genes = ['CD44','TMSB10','SPECC1','ITGB1','MMP7',"SLC4A11", "CDH1", "FGFR2", "PAH", "PIGR"]
sc.settings.figdir = "/data/ep924610/kidney_biopsy/trajectory_pec"
dp = sc.pl.dotplot(adata, marking_genes, group_key, standard_scale="var", save="small_dotplot.pdf")

# ==========================
# Cross-quintile heatmap of top-K pathways per quintile
#  - Rows: union of top-K per quintile (deduped, ordered by first appearance)
#  - Columns: quintiles (ordered_groups)
#  - Stars: significant cells (FDR / padj <= alpha)
#  - Uses ORA UNFILTERED table (ora) for values to avoid blanks
#  - Saves PNG + PDF
# ==========================
import os
import re
import textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patheffects as pe

def _clean_pathway_label(s):
    s = str(s)
    s = re.sub(r'^(REACTOME[_\s]+)', '', s, flags=re.IGNORECASE)
    return s.replace('_', ' ').strip()

def _wrap(s, width=38):
    return "\n".join(textwrap.wrap(str(s), width=width))

def _infer_group_col(df: pd.DataFrame):
    for cand in ["group", "quintile", "bin_group", "cluster", "Q", "bin"]:
        if cand in df.columns:
            return cand
    return None

def _infer_sig_col(df: pd.DataFrame):
    # Prefer adjusted P/FDR; fall back to raw p if needed
    for cand in ["padj", "FDR", "fdr", "qval", "q_value", "q", "BH_FDR", "p_fdr", "p_adj", "p_adjusted", "adj_p"]:
        if cand in df.columns:
            return cand
    for cand in ["pval", "pvalue", "p_val", "p"]:
        if cand in df.columns:
            return cand
    return None

def save_cross_quintile_topk_heatmap(
    ora: pd.DataFrame,
    ordered_groups: list,
    top_paths_by_group: dict,
    score_col: str,
    top_k: int = 5,
    alpha: float = 0.05,
    wrap_chars: int = 40,
    cmap_name: str = "viridis",
    out_base: str = "cross_quintile_pathways"
):
    """
    Build a cross-quintile heatmap:
      - Rows: union of top-K pathways per quintile (from top_paths_by_group)
      - Cols: ordered_groups
      - Values: score_col from ora (unfiltered)
      - Stars: cells with significance <= alpha (padj preferred)

    Parameters
    ----------
    ora : long-format DataFrame with columns including ['pathway', score_col, 'group', and preferably 'padj']
    ordered_groups : list of ordered quintiles, e.g., ['Q1','Q2','Q3','Q4','Q5']
    top_paths_by_group : dict[group] -> DataFrame with columns ['pathway', score_col, 'padj' (optional)]
    score_col : str, name of the score column to plot (e.g., 'odds_ratio' or 'score')
    top_k : int, number of pathways to take per group for the union of rows
    alpha : float, FDR/adj-P threshold for stars
    wrap_chars : int, line wrap width for pathway labels
    cmap_name : str, matplotlib/seaborn colormap name for non-diverging scores
    out_base : str, base path (no extension) for saving PNG/PDF
    """
    # Infer essential columns from ora
    if not isinstance(ora, pd.DataFrame) or ora.empty:
        print("ERROR: 'ora' must be a non-empty DataFrame.")
        return

    group_col = _infer_group_col(ora)
    if group_col is None:
        if "group" in ora.columns:
            group_col = "group"
        else:
            raise KeyError("Could not infer group column in 'ora'. Expected a column like 'group'.")

    if "pathway" not in ora.columns:
        raise KeyError("Expected 'pathway' column in 'ora'.")

    if score_col not in ora.columns:
        raise KeyError(f"score_col='{score_col}' not found in 'ora' columns.")

    sig_col_name = _infer_sig_col(ora)  # padj preferred
    # If sig_col is raw p-values, stars will correspond to raw p <= alpha (last resort).
    # Best is to have 'padj' computed beforehand.

    # Ensure output directory exists
    out_dir = os.path.dirname(out_base)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # 1) Collect top-K per group from top_paths_by_group
    x_groups = list(ordered_groups)
    top_by_group = {}
    for g in x_groups:
        gdf = top_paths_by_group.get(g, None)
        if gdf is None or getattr(gdf, "empty", True) or "pathway" not in gdf.columns:
            top_by_group[g] = []
            continue
        # Sort by score_col if present
        gdf2 = gdf.copy()
        if score_col in gdf2.columns:
            gdf2 = gdf2.sort_values(score_col, ascending=False)
        top_by_group[g] = gdf2["pathway"].head(top_k).tolist()

    # 2) Deduplicated row order by first appearance across groups
    row_paths = []
    for g in x_groups:
        for p in top_by_group.get(g, []):
            if p not in row_paths:
                row_paths.append(p)

    if len(row_paths) == 0:
        print("No pathways available to build the cross-quintile heatmap (top_k may be too small or inputs empty).")
        return

    # 3) Build score and significance maps from UNFILTERED ora
    scores = {}
    sigs = {}

    # Use only columns we need
    cols = [group_col, "pathway", score_col] + ([sig_col_name] if sig_col_name else [])
    sub = ora.loc[:, [c for c in cols if c in ora.columns]].copy()

    # Collapse duplicates: score=max (best), significance=min (strongest)
    agg = {score_col: "max"}
    if sig_col_name and sig_col_name in sub.columns:
        agg[sig_col_name] = "min"
    tmp = sub.groupby([group_col, "pathway"], as_index=False).agg(agg)

    for _, row in tmp.iterrows():
        g = row[group_col]
        p = row["pathway"]
        s = row[score_col]
        try:
            s = float(s)
        except Exception:
            s = np.nan
        scores[(g, p)] = s
        if sig_col_name and sig_col_name in row:
            try:
                sig_val = float(row[sig_col_name])
                if np.isfinite(sig_val):
                    sigs[(g, p)] = sig_val
            except Exception:
                pass

    # 4) Build matrix and star mask
    n_rows = len(row_paths)
    n_cols = len(x_groups)
    mat = np.full((n_rows, n_cols), np.nan, dtype=float)
    star = np.zeros((n_rows, n_cols), dtype=bool)

    for i, p in enumerate(row_paths):
        for j, g in enumerate(x_groups):
            s = scores.get((g, p), np.nan)
            mat[i, j] = s
            if (g, p) in sigs:
                try:
                    star[i, j] = (float(sigs[(g, p)]) <= float(alpha))
                except Exception:
                    star[i, j] = False

    # 5) Labels (clean + wrap), dynamic sizing
    labels_clean = [_clean_pathway_label(p) for p in row_paths]
    labels_disp = [_wrap(s, width=wrap_chars) for s in labels_clean]
    line_counts = [lab.count("\n") + 1 for lab in labels_disp]
    total_lines = sum(line_counts)

    fig_h = float(np.clip(1.8 + 0.45 * total_lines, 6, 60))  # grows with number/length of labels
    fig_w = float(np.clip(4.5 + 0.8 * n_cols, 8, 26))        # grows with number of groups

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # 6) Choose colormap; handle negatives if present
    valid = mat[np.isfinite(mat)]
    has_neg = (valid.size > 0 and np.nanmin(valid) < 0)

    if has_neg:
        # Diverging around 0 if negative values present
        cmap = sns.diverging_palette(240, 10, as_cmap=True)
        from matplotlib.colors import TwoSlopeNorm
        vmin = np.nanpercentile(valid, 2) if valid.size else -1.0
        vmax = np.nanpercentile(valid, 98) if valid.size else 1.0
        if not np.isfinite(vmin): vmin = -1.0
        if not np.isfinite(vmax): vmax = 1.0
        norm = TwoSlopeNorm(vcenter=0.0, vmin=vmin, vmax=vmax)
        im = ax.imshow(np.ma.masked_invalid(mat), aspect="auto", interpolation="nearest",
                       cmap=cmap, norm=norm, origin="upper")
    else:
        cmap = sns.color_palette(cmap_name, as_cmap=True)
        try:
            cmap = cmap.copy()
            cmap.set_bad(color="#f0f0f0")
        except Exception:
            pass
        if valid.size:
            vmin = max(0.0, float(np.nanmin(valid)))
            vmax = float(np.nanpercentile(valid, 98))
            if not np.isfinite(vmin): vmin = 0.0
            if not np.isfinite(vmax) or vmax <= vmin:
                vmax = vmin + 1e-6
        else:
            vmin, vmax = 0.0, 1.0
        im = ax.imshow(np.ma.masked_invalid(mat), aspect="auto", interpolation="nearest",
                       cmap=cmap, vmin=vmin, vmax=vmax, origin="upper")

    # 7) Axes ticks/labels
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(x_groups, rotation=0, fontsize=12)
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(labels_disp, fontsize=11)
    ax.set_xlabel("Quintile", fontsize=13)
    ax.set_ylabel("Pathway", fontsize=13)

    # Light grid between cells
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.0, alpha=0.6)
    ax.tick_params(which="both", bottom=False, left=False)

    # 8) Stars for significant cells (if significance data exists)
    if len(sigs) > 0:
        for i in range(n_rows):
            for j in range(n_cols):
                if star[i, j]:
                    ax.text(
                        j, i, "★", ha="center", va="center",
                        fontsize=12, color="black",
                        path_effects=[pe.withStroke(linewidth=3, foreground="white")]
                    )

    # 9) Colorbar label from score_col
    cbar_label = "Odds ratio" if score_col.lower() in ["odds_ratio", "or"] else "Enrichment score"
    cb = fig.colorbar(im, ax=ax, pad=0.01, fraction=0.03)
    cb.set_label(cbar_label, fontsize=12)
    cb.ax.tick_params(labelsize=11)

    # 10) Save
    out_png = f"{out_base}_top{top_k}_pathways_cross_quintile_heatmap.png"
    out_pdf = f"{out_base}_top{top_k}_pathways_cross_quintile_heatmap.pdf"
    fig.tight_layout()
    fig.savefig(out_png, dpi=600, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=600, bbox_inches="tight")
    plt.close(fig)

    # Optional: also export the numeric heatmap matrix
    heat_df = pd.DataFrame(mat, index=row_paths, columns=x_groups)
    heat_csv = f"{out_base}_top{top_k}_pathways_cross_quintile_heatmap_values.csv"
    heat_df.to_csv(heat_csv, index=True)

    print("Saved cross-quintile pathways heatmap:")
    print(f"  {out_png}")
    print(f"  {out_pdf}")
    print(f"Saved underlying values CSV:")
    print(f"  {heat_csv}")

# --------------------------
# Example call (uses your existing variables)
# --------------------------
# Make sure you already have:
#   - ora (unfiltered ORA long table with 'pathway', score_col, 'group', and 'padj' if available)
#   - ordered_groups (['Q1','Q2','Q3','Q4','Q5'])
#   - top_paths_by_group (built after BH/FDR filtering)
#   - score_col (chosen by your choose_score_column)
#   - alpha (FDR threshold)
#   - base_no_ext (derived from your output_path)
#
save_cross_quintile_topk_heatmap(
    ora=ora,
    ordered_groups=ordered_groups,
    top_paths_by_group=top_paths_by_group,
    score_col=score_col,
    top_k=top_k_pathways,   # e.g., 5
    alpha=alpha,            # e.g., 0.05
    wrap_chars=40,
    cmap_name="viridis",
    out_base=base_no_ext
)



# Hallmark

# ==========================
# ORA on peaking genes per quintile
# ==========================
import os
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.sparse as sp

# Where to store ORA-specific combined outputs
ora_root = os.path.join(results_root, "ora_peak_genes_hallmark")
os.makedirs(ora_root, exist_ok=True)
group_key = "pseudotime_quintile"
# Ensure quintiles exist and are ordered
if group_key not in adata.obs:
    adata.obs[group_key] = pd.qcut(
        adata.obs[pseudotime_key], q=5,
        labels=[f"Q{i}" for i in range(1, 6)],
        duplicates="drop"
    )

if not pd.api.types.is_categorical_dtype(adata.obs[group_key]):
    adata.obs[group_key] = pd.Categorical(
        adata.obs[group_key],
        categories=[f"Q{i}" for i in range(1, 6) if f"Q{i}" in adata.obs[group_key].astype(str).unique()],
        ordered=True
    )

ordered_groups = list(adata.obs[group_key].cat.categories)

def _detect_ora_keys(ad):
    """
    Detect score and p-value keys for ORA in ad.obsm.
    Returns (score_key, pval_key).
    """
    avail = list(ad.obsm_keys())
    lower_map = {k: k.lower() for k in avail}
    ora_keys = [k for k in avail if "ora" in lower_map[k]]
    if not ora_keys:
        raise KeyError(f"No ORA outputs found in .obsm. Available keys: {avail}")

    def is_pval_key(k):
        lk = lower_map[k]
        return any(s in lk for s in ("pval", "pvals", "pvalue", "padj", "fdr"))

    score_candidates = [k for k in ora_keys if not is_pval_key(k)]
    pval_candidates = [k for k in ora_keys if is_pval_key(k)]

    if not score_candidates:
        raise KeyError(f"Could not find an ORA score matrix in .obsm. ORA-related keys: {ora_keys}")
    if not pval_candidates:
        raise KeyError(f"Could not find an ORA p-value matrix in .obsm. ORA-related keys: {ora_keys}")

    # Heuristic preferences
    prefer_score = ["ora_z", "ora_estimate", "ora_odds", "ora_scores", "ora"]
    for pref in prefer_score:
        for k in score_candidates:
            if k.lower() == pref:
                score_key = k
                break
        else:
            continue
        break
    else:
        score_key = score_candidates[0]

    prefer_p = ["ora_pvals", "ora_pval", "ora_pvalues", "ora_fdr", "ora_padj"]
    for pref in prefer_p:
        for k in pval_candidates:
            if k.lower() == pref:
                pval_key = k
                break
        else:
            continue
        break
    else:
        pval_key = pval_candidates[0]

    return score_key, pval_key

# --------------------------
# Identify peaking genes per quintile
# Definition: a gene peaks in Qk if Qk attains both:
# - the highest mean expression across quintiles (ties allowed), and
# - the highest single-cell expression across quintiles (ties allowed).
# Genes with zero expression across all cells are excluded.
# --------------------------
group_key = "pseudotime_quintile"
if group_key not in adata.obs:
    raise KeyError(f"Missing {group_key} in adata.obs.")

groups = list(adata.obs[group_key].cat.categories) if hasattr(adata.obs[group_key], "cat") else sorted(adata.obs[group_key].unique())
var_names = adata.var_names.astype(str).to_numpy()
n_genes = adata.shape[1]
G = len(groups)

import scipy.sparse as sp
import numpy as np

def axis_stat_1d(X, stat="max", axis=0):
    # Ensure consistent sparse type
    if sp.issparse(X):
        X = X.tocsr()

    # Compute reduction
    if stat == "max":
        out = X.max(axis=axis)
    elif stat == "min":
        out = X.min(axis=axis)
    elif stat == "mean":
        out = X.mean(axis=axis)
    elif stat == "sum":
        out = X.sum(axis=axis)
    else:
        raise ValueError(f"Unsupported stat: {stat}")

    # Convert to 1D numeric ndarray
    if sp.issparse(out):
        vec = out.toarray().ravel()
    elif hasattr(out, "A1"):  # numpy.matrix path
        vec = out.A1
    else:
        vec = np.asarray(out).ravel()
    return vec

mean_mat = np.zeros((G, n_genes), dtype=np.float64)
max_mat  = np.zeros((G, n_genes), dtype=np.float64)

for gi, g in enumerate(groups):
    idx = np.where(adata.obs[group_key].to_numpy() == g)[0]
    if idx.size == 0:
        continue
    Xg = adata.X[idx, :]
    if sp.issparse(Xg):
        mean_mat[gi, :] = axis_stat_1d(Xg, stat="mean", axis=0)
        max_mat[gi, :] = axis_stat_1d(Xg, stat="max", axis=0)
    else:
        mean_mat[gi, :] = Xg.mean(axis=0)
        max_mat[gi, :]  = Xg.max(axis=0)

mean_max = mean_mat.max(axis=0)
max_max  = max_mat.max(axis=0)

# Exclude genes that are zero across all cells (never expressed)
nonzero_mask = (mean_max > 0) | (max_max > 0)
eps = 1e-12

peak_genes_by_group = {}
peak_counts = []

for gi, g in enumerate(groups):
    peak_mask = (
        (mean_mat[gi, :] >= (mean_max - eps)) &
        (max_mat[gi, :]  >= (max_max  - eps)) &
        nonzero_mask
    )
    peak_genes = var_names[peak_mask]
    peak_genes_by_group[g] = peak_genes.tolist()
    peak_counts.append({"group": g, "n_peak_genes": int(peak_mask.sum())})

# Save peak gene lists per quintile and a summary table
pd.DataFrame(peak_counts).to_csv(os.path.join(ora_root, "peak_genes_counts_per_quintile.csv"), index=False)
for g in groups:
    gdir = os.path.join(results_root, g)
    os.makedirs(gdir, exist_ok=True)
    with open(os.path.join(gdir, f"{g}_peak_genes.txt"), "w") as fh:
        fh.write("\n".join(peak_genes_by_group[g]))

# --------------------------
# Run ORA per quintile on peak gene lists
# Using dc.mt.ora with a one-sample AnnData (binary indicator vector).
# --------------------------
ora_long_list = []
n_top_ora = 10  # how many top pathways to display per group


# --------------------------
# Load MSigDB Reactome and filter
# --------------------------
hallmark = dc.op.hallmark()

# Intersect with genes present in data
genes_present = set(adata.var_names.astype(str))
hallmark = hallmark[hallmark["target"].astype(str).isin(genes_present)]

for g in groups:
    gdir = os.path.join(results_root, g)
    os.makedirs(gdir, exist_ok=True)

    peak_genes = set(peak_genes_by_group.get(g, []))
    if len(peak_genes) == 0:
        # Save empty artifacts and continue
        empty_df = pd.DataFrame(columns=["pathway", "score", "pval", "group"])
        empty_df.to_csv(os.path.join(gdir, f"{g}_ora_summary.csv"), index=False)
        continue

    indicator = np.isin(var_names, list(peak_genes)).astype(float).reshape(1, -1)
    ad_ora = sc.AnnData(
        X=indicator,
        var=pd.DataFrame(index=var_names),
        obs=pd.DataFrame(index=[f"{g}_peaks"])
    )

    # Run ORA; same API pattern as GSEA
    dc.mt.ora(ad_ora, net=hallmark, verbose=False)

    # Detect returned matrices
    available_obsm = list(ad_ora.obsm_keys())
    try:
        ora_score_key, ora_pval_key = _detect_ora_keys(ad_ora)
    except KeyError as e:
        raise KeyError(
            f"{e}\nAvailable `.obsm` keys after ORA for group {g}: {available_obsm}"
        )

    # Normalize to DataFrames
    ora_scores_wide = _ensure_df(ad_ora.obsm[ora_score_key], ad_ora.obs_names, all_cols=None)
    ora_pvals_wide  = _ensure_df(ad_ora.obsm[ora_pval_key],  ad_ora.obs_names, all_cols=None)

    # Attach pathway names if necessary
    if (ora_scores_wide.shape[1] == len(all_pathways) and
        (isinstance(ora_scores_wide.columns, pd.RangeIndex) or
         all(isinstance(c, (int, np.integer)) for c in ora_scores_wide.columns))):
        ora_scores_wide.columns = all_pathways
    if (ora_pvals_wide.shape[1] == len(all_pathways) and
        (isinstance(ora_pvals_wide.columns, pd.RangeIndex) or
         all(isinstance(c, (int, np.integer)) for c in ora_pvals_wide.columns))):
        ora_pvals_wide.columns = all_pathways

    # Save raw ORA matrices
    ora_scores_wide.to_parquet(os.path.join(gdir, f"{g}_ora_scores.parquet"))
    ora_pvals_wide.to_parquet(os.path.join(gdir, f"{g}_ora_pvals.parquet"))
    ora_scores_wide.to_csv(os.path.join(gdir, f"{g}_ora_scores.csv"))
    ora_pvals_wide.to_csv(os.path.join(gdir, f"{g}_ora_pvals.csv"))

    # Long format (one sample row -> straightforward melt)
    scores_long = (
        ora_scores_wide
        .reset_index(names="sample")
        .melt(id_vars="sample", var_name="pathway", value_name="score")
    )
    pvals_long = (
        ora_pvals_wide
        .reset_index(names="sample")
        .melt(id_vars="sample", var_name="pathway", value_name="pval")
    )

    merged = scores_long.merge(pvals_long, on=["sample", "pathway"], how="left")
    merged["group"] = g
    merged["neg_log10_p"] = -np.log10(pd.to_numeric(merged["pval"], errors="coerce").clip(lower=1e-300))

    # Save per-group ORA summary
    merged.to_parquet(os.path.join(gdir, f"{g}_ora_summary.parquet"), index=False)
    merged.to_csv(os.path.join(gdir, f"{g}_ora_summary.csv"), index=False)

    # Per-group plot: top enrichment scores (score column)
    top_g = merged.sort_values("score", ascending=False).head(n_top_ora).copy()
    if not top_g.empty:
        fig = plt.figure(figsize=(6.5, max(3, 0.4 * len(top_g))))
        sns.barplot(
            data=top_g,
            y="pathway", x="score",
            color="#4477AA", edgecolor="black"
        )
        plt.title(f"ORA on peaking genes in {g}: top enrichment scores")
        plt.xlabel("ORA enrichment score")
        plt.ylabel("")
        plt.tight_layout()
        plt.savefig(os.path.join(gdir, f"{g}_ora_top{n_top_ora}_scores_bar.png"), dpi=300)
        plt.close(fig)

    ora_long_list.append(merged[["pathway", "score", "pval", "neg_log10_p", "group"]])

# --------------------------
# Combine ORA across quintiles and plot comparison heatmap
# --------------------------
if len(ora_long_list) > 0:
    ora_df = pd.concat(ora_long_list, ignore_index=True)
    ora_df.to_parquet(os.path.join(ora_root, "ora_all_quintiles.parquet"), index=False)
    ora_df.to_csv(os.path.join(ora_root, "ora_all_quintiles.csv"), index=False)

    # Choose pathways to show: union of top n per group by score
    top_union = (
        ora_df.sort_values(["group", "score"], ascending=[True, False])
              .groupby("group", as_index=False)
              .head(n_top_ora)["pathway"]
              .unique()
    )
    # Preserve the order of pseudotime quintiles for plotting
    group_key = "pseudotime_quintile"

    # Ensure it's categorical and ordered
    if not pd.api.types.is_categorical_dtype(adata.obs[group_key]):
        # If you created Q1..Q5 labels above, make them ordered here
        adata.obs[group_key] = pd.Categorical(adata.obs[group_key], 
                                            categories=[f"Q{i}" for i in range(1, 6) if f"Q{i}" in adata.obs[group_key].astype(str).unique()],
                                            ordered=True)

    ordered_groups = list(adata.obs[group_key].cat.categories)

    # Use the same order for any iteration per group
    groups = ordered_groups
    heat_scores = (
        ora_df[ora_df["pathway"].isin(top_union)]
        .pivot(index="pathway", columns="group", values="score")
        .reindex(columns=ordered_groups)
    )
    # Order rows by maximum score
    heat_scores = heat_scores.loc[heat_scores.max(axis=1).sort_values(ascending=False).index]
    import textwrap

    def wrap_pathway_label(s, width=22):
        s = str(s)
        if s.startswith("REACTOME_"):
            s = s[len("REACTOME_"):]
        s = s.replace("_", " ")
        return textwrap.fill(s, width=width, break_long_words=False)

    # Example:
    # heat_scores is your DataFrame used for the heatmap (pathways as index)
    yticklabels_wrapped = [wrap_pathway_label(p, width=22) for p in heat_scores.index]
    sns.heatmap(
        heat_scores,
        yticklabels=yticklabels_wrapped,
        cmap="mako",
        linewidths=0.5,
        linecolor="white",
        annot=True,
        fmt=".2f",
        cbar_kws={"label": "ORA enrichment score"}
    )
    plt.title("ORA of peaking genes across pseudotime quintiles (score)")
    plt.xlabel("Pseudotime quintile")
    plt.ylabel("")
    plt.tight_layout()
    plt.savefig(os.path.join(ora_root, f"ora_peak_genes_heatmap_top{n_top_ora}_by_group_scores.png"), dpi=300)
    plt.close(fig)
else:
    # No ORA results (e.g., empty peak lists everywhere)
    with open(os.path.join(ora_root, "NO_ORA_RESULTS.txt"), "w") as fh:
        fh.write("No peaking genes detected for any quintile; ORA was skipped.")

# --------------------------
# Inputs
# --------------------------
import os
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.ndimage import gaussian_filter1d
import textwrap
from typing import List, Optional

# --------------------------
# User inputs
# --------------------------
adata_path = "/data/ep924610/kidney_biopsy/trajectory_pec/pec_with_pseudotime.h5ad"
pseudotime_key = "pseudotime_via"
group_key = "pseudotime_quintile"

ora_csv = "/data/ep924610/kidney_biopsy/trajectory_pec/reactome_results/ora_peak_genes_hallmark/ora_all_quintiles.csv"
output_path = "/data/ep924610/kidney_biopsy/trajectory_pec/reactome_results/ora_peak_genes_hallmark/pseudotime_heatmap_plus_ORA.png"

# Provide peaking genes either as a CSV/TSV list or from adata.var
# Option A: a CSV/TSV with one column (e.g., "gene") of gene names
peaking_genes_file: Optional[str] = None  # e.g., "/path/to/peaking_genes.csv"
peaking_genes_column_in_file: Optional[str] = None  # column name that holds gene names; if None, use first column

# Option B: a boolean column in adata.var that flags peaking genes
peaking_gene_flag_in_var: Optional[str] = None  # e.g., "is_peaking"

# Option C: a rank column in adata.var; take top N by rank (ascending = more peaky)
peaking_gene_rank_in_var: Optional[str] = None  # e.g., "peak_rank"
top_n_from_rank: int = 300

# Plot and analysis parameters
n_bins = 120                    # pseudotime bins (y-axis)
smooth_sigma_for_order = 1.2    # smoothing along pseudotime to determine peak position
smooth_sigma_for_display = 0.8  # smoothing for display (heatmap), set 0 to disable
cap_z = 2.5                     # clamp z-scores in heatmap to [-cap_z, cap_z]
max_genes_to_plot = 100000        # optional cap to avoid huge plots
alpha = 0.01                    # FDR cutoff (BH within each quintile)
top_k_pathways = 5              # per quintile
cmap_name = "vlag"              # heatmap colormap

# --------------------------
# Helpers
# --------------------------
def bh_adjust(pvals: pd.Series) -> pd.Series:
    p = np.asarray(pvals, dtype=float)
    m = np.isfinite(p).sum()
    out = np.full_like(p, np.nan, dtype=float)
    if m == 0:
        return pd.Series(out, index=pvals.index)
    idx = np.where(np.isfinite(p))[0]
    p_nonan = p[idx]
    order = np.argsort(p_nonan)
    ranks = np.arange(1, order.size + 1)
    q = p_nonan[order] * m / ranks
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.minimum(q, 1.0)
    tmp = np.empty_like(q)
    tmp[order] = q
    out[idx] = tmp
    return pd.Series(out, index=pvals.index)

def wrap_pathway(s: str, width: int = 22) -> str:
    s = str(s)
    if s.upper().startswith("REACTOME_"):
        s = s[len("REACTOME_"):]
    s = s.replace("_", " ")
    return textwrap.fill(s, width=width, break_long_words=False)

def choose_score_column(df: pd.DataFrame) -> str:
    candidates = ["odds_ratio", "OR", "or", "score"]
    for c in candidates:
        if c in df.columns:
            return c
    # Fallback: if neither present, synthesize from neg_log10_p (not an OR, but preserves ranking)
    if "neg_log10_p" in df.columns:
        df["_synthetic_score"] = df["neg_log10_p"].astype(float)
        return "_synthetic_score"
    raise KeyError("No odds ratio/score column found. Expected one of: odds_ratio, OR, score, or neg_log10_p.")

def get_peaking_genes(
    adata,
    peaking_genes_file: Optional[str],
    peaking_genes_column_in_file: Optional[str],
    peaking_gene_flag_in_var: Optional[str],
    peaking_gene_rank_in_var: Optional[str],
    top_n_from_rank: int,
    max_genes_to_plot: int
) -> List[str]:
    var_names = pd.Index(adata.var_names.astype(str))

    # A) File
    if peaking_genes_file is not None:
        ext = os.path.splitext(peaking_genes_file)[1].lower()
        if ext in [".csv"]:
            df = pd.read_csv(peaking_genes_file)
        elif ext in [".tsv", ".txt"]:
            df = pd.read_csv(peaking_genes_file, sep="\t")
        else:
            df = pd.read_csv(peaking_genes_file)
        if peaking_genes_column_in_file is None:
            gene_col = df.columns[0]
        else:
            gene_col = peaking_genes_column_in_file
        genes = df[gene_col].astype(str).tolist()
        genes = [g for g in genes if g in var_names]
        if not genes:
            raise ValueError("No peaking genes from file match adata.var_names.")
        return genes[:max_genes_to_plot]

    # B) Boolean flag in var
    if peaking_gene_flag_in_var is not None and peaking_gene_flag_in_var in adata.var.columns:
        mask = adata.var[peaking_gene_flag_in_var].astype(bool).to_numpy()
        genes = var_names[mask].tolist()
        if not genes:
            raise ValueError(f"No genes with {peaking_gene_flag_in_var}=True.")
        return genes[:max_genes_to_plot]

    # C) Rank column in var
    if peaking_gene_rank_in_var is not None and peaking_gene_rank_in_var in adata.var.columns:
        ser = pd.to_numeric(adata.var[peaking_gene_rank_in_var], errors="coerce")
        genes = ser.sort_values(ascending=True).index.tolist()[:top_n_from_rank]
        genes = [g for g in genes if g in var_names]
        if not genes:
            raise ValueError(f"No genes selected from rank column {peaking_gene_rank_in_var}.")
        return genes[:max_genes_to_plot]

    raise ValueError(
        "Please provide peaking genes via one of:\n"
        "- peaking_genes_file (CSV/TSV) + column name, or\n"
        "- peaking_gene_flag_in_var (boolean column in adata.var), or\n"
        "- peaking_gene_rank_in_var (numeric rank in adata.var)."
    )

# --------------------------
# Load AnnData and prepare
# --------------------------
adata = sc.read_h5ad(adata_path)
if pseudotime_key not in adata.obs:
    raise KeyError(f"Missing {pseudotime_key} in adata.obs.")
adata = adata[~adata.obs[pseudotime_key].isna()].copy()
adata = adata[adata.obs.sort_values(pseudotime_key).index].copy()

# Ensure quintiles exist and are ordered
if group_key not in adata.obs:
    adata.obs[group_key] = pd.qcut(
        adata.obs[pseudotime_key],
        q=5,
        labels=[f"Q{i}" for i in range(1, 6)],
        duplicates="drop",
    )
if not pd.api.types.is_categorical_dtype(adata.obs[group_key]):
    adata.obs[group_key] = pd.Categorical(
        adata.obs[group_key],
        categories=[f"Q{i}" for i in range(1, 6) if f"Q{i}" in adata.obs[group_key].astype(str).unique()],
        ordered=True,
    )
ordered_groups = list(adata.obs[group_key].cat.categories)

from pathlib import Path
from pandas.api.types import CategoricalDtype  # optional: for deprecation fix

# --------------------------
# Peaking genes: use the ones computed for ORA (union of Q1..Q5 lists)
# --------------------------
results_root_guess = Path(ora_csv).resolve().parent.parent  # .../reactome_results
peak_files = [(g, results_root_guess / g / f"{g}_peak_genes.txt") for g in ordered_groups]
use_union_from_files = any(p.is_file() for _, p in peak_files)

if use_union_from_files:
    peak_union = []
    for g, p in peak_files:
        if p.is_file():
            with open(p) as fh:
                peak_union.extend([ln.strip() for ln in fh if ln.strip()])
    # keep only genes present in adata
    peak_union = pd.Index(peak_union).intersection(adata.var_names.astype(str)).unique().tolist()
    if len(peak_union) == 0:
        raise RuntimeError(
            f"Found peak-gene files under {results_root_guess}, but none matched adata.var_names."
        )
    peaking_genes = peak_union[:max_genes_to_plot]
else:
    # Fall back to user-provided options (file/flag/rank)
    peaking_genes = get_peaking_genes(
        adata=adata,
        peaking_genes_file=peaking_genes_file,
        peaking_genes_column_in_file=peaking_genes_column_in_file,
        peaking_gene_flag_in_var=peaking_gene_flag_in_var,
        peaking_gene_rank_in_var=peaking_gene_rank_in_var,
        top_n_from_rank=top_n_from_rank,
        max_genes_to_plot=max_genes_to_plot,
    )

if len(peaking_genes) == 0:
    raise ValueError("Peaking gene list is empty after matching to adata.var_names.")
gene_idx = pd.Index(adata.var_names.astype(str)).get_indexer(peaking_genes)
if np.any(gene_idx < 0):
    raise ValueError("Some peaking genes were not found in adata.var_names.")

# --------------------------
# Bin cells along pseudotime with equal-count bins (no empty bins)
# --------------------------
pt = adata.obs[pseudotime_key].to_numpy()
n_cells = pt.size
# adata is already sorted by pseudotime above
cell_order = np.arange(n_cells)  # already sorted
bin_ids = (cell_order * n_bins // n_cells).astype(int)  # 0..n_bins-1, all bins populated

# Compute per-bin mean expression for selected genes (shape: n_bins x n_genes)
X = adata.X
n_plot_genes = len(peaking_genes)
expr_bins = np.zeros((n_bins, n_plot_genes), dtype=np.float32)

for b in range(n_bins):
    idx = np.where(bin_ids == b)[0]
    Xb = X[idx, :][:, gene_idx]
    if sp.issparse(Xb):
        expr_bins[b, :] = np.asarray(Xb.mean(axis=0)).ravel()
    else:
        expr_bins[b, :] = Xb.mean(axis=0)

# Z-score per gene across bins (along y-axis)
with np.errstate(invalid="ignore", divide="ignore"):
    mu = expr_bins.mean(axis=0, keepdims=True)
    sd = expr_bins.std(axis=0, keepdims=True)
    sd[sd == 0] = 1.0
    z_expr = (expr_bins - mu) / sd

# Smooth for ordering (optional)
if smooth_sigma_for_order and smooth_sigma_for_order > 0:
    z_for_order = gaussian_filter1d(z_expr, sigma=smooth_sigma_for_order, axis=0, mode="nearest")
else:
    z_for_order = z_expr

# Determine peak bin for each gene and order columns
peak_bins = np.argmax(z_for_order, axis=0)
order = np.argsort(peak_bins)
z_ordered = z_expr[:, order]
genes_ordered = [peaking_genes[i] for i in order]

# Optional smoothing for display
if smooth_sigma_for_display and smooth_sigma_for_display > 0:
    z_disp = gaussian_filter1d(z_ordered, sigma=smooth_sigma_for_display, axis=0, mode="nearest")
else:
    z_disp = z_ordered

# Clamp extremes
z_disp = np.clip(z_disp, -cap_z, cap_z)

# --------------------------
# Quintile bands aligned to bins (by majority label per bin)
# --------------------------
bin_groups = []
obs_groups = adata.obs[group_key].to_numpy()
for b in range(n_bins):
    idx = np.where(bin_ids == b)[0]
    if idx.size == 0:
        bin_groups.append(None)
    else:
        # majority vote of the quintile within the bin
        g = pd.Series(obs_groups[idx]).mode().iat[0]
        bin_groups.append(g)
bin_groups = np.array(bin_groups, dtype=object)

# Count bins per quintile and boundary positions
bins_per_group = (
    pd.Series(bin_groups, dtype="category")
      .astype(pd.CategoricalDtype(categories=ordered_groups, ordered=True))
      .value_counts(sort=False)
      .reindex(ordered_groups, fill_value=0)
      .to_numpy()
)

boundary_bins = np.where(pd.Series(bin_groups).shift(-1) != pd.Series(bin_groups))[0][:-1]  # edges between groups
q_centers = []
for g in ordered_groups:
    pos = np.where(bin_groups == g)[0]
    q_centers.append(pos.mean() if pos.size else np.nan)
# --------------------------
# Load ORA results, FDR filter, select top-K per group
# --------------------------
ora = pd.read_csv(ora_csv)
if "group" not in ora.columns or "pval" not in ora.columns:
    raise KeyError("ORA CSV must contain at least 'group' and 'pval' columns.")
ora = ora[ora["group"].isin(ordered_groups)].copy()
# Compute BH within each group
ora["padj"] = ora.groupby("group", group_keys=False)["pval"].apply(bh_adjust)

score_col = choose_score_column(ora)

top_paths_by_group = {}
for g in ordered_groups:
    gdf = ora[(ora["group"] == g) & (ora["padj"] < alpha)].copy()
    if gdf.empty:
        top_paths_by_group[g] = pd.DataFrame(columns=["pathway", score_col, "padj"])
        continue
    gdf = gdf.sort_values(score_col, ascending=False).head(top_k_pathways)
    top_paths_by_group[g] = gdf[["pathway", score_col, "padj"]].copy()

# Global x-limit for barplots
all_disp = pd.concat(top_paths_by_group.values(), ignore_index=True) if len(top_paths_by_group) else pd.DataFrame()
global_max_score = float(all_disp[score_col].max()) if not all_disp.empty else 1.0
global_max_score = max(global_max_score, 1.0)

import numpy as np
from matplotlib.lines import Line2D

# ==========================
# Separate heatmap (arrows only) + individual barplots
# ==========================
import os
import re
import textwrap
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D

sns.set_theme(context="talk", style="white")
expr_cmap = sns.color_palette("inferno", as_cmap=True)

# Colors
marker_up_color = "#D55E00"   # orange
marker_dn_color = "#0072B2"   # blue

# Color per group (adjust palette if you want unique colors)
okabe_ito = ["#033D8F", "#E69F00", "#009E73", "#CC79A7", "#56B4E9"]
colors_by_group = {g: okabe_ito[i % len(okabe_ito)] for i, g in enumerate(ordered_groups)}
marker_genes_up = ['CD44','TMSB10','SPECC1','ITGB1','MMP7','PDLIM7','EPHB2','ITGA11','RUNX1','SOX4']
marker_genes_down = ["SLC4A11", "CDH1", "FGFR2", "PAH", "PIGR"]
# Reverse order so top barplot would correspond to top quintile if you later recombine
groups_plot = list(ordered_groups)[::-1]

# Global score max for consistent x-limits across per-group barplots
if hasattr(all_disp, "empty") and not all_disp.empty:
    global_max_score = float(all_disp[score_col].max())
else:
    # fall back to scanning the dict if needed
    vals = []
    for g in groups_plot:
        gdf = top_paths_by_group.get(g, None)
        if gdf is not None and hasattr(gdf, "empty") and not gdf.empty and score_col in gdf:
            vals.append(float(gdf[score_col].max()))
    global_max_score = max(vals) if len(vals) else 1.0
global_max_score = max(global_max_score, 1.0)
label_for_x = "Odds ratio" if score_col.lower() in ["odds_ratio", "or"] else "Enrichment score"

# Output naming
base_no_ext, _ = os.path.splitext(output_path)
outdir = os.path.dirname(base_no_ext)
if outdir == "":
    outdir = "."
os.makedirs(outdir, exist_ok=True)

def slugify(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", str(s)).strip("_")

def clean_pathway_label(s):
    s = str(s)
    s = re.sub(r'^(REACTOME[_\s]+)', '', s, flags=re.IGNORECASE)
    return s.replace("_", " ").strip()

def wrap_pathway(s, width=30):
    return "\n".join(textwrap.wrap(str(s), width=width))

# --------------------------
# HEATMAP (arrows only)
# --------------------------
def save_heatmap_arrows_only():
    # Figure size: wide; height scales mildly with number of bins
    fig_w = 10
    fig_h = float(np.clip(5 + 0.01 * n_bins, 4, 10))
    fig, ax_heat = plt.subplots(figsize=(fig_w, fig_h), constrained_layout=True)

    im = ax_heat.imshow(
        z_disp,
        origin="lower",
        aspect="auto",
        interpolation="nearest",
        cmap=expr_cmap,
        vmin=-cap_z, vmax=cap_z,
    )

    # Y ticks as % pseudotime
    n_yticks = 6
    yticks = np.linspace(0, n_bins - 1, n_yticks)
    ax_heat.set_yticks(yticks)
    ax_heat.set_yticklabels([f"{int(v)}%" for v in np.linspace(0, 100, n_yticks)], fontsize=11)
    ax_heat.set_ylabel("Pseudotime (binned)")

    # Remove any x tick labels entirely (we keep only arrows)
    ax_heat.set_xticks([])
    ax_heat.set_xlabel("")

    # Quintile boundaries and labels
    for yb in boundary_bins:
        ax_heat.axhline(y=yb + 0.5, color="k", lw=1.2, ls="--", alpha=0.7)
    for y, g in zip(q_centers, ordered_groups):
        if np.isfinite(y):
            ax_heat.text(-0.8, y, g, va="center", ha="right", fontsize=12, clip_on=False)

    # Colorbar
    cb = fig.colorbar(im, ax=ax_heat, orientation="vertical", pad=0.01, fraction=0.03)
    cb.set_label("Z-scored mean expression per pseudotime bin", fontsize=12)
    cb.ax.tick_params(labelsize=11)

    # Markers (triangles) at top/bottom edges
    pos_map = {g: i for i, g in enumerate(genes_ordered)}
    marker_genes_up_local = list(marker_genes_up)
    marker_genes_down_local = list(marker_genes_down)

    x_up = [pos_map[g] for g in marker_genes_up_local if g in pos_map]
    x_dn = [pos_map[g] for g in marker_genes_down_local if g in pos_map]

    if len(x_up):
        ax_heat.scatter(
            x_up, np.full(len(x_up), n_bins - 1 - 0.2),
            marker="v", s=70, c=marker_up_color,
            edgecolors="white", linewidths=0.6, zorder=6
        )
    if len(x_dn):
        ax_heat.scatter(
            x_dn, np.full(len(x_dn), 0 + 0.2),
            marker="^", s=70, c=marker_dn_color,
            edgecolors="white", linewidths=0.6, zorder=6
        )

    legend_handles = []
    if len(x_up):
        legend_handles.append(Line2D([], [], marker="v", linestyle="None",
                                     markerfacecolor=marker_up_color, markeredgecolor="white",
                                     markersize=8, label="Markers (up)"))
    if len(x_dn):
        legend_handles.append(Line2D([], [], marker="^", linestyle="None",
                                     markerfacecolor=marker_dn_color, markeredgecolor="white",
                                     markersize=8, label="Markers (down)"))
    if legend_handles:
        ax_heat.legend(handles=legend_handles, loc="upper left", frameon=True,
                       framealpha=0.9, fontsize=10, title="Genes of interest", title_fontsize=11)

    # Report marker order (left->right along x), plus full gene order
    marker_set = set(marker_genes_up_local) | set(marker_genes_down_local)
    markers_in_order = [g for g in genes_ordered if g in marker_set]

    # Print to console
    print("\nArrowed marker genes (left → right):")
    for i, g in enumerate(markers_in_order, 1):
        direction = "up" if g in marker_genes_up_local else ("down" if g in marker_genes_down_local else "?")
        print(f"{i:>3}. {g}\t[{direction}]")

    # Write TSVs
    arrow_order_path = f"{base_no_ext}_marker_arrow_order.tsv"
    with open(arrow_order_path, "w", encoding="utf-8") as fh:
        fh.write("rank\tgene\tdirection\tx_index\n")
        for i, g in enumerate(markers_in_order, 1):
            direction = "up" if g in marker_genes_up_local else ("down" if g in marker_genes_down_local else "")
            fh.write(f"{i}\t{g}\t{direction}\t{pos_map.get(g, '')}\n")

    full_gene_order_path = f"{base_no_ext}_full_gene_order.tsv"
    with open(full_gene_order_path, "w", encoding="utf-8") as fh:
        fh.write("rank\tgene\n")
        for i, g in enumerate(genes_ordered, 1):
            fh.write(f"{i}\t{g}\n")

    missing_up = [g for g in marker_genes_up_local if g not in pos_map]
    missing_dn = [g for g in marker_genes_down_local if g not in pos_map]
    if missing_up or missing_dn:
        print("\nMarkers not found among genes_ordered:")
        print({"up": missing_up, "down": missing_dn})

    # Save
    heat_png = f"{base_no_ext}_heatmap_arrows_only.png"
    heat_pdf = f"{base_no_ext}_heatmap_arrows_only.pdf"
    fig.savefig(heat_png, dpi=600, bbox_inches="tight")
    fig.savefig(heat_pdf, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved heatmap: {heat_png}\nSaved heatmap: {heat_pdf}")
    print(f"Saved marker order TSV: {arrow_order_path}")
    print(f"Saved full gene order TSV: {full_gene_order_path}")

# --------------------------
# Individual BARPLOTS per group
# --------------------------
def save_barplot_for_group(g):
    gdf = top_paths_by_group.get(g, None)

    # Build labels and sizing
    if gdf is None or (hasattr(gdf, "empty") and gdf.empty):
        labels = []
        n = 0
    else:
        gdf = gdf.sort_values(score_col, ascending=True)
        labels_clean = [clean_pathway_label(p) for p in gdf["pathway"]]
        labels = [wrap_pathway(p, width=30) for p in labels_clean]
        n = len(labels)

    # Dynamic figure height based on number of bars
    # Tweak these two numbers to get more/less vertical real estate per bar
    fig_w = 12
    fig_h = float(np.clip(1.6 + n * 0.55, 2.8, 18))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), constrained_layout=True)

    # If no data for the group, just write a message
    if n == 0:
        ax.text(0.5, 0.5, "No enriched\npathways\n(FDR < 0.05)", ha="center", va="center", fontsize=14)
        ax.set_axis_off()
    else:
        # bar geometry
        bar_thickness = 1.40
        bar_step = 1.58
        y = np.arange(n)
        y_pos = y * bar_step + bar_thickness / 2

        ax.barh(
            y_pos,
            gdf[score_col],
            color=colors_by_group.get(g, "#033D8F"),
            edgecolor="black",
            height=bar_thickness,
            alpha=0.95,
        )

        ax.set_xlim(0, global_max_score * 1.20)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=13)
        top = y_pos[-1] + bar_thickness * 0.8
        ax.set_ylim(-bar_thickness * 0.6, top)

        ax.grid(axis="x", color="0.9", linewidth=0.8)
        ax.tick_params(axis="x", labelsize=12)
        ax.set_xlabel(label_for_x, fontsize=13)

    ax.set_title(str(g), loc="left", fontsize=16, weight="bold")
    ax.margins(y=0.01)

    slug = slugify(g)
    png = f"{base_no_ext}_bar_{slug}.png"
    pdf = f"{base_no_ext}_bar_{slug}.pdf"
    fig.savefig(png, dpi=600, bbox_inches="tight")
    fig.savefig(pdf, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved barplot for {g}: {png}")
    print(f"Saved barplot for {g}: {pdf}")

# --------------------------
# Run the separate saves
# --------------------------
save_heatmap_arrows_only()

for g in groups_plot:
    save_barplot_for_group(g)

sc.pp.normalize_total(adata)
# Logarithmize the data
sc.pp.log1p(adata)
sc.tl.rank_genes_groups(adata, groupby=group_key, method="wilcoxon")
dp = sc.pl.rank_genes_groups_dotplot(adata,groupby=group_key, n_genes=5, dendrogram=False, return_fig = True)
dp.savefig("/data/ep924610/kidney_biopsy/trajectory_pec/marker_genes.pdf")

marking_genes = ['CD44','TMSB10','SPECC1','ITGB1','MMP7',"SLC4A11", "CDH1", "FGFR2", "PAH", "PIGR"]
sc.settings.figdir = "/data/ep924610/kidney_biopsy/trajectory_pec"
dp = sc.pl.dotplot(adata, marking_genes, group_key, standard_scale="var", save="small_dotplot.pdf")

# ==========================
# Cross-quintile heatmap of top-K pathways per quintile
#  - Rows: union of top-K per quintile (deduped, ordered by first appearance)
#  - Columns: quintiles (ordered_groups)
#  - Stars: significant cells (FDR / padj <= alpha)
#  - Uses ORA UNFILTERED table (ora) for values to avoid blanks
#  - Saves PNG + PDF
# ==========================
import os
import re
import textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patheffects as pe

def _clean_pathway_label(s):
    s = str(s)
    s = re.sub(r'^(REACTOME[_\s]+)', '', s, flags=re.IGNORECASE)
    return s.replace('_', ' ').strip()

def _wrap(s, width=38):
    return "\n".join(textwrap.wrap(str(s), width=width))

def _infer_group_col(df: pd.DataFrame):
    for cand in ["group", "quintile", "bin_group", "cluster", "Q", "bin"]:
        if cand in df.columns:
            return cand
    return None

def _infer_sig_col(df: pd.DataFrame):
    # Prefer adjusted P/FDR; fall back to raw p if needed
    for cand in ["padj", "FDR", "fdr", "qval", "q_value", "q", "BH_FDR", "p_fdr", "p_adj", "p_adjusted", "adj_p"]:
        if cand in df.columns:
            return cand
    for cand in ["pval", "pvalue", "p_val", "p"]:
        if cand in df.columns:
            return cand
    return None

def save_cross_quintile_topk_heatmap(
    ora: pd.DataFrame,
    ordered_groups: list,
    top_paths_by_group: dict,
    score_col: str,
    top_k: int = 5,
    alpha: float = 0.05,
    wrap_chars: int = 40,
    cmap_name: str = "viridis",
    out_base: str = "cross_quintile_pathways"
):
    """
    Build a cross-quintile heatmap:
      - Rows: union of top-K pathways per quintile (from top_paths_by_group)
      - Cols: ordered_groups
      - Values: score_col from ora (unfiltered)
      - Stars: cells with significance <= alpha (padj preferred)

    Parameters
    ----------
    ora : long-format DataFrame with columns including ['pathway', score_col, 'group', and preferably 'padj']
    ordered_groups : list of ordered quintiles, e.g., ['Q1','Q2','Q3','Q4','Q5']
    top_paths_by_group : dict[group] -> DataFrame with columns ['pathway', score_col, 'padj' (optional)]
    score_col : str, name of the score column to plot (e.g., 'odds_ratio' or 'score')
    top_k : int, number of pathways to take per group for the union of rows
    alpha : float, FDR/adj-P threshold for stars
    wrap_chars : int, line wrap width for pathway labels
    cmap_name : str, matplotlib/seaborn colormap name for non-diverging scores
    out_base : str, base path (no extension) for saving PNG/PDF
    """
    # Infer essential columns from ora
    if not isinstance(ora, pd.DataFrame) or ora.empty:
        print("ERROR: 'ora' must be a non-empty DataFrame.")
        return

    group_col = _infer_group_col(ora)
    if group_col is None:
        if "group" in ora.columns:
            group_col = "group"
        else:
            raise KeyError("Could not infer group column in 'ora'. Expected a column like 'group'.")

    if "pathway" not in ora.columns:
        raise KeyError("Expected 'pathway' column in 'ora'.")

    if score_col not in ora.columns:
        raise KeyError(f"score_col='{score_col}' not found in 'ora' columns.")

    sig_col_name = _infer_sig_col(ora)  # padj preferred
    # If sig_col is raw p-values, stars will correspond to raw p <= alpha (last resort).
    # Best is to have 'padj' computed beforehand.

    # Ensure output directory exists
    out_dir = os.path.dirname(out_base)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # 1) Collect top-K per group from top_paths_by_group
    x_groups = list(ordered_groups)
    top_by_group = {}
    for g in x_groups:
        gdf = top_paths_by_group.get(g, None)
        if gdf is None or getattr(gdf, "empty", True) or "pathway" not in gdf.columns:
            top_by_group[g] = []
            continue
        # Sort by score_col if present
        gdf2 = gdf.copy()
        if score_col in gdf2.columns:
            gdf2 = gdf2.sort_values(score_col, ascending=False)
        top_by_group[g] = gdf2["pathway"].head(top_k).tolist()

    # 2) Deduplicated row order by first appearance across groups
    row_paths = []
    for g in x_groups:
        for p in top_by_group.get(g, []):
            if p not in row_paths:
                row_paths.append(p)

    if len(row_paths) == 0:
        print("No pathways available to build the cross-quintile heatmap (top_k may be too small or inputs empty).")
        return

    # 3) Build score and significance maps from UNFILTERED ora
    scores = {}
    sigs = {}

    # Use only columns we need
    cols = [group_col, "pathway", score_col] + ([sig_col_name] if sig_col_name else [])
    sub = ora.loc[:, [c for c in cols if c in ora.columns]].copy()

    # Collapse duplicates: score=max (best), significance=min (strongest)
    agg = {score_col: "max"}
    if sig_col_name and sig_col_name in sub.columns:
        agg[sig_col_name] = "min"
    tmp = sub.groupby([group_col, "pathway"], as_index=False).agg(agg)

    for _, row in tmp.iterrows():
        g = row[group_col]
        p = row["pathway"]
        s = row[score_col]
        try:
            s = float(s)
        except Exception:
            s = np.nan
        scores[(g, p)] = s
        if sig_col_name and sig_col_name in row:
            try:
                sig_val = float(row[sig_col_name])
                if np.isfinite(sig_val):
                    sigs[(g, p)] = sig_val
            except Exception:
                pass

    # 4) Build matrix and star mask
    n_rows = len(row_paths)
    n_cols = len(x_groups)
    mat = np.full((n_rows, n_cols), np.nan, dtype=float)
    star = np.zeros((n_rows, n_cols), dtype=bool)

    for i, p in enumerate(row_paths):
        for j, g in enumerate(x_groups):
            s = scores.get((g, p), np.nan)
            mat[i, j] = s
            if (g, p) in sigs:
                try:
                    star[i, j] = (float(sigs[(g, p)]) <= float(alpha))
                except Exception:
                    star[i, j] = False

    # 5) Labels (clean + wrap), dynamic sizing
    labels_clean = [_clean_pathway_label(p) for p in row_paths]
    labels_disp = [_wrap(s, width=wrap_chars) for s in labels_clean]
    line_counts = [lab.count("\n") + 1 for lab in labels_disp]
    total_lines = sum(line_counts)

    fig_h = float(np.clip(1.8 + 0.45 * total_lines, 6, 60))  # grows with number/length of labels
    fig_w = float(np.clip(4.5 + 0.8 * n_cols, 8, 26))        # grows with number of groups

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # 6) Choose colormap; handle negatives if present
    valid = mat[np.isfinite(mat)]
    has_neg = (valid.size > 0 and np.nanmin(valid) < 0)

    if has_neg:
        # Diverging around 0 if negative values present
        cmap = sns.diverging_palette(240, 10, as_cmap=True)
        from matplotlib.colors import TwoSlopeNorm
        vmin = np.nanpercentile(valid, 2) if valid.size else -1.0
        vmax = np.nanpercentile(valid, 98) if valid.size else 1.0
        if not np.isfinite(vmin): vmin = -1.0
        if not np.isfinite(vmax): vmax = 1.0
        norm = TwoSlopeNorm(vcenter=0.0, vmin=vmin, vmax=vmax)
        im = ax.imshow(np.ma.masked_invalid(mat), aspect="auto", interpolation="nearest",
                       cmap=cmap, norm=norm, origin="upper")
    else:
        cmap = sns.color_palette(cmap_name, as_cmap=True)
        try:
            cmap = cmap.copy()
            cmap.set_bad(color="#f0f0f0")
        except Exception:
            pass
        if valid.size:
            vmin = max(0.0, float(np.nanmin(valid)))
            vmax = float(np.nanpercentile(valid, 98))
            if not np.isfinite(vmin): vmin = 0.0
            if not np.isfinite(vmax) or vmax <= vmin:
                vmax = vmin + 1e-6
        else:
            vmin, vmax = 0.0, 1.0
        im = ax.imshow(np.ma.masked_invalid(mat), aspect="auto", interpolation="nearest",
                       cmap=cmap, vmin=vmin, vmax=vmax, origin="upper")

    # 7) Axes ticks/labels
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(x_groups, rotation=0, fontsize=12)
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(labels_disp, fontsize=11)
    ax.set_xlabel("Quintile", fontsize=13)
    ax.set_ylabel("Pathway", fontsize=13)

    # Light grid between cells
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.0, alpha=0.6)
    ax.tick_params(which="both", bottom=False, left=False)

    # 8) Stars for significant cells (if significance data exists)
    if len(sigs) > 0:
        for i in range(n_rows):
            for j in range(n_cols):
                if star[i, j]:
                    ax.text(
                        j, i, "★", ha="center", va="center",
                        fontsize=12, color="black",
                        path_effects=[pe.withStroke(linewidth=3, foreground="white")]
                    )

    # 9) Colorbar label from score_col
    cbar_label = "Odds ratio" if score_col.lower() in ["odds_ratio", "or"] else "Enrichment score"
    cb = fig.colorbar(im, ax=ax, pad=0.01, fraction=0.03)
    cb.set_label(cbar_label, fontsize=12)
    cb.ax.tick_params(labelsize=11)

    # 10) Save
    out_png = f"{out_base}_top{top_k}_pathways_cross_quintile_heatmap.png"
    out_pdf = f"{out_base}_top{top_k}_pathways_cross_quintile_heatmap.pdf"
    fig.tight_layout()
    fig.savefig(out_png, dpi=600, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=600, bbox_inches="tight")
    plt.close(fig)

    # Optional: also export the numeric heatmap matrix
    heat_df = pd.DataFrame(mat, index=row_paths, columns=x_groups)
    heat_csv = f"{out_base}_top{top_k}_pathways_cross_quintile_heatmap_values.csv"
    heat_df.to_csv(heat_csv, index=True)

    print("Saved cross-quintile pathways heatmap:")
    print(f"  {out_png}")
    print(f"  {out_pdf}")
    print(f"Saved underlying values CSV:")
    print(f"  {heat_csv}")

# --------------------------
# Example call (uses your existing variables)
# --------------------------
# Make sure you already have:
#   - ora (unfiltered ORA long table with 'pathway', score_col, 'group', and 'padj' if available)
#   - ordered_groups (['Q1','Q2','Q3','Q4','Q5'])
#   - top_paths_by_group (built after BH/FDR filtering)
#   - score_col (chosen by your choose_score_column)
#   - alpha (FDR threshold)
#   - base_no_ext (derived from your output_path)
#
save_cross_quintile_topk_heatmap(
    ora=ora,
    ordered_groups=ordered_groups,
    top_paths_by_group=top_paths_by_group,
    score_col=score_col,
    top_k=top_k_pathways,   # e.g., 5
    alpha=alpha,            # e.g., 0.05
    wrap_chars=40,
    cmap_name="viridis",
    out_base=base_no_ext
)

"""Per-core expression of single genes in two cell types of the Xenium data, with swarm plot and test (Fig. 3e).

Functions of Sidrah Maryam's notebook `Fig3e.ipynb` (cells "normalised", "with stats anova test" and the second
`create_swarm_plot`), copied unchanged except that the name of the core column is an argument of
`get_sample_gene_expression` (`core_col`, default "TMA_ID" as in the notebook), so that the de-identified object, where
cores are named by `xenium_core_id`, can be used.

- get_sample_gene_expression: per cell log1p(count / library size x 10,000), library size = sum of the count layer over
  all panel genes; mean per core x condition x cell type (`level2_refined`); long table.
- anova_and_pairs: one-way ANOVA over the groups and Welch t-test for the requested pairs on the per-core means
  (Benjamini-Hochberg over the pairs). Every core x cell type value is one observation: values of two cell types from
  the same core are treated as independent.
- create_swarm_plot: one point per core x cell type coloured by condition, diamonds = mean per condition, brackets with
  stars from anova_and_pairs (**** p < 1e-4, *** < 1e-3, ** < 1e-2, * < 0.05).
"""
from collections import defaultdict

import numpy as np
import pandas as pd

def get_sample_gene_expression(
    adata,
    genes,
    normalize=True,
    target_sum=1e4,
    use_layer=None,
    use_raw=False,
    core_col="TMA_ID"
):
    # 1) Choose the matrix
    A = adata
    if use_raw and (adata.raw is not None):
        A = adata.raw.to_adata()

    # keep only requested genes that exist in chosen matrix
    available_genes = [g for g in genes if g in A.var_names]
    if len(available_genes) < len(genes):
        missing = set(genes) - set(available_genes)
        print(f"Warning: {len(missing)} genes not found in data: {missing}")

    # pick matrix
    if use_layer and (use_layer in A.layers):
        X_all = A.layers[use_layer]
    else:
        X_all = A.X

    # 2) Compute per-cell library size

    if hasattr(X_all, "toarray"):   # sparse
        libsize = np.asarray(X_all.sum(axis=1)).ravel()
    else:
        libsize = X_all.sum(axis=1)

    libsize = libsize.astype(float)
    libsize[libsize == 0] = 1.0

    # 3) Pull the selected genes matrix
    Xg = A[:, available_genes].X
    if hasattr(Xg, "toarray"):
        Xg = Xg.toarray()
    Xg = Xg.astype(float)

    # 4) CP10K + log1p
    if normalize:
        Xg = np.log1p((Xg / libsize[:, None]) * float(target_sum))

    # 5) Tidy dataframe
    expr_df = pd.DataFrame(Xg, columns=available_genes, index=A.obs.index)
    expr_df[core_col] = adata.obs[core_col].values
    expr_df['condition'] = adata.obs['condition'].values
    expr_df['level2_refined'] = adata.obs['level2_refined'].values

    # 6) Average per core × condition × cell_type
    sample_avg = (expr_df
                  .groupby([core_col, 'condition', 'level2_refined'])[available_genes]
                  .mean()
                  .reset_index())

    # 7) Long format
    melted = sample_avg.melt(
        id_vars=[core_col, 'condition', 'level2_refined'],
        value_vars=available_genes,
        var_name='gene',
        value_name='expression'
    )
    return melted


def anova_and_pairs(
    expr_data,
    gene_name,
    group_col="cell_type_manual",
    value_col="expression",
    pairs=(("PT1","aPT")),
    use_welch=True,
    fdr=True
):
    """Return ANOVA p-value and pairwise test (p, q, stars) for requested pairs."""
    import numpy as np
    import pandas as pd
    from scipy import stats

    df = expr_data[expr_data['gene'] == gene_name].copy()
    if df.empty:
        raise ValueError(f"No rows for gene '{gene_name}'")

    # --- One-way ANOVA across groups (pooled over conditions) ---
    groups = [g[value_col].values for _, g in df.groupby(group_col, observed=True)]
    anova_p = np.nan
    if len(groups) >= 2 and all(len(g) >= 2 for g in groups):
        try:
            F, anova_p = stats.f_oneway(*groups)
        except Exception:
            anova_p = np.nan

    # --- Pairwise tests (Welch t-test by default; fallback to Mann-Whitney if tiny n) ---
    results = []
    for a, b in pairs:
        ga = df.loc[df[group_col] == a, value_col].dropna().values
        gb = df.loc[df[group_col] == b, value_col].dropna().values
        if len(ga) < 2 or len(gb) < 2:
            p = np.nan
        else:
            # Prefer Welch t-test; if one of the groups has all equal values, fallback to MWU
            try:
                if use_welch:
                    stat, p = stats.ttest_ind(ga, gb, equal_var=False)
                else:
                    stat, p = stats.mannwhitneyu(ga, gb, alternative="two-sided")
            except Exception:
                try:
                    stat, p = stats.mannwhitneyu(ga, gb, alternative="two-sided")
                except Exception:
                    p = np.nan
        results.append({"pair": (a, b), "p": p})

    # FDR correction (Benjamini–Hochberg)
    qvals = [np.nan]*len(results)
    if fdr:
        try:
            from statsmodels.stats.multitest import multipletests
            p_list = [r["p"] if np.isfinite(r["p"]) else 1.0 for r in results]
            reject, q_list, _, _ = multipletests(p_list, method="fdr_bh")
            qvals = q_list
        except Exception:
            pass

    # Add significance
    def stars(p):
        if p is None or not np.isfinite(p): return "n.s."
        return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "n.s."

    out = []
    for (r, q) in zip(results, qvals):
        out.append({
            "pair": r["pair"],
            "p": r["p"],
            "q": q if fdr and np.isfinite(q) else r["p"],
            "stars": stars(q if fdr and np.isfinite(q) else r["p"])
        })
    return anova_p, pd.DataFrame(out)


def add_sig_brackets(
    ax,
    pair_results_df,
    x_order,
    y_data_per_group,
    line_height=0.02,
    text_offset=0.01,
    lw=1.0,
):

    import numpy as np
    import matplotlib.pyplot as plt


    current_tops = defaultdict(float)

    for _, row in pair_results_df.iterrows():
        a, b = row["pair"]
        if a not in x_order or b not in x_order:
            continue
        i, j = x_order.index(a), x_order.index(b)
        x1, x2 = min(i, j), max(i, j)

        y1 = y_data_per_group.get(a, 0)
        y2 = y_data_per_group.get(b, 0)
        base = max(y1, y2)


        key = (x1, x2)
        current = current_tops.get(key, base)
        h = current + line_height * (1.0 + 2.0*len([k for k in current_tops if k == key]))
        current_tops[key] = h


        ax.plot([x1, x1, x2, x2], [h, h+line_height, h+line_height, h], lw=lw, c="black")
        # add the text
        label = row["stars"]
        ax.text((x1 + x2) / 2.0, h + line_height + text_offset, label,
                ha='center', va='bottom', fontsize=10)


def create_swarm_plot(
    expr_data,
    gene_name,
    value_col="expression",
    figsize=(9, 6),
    show_means=True,
    mean_marker="D",
    mean_size=45,
    mean_alpha=0.95,
    sig_pairs=(("PT1","aPT")),
    annotate_sig=True,
    fdr=True
):
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from matplotlib.lines import Line2D

    gene_data = expr_data[expr_data['gene'] == gene_name].copy()
    if gene_data.empty:
        raise ValueError(f"No rows for gene '{gene_name}' in expr_data")

    # Order axes
    celltype_order = sorted(gene_data['level2_refined'].unique())
    cond_order = sorted(gene_data['condition'].unique())

    # Categorical
    gene_data['level2_refined'] = pd.Categorical(
        gene_data['level2_refined'],
        categories=celltype_order,
        ordered=True
    )

    # Colors
    palette = sns.color_palette("Set1", len(cond_order))
    condition_colors = dict(zip(cond_order, palette))

    fig, ax = plt.subplots(figsize=figsize)

    # Swarm
    sns.swarmplot(
        data=gene_data,
        x='level2_refined',
        y=value_col,
        hue='condition',
        hue_order=cond_order,
        palette=condition_colors,
        size=5,
        alpha=0.8,
        ax=ax
    )

    # Means
    if show_means:
        offsets = np.linspace(-0.22, 0.22, len(cond_order)) if len(cond_order) > 1 else [0.0]
        x_base = np.arange(len(celltype_order))
        means = (gene_data
                 .groupby(['level2_refined', 'condition'], observed=True)[value_col]
                 .mean()
                 .reset_index())
        xtick_index = {ct: i for i, ct in enumerate(celltype_order)}
        for ci, cond in enumerate(cond_order):
            m_sub = means[means['condition'] == cond]
            xs = [xtick_index[ct] + offsets[ci] for ct in m_sub['level2_refined']]
            ax.scatter(
                xs,
                m_sub[value_col].values,
                s=mean_size,
                marker=mean_marker,
                alpha=mean_alpha,
                edgecolor='black',
                linewidths=0.4,
                color=condition_colors[cond],
                zorder=5,
                label=f"{cond} mean"
            )

    # Cosmetics
    ax.set_title(f'Gene: {gene_name}', fontsize=14)
    ax.set_xlabel('Cell Type')
    ax.set_ylabel('Average Expression')
    ax.grid(True, alpha=0.3)
    ax.set_xticklabels(celltype_order, rotation=90)

    # Legend de-dup
    handles, labels = ax.get_legend_handles_labels()
    seen = {}
    for h, l in zip(handles, labels):
        if l not in seen:
            seen[l] = h
    ax.legend(seen.values(), seen.keys(), title=None,
              bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)

    if annotate_sig and len(sig_pairs) > 0:
        try:
            anova_p, pair_df = anova_and_pairs(
                expr_data, gene_name,
                group_col="level2_refined",
                value_col=value_col,
                pairs=sig_pairs,
                use_welch=True,
                fdr=fdr
            )
            # compute top y per x group to place brackets
            y_tops = (gene_data.groupby('level2_refined', observed=True)[value_col]
                      .max().reindex(celltype_order).fillna(0).to_dict())

            add_sig_brackets(
                ax,
                pair_df,
                x_order=celltype_order,
                y_data_per_group=y_tops,
                line_height=0.05 * (gene_data[value_col].max() - gene_data[value_col].min()),
                text_offset=0.02 * (gene_data[value_col].max() - gene_data[value_col].min()),
                lw=1.0,
            )

            # Optional: show ANOVA p-value in the title
            if anova_p is not None and np.isfinite(anova_p):
                ax.set_title(f'Gene: {gene_name} (ANOVA p={anova_p:.2e})', fontsize=14)

        except Exception as e:
            print(f"[warn] Could not annotate significance for {gene_name}: {e}")

    plt.tight_layout()
    return fig

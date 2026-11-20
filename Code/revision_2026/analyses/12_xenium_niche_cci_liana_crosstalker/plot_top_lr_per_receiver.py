"""Dot plot of the top ligand-receptor pair per receiver for one sender cell type (Fig. 6e).

Plotting function of Sidrah Maryam's notebook `Fig_6e.ipynb`, unchanged except for the `dpi` argument of `savefig`
(the notebook has the typo `dpi=s200`; 200 is used) and a PDF-friendly font setting left to the caller.
Input table: one row per significant interaction with the columns source, target, gene_A (`<ligand>|L`),
gene_B (`<receptor>|R`) and LRScore (`<run>_lr_table.csv` of 03_export_tables.R, or CrossTalkeR's LR table).
"""
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap


def plot_lr_dotplot_top_per_receiver(
    df,
    sender,
    score_col="LRScore",
    source_col="source",
    target_col="target",
    ligand_col="gene_A",
    receptor_col="gene_B",
    title=None,
    subtitle="as sender",
    dot_size=450,
    figsize=None,
    savepath=None,
    x_pad=0.7,
    y_pad=0.7,
):
    """
    For a given sender, find the top-1 LR pair for EACH receiver
    """
    d = df.copy()
    d[ligand_col] = d[ligand_col].astype(str).str.replace(r"\|L$", "", regex=True)
    d[receptor_col] = d[receptor_col].astype(str).str.replace(r"\|R$", "", regex=True)
    d["pair"] = d[ligand_col] + "-" + d[receptor_col]

    d = d[d[source_col] == sender]

    # top-1 pair per receiver
    top_per_target = (
        d.sort_values(score_col, ascending=False)
        .groupby(target_col, as_index=False)
        .head(1)
    )

    # winning pairs, ordered by score
    pair_best_score = (
        d[d["pair"].isin(top_per_target["pair"])]
        .groupby("pair")[score_col]
        .max()
        .sort_values(ascending=False)
    )
    pairs = list(pair_best_score.index)

    d = d[d["pair"].isin(pairs)]
    targets = sorted(d[target_col].unique())

    d_indexed = d.set_index(["pair", target_col])[score_col]

    fig, ax = plt.subplots(figsize=figsize or (1.6 + 1.1 * len(targets), 1.2 + 0.55 * len(pairs)))
    vmin, vmax = d[score_col].min(), d[score_col].max()
    cmap = LinearSegmentedColormap.from_list("lr_score", ["#eceaf2", "#8fb8dd", "#1a5c4a"])

    for yi, pair in enumerate(pairs):
        for xi, tgt in enumerate(targets):
            if (pair, tgt) in d_indexed.index:
                val = d_indexed.loc[(pair, tgt)]
                if isinstance(val, pd.Series):
                    val = val.max()
                color = cmap((val - vmin) / (vmax - vmin) if vmax > vmin else 0.5)
                ax.scatter(xi, yi, s=dot_size, color=color, edgecolor="black",
                           linewidth=1.2, zorder=3, clip_on=True)

    for xi in range(len(targets)):
        ax.axvline(xi, color="lightgray", linestyle=":", linewidth=1, zorder=1)

    ax.set_xlim(-x_pad, len(targets) - 1 + x_pad)
    ax.set_ylim(-y_pad, len(pairs) - 1 + y_pad)
    ax.set_xticks(range(len(targets)))
    ax.set_xticklabels(targets, rotation=45, ha="right", fontsize=13)
    ax.set_yticks(range(len(pairs)))
    ax.set_yticklabels(pairs, fontsize=13)
    ax.invert_yaxis()
    ax.set_xlabel("Receiver", fontsize=14, labelpad=10)

    for spine in ax.spines.values():
        spine.set_linewidth(1.3)

    if title:
        fig.suptitle(title, color="#b5432c", fontsize=26, y=1.05)
    ax.set_title(f"{sender} {subtitle}", fontsize=17, pad=10)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.06, pad=0.04)
    cbar.set_label(score_col, fontsize=13)

    fig.tight_layout()
    if savepath:
        fig.savefig(savepath, dpi=200, bbox_inches="tight")

    return fig, ax

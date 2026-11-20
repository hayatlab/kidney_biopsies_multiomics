"""Per-cell quality metrics of the Xenium data, shown as one box per core and run (Supplementary Fig. 11a-c).

Method of Sidrah Maryam's notebook `Supp_fig_4a-c.ipynb` (cells "collect QC metrics" / "one figure per page"):
`scanpy.pp.calculate_qc_metrics` on the raw counts gives, per cell, the total transcripts (`total_counts`, sum over the
468 panel genes) and the unique transcripts (`n_genes_by_counts`, genes with at least one count); the cell area is the
`cell_area` column of the Xenium cell table (square micrometres). Each metric is drawn with `seaborn.boxplot`, one box
per core (box = quartiles, whiskers = most extreme cell within 1.5 x IQR, circles = cells beyond).

Differences from the notebook as received (20 Nov 2026, Claude Fable 5.2 for Leon Frommer):
- the notebook loops over three runs read with spatialdata-io and writes one PDF per run; the published figure shows the
  five runs of the analysed object, so the metrics are computed once on that object and drawn run by run;
- cores are named by public ids (sample id, condition, `xenium_core_id`) instead of the internal core label, and every
  tick label is taken from the same table row as its box.
"""
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns

METRICS = [                                            # panel, column, title, y label (titles and labels of the notebook)
    ("a", "Total_Transcripts", "Total Transcripts per Cell", "Counts"),
    ("b", "Unique_Transcripts", "Unique Transcripts per Cell", "Counts"),
    ("c", "Cell_Area", "Cell Area", "Area"),
]


def qc_table(adata, core_ids, sample_ids, run_col="sample", condition_col="condition"):
    """One row per cell: run, public core and sample id, condition, label and the three metrics.
    `core_ids` / `sample_ids`: public ids per cell (xenium_core_id, K number)."""
    a = adata.copy()
    sc.pp.calculate_qc_metrics(a, percent_top=None, inplace=True)
    qc = pd.DataFrame({
        "run": a.obs[run_col].astype(str).to_numpy(),
        "xenium_core_id": list(core_ids),
        "sample_id": list(sample_ids),
        "condition": a.obs[condition_col].astype(str).to_numpy(),
        "Total_Transcripts": a.obs["total_counts"].to_numpy(),
        "Unique_Transcripts": a.obs["n_genes_by_counts"].to_numpy(),
        "Cell_Area": a.obs["cell_area"].to_numpy(),
    })
    qc["core_label"] = qc["sample_id"] + " " + qc["condition"] + " " + qc["xenium_core_id"]
    return qc


def plot_metric(qc, metric, title, ylabel, runs, width=10, row_height=4.2):
    """One figure for one metric: one row per run, one box per core (ordered by public core id)."""
    fig, axes = plt.subplots(len(runs), 1, figsize=(width, row_height * len(runs)))
    for ax, run in zip(axes, runs):
        d = qc[qc["run"] == run]
        order = d.drop_duplicates("xenium_core_id").sort_values("xenium_core_id")["core_label"].tolist()
        sns.boxplot(data=d, x="core_label", y=metric, order=order, ax=ax)
        ax.set_title(f"{title} - Sample: {run}", fontsize=11, fontweight="bold")
        ax.set_ylabel(ylabel)
        ax.set_xlabel("")
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order, rotation=90, fontsize=7)
    fig.tight_layout()
    return fig

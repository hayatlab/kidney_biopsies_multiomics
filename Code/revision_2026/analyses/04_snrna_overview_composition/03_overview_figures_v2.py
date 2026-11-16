#!/usr/bin/env python
"""Dataset-overview figures, version 2 (2026-09-17): redesigned, every figure as PNG + PDF + SVG.

Supersedes 02_overview_figures.py (kept for reference; its outputs are moved to figures/archive_v1_20260902/).
Changes, as asked by Leon ("this is not nice; all figures also as svg and pdf"):
  * metadata from the corrected sample table (vivien_takeover/outputs/metadata/sample_metadata_corrected.csv), so
    K68 / K95 / K98 carry their disease and sex instead of "nan"; disease names as in the explorer (IgA -> IgAN)
  * project palettes (kidney-biopsy-explorer public/data/meta/labels.json) for disease, level1 and level2
  * UMAP: 2 x 2 grid (cell type with labels on the clusters, disease, sex, sample), points drawn in random order
    (seed 42) so no group is systematically painted over another, small points, legends beside each panel, small
    UMAP axis arrows instead of heavy frames; the four panels are also written as single-panel figures
  * per-sample plots: samples grouped by disease (then K-number), bars / violins coloured by disease, log10 axis for
    transcripts and genes per nucleus, composition bars with a disease strip on top
  * scatter layers rasterised at 300 dpi inside PDF / SVG (363,656 points), all text and axes stay vector
Outputs: figures/{png,pdf,svg}/<figure>.{png,pdf,svg}; figures/tables/ (numbers behind each figure); figures/run_meta_v2.json
Env: agent_env (matplotlib 3.10, seaborn 0.13).  Human GRCh38, HGNC symbols.  Seed 42.
"""
from kbpaths import KB_ANALYSIS, KB_PROJECT  # repository path config (config/kbpaths.py)
import json
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

HERE = Path(f"{KB_ANALYSIS}/snrna_overview_20260831")
META = Path(f"{KB_ANALYSIS}/vivien_takeover/outputs/metadata/sample_metadata_corrected.csv")
LABELS = Path(f"{KB_PROJECT}/kidney-biopsy-explorer/public/data/meta/labels.json")
FIG = HERE / "figures"
SEED = 42
DISEASE_MAP = {"control": "Control", "control_AKI": "control_AKI", "DN": "DN", "FSGS": "FSGS", "IgA": "IgAN", "MN": "MN", "TIN": "TIN"}
SEX_MAP = {"f": "female", "m": "male"}
SEX_COLORS = {"female": "#CC6677", "male": "#4477AA"}

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 11, "axes.labelsize": 9,
                     "legend.fontsize": 8, "legend.title_fontsize": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "svg.fonttype": "none", "pdf.fonttype": 42})


def save(fig, name):
    for ext, kw in [("png", {"dpi": 300}), ("pdf", {"dpi": 300}), ("svg", {"dpi": 300})]:
        (FIG / ext).mkdir(parents=True, exist_ok=True)
        fig.savefig(FIG / ext / f"{name}.{ext}", bbox_inches="tight", facecolor="white", **kw)
    plt.close(fig)
    print(f"[save] {name}", flush=True)


def natural_k(s):
    return int(str(s).lstrip("K"))


def sample_palette(samples):
    base = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors) + list(plt.get_cmap("tab20c").colors) \
        + list(plt.get_cmap("Set1").colors)
    return {s: base[i % len(base)] for i, s in enumerate(samples)}


def umap_arrows(ax):
    """Small L-shaped UMAP 1 / UMAP 2 indicator in the lower-left corner instead of a full frame."""
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    kw = dict(xycoords="axes fraction", textcoords="axes fraction", arrowprops=dict(arrowstyle="-|>", lw=0.8, color="#444444"))
    ax.annotate("", xy=(0.14, 0.02), xytext=(0.02, 0.02), **kw)
    ax.annotate("", xy=(0.02, 0.14), xytext=(0.02, 0.02), **kw)
    ax.text(0.08, -0.01, "UMAP 1", transform=ax.transAxes, ha="center", va="top", fontsize=7, color="#444444")
    ax.text(-0.01, 0.08, "UMAP 2", transform=ax.transAxes, ha="right", va="center", rotation=90, fontsize=7, color="#444444")


def scatter_cat(ax, df, col, colors, order, title, legend=True, legend_cols=1, counts=True, labels_on_plot=False, s=0.25):
    ax.scatter(df.UMAP1, df.UMAP2, c=[colors[v] for v in df[col]], s=s, linewidths=0, rasterized=True)
    ax.set_title(title, loc="left", fontweight="bold")
    umap_arrows(ax)
    pad_x, pad_y = 0.02 * np.ptp(df.UMAP1), 0.02 * np.ptp(df.UMAP2)
    ax.set_xlim(df.UMAP1.min() - pad_x, df.UMAP1.max() + pad_x)
    ax.set_ylim(df.UMAP2.min() - pad_y, df.UMAP2.max() + pad_y)
    ax.set_aspect("equal", adjustable="box", anchor="NW")
    if labels_on_plot:
        med = df.groupby(col, observed=True)[["UMAP1", "UMAP2"]].median()
        for lab, r in med.iterrows():
            ax.text(r.UMAP1, r.UMAP2, lab, fontsize=7.5, fontweight="bold", ha="center", va="center", color="black",
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75))
    if legend:
        n = df[col].value_counts()
        handles = [Line2D([], [], marker="o", ls="", markersize=5, markerfacecolor=colors[v], markeredgewidth=0,
                          label=f"{v} ({n.get(v, 0):,})" if counts else str(v)) for v in order if n.get(v, 0) > 0]
        ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False, ncol=legend_cols,
                  handletextpad=0.2, columnspacing=0.8, borderaxespad=0.2)


# Display names (13 Nov 2026): the co-authors renamed the level-2 label dPT to aPT2 in all figures. The object keeps
# "dPT"; LEVEL2_DISPLAY = {"dPT": "aPT2"} (set by the Supp. Fig. 2 entry point) prints the new name with the same colour.
LEVEL2_DISPLAY = {}


def main():
    rng = np.random.default_rng(SEED)
    labels = json.loads(LABELS.read_text())
    if LEVEL2_DISPLAY:
        l2 = labels["level2"]
        l2["order"] = [LEVEL2_DISPLAY.get(c, c) for c in l2["order"]]
        l2["colors"] = {LEVEL2_DISPLAY.get(c, c): v for c, v in l2["colors"].items()}
    df = pd.read_parquet(HERE / "umap_and_obs.parquet")
    meta = pd.read_csv(META).set_index("sample")
    df["sample"] = df["sample"].astype(str)
    assert df["sample"].isin(meta.index).all()
    before = {"sex_nan": int((df["sex_clinical"].astype(str) == "nan").sum()), "disease_nan": int((df["disease"].astype(str) == "nan").sum())}
    df["disease_c"] = df["sample"].map(meta["disease"]).map(DISEASE_MAP)
    df["sex_c"] = df["sample"].map(meta["sex_clinical"]).map(SEX_MAP)
    assert df["disease_c"].notna().all() and df["sex_c"].notna().all(), "unmapped disease / sex after correction"
    df = df.iloc[rng.permutation(len(df))].reset_index(drop=True)      # random drawing order

    dis_order = labels["disease"]["order"]
    dis_col = labels["disease"]["colors"]
    l2_order, l2_col = labels["level2"]["order"], labels["level2"]["colors"]
    l1 = labels.get("level1", {})
    l1_order = l1.get("order") or sorted(df["level1"].astype(str).unique())
    l1_col = l1.get("colors") or dict(zip(l1_order, plt.get_cmap("tab10").colors))
    df["level1"], df["level2"] = df["level1"].astype(str), df["level2"].astype(str)
    if LEVEL2_DISPLAY:
        df["level2"] = df["level2"].replace(LEVEL2_DISPLAY)
    assert set(df.level2) <= set(l2_order) and set(df.level1) <= set(l1_col), "labels missing from the palette"

    per_sample = df.drop_duplicates("sample").set_index("sample")["disease_c"]
    samples = sorted(per_sample.index, key=lambda s_: (dis_order.index(per_sample[s_]), natural_k(s_)))
    s_col = sample_palette(sorted(per_sample.index, key=natural_k))

    # archive version-1 outputs once
    arch = FIG / "archive_v1_20260902"
    if not arch.exists():
        arch.mkdir(parents=True)
        for f in ["umap_plots.pdf", "sample_plots.pdf", "run_meta.json"]:
            if (FIG / f).exists():
                shutil.move(str(FIG / f), arch / f)
        if (FIG / "png").exists():
            shutil.move(str(FIG / "png"), arch / "png")

    # ---- UMAPs -----------------------------------------------------------------------------------------------------
    panels = [("level2", l2_col, l2_order, "Cell type (level 2)", dict(labels_on_plot=True, legend_cols=1)),
              ("disease_c", dis_col, dis_order, "Disease", dict()),
              ("sex_c", SEX_COLORS, ["female", "male"], "Sex", dict()),
              ("sample", s_col, sorted(per_sample.index, key=natural_k), f"Sample ({per_sample.size} biopsies)",
               dict(legend_cols=4, counts=False))]
    fig, axes = plt.subplots(2, 2, figsize=(15, 10.2), gridspec_kw={"wspace": 0.62, "hspace": 0.1, "top": 0.93, "bottom": 0.03,
                                                                      "left": 0.03, "right": 0.86})
    for ax, (col, colors, order, title, kw) in zip(axes.ravel(), panels):
        scatter_cat(ax, df, col, colors, order, title, **kw)
    fig.suptitle(f"snRNA-seq: {len(df):,} nuclei from {per_sample.size} kidney biopsies", x=0.03, y=0.985, ha="left", fontsize=13, fontweight="bold")
    save(fig, "umap_overview")
    for col, colors, order, title, kw in panels:
        fig, ax = plt.subplots(figsize=(8.6 if col == "sample" else 7.0, 4.9))
        scatter_cat(ax, df, col, colors, order, title, **kw)
        save(fig, f"umap_{col.replace('_c', '')}")

    # ---- per-sample helpers -----------------------------------------------------------------------------------------
    xs = np.arange(len(samples))

    def disease_strip(ax_strip):
        ax_strip.bar(xs, 1, width=1.0, color=[dis_col[per_sample[s_]] for s_ in samples], linewidth=0)
        ax_strip.set_xlim(-0.5, len(samples) - 0.5); ax_strip.set_ylim(0, 1)
        ax_strip.axis("off")

    def disease_legend(ax, loc="upper left", anchor=(1.0, 1.0)):
        present = [d for d in dis_order if d in set(per_sample)]
        n = per_sample.value_counts()
        return ax.legend(handles=[Patch(color=dis_col[d], label=f"{d} ({n[d]})") for d in present], title="Disease (samples)",
                         frameon=False, loc=loc, bbox_to_anchor=anchor)

    def sample_axis(ax):
        ax.set_xticks(xs, samples, rotation=90, fontsize=7)
        ax.set_xlim(-0.6, len(samples) - 0.4)
        ax.set_xlabel("Sample (grouped by disease)")

    tables = FIG / "tables"
    tables.mkdir(parents=True, exist_ok=True)

    # ---- QC violins ---------------------------------------------------------------------------------------------------
    qc = [("nCount_RNA", "Transcripts (UMIs) per nucleus", True), ("nFeature_RNA", "Genes detected per nucleus", True),
          ("pct_counts_mt", "Mitochondrial transcripts (% of UMIs)", False)]
    for metric, ylab, logy in qc:
        fig, ax = plt.subplots(figsize=(15, 4.6))
        d = df[["sample", metric, "disease_c"]].copy()
        if logy:
            d = d[d[metric] > 0]
        sns.violinplot(data=d, x="sample", y=metric, order=samples, hue="sample", hue_order=samples, legend=False,
                       palette={s_: dis_col[per_sample[s_]] for s_ in samples}, inner="quart", cut=0, linewidth=0.5,
                       density_norm="width", log_scale=logy, ax=ax)
        med = d.groupby("sample")[metric].median().reindex(samples)
        ax.set_ylabel(ylab + (" (log scale)" if logy else ""))
        sample_axis(ax)
        ax.grid(axis="y", color="#e5e5e5", lw=0.6); ax.set_axisbelow(True)
        disease_legend(ax)
        ax.set_title(f"{ylab}: distribution per sample (lines: quartiles); median of sample medians {med.median():,.2f}" if not logy else
                     f"{ylab}: distribution per sample (lines: quartiles); median of sample medians {med.median():,.0f}", loc="left")
        save(fig, f"violin_{metric}")
        d.groupby("sample")[metric].describe().reindex(samples).to_csv(tables / f"per_sample_{metric}.csv")

    # ---- nuclei per sample ----------------------------------------------------------------------------------------------
    counts = df["sample"].value_counts().reindex(samples)
    fig, ax = plt.subplots(figsize=(15, 4.2))
    ax.bar(xs, counts.to_numpy(), color=[dis_col[per_sample[s_]] for s_ in samples], width=0.8, linewidth=0)
    med_line = ax.axhline(counts.median(), color="#555555", lw=0.8, ls="--", label=f"median {counts.median():,.0f} nuclei")
    ax.set_ylabel("Nuclei per sample")
    sample_axis(ax)
    ax.grid(axis="y", color="#e5e5e5", lw=0.6); ax.set_axisbelow(True)
    ax.add_artist(disease_legend(ax))
    ax.legend(handles=[med_line], frameon=False, loc="lower left", bbox_to_anchor=(1.0, 0.0))
    ax.set_title(f"Nuclei per sample ({len(df):,} nuclei, {len(samples)} samples)", loc="left")
    save(fig, "bar_nuclei_per_sample")
    counts.rename("n_nuclei").to_frame().assign(disease=per_sample.reindex(samples)).to_csv(tables / "nuclei_per_sample.csv")

    # ---- composition ----------------------------------------------------------------------------------------------------
    for level, order, colors, title in [("level1", l1_order, l1_col, "Cell-type composition per sample (level 1)"),
                                        ("level2", l2_order, l2_col, "Cell-type composition per sample (level 2)")]:
        ct = pd.crosstab(df["sample"], df[level], normalize="index").reindex(index=samples)
        ct = ct[[c for c in order if c in ct.columns]]
        fig = plt.figure(figsize=(15, 5.4))
        gs = fig.add_gridspec(2, 1, height_ratios=[0.035, 1], hspace=0.02)
        ax_s, ax = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
        disease_strip(ax_s)
        bottom = np.zeros(len(samples))
        for c in ct.columns:
            ax.bar(xs, ct[c].to_numpy(), bottom=bottom, width=0.85, color=colors[c], linewidth=0, label=c)
            bottom += ct[c].to_numpy()
        ax.set_ylim(0, 1); ax.set_ylabel("Fraction of nuclei")
        sample_axis(ax)
        leg1 = ax.legend(handles=[Patch(color=colors[c], label=c) for c in ct.columns][::-1], title="Cell type", frameon=False,
                         loc="upper left", bbox_to_anchor=(1.0, 1.02))
        ax.add_artist(leg1)
        disease_legend(ax, anchor=(1.1, 1.02), loc="upper left")
        ax_s.set_title(title + " (top strip: disease)", loc="left")
        save(fig, f"stacked_{level}")
        ct.to_csv(tables / f"composition_{level}_per_sample.csv")

    (FIG / "run_meta_v2.json").write_text(json.dumps({
        "script": "03_overview_figures_v2.py", "date": "2026-09-17", "n_nuclei": int(len(df)), "n_samples": int(per_sample.size),
        "metadata": str(META), "nan_in_object_replaced_by_corrected_table": before, "palettes": str(LABELS),
        "formats": ["png (300 dpi)", "pdf (vector text, rasterised points)", "svg (vector text, rasterised points)"],
        "seed": SEED, "versions": {"matplotlib": matplotlib.__version__, "seaborn": sns.__version__, "pandas": pd.__version__}},
        indent=2))
    print("[done]", flush=True)


if __name__ == "__main__":
    main()

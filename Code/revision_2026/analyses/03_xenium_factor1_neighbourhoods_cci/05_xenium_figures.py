#!/usr/bin/env python
"""Stage D -- the Xenium figure panels (Vivien takeover).

  fig3d   mean Factor1 per TMA core vs eGFR, spatial thumbnail per core, OLS + 95% CI
          <- plot_egfr_f2_on_spatial/4_scatterplot_of_cores.py
  fig3g   distance-bounded Voronoi graph of one core + 2-hop neighbourhood insets
          <- nh_analysis/1_compute_nh.py
  fig3h   cell-type composition of 2-hop neighbourhoods vs Factor1 (stacked spline ribbon)
          <- nh_analysis/2_plot_nh_compo_ct.py
  ed3b    2-hop neighbourhood cell-type x cell-type proportion heatmap
          <- nh_analysis/2_plot_nh_compo_ct.py
  ed3f    niche composition of 2-hop neighbourhoods vs Factor1 (stacked spline ribbon)
          <- nh_analysis/3_plot_nh_compo_niche.py
  suppl7  per-core cell-boundary maps: cell type, Factor1, niches
          <- cell_boundaries/2_plot_cell_boundaries.py + 3_plot_niches_cell_boundaries.py
  (`<-` names the first-submission script that a panel function replaces; those scripts are no longer in this
  repository.)

Ribbon method is the original's: per level, LSQUnivariateSpline(k=3) with 8 interior knots
over the crop range, clip to [0,1], stack with fill_between. The original HARD-CODED the
crop to Factor 2's range ([-4, 8] / [-2, 5]); Factor 1 is on a different scale, so the crop
is taken from the 1st/99th percentile of the observed Factor1 and the range is printed and
recorded rather than copied across.

Labels: level2_refined (16 cell types) and latent_leiden_0.4 (10 niches), per Leon 2026-09-01.
Human GRCh38, HGNC symbols. Seed 42.
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
from cores import core_label  # public core id -> internal label (config/cores.py)
import argparse, json, re, warnings
from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scipy.stats as ss
import seaborn as sns
import statsmodels.api as sm
from matplotlib.colors import ListedColormap
from scipy.interpolate import LSQUnivariateSpline

HERE = Path(f"{KB_ANALYSIS}/vivien_takeover")
OUT = HERE / "outputs"
FIG = HERE / "figures"
SEED = 42
FACTOR = "Factor1_ulm"
CT_COL, NICHE_COL = "level2_refined", "latent_leiden_0.4"
COND_COLORS = {"FSGS": "#1f77b4", "IgAN": "#ff7f0e", "MN": "#2ca02c",
               "TIN": "#d62728", "DN": "#9467bd", "CTRL": "#8c564b"}

# Fig 3g: the core Vivien portrayed in the published panel (Leon, 2026-09-07).
FIG3G_CORE = core_label("XEN085")  # Utrecht1 E4, K16 TIN (eGFR 5)
# Suppl 7: the two cores Leon asked for (2026-09-07), replacing the automatic
# lowest-/highest-eGFR pick.
SUPPL7_CORES = [core_label("XEN116"), core_label("XEN008")]  # K18 TIN (Utrecht2 F2) and the Apollo d0 B5 control core
# Fig 3h / ED 3f ribbon crop. Vivien hard-coded Factor 2's range with no criterion;
# Factor 1 is on a different scale, so the crop is the largest contiguous stretch of
# 0.25-unit bins each holding >= 500 2-hop neighbourhoods -- 99.07 % of all 258,375
# neighbourhoods, and the weakest bin inside still draws on 58 of the 119 cores.
# Evidence: scripts/13_crop_justification.py -> figures/xenium/crop_justification_factor1.pdf
# Set to None to fall back to the 1st-99th percentile.
RIBBON_CROP = (-7.25, 5.25)


def palette(levels, name="tab20"):
    base = plt.get_cmap(name)(range(20))
    extra = sns.color_palette("husl", max(len(levels) - 20, 0))
    cols = list(base) + list(extra)
    return {lv: cols[i] for i, lv in enumerate(levels)}, ListedColormap(cols[:len(levels)])


def corrected_core_egfr(obs):
    """Per-core eGFR with the object's broadcast values replaced by the source workbooks.

    xenium_labeled.h5ad's `EGFR` substitutes one group value for many donors (22 CTRL cores
    all read 67.0 = the first Appollo control; 15 FSGS all read 58.0). See
    scripts/10_build_authoritative_egfr.py and outputs/clinical_column_audit.csv.
    Metadata wins where it covers the core; the object value is kept only for the 32 cores
    from the two NEW Xenium runs, which the workbooks predate.
    """
    obj = (obs.groupby("TMA_ID", observed=True)["EGFR"].first()
              .pipe(pd.to_numeric, errors="coerce"))
    auth = pd.read_csv(OUT / "authoritative_core_egfr.csv", index_col=0)["eGFR"]
    key = pd.Index([re.sub(r"\s+", " ", str(i)).strip() for i in obj.index])
    fixed = auth.reindex(key)
    out = pd.Series(np.where(fixed.notna(), fixed.to_numpy(), obj.to_numpy()),
                    index=obj.index, name="eGFR")
    n_auth = int(fixed.notna().sum())
    n_changed = int((out != obj).sum())
    print(f"[egfr] {n_auth}/{len(out)} cores from source metadata, {n_changed} values "
          f"corrected; distinct {out.nunique()} (object had {obj.nunique()}), "
          f"max cores sharing one value {int(out.value_counts().max())} "
          f"(was {int(obj.value_counts().max())})", flush=True)
    return out, obj, n_auth, n_changed


# ---------------------------------------------------------------- fig 3d
def fig3d(obs, outdir, stem="fig3d_factor1_vs_egfr", write_table=True):
    core = (obs.groupby("TMA_ID", observed=True)
               .agg(F=(FACTOR, "mean"), eGFR_obj=("EGFR", "first"),
                    cond=("condition", "first"), n=("cell_id", "size")))
    core["eGFR_obj"] = pd.to_numeric(core["eGFR_obj"], errors="coerce")

    # The object's EGFR column is CORRUPTED: one group value broadcast across many donors
    # (22 CTRL cores all read 67.0, 15 FSGS all read 58.0; Vivien's object has 19 distinct
    # values for those same 22 CTRL cores). Use the source metadata workbooks instead
    # (scripts/10_build_authoritative_egfr.py), keeping the object value only for cores the
    # workbooks do not cover -- the 32 from the two NEW Xenium runs, whose values are
    # distinct and plausible. TMA_ID whitespace is normalised before joining.
    fixed, _, n_auth, n_changed = corrected_core_egfr(obs)
    core["eGFR"] = fixed.reindex(core.index).to_numpy()
    core = core.dropna(subset=["eGFR"])
    x, y = core["eGFR"].to_numpy(float), core["F"].to_numpy(float)
    r, p = ss.pearsonr(x, y)
    m = sm.OLS(y, sm.add_constant(x)).fit()
    ci = m.conf_int(0.05)[1]
    bp = sm.stats.diagnostic.het_breuschpagan(m.resid, m.model.exog)[1]
    sw = ss.shapiro(m.resid)[1]
    print(f"[fig3d] n={len(core)} cores | pearson r={r:+.3f} (p={p:.3g}) | "
          f"slope={m.params[1]:+.4f} 95% CI [{ci[0]:+.4f}, {ci[1]:+.4f}] | R2={m.rsquared:.3f} | "
          f"Breusch-Pagan p={bp:.3g}, Shapiro-Wilk p={sw:.3g}", flush=True)

    xs = np.linspace(x.min(), x.max(), 100)
    pred = m.get_prediction(sm.add_constant(xs)).summary_frame(alpha=0.05)

    # Square canvas (Leon, 2026-09-07). The per-core marker is NOT drawn: each core is
    # represented by its spatial thumbnail, and a dot underneath it was only ever visible
    # as a speck in the middle of the thumbnail. Condition is carried by the thumbnail
    # frame colour, so the legend is built from proxy handles that match those frames.
    fig, ax = plt.subplots(figsize=(16, 16))
    ax.plot(xs, pred["mean"], color="black", lw=1.5, zorder=3)
    ax.fill_between(xs, pred["mean_ci_lower"], pred["mean_ci_upper"],
                    color="grey", alpha=0.3, lw=0, zorder=2)
    # limits must be set from the data, not autoscaled from a scatter that no longer
    # exists; the padding also keeps edge thumbnails from being clipped
    xpad, ypad = 0.07 * np.ptp(x), 0.07 * np.ptp(y)
    ax.set_xlim(x.min() - xpad, x.max() + xpad)
    ax.set_ylim(y.min() - ypad, y.max() + ypad)
    conds = [c for c in core["cond"].astype(str).unique()]
    handles = [plt.Line2D([], [], marker="s", linestyle="none", markersize=9,
                          markerfacecolor="none", markeredgewidth=1.6,
                          markeredgecolor=COND_COLORS.get(c, "grey"), label=c)
               for c in sorted(conds)]

    ax.set_xlabel("eGFR (ml/min/1.73m$^2$)", fontsize=13)
    ax.set_ylabel(f"Mean {FACTOR} in all cells per core", fontsize=13)
    ax.set_title(f"correlation: {r:.2f}\nslope: {m.params[1]:.3f}  (95% CI {ci[0]:.3f} to {ci[1]:.3f})\n"
                 f"p-value: {p:.2e}", loc="right", fontsize=10)
    ax.legend(handles=handles, frameon=False, fontsize=11, loc="upper right",
              title="core outline")
    ax.set_aspect("auto")
    sns.despine(ax=ax)
    fig.tight_layout()

    # --- spatial thumbnail per core, at its scatter position ------------------
    # Placed in FIGURE coordinates via fig.add_axes, NOT inset_axes with
    # bbox_transform=ax.transData: that transform is re-evaluated during the two-pass
    # render that bbox_inches="tight" triggers, which left every inset empty in the PDF
    # and collapsed their contents into one corner blob. Layout is frozen first
    # (tight_layout above, explicit limits below) and the figure is saved with a fixed
    # bbox so the positions computed here are the ones that get drawn.
    fig.canvas.draw()          # limits were set explicitly above; freeze the layout
    vmin, vmax = obs[FACTOR].quantile([0.02, 0.98])
    side = 0.050                      # thumbnail size as a fraction of the figure
    n_drawn = 0
    for tma, row in core.iterrows():
        pts = obs.loc[obs["TMA_ID"] == tma, ["spatial_x", "spatial_y", FACTOR]]
        if pts.empty:
            continue
        # data -> display -> figure fraction
        xd, yd = ax.transData.transform((row["eGFR"], row["F"]))
        xf, yf = fig.transFigure.inverted().transform((xd, yd))
        iax = fig.add_axes([xf - side / 2, yf - side / 2, side, side], zorder=5)
        iax.scatter(pts["spatial_x"], pts["spatial_y"], c=pts[FACTOR], cmap="viridis",
                    s=0.14, vmin=vmin, vmax=vmax, linewidths=0, rasterized=True)
        iax.set_xticks([]); iax.set_yticks([]); iax.set_aspect("equal")
        iax.invert_yaxis(); iax.patch.set_alpha(0.0)
        for sp in iax.spines.values():
            sp.set_color(COND_COLORS.get(str(row["cond"]), "grey"))
            sp.set_linewidth(1.1)
        n_drawn += 1
    print(f"[fig3d] {n_drawn} spatial thumbnails placed in figure coordinates", flush=True)

    for ext in ("pdf", "png"):
        fig.savefig(outdir / f"{stem}.{ext}", dpi=300)
    plt.close(fig)
    core["donor"] = obs.groupby("TMA_ID", observed=True)["H.number"].first().reindex(core.index)
    if write_table:   # rows are internal core labels, with donor numbers
        core.to_csv(outdir / "fig3d_per_core.csv")

    # Cores are NOT independent: a TMA carries several punches per donor (up to 22 share one
    # eGFR here), so the per-core test overstates the effective N. H.number is populated for
    # only 31 donors, so cores lacking it are grouped by (condition, eGFR) -- an approximation.
    k = core["donor"].astype(str)
    miss = k.isin(["nan", "", "None"])
    k = k.mask(miss, core["cond"].astype(str) + "_eGFR" + core["eGFR"].astype(str))
    don = core.assign(k=k).groupby("k").agg(F=("F", "mean"), eGFR=("eGFR", "first"))
    dr, dp = ss.pearsonr(don["eGFR"].to_numpy(float), don["F"].to_numpy(float))
    dm = sm.OLS(don["F"].to_numpy(float), sm.add_constant(don["eGFR"].to_numpy(float))).fit()
    dci = dm.conf_int(0.05)[1]
    print(f"[fig3d] donor-level sensitivity ({len(don)} donors from {len(core)} cores, "
          f"{int(miss.sum())} without H.number): r={dr:+.3f} (p={dp:.3g}) "
          f"slope={dm.params[1]:+.4f} 95% CI [{dci[0]:+.4f}, {dci[1]:+.4f}]", flush=True)

    return {"n_cores": int(len(core)), "pearson_r": r, "pearson_p": p,
            "egfr_source": {"from_metadata": n_auth, "from_object": int(len(core) - n_auth),
                            "values_corrected": n_changed},
            "donor_level": {"n_donors": int(len(don)), "pearson_r": float(dr),
                            "pearson_p": float(dp), "slope": float(dm.params[1]),
                            "slope_ci95": [float(dci[0]), float(dci[1])],
                            "note": "cores are clustered within donors; H.number covers only "
                                    "31 donors, the rest grouped by (condition, eGFR)"},
            "slope": float(m.params[1]), "slope_ci95": [float(ci[0]), float(ci[1])],
            "r2": float(m.rsquared), "breusch_pagan_p": float(bp), "shapiro_p": float(sw)}


# ---------------------------------------------------------------- fig 3g
def pick_neighbourhood(gc, g, root_type="aPT", lo=12, hi=25):
    """Deterministic pick of the 2-hop neighbourhood to inset.

    The published panel shows an aPT-dominant neighbourhood that also contains
    Fib/MyoFib, Immune and PC -- i.e. an aPT root sitting on a cell-type interface.
    Vivien's own script simply plotted the first ten cells of the core, so his exact
    root cannot be recovered from a different cell set. The rule used here reproduces
    that character explicitly: among roots of `root_type` whose 2-hop neighbourhood
    holds between `lo` and `hi` cells, take the one with the most distinct cell types,
    breaking ties on the root type's own fraction (so it stays aPT-dominant) and then
    on cell id, so the choice is reproducible.
    """
    cand = gc.index[gc[CT_COL].astype(str) == root_type]
    if not len(cand):
        cand = gc.index
    best = None
    for i in cand:
        if i not in g:
            continue
        mem = list(nx.single_source_shortest_path_length(g, i, cutoff=2).keys())
        if not (lo <= len(mem) <= hi):
            continue
        types = gc.loc[mem, CT_COL].astype(str)
        key = (types.nunique(), float((types == root_type).mean()), str(i))
        if best is None or key[:2] > best[0][:2]:
            best = (key, i, mem)
    if best is None:                       # no root inside the size window: fall back
        i = max(cand, key=lambda i: g.degree(i) if i in g else -1)
        mem = list(nx.single_source_shortest_path_length(g, i, cutoff=2).keys())
        return i, mem
    return best[1], best[2]


def draw_net(ax, subgraph, pos, root=None, lw=1.0, s=10, s_root=80, color="white"):
    """Vivien's 'line net': white Voronoi edges + centroids, root centroid enlarged.

    His neighbourhood was compact, so plain white worked. Here the 2-hop set includes cells
    separated by gaps, and a white line crossing the white figure background disappears --
    so every white stroke carries a thin dark outline (path effect) and stays visible both
    over a coloured polygon and over the background.
    """
    halo = [pe.Stroke(linewidth=lw + 0.9, foreground="0.35"), pe.Normal()]
    for u, v in subgraph.edges():
        ax.plot(*zip(pos[u], pos[v]), color=color, lw=lw, zorder=6, path_effects=halo,
                solid_capstyle="round")
    xs = [pos[i][0] for i in subgraph.nodes()]; ys = [pos[i][1] for i in subgraph.nodes()]
    ax.scatter(xs, ys, color=color, s=s, zorder=7, linewidths=0.5, edgecolors="0.35")
    if root is not None and root in pos:
        ax.scatter(*pos[root], color=color, s=s_root, zorder=8, linewidths=0.8,
                   edgecolors="0.2")


def fig3g(gdf, graphs, outdir, core=None, stem_suffix="", stem="fig3g_voronoi_2hop", core_name=None,
          factor_cbar_kwds=None):
    if core is None:
        core = FIG3G_CORE if FIG3G_CORE in graphs else gdf["TMA_ID"].value_counts().idxmax()
    name = core if core_name is None else core_name   # shown in the title and the log
    if core not in graphs:
        raise SystemExit(f"core {name!r} is not in the Voronoi graphs "
                         f"({len(graphs)} available)")
    gc = gdf[gdf["TMA_ID"] == core]
    g = graphs[core]
    levels = list(pd.Categorical(gdf[CT_COL]).categories)
    pal, cmap = palette(levels)

    cent = gc.geometry.centroid
    pos = {i: (cent.x[i], cent.y[i]) for i in gc.index}

    seed_cell, members = pick_neighbourhood(gc, g)
    sub = gc.loc[members]
    subgraph = g.subgraph(members)
    compo = sub[CT_COL].astype(str).value_counts()

    fig, axes = plt.subplots(1, 3, figsize=(26, 9))
    gc.plot(column=CT_COL, categorical=True, legend=True, cmap=cmap, linewidth=0.1,
            edgecolor="none", ax=axes[0],
            legend_kwds={"bbox_to_anchor": (1.02, 0.5), "loc": "center left", "fontsize": 8})
    # whole-core graph: faint so the cell-type map stays readable underneath (Leon, 2026-09-07)
    for u, v in g.edges():
        axes[0].plot(*zip(pos[u], pos[v]), color="white", lw=0.25, alpha=0.45, zorder=3)
    # mark where the inset neighbourhood sits, as the published panel does
    xmin, ymin, xmax, ymax = sub.total_bounds
    mx, my = 0.35 * (xmax - xmin) + 8, 0.35 * (ymax - ymin) + 8
    axes[0].add_patch(plt.Rectangle((xmin - mx, ymin - my),
                                    (xmax - xmin) + 2 * mx, (ymax - ymin) + 2 * my,
                                    fill=False, edgecolor="black", lw=1.4, zorder=9))
    axes[0].set_title(f"Distance-bounded Voronoi Graph\n{name}", fontsize=13)

    sub.plot(column=CT_COL, categorical=True, legend=False, cmap=cmap, linewidth=0.5,
             edgecolor="black", ax=axes[1])
    draw_net(axes[1], subgraph, pos, seed_cell)
    # only the cell types actually in the neighbourhood, in descending abundance
    axes[1].legend(handles=[plt.Line2D([], [], marker="s", linestyle="none", markersize=9,
                                       color=pal[t], label=f"{t} ({n})")
                            for t, n in compo.items()],
                   loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9, frameon=False)
    axes[1].set_title(f"2-hop neighbourhood ({len(members)} cells)", fontsize=13)

    sub.plot(column=FACTOR, cmap="viridis", legend=True, linewidth=0.5,
             edgecolor="black", ax=axes[2],
             legend_kwds=({"shrink": 0.55, "pad": 0.01, "aspect": 30} if factor_cbar_kwds is None
                          else factor_cbar_kwds))
    draw_net(axes[2], subgraph, pos, seed_cell)
    axes[2].set_title(FACTOR, fontsize=13)

    for ax in axes:
        ax.set_aspect("equal"); ax.axis("off")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(outdir / f"{stem}{stem_suffix}.{ext}", dpi=300,
                    bbox_inches="tight")
    plt.close(fig)
    print(f"[fig3g] core '{name}' ({len(gc)} cells, {g.number_of_edges()} edges); "
          f"inset root {str(seed_cell).replace(core, name)} ({gc.loc[seed_cell, CT_COL]}) -> {len(members)} cells, "
          f"{compo.to_dict()}", flush=True)
    return {"core": core, "n_cells": int(len(gc)), "n_edges": int(g.number_of_edges()),
            "inset_seed": str(seed_cell), "inset_size": len(members),
            "inset_root_type": str(gc.loc[seed_cell, CT_COL]),
            "inset_composition": compo.to_dict(),
            "inset_selection_rule": "aPT root, 12-25 cells in the 2-hop neighbourhood, "
                                    "most distinct cell types, ties on aPT fraction then cell id"}


# ------------------------------------------------------- ribbon (3h, ed3f)
def ribbon(compo, levels, root_col, factor, title, ylab, outdir, stem, cmap_name="tab20",
           crop=None, xlabel=None, crop_note=None):
    """Stacked spline ribbon of neighbourhood composition vs the factor.

    `crop`/`crop_note` let another score (27_ed3f_nfkb.py) bring its own crop; left None,
    the Factor-1 panels use RIBBON_CROP exactly as before.
    """
    prop = compo[levels].div(compo[levels].sum(axis=1), axis=0)
    assert np.allclose(prop.sum(axis=1), 1)
    df = prop.join(compo[[factor]]).dropna(subset=[factor]).sort_values(factor)

    q_lo, q_hi = df[factor].quantile(0.01), df[factor].quantile(0.99)
    if crop is None:
        crop = RIBBON_CROP if RIBBON_CROP is not None else (q_lo, q_hi)
        crop_note = ("largest contiguous run of 0.25-unit bins with >= 500 "
                     "neighbourhoods each (scripts/13_crop_justification.py)"
                     if RIBBON_CROP is not None else "1st-99th percentile")
    lo, hi = crop
    cov = float(((df[factor] >= lo) & (df[factor] <= hi)).mean())
    print(f"[{stem}] crop [{lo:.2f}, {hi:.2f}] covers {cov*100:.2f}% of neighbourhoods "
          f"(1st-99th pct would be [{q_lo:.2f}, {q_hi:.2f}]); rule: {crop_note}", flush=True)
    x_seq = np.linspace(lo, hi, 100)
    knots = np.linspace(lo, hi, num=10)[1:-1]
    x = df[factor].to_numpy(float)
    inside = (x >= lo) & (x <= hi)

    fitted = {}
    for lv in levels:
        sp = LSQUnivariateSpline(x[inside], df[lv].to_numpy(float)[inside], t=knots, k=3)
        fitted[lv] = np.clip(sp(x_seq), 0, 1)
    fit = pd.DataFrame(fitted, index=x_seq)
    dev = float(np.abs(fit.sum(axis=1) - 1).max())
    fit = fit.div(fit.sum(axis=1), axis=0)          # renormalise so the stack closes at 1
    print(f"[{stem}] max |sum(spline)-1| before renormalising = {dev:.4f}", flush=True)

    pal, _ = palette(levels, cmap_name)
    cum = fit.cumsum(axis=1)
    fig, ax = plt.subplots(figsize=(12, 8))
    lower = np.zeros(len(x_seq))
    for lv in levels:
        ax.fill_between(x_seq, lower, cum[lv].to_numpy(), color=pal[lv], label=str(lv),
                        alpha=0.8, linewidth=0)
        lower = cum[lv].to_numpy()
    ax.set_xlim(lo, hi); ax.set_ylim(0, 1)
    ax.set_xlabel(f"{xlabel or factor} (cropped to [{lo:.1f}, {hi:.1f}])", fontsize=13)
    ax.set_ylabel(ylab, fontsize=13)
    ax.set_title(title, fontsize=14)
    ax.legend(bbox_to_anchor=(1.005, 1), loc="upper left", fontsize=9, title=root_col)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(outdir / f"{stem}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    fit.to_csv(outdir / f"{stem}_fitted.csv")
    return {"crop": [float(lo), float(hi)], "crop_coverage": cov,
            "crop_rule": crop_note,
            "quantile_crop_1_99": [float(q_lo), float(q_hi)],
            "max_abs_sum_dev": dev, "n_cells": int(len(df))}


# ---------------------------------------------------------------- ed3b
def ed3b(compo, levels, outdir, stem="ed3b_2hop_celltype_heatmap", matrix_stem="ed3b_2hop_celltype_matrix"):
    prop = compo[levels].div(compo[levels].sum(axis=1), axis=0)
    mat = prop.join(compo[["root_cell_type"]]).groupby("root_cell_type",
                                                       observed=True).mean()
    mat = mat.reindex(index=levels, columns=levels)
    mat.to_csv(outdir / f"{matrix_stem}.csv")

    masked = mat.copy()
    np.fill_diagonal(masked.values, np.nan)   # self-proportion dominates; masked as published
    fig, ax = plt.subplots(figsize=(10, 8.5))
    # thin grid lines so the masked (white) diagonal stays visible at the corners,
    # where it would otherwise merge with the figure background
    sns.heatmap(masked, cmap="viridis", cbar_kws={"label": "Proportions"}, ax=ax,
                linewidths=0.5, linecolor="lightgrey", square=True)
    ax.set_xlabel(""); ax.set_ylabel("")
    ax.set_title("Cell-type composition of 2-hop neighbourhoods", fontsize=13)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(outdir / f"{stem}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    d = np.diag(mat.to_numpy())
    print(f"[ed3b] {mat.shape[0]}x{mat.shape[1]} matrix; self-proportion (masked diagonal) "
          f"median {np.nanmedian(d):.3f}, range [{np.nanmin(d):.3f}, {np.nanmax(d):.3f}]", flush=True)
    return {"shape": list(mat.shape), "diag_median": float(np.nanmedian(d))}


# ---------------------------------------------------------------- suppl 7
def suppl7(gdf, outdir, cores, egfr, names=None, prefix="suppl7_"):
    levels = list(pd.Categorical(gdf[CT_COL]).categories)
    _, cmap = palette(levels)
    nl = sorted(gdf[NICHE_COL].astype(str).unique(), key=int)
    _, ncmap = palette(nl, "tab10")
    made = []
    for core in cores:
        gc = gdf[gdf["TMA_ID"] == core]
        name = core if names is None else names[core]   # shown in the title, file name and log
        fig, axes = plt.subplots(1, 3, figsize=(30, 10))
        gc.plot(column=CT_COL, categorical=True, legend=True, cmap=cmap, linewidth=0.05,
                edgecolor="none", ax=axes[0],
                legend_kwds={"bbox_to_anchor": (1.02, 0.5), "loc": "center left", "fontsize": 8})
        axes[0].set_title("Cell type", fontsize=14)
        gc.plot(column=FACTOR, cmap="viridis", legend=True, linewidth=0.05,
                edgecolor="none", ax=axes[1],
                legend_kwds={"shrink": 0.55, "pad": 0.01, "aspect": 30})
        axes[1].set_title(FACTOR, fontsize=14)
        gc.assign(**{NICHE_COL: gc[NICHE_COL].astype(str)}).plot(
            column=NICHE_COL, categorical=True, legend=True, cmap=ncmap, linewidth=0.05,
            edgecolor="none", ax=axes[2],
            legend_kwds={"bbox_to_anchor": (1.02, 0.5), "loc": "center left", "fontsize": 9})
        axes[2].set_title("Niches", fontsize=14)
        for ax in axes:
            ax.set_aspect("equal"); ax.axis("off")
        cond = str(gc["condition"].iloc[0]); eg = egfr.get(core, float("nan"))
        fig.suptitle(f"{name}   |   {cond}   |   eGFR {eg}", fontsize=15)
        fig.tight_layout()
        safe = name.replace(" ", "_").replace("/", "-")
        for ext in ("pdf", "png"):
            fig.savefig(outdir / f"{prefix}{safe}.{ext}", dpi=200, bbox_inches="tight")
        plt.close(fig)
        made.append(core)
        print(f"[suppl7] {name} ({len(gc)} cells, {cond}, eGFR {eg})", flush=True)
    return {"cores": made}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels", default="fig3d,fig3g,fig3h,ed3b,ed3f,suppl7")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--fig3g-core", default=None,
                    help="override the Fig 3g core (default: FIG3G_CORE)")
    ap.add_argument("--outdir", default=None, help="override the output directory")
    ap.add_argument("--stem-suffix", default="", help="appended to output file stems")
    args = ap.parse_args()
    tag = "_smoke" if args.smoke else ""
    want = [p.strip() for p in args.panels.split(",")]
    np.random.seed(SEED)
    warnings.filterwarnings("ignore", message=".*geographic CRS.*")
    sns.set_style("white")

    outdir = Path(args.outdir) if args.outdir else FIG / f"xenium{tag}"
    outdir.mkdir(parents=True, exist_ok=True)
    obs = pd.read_parquet(OUT / "xenium_factor_scores.parquet")
    print(f"[load] obs {obs.shape}, {obs.TMA_ID.nunique()} cores", flush=True)
    meta = {}

    if "fig3d" in want:
        meta["fig3d"] = fig3d(obs, outdir)

    need_gdf = {"fig3g", "suppl7"} & set(want)
    if need_gdf:
        gdf = gpd.read_parquet(OUT / "geodataframes" / "cell_boundaries.parquet")
        gdf.index = gdf["sample_cell_id"]
        gdf[FACTOR] = obs.set_index("sample_cell_id")[FACTOR].reindex(gdf.index)
        if "fig3g" in want:
            import pickle
            with open(OUT / "neighborhoods" / "voronoi_graphs_10.pkl", "rb") as f:
                graphs = pickle.load(f)
            meta["fig3g"] = fig3g(gdf, graphs, outdir, core=args.fig3g_core,
                                  stem_suffix=args.stem_suffix)
        if "suppl7" in want:
            eg, _, _, _ = corrected_core_egfr(obs)        # NOT the corrupted EGFR column
            missing = [c for c in SUPPL7_CORES if c not in set(gdf["TMA_ID"])]
            if missing:
                raise SystemExit(f"Suppl 7 cores not in the object: {missing}")
            meta["suppl7"] = suppl7(gdf, outdir, SUPPL7_CORES, eg)

    if {"fig3h", "ed3b"} & set(want):
        compo = pd.read_csv(OUT / "neighborhoods" / "nh_2hop_compo_10.csv", index_col=0)
        levels = [c for c in compo.columns
                  if c not in ("TMA_ID", "root_cell_type", "disease")]
        compo[FACTOR] = obs.set_index("sample_cell_id")[FACTOR].reindex(compo.index)
        if "ed3b" in want:
            meta["ed3b"] = ed3b(compo, levels, outdir)
        if "fig3h" in want:
            meta["fig3h"] = ribbon(compo, levels, "Cell type", FACTOR,
                                   "Cell-type composition of 2-hop neighbourhoods vs Factor 1",
                                   "Proportions in 2-hop neighbourhood", outdir,
                                   "fig3h_celltype_vs_factor1")

    if "ed3f" in want:
        nz = pd.read_csv(OUT / "neighborhoods" / "nh_2hop_niche_compo_10.csv", index_col=0)
        nlv = [c for c in nz.columns if c not in ("TMA_ID", "root_niche", "disease")]
        nz[FACTOR] = obs.set_index("sample_cell_id")[FACTOR].reindex(nz.index)
        meta["ed3f"] = ribbon(nz, nlv, "Niche", FACTOR,
                              "Niche composition of 2-hop neighbourhoods vs Factor 1",
                              "Proportions in 2-hop neighbourhood", outdir,
                              "ed3f_niche_vs_factor1", cmap_name="tab10")

    (outdir / f"05_figures_meta{tag}.json").write_text(json.dumps({
        "factor": FACTOR, "cell_type_col": CT_COL, "niche_col": NICHE_COL,
        "panels": meta, "seed": SEED,
        "versions": {"geopandas": gpd.__version__, "matplotlib": matplotlib.__version__,
                     "seaborn": sns.__version__, "numpy": np.__version__,
                     "pandas": pd.__version__, "statsmodels": sm.__version__},
    }, indent=2, default=str))
    print(f"[done] {outdir}", flush=True)


if __name__ == "__main__":
    main()

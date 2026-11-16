#!/usr/bin/env python
"""NELL1 and APP: expression across cell types and association with eGFR, snRNA + Xenium.

Leon (2026-09-10): "show me the expression levels of NELL1 and APP across all cell types and
their correlation with eGFR? Single Nuc and Xenium (if NELL1 is included in Xenium)".

NELL1 is NOT on the 468-gene Xenium panel (checked against var_names, not assumed; NELL2 is
not either), so Xenium covers APP only.

Data
  snRNA   adata_20260717_with_clinical.h5ad (363,656 nuclei x 33,000; ambient-corrected
          non-integer counts in .X). Cell types = level2 (16). Per-sample eGFR, patient_ID and
          disease come from the corrected per-sample table
          (vivien_takeover/outputs/metadata/sample_metadata_corrected.csv, 2026-09-10), which
          fills K68 / K95 / K98 where the object has nothing.
  Xenium  xenium_labeled.h5ad (258,375 cells x 468; raw counts in layers["counts"]).
          Cell types = level2_refined (the same 16 names). Per-core eGFR from
          vivien_takeover 05_xenium_figures.corrected_core_egfr(): the object's EGFR column is
          corrupted (one group value broadcast across donors; vivien_takeover/PLAN.md).

Expression across cell types: per cell log1p(counts / cell total x 1e4); per cell type the mean
of that and the % of cells with a nonzero count.

Association with eGFR -- the unit is the DONOR, never the cell: all cells of a donor share one
eGFR, so a cell-level correlation would count thousands of nuclei as independent replicates.
Per donor x cell type: pseudobulk = log1p(sum of gene counts / sum of cell totals x 1e4),
kept when the donor has >= MIN_CELLS cells of that type; plus an "all cells" pseudobulk.
  snRNA  donor = patient_ID (61 donors; the donors of K94/K97 and K95/K98 give 2 samples each, with
         the same eGFR) -- a donor's samples are pooled.
  Xenium several TMA punches can come from one donor. Key = T-number in TMA_ID (55 cores, one
         core each) -> else H.number (31 cores) -> else (condition, eGFR) for the 33 Appollo
         controls, which carry neither (an approximation, as in fig3d). Core-level results are
         written alongside.
Statistics per gene x cell type: Spearman rho with 95 % CI (Fisher z, Bonett-Wright SE); p and
BH-FDR across the cell types within gene x modality x unit; OLS slope per 10 ml/min/1.73 m2
with 95 % CI and Breusch-Pagan / Shapiro-Wilk p on its residuals; partial Spearman adjusted
for disease (ranks residualised on disease dummies) -- eGFR differs systematically between
diseases, so a raw correlation can be a disease contrast in disguise.
Human GRCh38, HGNC symbols. Seed 42.
"""
from kbpaths import KB_ANALYSIS, KB_DATA  # repository path config (config/kbpaths.py)
import argparse, importlib.util, json, sys
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as ss
import seaborn as sns
import statsmodels.api as sm
from statsmodels.stats.diagnostic import het_breuschpagan
from statsmodels.stats.multitest import multipletests

HERE = Path(f"{KB_ANALYSIS}/nell1_app_egfr")
VT = Path(f"{KB_ANALYSIS}/vivien_takeover")
SNRNA = f"{KB_DATA}/charlotte_new_data/adata_20260717_with_clinical.h5ad"
XENIUM = f"{KB_ANALYSIS}/celltype_concordance/outputs/xenium_labeled.h5ad"
META = VT / "outputs" / "metadata" / "sample_metadata_corrected.csv"
GENES = ["NELL1", "APP"]
SEED, MIN_CELLS, MIN_UNITS, SMOKE_PER_GROUP = 42, 20, 8, 40
ORDER = ["PT", "aPT", "dPT", "TAL", "DCT", "CNT", "PC", "IC-A", "IC-B",
         "POD", "PEC", "EC", "MES", "VSM-P", "FIB", "Immune"]
ALL = "all cells"
DCOL = {"CTRL": "#8c564b", "control": "#8c564b", "control_AKI": "#e377c2", "DN": "#9467bd",
        "FSGS": "#1f77b4", "IgAN": "#ff7f0e", "IgA": "#ff7f0e", "MN": "#2ca02c", "TIN": "#d62728"}

_spec = importlib.util.spec_from_file_location("xenium_figures",
                                               Path(__file__).resolve().parents[1] / "03_xenium_factor1_neighbourhoods_cci" / "05_xenium_figures.py")
xf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(xf)


def stratify(obs, keys, k):
    return (obs.groupby(keys, observed=True)
               .apply(lambda d: d.sample(min(len(d), k), random_state=SEED), include_groups=False)
               .index.get_level_values(-1))


def extract(X, var_names):
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
    lib = np.asarray(X.sum(axis=1)).ravel()
    genes = [g for g in GENES if g in var_names]
    G = X[:, [var_names.get_loc(g) for g in genes]].toarray()
    neg = int((X.data < 0).sum())
    return G, lib, genes, neg


# ------------------------------------------------------------------- loaders
def load_snrna(smoke):
    a = ad.read_h5ad(SNRNA)
    print(f"[snRNA] {a.shape}", flush=True)
    if smoke:
        a = a[stratify(a.obs, ["sample", "level2"], SMOKE_PER_GROUP)].copy()
        print(f"[snRNA][smoke] stratified to {a.shape}", flush=True)
    G, lib, genes, neg = extract(a.X, a.var_names)
    m = pd.read_csv(META).set_index("sample")
    s = a.obs["sample"].astype(str)
    miss = sorted(set(s) - set(m.index))
    if miss:
        raise SystemExit(f"[snRNA] samples absent from {META}: {miss}")
    obs = pd.DataFrame({"cell_type": a.obs["level2"].astype(str).to_numpy(), "sample": s.to_numpy(),
                        "donor": s.map(m["patient_ID"]).astype(str).to_numpy(),
                        "egfr": s.map(m["egfr_at_time_of_biopsy"]).astype(float).to_numpy(),
                        "disease": s.map(m["disease"]).astype(str).to_numpy()}, index=a.obs_names)
    info = {"n_cells": int(a.n_obs), "genes": genes, "negative_values": neg,
            "cells_zero_total": int((lib <= 0).sum()),
            "cells_without_egfr_pct": float(np.isnan(obs["egfr"]).mean() * 100)}
    print(f"[snRNA] genes {genes}; {info}", flush=True)
    return obs, G, lib, genes, info


def load_xenium(smoke):
    a = ad.read_h5ad(XENIUM)
    print(f"[Xenium] {a.shape}", flush=True)
    if smoke:
        a = a[stratify(a.obs, ["TMA_ID", "level2_refined"], SMOKE_PER_GROUP)].copy()
        print(f"[Xenium][smoke] stratified to {a.shape}", flush=True)
    G, lib, genes, neg = extract(a.layers["counts"], a.var_names)
    o = a.obs
    egfr_core, _, n_auth, n_changed = xf.corrected_core_egfr(o[["TMA_ID", "EGFR"]])
    egfr_core.index = egfr_core.index.astype(str)
    core = o["TMA_ID"].astype(str)
    obs = pd.DataFrame({"cell_type": o["level2_refined"].astype(str).to_numpy(),
                        "core": core.to_numpy(), "disease": o["condition"].astype(str).to_numpy()},
                       index=a.obs_names)
    obs["egfr"] = obs["core"].map(egfr_core).astype(float)
    tnum = obs["core"].str.extract(r"(T\d{2}-\d{4,5})")[0]
    hnum = pd.Series(o["H.number"].astype(str).to_numpy(), index=obs.index).replace(
        {"nan": np.nan, "None": np.nan, "": np.nan})
    fallback = obs["disease"] + "_eGFR" + obs["egfr"].astype(str)
    obs["donor"] = tnum.fillna(hnum).fillna(fallback)
    src = np.where(tnum.notna(), "T-number", np.where(hnum.notna(), "H.number", "condition+eGFR"))
    per_core = pd.DataFrame({"core": obs["core"], "src": src}).drop_duplicates("core")
    info = {"n_cells": int(a.n_obs), "genes": genes, "not_on_panel": [g for g in GENES if g not in genes],
            "negative_values": neg, "cells_zero_total": int((lib <= 0).sum()),
            "cores": int(obs["core"].nunique()), "donors": int(obs["donor"].nunique()),
            "donor_key_source_cores": per_core["src"].value_counts().to_dict(),
            "cores_without_egfr": int(obs.groupby("core")["egfr"].first().isna().sum()),
            "egfr_cores_from_metadata": int(n_auth), "egfr_values_corrected": int(n_changed)}
    print(f"[Xenium] genes {genes}; {info}", flush=True)
    return obs, G, lib, genes, info


# ------------------------------------------------------------------- summaries
def celltype_summary(obs, G, lib, genes, modality):
    ok = lib > 0
    d = pd.DataFrame({"cell_type": obs["cell_type"].to_numpy()})
    for i, g in enumerate(genes):
        d[f"{g}|mean"] = np.where(ok, np.log1p(G[:, i] / np.where(ok, lib, 1.0) * 1e4), np.nan)
        d[f"{g}|pct"] = (G[:, i] > 0) * 100.0
    agg = d.groupby("cell_type").agg(["mean"]).droplevel(1, axis=1)
    n = d.groupby("cell_type").size()
    rows = [{"modality": modality, "gene": g, "cell_type": ct, "n_cells": int(n[ct]),
             "mean_log1p_cp10k": float(agg.at[ct, f"{g}|mean"]),
             "pct_expressing": float(agg.at[ct, f"{g}|pct"])}
            for g in genes for ct in agg.index]
    return pd.DataFrame(rows)


def pseudobulk(obs, G, lib, genes, unit, modality):
    d = obs[[unit, "cell_type"]].copy()
    for i, g in enumerate(genes):
        d[g] = G[:, i]
    d["lib"], d["n_cells"] = lib, 1
    per = d.groupby([unit, "cell_type"]).sum(numeric_only=True)
    allc = d.drop(columns="cell_type").groupby(unit).sum(numeric_only=True)
    allc.index = pd.MultiIndex.from_arrays([allc.index, [ALL] * len(allc)], names=[unit, "cell_type"])
    pb = pd.concat([per, allc]).reset_index()
    cov = obs.groupby(unit).agg(egfr=("egfr", "mean"), egfr_n=("egfr", "nunique"),
                                disease=("disease", "first"), disease_n=("disease", "nunique"))
    conflicts = {"egfr": int((cov["egfr_n"] > 1).sum()), "disease": int((cov["disease_n"] > 1).sum())}
    pb = pb.join(cov[["egfr", "disease"]], on=unit)
    n_before = len(pb)
    pb = pb[pb["n_cells"] >= MIN_CELLS]
    long = []
    for g in genes:
        t = pb[[unit, "cell_type", "n_cells", "egfr", "disease"]].copy()
        t["gene"], t["counts"] = g, pb[g].to_numpy()
        t["pseudobulk"] = np.log1p(pb[g] / pb["lib"] * 1e4).to_numpy()
        long.append(t)
    out = pd.concat(long).rename(columns={unit: "unit_id"})
    out.insert(0, "unit", unit); out.insert(0, "modality", modality)
    print(f"[{modality}/{unit}] {pb[unit].nunique()} units; {n_before - len(pb)} of {n_before} "
          f"unit x cell-type groups below {MIN_CELLS} cells dropped; conflicts within unit: "
          f"{conflicts}", flush=True)
    return out, conflicts


def fisher_ci(rho, n):
    se = np.sqrt((1 + rho ** 2 / 2) / (n - 3))
    z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    return float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se))


def association(d):
    d = d.dropna(subset=["egfr", "pseudobulk"])
    n = len(d)
    res = {"n": n, "pct_units_detected": float((d["counts"] > 0).mean() * 100) if n else np.nan}
    if n < MIN_UNITS or d["pseudobulk"].nunique() < 3:
        return res
    rho, p = ss.spearmanr(d["egfr"], d["pseudobulk"])
    lo, hi = fisher_ci(rho, n)
    x = sm.add_constant(d["egfr"].to_numpy(float) / 10)
    m = sm.OLS(d["pseudobulk"].to_numpy(float), x).fit()
    ci = m.conf_int(0.05)[1]
    res.update({"rho": float(rho), "rho_ci_lo": lo, "rho_ci_hi": hi, "p": float(p),
                "slope_per_10_egfr": float(m.params[1]), "slope_ci_lo": float(ci[0]),
                "slope_ci_hi": float(ci[1]),
                "breusch_pagan_p": float(het_breuschpagan(m.resid, m.model.exog)[1]),
                "shapiro_p": float(ss.shapiro(m.resid)[1])})
    dums = pd.get_dummies(d["disease"], drop_first=True, dtype=float)
    k = dums.shape[1]
    if k and n - 2 - k > 2:
        z = sm.add_constant(dums.to_numpy())
        rx = sm.OLS(ss.rankdata(d["egfr"]), z).fit().resid
        ry = sm.OLS(ss.rankdata(d["pseudobulk"]), z).fit().resid
        pr = float(np.corrcoef(rx, ry)[0, 1]); df = n - 2 - k
        tt = pr * np.sqrt(df / max(1 - pr ** 2, 1e-12))
        res.update({"partial_rho_disease": pr, "partial_p": float(2 * ss.t.sf(abs(tt), df)),
                    "n_disease_groups": int(k + 1)})
    return res


def associations(pb):
    rows = []
    for (mod, unit, g, ct), d in pb.groupby(["modality", "unit", "gene", "cell_type"]):
        rows.append({"modality": mod, "unit": unit, "gene": g, "cell_type": ct, **association(d)})
    res = pd.DataFrame(rows)
    res["q"] = np.nan
    for _, idx in res.groupby(["modality", "unit", "gene"]).groups.items():
        sub = res.loc[idx, "p"].dropna()
        if len(sub):
            res.loc[sub.index, "q"] = multipletests(sub, method="fdr_bh")[1]
    return res


# ------------------------------------------------------------------- figures
def rows_present(labels):
    extra = sorted(set(labels) - set(ORDER) - {ALL})
    return [c for c in ORDER if c in set(labels)] + extra


def dotplot(summ, out, tag, xinfo):
    cols = [(m, g) for m, g in [("snRNA", "NELL1"), ("snRNA", "APP"), ("Xenium", "APP")]
            if ((summ.modality == m) & (summ.gene == g)).any()]
    rows = rows_present(summ["cell_type"])
    y = {c: len(rows) - 1 - i for i, c in enumerate(rows)}
    fig, axes = plt.subplots(1, len(cols), figsize=(2.3 * len(cols) + 2.2, 0.42 * len(rows) + 2.6),
                             sharey=True)
    axes = np.atleast_1d(axes)
    size = lambda pct: 8 + np.asarray(pct) * 3.2
    for ax, (m, g) in zip(axes, cols):
        s = summ[(summ.modality == m) & (summ.gene == g)].set_index("cell_type").reindex(rows)
        norm = plt.Normalize(0, float(np.nanmax(s["mean_log1p_cp10k"])))
        pts = ax.scatter(np.zeros(len(s)), [y[c] for c in s.index], s=size(s["pct_expressing"].fillna(0)),
                         c=s["mean_log1p_cp10k"], cmap="Reds", norm=norm, edgecolors="grey",
                         linewidths=0.4)
        for c, r in s.iterrows():
            ax.text(0.34, y[c], f"{r['pct_expressing']:.0f}%", va="center", fontsize=7, color="dimgrey")
        ax.set_xlim(-0.5, 0.75); ax.set_xticks([])
        ax.set_title(f"{m}\n{g}", fontsize=12)
        cb = fig.colorbar(pts, ax=ax, orientation="horizontal", pad=0.03, fraction=0.045, aspect=12)
        cb.set_label("mean log1p(CP10k)", fontsize=8); cb.ax.tick_params(labelsize=7)
        sns.despine(ax=ax, left=True, bottom=True)
    axes[0].set_yticks(range(len(rows))); axes[0].set_yticklabels(rows[::-1])
    handles = [plt.scatter([], [], s=size(p), color="lightgrey", edgecolors="grey", label=f"{p}%")
               for p in (10, 25, 50, 100)]
    axes[-1].legend(handles=handles, title="% cells\nexpressing", frameon=False, fontsize=8,
                    title_fontsize=8, loc="upper left", bbox_to_anchor=(1.02, 1.0), labelspacing=1.3)
    note = (f"NELL1 is not on the 468-gene Xenium panel" if "NELL1" in xinfo.get("not_on_panel", [])
            else "")
    fig.suptitle("NELL1 and APP expression per cell type", fontsize=13, y=0.995)
    if note:
        fig.text(0.5, 0.005, note, ha="center", fontsize=8, color="dimgrey")
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))
    for e in ("pdf", "png"):
        fig.savefig(out / f"nell1_app_expression_by_celltype{tag}.{e}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def forest(res, out, tag):
    panels = [(m, "donor", g) for m, g in [("snRNA", "NELL1"), ("snRNA", "APP"), ("Xenium", "APP")]
              if ((res.modality == m) & (res.gene == g)).any()]
    rows = rows_present(res["cell_type"]) + [ALL]
    yy = np.arange(len(rows))[::-1]
    fig, axes = plt.subplots(1, len(panels), figsize=(4.6 * len(panels), 0.36 * len(rows) + 2.2),
                             sharey=True)
    axes = np.atleast_1d(axes)
    for ax, (m, u, g) in zip(axes, panels):
        r = res[(res.modality == m) & (res.unit == u) & (res.gene == g)].set_index("cell_type").reindex(rows)
        for yi, (ct, row) in zip(yy, r.iterrows()):
            if pd.isna(row.get("rho")):
                ax.text(0, yi, "too few donors / not detected", ha="center", va="center",
                        fontsize=7, color="grey")
                continue
            sig = row["q"] < 0.05
            ax.plot([row["rho_ci_lo"], row["rho_ci_hi"]], [yi, yi], color="black", lw=1.1)
            ax.scatter(row["rho"], yi, s=34, facecolor="black" if sig else "white",
                       edgecolor="black", zorder=3)
            if not pd.isna(row.get("partial_rho_disease")):
                ax.scatter(row["partial_rho_disease"], yi, marker="D", s=16, color="#d62728", zorder=4)
            ax.text(1.02, yi, f"n={int(row['n'])}", transform=ax.get_yaxis_transform(),
                    fontsize=7, va="center")
        ax.axvline(0, color="grey", ls="--", lw=0.8)
        ax.axhline(yy[-1] + 0.5, color="lightgrey", lw=0.8)
        ax.set_xlim(-1, 1)
        ax.set_title(f"{m}: {g}", fontsize=12)
        ax.set_xlabel("Spearman ρ with eGFR (per donor)", fontsize=10)
        sns.despine(ax=ax)
    axes[0].set_yticks(yy); axes[0].set_yticklabels(rows)
    handles = [plt.scatter([], [], s=34, facecolor="black", edgecolor="black", label="BH q < 0.05"),
               plt.scatter([], [], s=34, facecolor="white", edgecolor="black", label="q ≥ 0.05"),
               plt.scatter([], [], marker="D", s=16, color="#d62728", label="partial ρ | disease")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=9,
               bbox_to_anchor=(0.5, -0.01))
    fig.suptitle("Pseudobulk expression vs eGFR per cell type (95 % CI)", fontsize=13)
    fig.tight_layout(rect=(0, 0.04, 0.97, 0.96))
    for e in ("pdf", "png"):
        fig.savefig(out / f"nell1_app_egfr_forest{tag}.{e}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def scatter_grid(pb, res, m, u, g, out, tag):
    rows = rows_present(pb.loc[pb.modality == m, "cell_type"]) + [ALL]
    nc = 5; nr = int(np.ceil(len(rows) / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(4.1 * nc, 3.4 * nr))
    axes = axes.ravel()
    for ax, ct in zip(axes, rows):
        d = pb[(pb.modality == m) & (pb.unit == u) & (pb.gene == g) & (pb.cell_type == ct)]
        d = d.dropna(subset=["egfr"])
        for dis, dd in d.groupby("disease"):
            ax.scatter(dd["egfr"], dd["pseudobulk"], s=16, color=DCOL.get(dis, "grey"),
                       edgecolors="none", alpha=0.85, label=dis)
        r = res[(res.modality == m) & (res.unit == u) & (res.gene == g) & (res.cell_type == ct)]
        if len(r) and not pd.isna(r.iloc[0].get("rho")):
            r = r.iloc[0]
            xs = np.linspace(d["egfr"].min(), d["egfr"].max(), 50)
            fit = sm.OLS(d["pseudobulk"].to_numpy(float), sm.add_constant(d["egfr"].to_numpy(float))).fit()
            pr = fit.get_prediction(sm.add_constant(xs)).summary_frame(alpha=0.05)
            ax.plot(xs, pr["mean"], color="black", lw=1)
            ax.fill_between(xs, pr["mean_ci_lower"], pr["mean_ci_upper"], color="grey", alpha=0.25, lw=0)
            ax.set_title(f"{ct}\nρ = {r['rho']:+.2f} [{r['rho_ci_lo']:+.2f}, {r['rho_ci_hi']:+.2f}], "
                         f"q = {r['q']:.2g}, n = {int(r['n'])}", fontsize=9)
        else:
            ax.set_title(f"{ct}\n(n = {len(d)}: too few donors)", fontsize=9)
        ax.tick_params(labelsize=7)
        sns.despine(ax=ax)
    for ax in axes[len(rows):]:
        ax.axis("off")
    h, l = [], []
    for ax in axes[:len(rows)]:
        for hh, ll in zip(*ax.get_legend_handles_labels()):
            if ll not in l:
                h.append(hh); l.append(ll)
    axes[-1].legend(h, l, title="disease", frameon=False, fontsize=9, loc="center")
    fig.supxlabel("eGFR at biopsy (ml/min/1.73 m²)", fontsize=12)
    fig.supylabel(f"{g} pseudobulk, log1p(CP10k)", fontsize=12)
    fig.suptitle(f"{m}: {g} vs eGFR, one point per {u} (>= {MIN_CELLS} cells of that type); "
                 f"OLS fit with 95 % CI", fontsize=13)
    fig.tight_layout(rect=(0.01, 0.01, 1, 0.97))
    for e in ("pdf", "png"):
        fig.savefig(out / f"{g}_{m}_{u}_egfr_scatter{tag}.{e}", dpi=250, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    tag = "_smoke" if args.smoke else ""
    np.random.seed(SEED)
    sns.set_style("white")
    fig_dir, tab_dir = HERE / "figures", HERE / "tables"
    fig_dir.mkdir(parents=True, exist_ok=True); tab_dir.mkdir(parents=True, exist_ok=True)

    so, sG, slib, sgenes, sinfo = load_snrna(args.smoke)
    xo, xG, xlib, xgenes, xinfo = load_xenium(args.smoke)

    summ = pd.concat([celltype_summary(so, sG, slib, sgenes, "snRNA"),
                      celltype_summary(xo, xG, xlib, xgenes, "Xenium")])
    summ.to_csv(tab_dir / f"expression_by_celltype{tag}.csv", index=False)
    print(f"\n[expression by cell type]\n"
          f"{summ.pivot_table(index='cell_type', columns=['modality','gene'], values=['mean_log1p_cp10k','pct_expressing']).round(3).to_string()}",
          flush=True)

    pbs, conflicts = [], {}
    for obs, G, lib, genes, mod, units in [(so, sG, slib, sgenes, "snRNA", ["donor"]),
                                           (xo, xG, xlib, xgenes, "Xenium", ["donor", "core"])]:
        for u in units:
            p, c = pseudobulk(obs, G, lib, genes, u, mod)
            pbs.append(p); conflicts[f"{mod}/{u}"] = c
    pb = pd.concat(pbs)
    pb.to_csv(tab_dir / f"pseudobulk_per_unit{tag}.csv", index=False)
    res = associations(pb)
    res.to_csv(tab_dir / f"egfr_association{tag}.csv", index=False)
    show = res[["modality", "unit", "gene", "cell_type", "n", "rho", "rho_ci_lo", "rho_ci_hi", "q",
                "partial_rho_disease", "partial_p", "slope_per_10_egfr", "pct_units_detected"]]
    with pd.option_context("display.width", 220, "display.max_rows", 200):
        print(f"\n[eGFR association]\n{show.round(3).to_string(index=False)}", flush=True)

    dotplot(summ, fig_dir, tag, xinfo)
    forest(res, fig_dir, tag)
    for m, u, g in [("snRNA", "donor", g) for g in sgenes] + [("Xenium", "donor", g) for g in xgenes]:
        scatter_grid(pb, res, m, u, g, fig_dir, tag)

    (tab_dir / f"01_meta{tag}.json").write_text(json.dumps({
        "snRNA": sinfo, "Xenium": xinfo, "metadata": str(META), "min_cells": MIN_CELLS,
        "min_units": MIN_UNITS, "unit_conflicts": conflicts, "smoke": bool(args.smoke),
        "seed": SEED, "species": "human GRCh38", "gene_id_format": "HGNC symbols",
        "versions": {"anndata": ad.__version__, "numpy": np.__version__, "pandas": pd.__version__,
                     "scipy": sp.__name__ and __import__("scipy").__version__,
                     "statsmodels": sm.__version__, "matplotlib": matplotlib.__version__,
                     "python": sys.version.split()[0]}}, indent=2, default=str))
    print(f"[done] {fig_dir}", flush=True)


if __name__ == "__main__":
    main()

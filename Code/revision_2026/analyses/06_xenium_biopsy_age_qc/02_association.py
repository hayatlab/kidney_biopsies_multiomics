#!/usr/bin/env python
"""Associations between biopsy age and per-core Xenium quality metrics. See PLAN.md.
Spearman + bootstrap CI (seed 42); OLS with condition (or TMA block) covariate, HC3 SE; BH across metrics.
Figures: year distribution, scatter grid, forest plot, residual QQ."""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
import statsmodels.api as sm
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests

D = Path(f"{KB_ANALYSIS}/xenium_biopsy_age_qc")
SEED, NBOOT = 42, 2000
rng = np.random.default_rng(SEED)
core = pd.read_csv(D / "outputs/core_qc_metrics.csv")
core["age10"] = core.biopsy_age_years / 10
COND_ORDER = ["CTRL", "FSGS", "DN", "IgAN", "MN", "TIN"]
PAL = dict(zip(COND_ORDER, ["#000000", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7"]))  # Okabe-Ito

# metric -> (label, transform) ; "log" = log scale (effect reported as % per 10 y), "raw" = linear (abs change per 10 y)
METRICS = {
    "n_cells": ("cells per core", "log"),
    "cell_density_per_mm2": ("cell density (cells/mm² hull)", "log"),
    "median_tx_per_cell": ("median transcripts per cell", "log"),
    "median_genes_per_cell": ("median genes per cell", "log"),
    "tx_per_100um2": ("transcripts per 100 µm² cell area", "log"),
    "genes_at_100": ("complexity: genes at 100 transcripts", "log"),
    "frac_cells_lt20": ("fraction of cells < 20 transcripts", "raw"),
    "negprobe_per_10k_tx": ("negative-probe counts per 10k transcripts", "raw"),
    "unassigned_per_10k_tx": ("unassigned codewords per 10k transcripts", "raw"),
    "median_cell_area": ("median cell area (µm²)", "log"),
    "median_nucleus_area": ("median nucleus area (µm²)", "log"),
    "median_conf_level2": ("median scANVI level2 confidence", "raw"),
    "frac_conf_lt05": ("fraction of cells with confidence < 0.5", "raw"),
    "PT_median_tx_per_cell": ("PT cells: median transcripts per cell", "log"),
    "PT_median_genes_per_cell": ("PT cells: median genes per cell", "log"),
    "TAL_median_tx_per_cell": ("TAL cells: median transcripts per cell", "log"),
    "TAL_median_genes_per_cell": ("TAL cells: median genes per cell", "log"),
    "Factor1_ulm_mean_all": ("Factor 1 score, all cells (mean)", "raw"),
    "Factor1_ulm_mean_PTview": ("Factor 1 score, PT view (mean)", "raw"),
    "Factor2_ulm_mean_all": ("Factor 2 score, all cells (mean)", "raw"),
    "Factor3_ulm_mean_all": ("Factor 3 score, all cells (mean)", "raw"),
    "Factor4_ulm_mean_all": ("Factor 4 score, all cells (mean)", "raw"),
    "Factor5_ulm_mean_all": ("Factor 5 score, all cells (mean)", "raw"),
}
QC_KEY = ["n_cells", "cell_density_per_mm2", "median_tx_per_cell", "median_genes_per_cell", "tx_per_100um2",
          "genes_at_100", "frac_cells_lt20", "negprobe_per_10k_tx", "median_cell_area", "median_conf_level2",
          "PT_median_tx_per_cell", "Factor1_ulm_mean_PTview"]
SUBSETS = {"all": core, "disease": core[core.condition != "CTRL"]}
COVARS = {"condition": "C(condition)", "tma_block": "C(tma_block)"}

def spearman_boot(x, y):
    ok = x.notna() & y.notna(); x, y = x[ok].values, y[ok].values
    rho = stats.spearmanr(x, y).statistic
    idx = rng.integers(0, len(x), size=(NBOOT, len(x)))
    bs = np.array([stats.spearmanr(x[i], y[i]).statistic for i in idx])
    lo, hi = np.nanpercentile(bs, [2.5, 97.5])
    return rho, lo, hi, int(ok.sum())

res, resid_store = [], {}
for sub_name, df in SUBSETS.items():
    for metric, (label, tf) in METRICS.items():
        d = df[["age10", "biopsy_age_years", "condition", "tma_block", metric]].dropna().copy()
        rho, lo, hi, n = spearman_boot(d.biopsy_age_years, d[metric])
        d["y"] = np.log(d[metric]) if tf == "log" else d[metric]
        row = dict(subset=sub_name, metric=metric, label=label, transform=tf, n=n, spearman_rho=rho, rho_ci_lo=lo, rho_ci_hi=hi)
        for cov_name, cov in COVARS.items():
            fit = smf.ols(f"y ~ age10 + {cov}", data=d).fit(cov_type="HC3")
            b, ci = fit.params["age10"], fit.conf_int().loc["age10"]
            if tf == "log":
                eff, elo, ehi = 100 * (np.exp(b) - 1), 100 * (np.exp(ci[0]) - 1), 100 * (np.exp(ci[1]) - 1)
            else:
                eff, elo, ehi = b, ci[0], ci[1]
            sw = stats.shapiro(fit.resid).pvalue if len(fit.resid) >= 3 else np.nan
            row.update({f"{cov_name}_effect_per10y": eff, f"{cov_name}_ci_lo": elo, f"{cov_name}_ci_hi": ehi,
                        f"{cov_name}_p": fit.pvalues["age10"], f"{cov_name}_r2": fit.rsquared, f"{cov_name}_shapiro_p": sw})
            if cov_name == "condition":
                resid_store[(sub_name, metric)] = fit.resid.values
                Xr = pd.get_dummies(d[["age10", "condition"]], columns=["condition"], drop_first=True, dtype=float)
                rfit = sm.RLM(d["y"], sm.add_constant(Xr), M=sm.robust.norms.HuberT()).fit()
                rb, rci = rfit.params["age10"], rfit.conf_int().loc["age10"]
                row.update(dict(huber_effect_per10y=100 * (np.exp(rb) - 1) if tf == "log" else rb,
                                huber_ci_lo=100 * (np.exp(rci[0]) - 1) if tf == "log" else rci[0],
                                huber_ci_hi=100 * (np.exp(rci[1]) - 1) if tf == "log" else rci[1]))
        # unadjusted model (for transparency)
        fit0 = smf.ols("y ~ age10", data=d).fit(cov_type="HC3")
        b0, ci0 = fit0.params["age10"], fit0.conf_int().loc["age10"]
        row.update(dict(unadj_effect_per10y=100 * (np.exp(b0) - 1) if tf == "log" else b0,
                        unadj_ci_lo=100 * (np.exp(ci0[0]) - 1) if tf == "log" else ci0[0],
                        unadj_ci_hi=100 * (np.exp(ci0[1]) - 1) if tf == "log" else ci0[1], unadj_p=fit0.pvalues["age10"]))
        res.append(row)
R = pd.DataFrame(res)
for sub_name in SUBSETS:
    for cov_name in list(COVARS) + ["unadj"]:
        m = R.subset == sub_name
        R.loc[m, f"{cov_name}_q_bh"] = multipletests(R.loc[m, f"{cov_name}_p"], method="fdr_bh")[1]
R.to_csv(D / "outputs/associations_biopsy_age.csv", index=False)

# age-bin descriptive table (disease cores + controls separately)
core["age_bin"] = pd.cut(core.biopsy_age_years, [0, 5, 10, 20, 30], labels=["<5 y", "5-10 y", "10-20 y", ">20 y"])
core["group"] = np.where(core.condition == "CTRL", "CTRL", "disease")
T = core.groupby(["group", "age_bin"], observed=True).agg(n=("core_id", "size"), **{k: (k, "median") for k in QC_KEY}).round(3)
T.to_csv(D / "outputs/age_bin_medians.csv")
print(T.to_string())

# ---------------- figures ----------------
plt.rcParams.update({"font.size": 7, "axes.spines.top": False, "axes.spines.right": False})
# F1 biopsy year by condition
fig, ax = plt.subplots(figsize=(4.2, 2.4))
for i, c in enumerate(COND_ORDER):
    y = core.loc[core.condition == c, "biopsy_year"]
    ax.scatter(y + rng.uniform(-0.25, 0.25, len(y)), np.full(len(y), i) + rng.uniform(-0.2, 0.2, len(y)), s=9, color=PAL[c], alpha=0.85, lw=0)
ax.set_yticks(range(len(COND_ORDER))); ax.set_yticklabels([f"{c} (n={int((core.condition == c).sum())})" for c in COND_ORDER])
ax.axvline(2020.3, ls="--", lw=0.7, color="0.5"); ax.text(2020.4, -0.6, "5 y before run", fontsize=6, color="0.4")
if (core.loc[core.condition == "CTRL", "biopsy_year"] == 2000).all():
    ax.set_ylim(-0.6, len(COND_ORDER) - 0.2); ax.text(1999.6, len(COND_ORDER) - 0.45, "controls dated 2000 per collaborator (clinical file: 2021–2024)", fontsize=5.5, color="0.4", ha="left", va="top")
ax.set_xlabel("biopsy year"); ax.set_title("Biopsy year of the 119 Xenium cores (run 2025-04-16)", fontsize=7)
fig.tight_layout(); fig.savefig(D / "figures/F1_biopsy_year_by_condition.pdf"); fig.savefig(D / "figures/F1_biopsy_year_by_condition.png", dpi=200); plt.close(fig)

# F2 scatter grid
fig, axes = plt.subplots(3, 4, figsize=(9, 6.2))
for ax, metric in zip(axes.ravel(), QC_KEY):
    label, tf = METRICS[metric]
    for c in COND_ORDER:
        s = core[core.condition == c]
        ax.scatter(s.biopsy_age_years, s[metric], s=8, color=PAL[c], alpha=0.8, lw=0, label=c)
    dz = core[(core.condition != "CTRL")].dropna(subset=[metric])
    xx = np.linspace(dz.biopsy_age_years.min(), dz.biopsy_age_years.max(), 50)
    if tf == "log":
        b = np.polyfit(dz.biopsy_age_years, np.log(dz[metric]), 1); ax.plot(xx, np.exp(np.polyval(b, xx)), color="0.3", lw=0.8); ax.set_yscale("log")
    else:
        b = np.polyfit(dz.biopsy_age_years, dz[metric], 1); ax.plot(xx, np.polyval(b, xx), color="0.3", lw=0.8)
    r = R[(R.subset == "disease") & (R.metric == metric)].iloc[0]
    ax.set_title(f"{label}\ndisease-only ρ={r.spearman_rho:.2f} [{r.rho_ci_lo:.2f}, {r.rho_ci_hi:.2f}]", fontsize=6)
    ax.set_xlabel("biopsy age at run (years)", fontsize=6); ax.tick_params(labelsize=6)
axes[0, 0].legend(fontsize=5, frameon=False, loc="lower left")
fig.suptitle("Per-core Xenium quality vs biopsy age (line: unadjusted fit on disease cores)", fontsize=8)
fig.tight_layout(); fig.savefig(D / "figures/F2_metrics_vs_biopsy_age.pdf"); fig.savefig(D / "figures/F2_metrics_vs_biopsy_age.png", dpi=200); plt.close(fig)

# F3 forest: condition-adjusted effect per 10 y, log-metrics (% change) and raw-metrics (absolute) in two panels
fig, axes = plt.subplots(1, 2, figsize=(8.5, 4.6), gridspec_kw={"width_ratios": [1, 1]})
for ax, tf, xlab in [(axes[0], "log", "% change per 10 years of biopsy age (95% CI)"), (axes[1], "raw", "absolute change per 10 years (95% CI), metric units")]:
    sub = R[(R["transform"] == tf) & (R.subset == "disease")].reset_index(drop=True)
    suba = R[(R["transform"] == tf) & (R.subset == "all")].set_index("metric").loc[sub.metric].reset_index()
    y = np.arange(len(sub))
    ax.errorbar(sub.condition_effect_per10y, y + 0.2, xerr=[sub.condition_effect_per10y - sub.condition_ci_lo, sub.condition_ci_hi - sub.condition_effect_per10y], fmt="o", ms=3, color="#D55E00", lw=0.8, label="disease cores, + condition")
    ax.errorbar(sub.huber_effect_per10y, y, xerr=[sub.huber_effect_per10y - sub.huber_ci_lo, sub.huber_ci_hi - sub.huber_effect_per10y], fmt="^", ms=3, color="#009E73", lw=0.8, label="disease cores, + condition, Huber robust")
    ax.errorbar(suba.condition_effect_per10y, y - 0.2, xerr=[suba.condition_effect_per10y - suba.condition_ci_lo, suba.condition_ci_hi - suba.condition_effect_per10y], fmt="s", ms=3, color="#0072B2", lw=0.8, label="all cores, + condition")
    ax.axvline(0, color="0.5", lw=0.7); ax.set_yticks(y); ax.set_yticklabels(sub.label, fontsize=6); ax.invert_yaxis(); ax.set_xlabel(xlab, fontsize=6)
    ax.tick_params(labelsize=6)
axes[0].legend(fontsize=6, frameon=False, loc="lower right")
fig.suptitle("OLS effect of biopsy age, adjusted for condition (HC3 SE)", fontsize=8)
fig.tight_layout(); fig.savefig(D / "figures/F3_forest_adjusted_effects.pdf"); fig.savefig(D / "figures/F3_forest_adjusted_effects.png", dpi=200); plt.close(fig)

# F4 residual QQ (disease, condition model) for the key metrics
fig, axes = plt.subplots(3, 4, figsize=(8, 6))
for ax, metric in zip(axes.ravel(), QC_KEY):
    r = resid_store[("disease", metric)]
    sm.qqplot(r, line="s", ax=ax, markersize=2, alpha=0.7)
    ax.set_title(f"{METRICS[metric][0]}\nShapiro p={stats.shapiro(r).pvalue:.2g}", fontsize=6); ax.tick_params(labelsize=5); ax.set_xlabel(""); ax.set_ylabel("")
fig.suptitle("Residual QQ plots, disease cores, log/raw metric ~ biopsy age + condition", fontsize=8)
fig.tight_layout(); fig.savefig(D / "figures/F4_residual_qq.pdf"); fig.savefig(D / "figures/F4_residual_qq.png", dpi=200); plt.close(fig)

cols = ["subset", "metric", "n", "spearman_rho", "rho_ci_lo", "rho_ci_hi", "condition_effect_per10y", "condition_ci_lo", "condition_ci_hi", "condition_q_bh", "huber_effect_per10y", "huber_ci_lo", "huber_ci_hi", "tma_block_effect_per10y", "tma_block_ci_lo", "tma_block_ci_hi", "condition_shapiro_p"]
pd.set_option("display.width", 250)
print(R[cols].round(3).to_string())
json.dump(dict(seed=SEED, nboot=NBOOT, n_all=int(len(core)), n_disease=int((core.condition != "CTRL").sum()), metrics=list(METRICS),
               versions={"statsmodels": sm.__version__, "scipy": __import__("scipy").__version__, "pandas": pd.__version__}),
          open(D / "outputs/02_meta.json", "w"), indent=1)

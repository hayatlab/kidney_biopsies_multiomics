#!/usr/bin/env python
"""Two simple figures for the reviewer response / Slack: scatter + box plots of key QC metrics vs biopsy age.
Numbers are read from outputs/associations_biopsy_age.csv (disease-only, condition-adjusted) — nothing recomputed."""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.api as sm

D = Path(f"{KB_ANALYSIS}/xenium_biopsy_age_qc")
core = pd.read_csv(D / "outputs/core_qc_metrics.csv")
R = pd.read_csv(D / "outputs/associations_biopsy_age.csv"); R = R[R.subset == "disease"].set_index("metric")
COND = ["CTRL", "FSGS", "DN", "IgAN", "MN", "TIN"]
PAL = dict(zip(COND, ["#000000", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7"]))
M = [("median_tx_per_cell", "transcripts per cell (median)", "log"),
     ("median_genes_per_cell", "genes per cell (median)", "log"),
     ("tx_per_100um2", "transcripts per 100 µm² cell area", "log"),
     ("genes_at_100", "complexity: genes at 100 transcripts", "log"),
     ("negprobe_per_10k_tx", "negative-probe counts per 10k transcripts", "raw"),
     ("Factor1_ulm_mean_PTview", "Factor 1 score, PT cells (mean)", "raw")]
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})

def effect_text(m, tf):
    r = R.loc[m]
    def e(pre):
        return (f"{r[pre+'effect_per10y']:+.1f}% [{r[pre+'ci_lo']:+.0f}, {r[pre+'ci_hi']:+.0f}] per 10 y" if tf == "log"
                else f"{r[pre+'effect_per10y']:+.2f} [{r[pre+'ci_lo']:+.2f}, {r[pre+'ci_hi']:+.2f}] per 10 y")
    return f"ρ = {r.spearman_rho:.2f} [{r.rho_ci_lo:.2f}, {r.rho_ci_hi:.2f}]\nslope (band): {e('unadj_')}\nadjusted for condition: {e('condition_')}"

# ---- Figure A: scatter (coloured by condition, and a plain version) ----
for coloured in (True, False):
  fig, axes = plt.subplots(2, 3, figsize=(9.5, 5.6))
  for ax, (m, lab, tf) in zip(axes.ravel(), M):
      if coloured:
          for c in COND:
              s = core[core.condition == c]
              ax.scatter(s.biopsy_age_years, s[m], s=14, color=PAL[c], alpha=0.8, lw=0, label=c)
      else:
          ctl = core.condition == "CTRL"
          ax.scatter(core.biopsy_age_years[~ctl], core[m][~ctl], s=14, color="0.35", alpha=0.75, lw=0, label="disease cores")
          ax.scatter(core.biopsy_age_years[ctl], core[m][ctl], s=14, facecolors="none", edgecolors="0.35", lw=0.7, alpha=0.9, label="control cores")
      dz = core[core.condition != "CTRL"].dropna(subset=[m])
      xx = np.linspace(dz.biopsy_age_years.min(), dz.biopsy_age_years.max(), 50)
      yv = np.log(dz[m]) if tf == "log" else dz[m]
      fit = sm.OLS(yv.values, sm.add_constant(dz.biopsy_age_years.values)).fit()
      pr = fit.get_prediction(sm.add_constant(xx)).summary_frame(alpha=0.05)
      back = np.exp if tf == "log" else (lambda v: v)
      ax.fill_between(xx, back(pr["mean_ci_lower"]), back(pr["mean_ci_upper"]), color="0.5", alpha=0.25, lw=0, label="95% CI of fit")
      ax.plot(xx, back(pr["mean"]), color="0.25", lw=1.2)
      if tf == "log": ax.set_yscale("log")
      ax.axvline(5, ls=":", lw=0.8, color="0.6")
      ax.set_title(lab, fontsize=8.5, fontweight="bold")
      ax.text(0.03, 0.03, effect_text(m, tf), transform=ax.transAxes, fontsize=6, va="bottom",
              bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.8", alpha=0.9))
      ax.set_xlabel("biopsy age at Xenium run (years)")
  h, l = axes[0, 0].get_legend_handles_labels()
  fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 0.905), fontsize=6.5, frameon=False, ncol=7)
  fig.suptitle("Xenium quality per TMA core vs biopsy age (119 cores, one run on 2025-04-16)\n"
               "Line and shaded band: linear fit through the disease cores (n = 86) with its 95% confidence band; controls (year 2000 per collaborator) shown but not fitted\n"
               "Box: Spearman ρ, slope of the band, and the slope adjusted for condition, each with 95% CI",
               fontsize=8)
  fig.tight_layout(rect=(0, 0, 1, 0.885))
  tag = "S1_simple_scatter" if coloured else "S1_simple_scatter_nocolour"
  fig.savefig(D / f"figures/{tag}.pdf"); fig.savefig(D / f"figures/{tag}.png", dpi=220); plt.close(fig)

# ---- Figure B: box plots by age bin ----
core["bin"] = np.where(core.condition == "CTRL", "CTRL\n(year 2000)",
               np.where(core.biopsy_age_years < 10, "disease\n< 10 y", np.where(core.biopsy_age_years < 20, "disease\n10–20 y", "disease\n> 20 y")))
order = ["CTRL\n(year 2000)", "disease\n< 10 y", "disease\n10–20 y", "disease\n> 20 y"]
rng = np.random.default_rng(42)
fig, axes = plt.subplots(2, 3, figsize=(9.5, 5.6))
for ax, (m, lab, tf) in zip(axes.ravel(), M):
    data = [core.loc[core.bin == b, m].dropna().values for b in order]
    ax.boxplot(data, positions=range(len(order)), widths=0.55, showfliers=False,
               medianprops=dict(color="black", lw=1.2), boxprops=dict(color="0.4"), whiskerprops=dict(color="0.4"), capprops=dict(color="0.4"))
    for i, b in enumerate(order):
        s = core[core.bin == b]
        ax.scatter(i + rng.uniform(-0.18, 0.18, len(s)), s[m], s=9, color=[PAL[c] for c in s.condition], alpha=0.7, lw=0, zorder=3)
    ax.set_xticks(range(len(order))); ax.set_xticklabels([f"{b}\n(n = {len(d)})" for b, d in zip(order, data)], fontsize=6.5)
    if tf == "log": ax.set_yscale("log")
    ax.set_title(lab, fontsize=8.5, fontweight="bold")
fig.suptitle("Same metrics by biopsy-age bin (box = median and IQR; points = cores coloured by condition as in the scatter figure)", fontsize=8)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(D / "figures/S2_simple_boxplots.pdf"); fig.savefig(D / "figures/S2_simple_boxplots.png", dpi=220); plt.close(fig)
print(core.groupby("bin")[["median_tx_per_cell", "median_genes_per_cell", "tx_per_100um2", "genes_at_100", "negprobe_per_10k_tx", "Factor1_ulm_mean_PTview"]].median().round(2).loc[order])

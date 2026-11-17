#!/usr/bin/env python
"""Supp. Fig. 10d — per-core Xenium quality metrics against biopsy age at the Xenium run: median transcripts and genes per
cell, transcripts per 100 um2 cell area, genes at 100 transcripts, negative-probe counts per 10k transcripts and the
mean Factor 1 score of PT cells; disease cores filled, control cores open; line and band = linear fit through the
disease cores (n = 86) with its 95% confidence band; each box gives Spearman rho, the slope of the band and the slope
adjusted for condition, each with its 95% CI (read from 02_association.py, not recomputed).
Panel        : All Supp Figures.pdf, Supplementary Fig. 10, panel d (all six subpanels)
Reads        : $KB_ANALYSIS/xenium_biopsy_age_qc/outputs/core_qc_metrics.csv
                 (analyses/06_xenium_biopsy_age_qc/01_core_qc_metrics.py)
               $KB_ANALYSIS/xenium_biopsy_age_qc/outputs/associations_biopsy_age.csv, rows subset == "disease"
                 (analyses/06_xenium_biopsy_age_qc/02_association.py)
Writes       : $KB_FIG_OUT/S10/S10d_qc_metrics_vs_biopsy_age.pdf (+ .png)
Plotting code: extracted verbatim from analyses/06_xenium_biopsy_age_qc/04_simple_plots.py, "Figure A" with
               coloured=False (that script has no functions and writes all its figures into the analysis workspace
               when run, so it is not imported) (original run: python 04_simple_plots.py, 2026-10-16 ->
               figures/S1_simple_scatter_nocolour)
Illustrator  : three-line suptitle removed; legend moved above the subpanels
Environment  : conda env agent_env
Verified     : 2026-11-12, PIXEL-IDENTICAL vs $KB_ANALYSIS/xenium_biopsy_age_qc/figures/S1_simple_scatter_nocolour.pdf
               (2026-10-16). The published panel is an EARLIER version of that file: there the 33 control cores sit at biopsy
               ages of 0-5 years (biopsy dates of the metadata workbook, run of 2026-10-15); in the file on disk and in this
               script they sit at 24.8 years (01_core_qc_metrics.py CTRL_YEAR_OVERRIDE = 2000, collaborator statement,
               2026-10-16). Disease cores, fits and all printed statistics are identical in both versions

Biopsy age of the control cores is taken as the year 2000 (collaborator information); controls are shown, not fitted.
"""
from pathlib import Path

import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.api as sm

from kbpaths import KB_ANALYSIS, fig_out  # repository path config (config/kbpaths.py)

# ---- verbatim from 04_simple_plots.py (module level) ----
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


def main():
    # ---- verbatim from 04_simple_plots.py, Figure A, the coloured=False pass ----
    coloured = False
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
    out = Path(fig_out("S10", "S10d_qc_metrics_vs_biopsy_age.pdf"))
    fig.savefig(out); fig.savefig(out.with_suffix(".png"), dpi=220); plt.close(fig)
    print(f"[S10d] wrote {out}", flush=True)


if __name__ == "__main__":
    main()

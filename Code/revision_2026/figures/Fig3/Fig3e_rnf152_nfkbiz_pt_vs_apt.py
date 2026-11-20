#!/usr/bin/env python
"""Fig. 3e — RNF152 and NFKBIZ expression in PT and aPT cells per Xenium core (swarm plots with condition means and a
Welch t-test PT vs aPT).
Panel        : Main Figures.pdf, Figure 3, panel e
Reads        : $KB_RESULTS/xenium/leon/new_data/overall_with_glom_04_level2/xenium_with_niches_overall_with_glom.h5ad
               (analyses/01_xenium_preprocessing_niches; 258,375 cells x 468 genes, layers["counts"], `level2_refined`,
               `condition`, core column `TMA_ID`, or `xenium_core_id` in the de-identified object)
Writes       : $KB_FIG_OUT/Fig3/Fig3e_rnf152_nfkbiz_pt_vs_apt.pdf (page 1 RNF152, page 2 NFKBIZ)
Plotting code: analyses/13_xenium_pt_apt_gene_expression/gene_expression_per_core.py (functions of the notebook `Fig3e.ipynb`
               of Sidrah Maryam, imported; calls as in the notebook). The notebook as received lists the genes PROM1,
               NFKBIZ, RNF212B, GRAMD1B, SLC13A3 and PDZD2; the published panel shows RNF152 and NFKBIZ, so the same
               calls are made for these two genes.
Illustrator  : the two plots stacked; gene name moved into the plot ("Gene: ... (ANOVA p=...)" title removed); legend
               entry CTRL relabelled CTRL_Xenium; fonts
Environment  : conda env nichecompass (scanpy, seaborn, scipy, statsmodels)
Verified     : 2026-11-20, values equal to the printed panel: the 24 condition means (diamonds) measured in the printed
               PDF agree within 0.004; both brackets "****" (Welch t-test on the per-core means, RNF152 p = 3.3e-63,
               NFKBIZ p = 6.6e-23; 117 cores with PT cells, 119 with aPT cells). Not pixel-compared (the notebook figure
               of the two genes was not saved as a file).
"""
import importlib.util
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "config"))
from kbpaths import KB_RESULTS, fig_out  # noqa: E402
import cores  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import scanpy as sc  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

spec = importlib.util.spec_from_file_location(
    "gene_expression_per_core", REPO / "analyses" / "13_xenium_pt_apt_gene_expression" / "gene_expression_per_core.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

OBJ = os.path.join(KB_RESULTS, "xenium", "leon", "new_data", "overall_with_glom_04_level2", "xenium_with_niches_overall_with_glom.h5ad")
GENES = ["RNF152", "NFKBIZ"]

adata = sc.read_h5ad(OBJ)
# as in the notebook
adata = adata[adata.obs["condition"].isin(["CTRL", "TIN", "FSGS", "MN", "DN", "IgAN"])]
cell_types = ["PT", "aPT"]
adata_sub = adata[adata.obs["level2_refined"].isin(cell_types)].copy()

# cores are grouped by their public id; the internal label (pathology number) is translated with the private key
if "xenium_core_id" not in adata_sub.obs:
    adata_sub.obs["xenium_core_id"] = [x for x, _ in cores.public_core_ids(adata_sub.obs["TMA_ID"].astype(str))]

expr_data = m.get_sample_gene_expression(adata_sub, GENES, normalize=True, target_sum=1e4, use_layer="counts",
                                         use_raw=False, core_col="xenium_core_id")
pairs_to_test = [("PT", "aPT")]

# the statistics drawn as brackets, checked against the published values
expected_p = {"RNF152": 3.30e-63, "NFKBIZ": 6.56e-23}
for g in GENES:
    anova_p, pair_df = m.anova_and_pairs(expr_data, g, group_col="level2_refined", value_col="expression",
                                         pairs=pairs_to_test, use_welch=True, fdr=True)
    d = expr_data[expr_data["gene"] == g].dropna(subset=["expression"])
    n = d.groupby("level2_refined", observed=True)["expression"].size()
    print(f"{g}: Welch p = {pair_df['p'].iloc[0]:.3e} ({pair_df['stars'].iloc[0]}), ANOVA p = {anova_p:.3e}, "
          f"cores PT {n['PT']} / aPT {n['aPT']}")
    assert np.isclose(pair_df["p"].iloc[0], expected_p[g], rtol=0.01) and pair_df["stars"].iloc[0] == "****"
    assert (n["PT"], n["aPT"]) == (117, 119)

out = fig_out("Fig3", "Fig3e_rnf152_nfkbiz_pt_vs_apt.pdf")
with PdfPages(out) as pdf:
    for g in GENES:
        fig = m.create_swarm_plot(expr_data, g, value_col="expression", show_means=True, sig_pairs=pairs_to_test,
                                  annotate_sig=True, fdr=True)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)
print("written:", out)

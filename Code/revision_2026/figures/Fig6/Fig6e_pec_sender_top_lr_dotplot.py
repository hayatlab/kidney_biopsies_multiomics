#!/usr/bin/env python
"""Fig. 6e — ligand-receptor pairs with PEC as sender in glomerular niche 20 (Xenium): for each receiver cell type the
pair with the highest LR score, and every significant PEC-sender interaction of these pairs (4 pairs, 16 dots).
Panel        : Main Figures.pdf, Figure 6, panel e
Reads        : $KB_CCI_EXPORT/niche20_lr_table.csv (analyses/12_xenium_niche_cci_liana_crosstalker/03_export_tables.R; default
               $KB_ANALYSIS/sidrah_scripts_20261120/outputs/cci/). The published panel was drawn from `LR_Report.csv` in
               $KB_RESULTS/xenium/sidrah/whole_data_onlyglom_20/CCI_CrossTalkeR/, which holds the same 1,863 interactions
               (checked below when that file is present).
Writes       : $KB_FIG_OUT/Fig6/Fig6e_pec_sender_top_lr_dotplot.pdf
Plotting code: analyses/12_xenium_niche_cci_liana_crosstalker/plot_top_lr_per_receiver.py::plot_lr_dotplot_top_per_receiver
               (the function of the notebook `Fig_6e.ipynb` of Sidrah Maryam, imported unchanged; call as in the notebook)
Illustrator  : colour bar moved below the panel and relabelled "LR Score"; fonts; "PEC as sender" kept as title
Environment  : conda env nichecompass (pandas, matplotlib)
Verified     : 2026-11-20, values equal to the printed panel (rows APP-CD74, COL4A2-ITGA3*ITGB1, COL4A2-ITGA1*ITGB1,
               SPARC-ENG from top to bottom; 16 dots; colour range 4.03-5.55); not pixel-compared (the notebook figure was
               not saved as a file).
"""
import importlib.util
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "config"))
from kbpaths import KB_ANALYSIS, KB_RESULTS, fig_out  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

spec = importlib.util.spec_from_file_location(
    "plot_top_lr_per_receiver", REPO / "analyses" / "12_xenium_niche_cci_liana_crosstalker" / "plot_top_lr_per_receiver.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

EXPORT = os.environ.get("KB_CCI_EXPORT", os.path.join(KB_ANALYSIS, "sidrah_scripts_20261120", "outputs", "cci"))
df = pd.read_csv(os.path.join(EXPORT, "niche20_lr_table.csv"))[["source", "target", "gene_A", "gene_B", "LRScore"]]

report = os.path.join(KB_RESULTS, "xenium", "sidrah", "whole_data_onlyglom_20", "CCI_CrossTalkeR", "LR_Report.csv")
if os.path.exists(report):                                   # the table of the published panel
    rep = pd.read_csv(report, index_col=0)
    key = ["source", "target", "gene_A", "gene_B"]
    a = df.sort_values(key).reset_index(drop=True); b = rep.sort_values(key).reset_index(drop=True)
    assert a[key].equals(b[key]) and np.allclose(a.LRScore, b.LRScore), "exported table differs from LR_Report.csv"

# as in the notebook: df2 = df[df['source']=='PEC']; dff = df2.sort_values(by="LRScore", ascending=False)
dff = df[df["source"] == "PEC"].sort_values(by="LRScore", ascending=False)
fig, ax = m.plot_lr_dotplot_top_per_receiver(dff, sender="PEC")
labels = [t.get_text() for t in ax.get_yticklabels()]
assert labels == ["APP-CD74", "COL4A2-ITGA3*ITGB1", "COL4A2-ITGA1*ITGB1", "SPARC-ENG"], labels
assert len(ax.collections) == 16
out = fig_out("Fig6", "Fig6e_pec_sender_top_lr_dotplot.pdf")
fig.savefig(out, bbox_inches="tight")
print("written:", out)

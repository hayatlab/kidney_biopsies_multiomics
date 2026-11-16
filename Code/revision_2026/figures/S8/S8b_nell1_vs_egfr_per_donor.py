#!/usr/bin/env python
"""Supp. Fig. 8b — snRNA-seq NELL1 pseudobulk expression vs eGFR at biopsy, one point per donor (log1p of the summed
NELL1 counts / summed library sizes x 1e4 over the donor's nuclei of that type, kept if >= 20 nuclei), coloured by
disease; OLS fit with its 95% confidence band; Spearman rho with 95% CI (Fisher z), BH q over the 17 cell-type rows
and n donors in each title.
Panel        : All Supp Figures.pdf, Supplementary Fig. 8, panel b (= the PT, POD, TAL and DCT subpanels of the
               17-subpanel grid drawn here)
Reads        : $KB_ANALYSIS/nell1_app_egfr/tables/pseudobulk_per_unit.csv and egfr_association.csv
                 (analyses/05_nell1_snrna/01_nell1_app_egfr.py via 01_nell1_app_egfr.sh)
               only when those are missing, or with --from-objects: $KB_DATA/charlotte_new_data/
                 adata_20260717_with_clinical.h5ad (analyses/02_xenium_celltyping/18_add_clinical_to_snRNA_ref.py) and
                 $KB_ANALYSIS/vivien_takeover/outputs/metadata/sample_metadata_corrected.csv (eGFR, donor, disease;
                 analyses/03_xenium_factor1_neighbourhoods_cci/24_apply_corrected_metadata.py), from which the snRNA
                 per-donor pseudobulk and its association with eGFR are recomputed
                 (01_nell1_app_egfr.py::load_snrna, ::pseudobulk, ::associations; light aggregation only)
Writes       : $KB_FIG_OUT/S8/S8b_NELL1_snRNA_donor_egfr_scatter.pdf (+ .png)
Plotting code: analyses/05_nell1_snrna/01_nell1_app_egfr.py::scatter_grid (imported unchanged; written under its
               fixed name, then renamed) (original run: 01_nell1_app_egfr.sh, SLURM job 329592, 2026-09-10)
Illustrator  : PT, POD, TAL and DCT subpanels kept and arranged 2 x 2; each subpanel title put on one line; legend
               moved into the PT subpanel; y label "NELL1 pseudobulk, log1p(CP10k)" -> "donor-level pseudobulk mean
               log1p(CP10K)"; suptitle shortened to "snRNA-seq: NELL1 vs eGFR, one point per donor; OLS fit with 95% CI"
Environment  : conda env nichecompass
Verified     : 2026-11-12, PIXEL-IDENTICAL vs $KB_ANALYSIS/nell1_app_egfr/figures/NELL1_snRNA_donor_egfr_scatter.pdf,
               both from the saved tables and with --from-objects
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from xenium_figs import load_script, rename  # noqa: E402
from kbpaths import fig_out  # noqa: E402

na = load_script("05_nell1_snrna/01_nell1_app_egfr.py")


def from_objects():
    """snRNA donor pseudobulk + eGFR association as in 01_nell1_app_egfr.py::main (BH within modality x unit x gene,
    so leaving out the Xenium rows does not change the snRNA q values)."""
    so, sG, slib, sgenes, _ = na.load_snrna(False)
    pb, _ = na.pseudobulk(so, sG, slib, sgenes, "donor", "snRNA")
    return pb, na.associations(pb)


def main():
    ap = argparse.ArgumentParser(description="Supp. Fig. 8b")
    ap.add_argument("--from-objects", action="store_true", help="recompute the pseudobulk even if the tables exist")
    args = ap.parse_args()
    tab = na.HERE / "tables"
    if not args.from_objects and (tab / "pseudobulk_per_unit.csv").exists() and (tab / "egfr_association.csv").exists():
        pb = pd.read_csv(tab / "pseudobulk_per_unit.csv")
        res = pd.read_csv(tab / "egfr_association.csv")
        print(f"[S8b] pseudobulk and associations read from {tab}", flush=True)
    else:
        pb, res = from_objects()
        print("[S8b] pseudobulk and associations computed from the snRNA object", flush=True)
    np.random.seed(na.SEED)                       # global state as in 01_nell1_app_egfr.py::main
    sns.set_style("white")
    outdir = Path(fig_out("S8", "S8b_NELL1_snRNA_donor_egfr_scatter.pdf")).parent
    na.scatter_grid(pb, res, "snRNA", "donor", "NELL1", outdir, "")
    for f in rename(outdir, {f"NELL1_snRNA_donor_egfr_scatter.{e}": f"S8b_NELL1_snRNA_donor_egfr_scatter.{e}"
                             for e in ("pdf", "png")}):
        print(f"[S8b] wrote {f}", flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Factor-1 top-positive PT gene list used by the eGRN views (45, 47, 48): rebuilds the columns those scripts read.

The eGRN scripts 45_egrn_paper_views.py, 47_egrn_factor1_tf_side.py and 48_egrn_tf_centric_panels.py read
  $KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/cis_regulatory_pt_apt/outputs/task3/C4_factor1_top200_positive_genes.csv
and use only its columns `variable` (gene, HGNC symbol) and `Loading` (45, 47: set_index("variable")["Loading"];
48: the `variable` column). In the analysis that file was written by cis_regulatory_pt_apt/scripts/07_task3_nfkb.py
(part C4) at the end of a long chain of reviewer analyses (pseudobulks, limma, peak annotation, motif matching, peak-gene
links) that only add further columns. This script replaces that chain for the two used columns with the same
selection code (07_task3_nfkb.py::part_c):

    MCFA loadings ($KB_ANALYSIS/loadings~hvg~5.csv; MCFA model with 5 factors on HVGs, the manuscript's factor model)
    -> view == "PT" (the PT / aPT / dPT view) and Factor == "Factor1"
    -> the 200 genes with the largest |Loading| (pandas nlargest, keep="first")
    -> of these, the genes with Loading > 0 (80 genes)
    -> sorted by Loading, descending; written with the gene as the index column `variable`.

Checked on 2026-11-12: the output equals the `variable` and `Loading` columns of the original file (80 rows, same order,
identical values); there are no duplicated genes in the PT view and no tie at the 200th |loading|.

Default output = the path the eGRN scripts read. If that file already exists it is NOT overwritten (the original holds
26 more columns used by the reviewer reports): the used columns are compared instead and the script exits with 1 if they
differ. --force overwrites; --out writes elsewhere.
Usage: python analyses/10_atac_pt_apt_egrn/make_factor1_pt_list.py [--out FILE] [--force]
Env: agent_env (pandas). Human, HGNC symbols. Deterministic (no random step).
"""
import argparse
import io
import sys
from pathlib import Path

import pandas as pd

from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py; `source config/env.sh` first)

LOADINGS = Path(f"{KB_ANALYSIS}/loadings~hvg~5.csv")
DEFAULT_OUT = Path(f"{KB_ANALYSIS}/aki_c5_exclusion_20261022/atac_rerun/cis_regulatory_pt_apt/outputs/task3"
                   "/C4_factor1_top200_positive_genes.csv")
USED = ["variable", "Loading"]


def build(loadings=LOADINGS):
    """The Factor-1 top-200-|loading| positive PT genes as a DataFrame indexed by `variable` with column `Loading`."""
    lo = pd.read_csv(loadings)
    f1 = lo[(lo.view == "PT") & (lo.Factor == "Factor1")].copy()
    if f1.variable.duplicated().any():
        raise ValueError("duplicated genes in the PT view of Factor 1")
    f1["abs"] = f1.Loading.abs()
    top200 = f1.nlargest(200, "abs")
    if len(f1) > 200 and (f1["abs"] == top200["abs"].iloc[-1]).sum() > 1:
        raise ValueError("tie at the 200th |loading|: the selection would depend on row order")
    positive = set(top200.loc[top200.Loading > 0, "variable"])
    lst = f1.set_index("variable")[["Loading"]]
    lst = lst[lst.index.isin(positive)]
    return lst.sort_values("Loading", ascending=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT, help=f"output CSV (default: {DEFAULT_OUT})")
    ap.add_argument("--force", action="store_true", help="overwrite an existing output file")
    args = ap.parse_args()
    lst = build()
    print(f"[factor1] {LOADINGS.name}: PT view, Factor1 -> top 200 |loading| -> {len(lst)} positive genes "
          f"(loading {lst.Loading.min():.3f} to {lst.Loading.max():.3f})", flush=True)
    if args.out.exists() and not args.force:
        old = pd.read_csv(args.out)
        new = pd.read_csv(io.StringIO(lst.to_csv()))   # compare as the eGRN scripts read it (same CSV parser)
        same = list(old.columns[:2]) == USED and old[USED].equals(new[USED])
        print(f"[factor1] {args.out} exists; not overwritten (--force to overwrite). Columns {USED} "
              f"{'IDENTICAL to' if same else 'DIFFER from'} the rebuilt list.", flush=True)
        sys.exit(0 if same else 1)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    lst.to_csv(args.out)
    print(f"[factor1] wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()

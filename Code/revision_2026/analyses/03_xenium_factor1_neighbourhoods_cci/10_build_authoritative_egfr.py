#!/usr/bin/env python
"""Rebuild the authoritative per-core eGFR table from the source metadata workbooks.

WHY THIS EXISTS. The `EGFR` column of the new Xenium object
(celltype_concordance/outputs/xenium_labeled.h5ad) is CORRUPTED for a large subset of cores:
a single group value has been broadcast across many donors. 22 CTRL cores all read 67.0 and
15 FSGS cores all read 58.0. Cross-checking Vivien's object shows those same 22 CTRL cores
carried 19 DISTINCT values spanning 35-138, and 49 of the 85 cores common to both objects
disagree. 67.0 is simply the first Appollo control's eGFR; 58.0 likewise for FSGS.

Sources (both under $KB_DATA/metadata/, the files Vivien's
1_preprocess_xenium.py read):
  shared_for_xenium_analysis.xlsx                  sheet 0  -> disease cores, 66 values
  shared_for_xenium_analysis_with_kpmp_egfr.xlsx   'Sheet1' -> the 36 Appollo controls,
                                                               which are absent from sheet 0
Both are keyed on `Xenium ID`, which is the Xenium `TMA_ID` (whitespace-normalised: some
TMA_IDs differ between objects only by a double space).

Output: outputs/authoritative_core_egfr.csv  (core, eGFR, source)

Needs openpyxl -> run under `agent_env`. Downstream figure code reads only the CSV.
"""
from kbpaths import KB_ANALYSIS, KB_DATA  # repository path config (config/kbpaths.py)
import re
from pathlib import Path

import pandas as pd

HERE = Path(f"{KB_ANALYSIS}/vivien_takeover")
MD = Path(f"{KB_DATA}/metadata")
norm = lambda s: re.sub(r"\s+", " ", str(s)).strip()


def main():
    main_sheet = pd.read_excel(MD / "shared_for_xenium_analysis.xlsx")
    main_sheet["core"] = main_sheet["Xenium ID"].map(norm)
    disease = (main_sheet.set_index("core")["egfr_at_time_of_biopsy"]
               .pipe(pd.to_numeric, errors="coerce").dropna())

    ctrl_sheet = pd.read_excel(MD / "shared_for_xenium_analysis_with_kpmp_egfr.xlsx",
                               sheet_name="Sheet1")
    ctrl_sheet["core"] = ctrl_sheet["Xenium ID"].map(norm)
    egcol = [c for c in ctrl_sheet.columns if str(c).strip().lower() == "egfr"][0]
    control = (ctrl_sheet.set_index("core")[egcol]
               .pipe(pd.to_numeric, errors="coerce").dropna())

    overlap = set(disease.index) & set(control.index)
    print(f"[md] disease cores {len(disease)} (distinct {disease.nunique()}); "
          f"control cores {len(control)} (distinct {control.nunique()}, "
          f"range {control.min()}-{control.max()}); overlap {len(overlap)}", flush=True)
    if overlap:
        raise SystemExit(f"unexpected overlap between the two sheets: {sorted(overlap)[:5]}")

    out = pd.concat([
        disease.rename("eGFR").to_frame().assign(source="shared_for_xenium_analysis.xlsx"),
        control.rename("eGFR").to_frame().assign(source="with_kpmp_egfr.xlsx/Sheet1"),
    ])
    out.index.name = "core"
    p = HERE / "outputs" / "authoritative_core_egfr.csv"
    out.to_csv(p)
    print(f"[save] {len(out)} cores, {out['eGFR'].nunique()} distinct -> {p}", flush=True)
    vc = out["eGFR"].value_counts()
    print(f"[check] max cores sharing one value: {vc.max()} "
          f"(the corrupted object column had 22)", flush=True)


if __name__ == "__main__":
    main()

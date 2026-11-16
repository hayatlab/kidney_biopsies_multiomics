#!/usr/bin/env python
"""Adopt the colleague's corrected per-sample metadata (Leon, 2026-09-10).

Leon: "my colleague went over the meta data and found some issues with them. Can you please
adjust it with the info from sample_metadata_milo_korrigiert ... and then rerun milo with
those updated progression classifications".

Input     : <takeover>/sample_metadata_milo_korrigiert - sample_metadata_milo_korrigiert.csv
            (the 63-row table from 22_sample_metadata_table.py, corrected by hand)
Reference : figures/milo_donor_composition/sample_metadata_milo.csv
            (the same table as read from adata_20260717_with_clinical.h5ad, 2026-09-08)
Output    : outputs/metadata/sample_metadata_corrected.csv  -- the per-sample table from now on
            outputs/metadata/metadata_corrections_log.csv   -- every differing cell, old vs new
            outputs/metadata/24_metadata_meta.json

Every difference is typed:
  value   a real change (e.g. a progression class)
  filled  missing in the object, supplied by the correction
  removed present in the object, blank in the correction
  format  same quantity, different notation (decimal comma; prop_gs fraction -> percent)
The Milo bookkeeping columns are recomputed downstream and are not compared; the stale Milo
counts (n_in_milo, n_milo_pos, ...) are dropped from the corrected table.

The shared object is NOT edited; scripts/sample_metadata.py applies this table at load time.
Human GRCh38, HGNC symbols. Seed 42.
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import json
from pathlib import Path

import numpy as np
import pandas as pd

from sample_metadata import sha256

HERE = Path(f"{KB_ANALYSIS}/vivien_takeover")
SRC = HERE / "sample_metadata_milo_korrigiert - sample_metadata_milo_korrigiert.csv"
REF = HERE / "figures" / "milo_donor_composition" / "sample_metadata_milo.csv"
OUT = HERE / "outputs" / "metadata"
BOOKKEEPING = ["milo_arm", "in_milo_contrast", "in_milo_PRIMARY", "excluded_because",
               "classification_basis"]
STALE = ["n_in_milo", "n_milo_pos", "pct_milo_pos", "pct_of_all_milo_pos"]
FP_BELOW, IMPROVED_ABOVE = -3.0, 3.0        # the rule written into classification_basis


def num(v):
    return pd.to_numeric(str(v).replace(",", ".").strip(), errors="coerce")


def blank(v):
    return pd.isna(v) or str(v).strip() == ""


def diff_type(col, o, n):
    """None if the object value `o` and corrected value `n` agree, else the difference type."""
    if blank(o) and blank(n):
        return None
    if blank(o):
        return "filled"
    if blank(n):
        return "removed"
    on, nn = num(o), num(n)
    if not (np.isnan(on) or np.isnan(nn)):
        if np.isclose(on, nn, atol=0.006):
            return "format" if "," in str(o) else None
        if col == "prop_gs" and np.isclose(on * 100, nn, atol=0.006):
            return "format"
        return "value"
    return None if str(o).strip().lower() == str(n).strip().lower() else "value"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    src = pd.read_csv(SRC, dtype=str)
    ref = pd.read_csv(REF, dtype=str)
    if src["sample"].duplicated().any() or set(src["sample"]) != set(ref["sample"]):
        raise SystemExit("sample sets differ or duplicate samples in the corrected sheet")
    src.index, ref.index = src["sample"], ref["sample"]
    ref = ref.loc[src.index]
    print(f"[load] corrected {src.shape}, object-derived {ref.shape}", flush=True)

    # --- every differing cell --------------------------------------------------
    compare = [c for c in src.columns if c in ref.columns and c != "sample"
               and c not in BOOKKEEPING + STALE]
    rows = [{"sample": s, "column": c, "object": ref.at[s, c], "corrected": src.at[s, c],
             "type": t}
            for c in compare for s in src.index
            if (t := diff_type(c, ref.at[s, c], src.at[s, c])) is not None]
    log = pd.DataFrame(rows)
    log.to_csv(OUT / "metadata_corrections_log.csv", index=False)
    by = log.groupby(["type", "column"]).size().rename("n_samples").reset_index()
    print(f"[diff] {len(log)} differing cells\n{by.to_string(index=False)}", flush=True)
    for t in ["value", "removed"]:
        sub = log[log.type == t]
        if len(sub):
            print(f"\n[diff] type={t}\n{sub.to_string(index=False)}", flush=True)
    filled = log[log.type == "filled"].groupby("sample")["column"].apply(len)
    print(f"\n[diff] filled per sample: {filled.to_dict()}", flush=True)

    # spreadsheet artefact: "03.07" for 3.07. Values agree with the object; notation only.
    lz = [(s, c, src.at[s, c]) for c in compare for s in src.index
          if isinstance(src.at[s, c], str) and len(src.at[s, c]) == 5
          and src.at[s, c][0] == "0" and src.at[s, c][2] == "."]
    lz_ok = all(diff_type(c, ref.at[s, c], v) in (None, "format", "filled") for s, c, v in lz)
    print(f"[notation] {len(lz)} leading-zero cells like '03.07' "
          f"(all numerically equal to the object: {lz_ok})", flush=True)

    # --- classification: before/after, and does it follow the stated rule? ------
    old_c, new_c = ref["classification"].fillna("(unclassified)"), \
        src["classification"].fillna("(unclassified)")
    trans = pd.crosstab(old_c.rename("object"), new_c.rename("corrected"))
    print(f"\n[class] transitions (samples)\n{trans}", flush=True)
    jm = src["jm_slope"].map(num)
    rule = pd.Series(np.select([jm < FP_BELOW, jm > IMPROVED_ABOVE, jm.notna()],
                               ["Fast Progressor", "improved", "stable"],
                               default="(unclassified)"), index=src.index)
    off = src.loc[rule != new_c, ["disease", "jm_slope", "classification",
                                  "classification_basis"]].assign(rule=rule[rule != new_c])
    print(f"\n[class] rule FP < {FP_BELOW} / stable / improved > {IMPROVED_ABOVE}: "
          f"{len(src) - len(off)}/{len(src)} samples follow it; exceptions:\n"
          f"{off.to_string()}", flush=True)
    tin_primary = src.index[(src["disease"] == "TIN")
                            & (src["in_milo_PRIMARY"].str.upper() == "TRUE")].tolist()
    print(f"\n[note] the corrected sheet marks TIN sample(s) {tin_primary} as in_milo_PRIMARY; "
          f"04_run_milo.py decides TIN with --exclude-tin, not with this column", flush=True)

    # --- corrected table --------------------------------------------------------
    corr = pd.read_csv(SRC).drop(columns=[c for c in STALE if c in src.columns])
    corr.to_csv(OUT / "sample_metadata_corrected.csv", index=False)
    arm = corr.assign(classification=corr["classification"].fillna("(unclassified)"))
    print(f"\n[save] {OUT / 'sample_metadata_corrected.csv'} {corr.shape}\n"
          f"{pd.crosstab(arm['disease'], arm['classification'], margins=True)}", flush=True)

    (OUT / "24_metadata_meta.json").write_text(json.dumps({
        "source": str(SRC), "source_sha256": sha256(SRC),
        "reference": str(REF), "reference_sha256": sha256(REF),
        "n_differing_cells": int(len(log)),
        "by_type": log["type"].value_counts().to_dict(),
        "by_type_and_column": {f"{r.type}:{r.column}": int(r.n_samples) for r in by.itertuples()},
        "classification_transitions": {f"{o} -> {n}": int(trans.at[o, n])
                                       for o in trans.index for n in trans.columns
                                       if trans.at[o, n] and o != n},
        "rule": f"Fast Progressor jm_slope < {FP_BELOW}; stable {FP_BELOW}..{IMPROVED_ABOVE}; "
                f"improved > {IMPROVED_ABOVE}",
        "rule_exceptions": off.reset_index().to_dict(orient="records"),
        "leading_zero_cells": len(lz), "leading_zero_cells_numerically_equal": bool(lz_ok),
        "tin_marked_in_milo_PRIMARY": tin_primary,
        "dropped_stale_columns": [c for c in STALE if c in src.columns],
        "prop_gs_note": "percent in the corrected table, fraction in the object",
    }, indent=2, default=str))
    print("[done]", flush=True)


if __name__ == "__main__":
    main()

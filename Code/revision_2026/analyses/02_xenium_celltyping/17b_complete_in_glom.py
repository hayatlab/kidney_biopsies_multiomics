"""Step 17b - complete the manual-glomerulus flag `in_glom` for the two Region_1 samples.

RECONSTRUCTION (written 2026-11-12). On 2026-07-09, after 17_recover_pec_from_apt.py, `in_glom` in
outputs/xenium_labeled.h5ad was completed by hand and no script was kept; this file re-implements that step.

What and why
------------
`in_glom` = the cell lies inside a manually drawn (QuPath) glomerulus polygon; it is NOT a NicheCompass niche. The
flag came into xenium_labeled.h5ad from the NicheCompass object, whose parent run had glomerulus polygons for the three
original Xenium samples only, so every cell of the two Region_1 samples added during review read False. Their flags
were computed later by analyses/01_xenium_preprocessing_niches/annotate_in_glom_new_samples.py
(-> in_glom_for_new_samples/in_glom_new_samples.parquet: sample, cell_id, in_glom) but never merged. This step merges
them. Cells are matched on (sample, cell_id), because Xenium re-uses cell ids across regions; only the cells listed
in the parquet change, every other cell keeps its value.

Check (2026-11-12)
------------------
The object between step 17 and this step was not kept, but step 17 does not touch `in_glom`, so the input flag equals
that of the archived pre-step-17 copy (outputs/xenium_labeled.PRE_PEC_RECOVERY_20260709.h5ad: 15,618 True, all in the
three original samples). Applied to it, complete_in_glom() returns exactly the `in_glom` of the final
xenium_labeled.h5ad for all 258,375 cells (dtype bool; 20,081 True = 15,618 + 2,587 + 1,876 in the two Region_1
samples; 0 disagreements). The final object differs from the pre-step-17 copy only in in_glom (4,463 cells),
level2_refined and subcluster_level2 (182 cells each, step 17).

Order matters
-------------
Step 17 used `in_glom` (subcluster in_glom fraction >= 0.25) BEFORE this completion, i.e. with all Region_1 cells
False; that is how the published labels were made. Running this step before step 17 can change the PEC recovery.

Default: updates outputs/xenium_labeled.h5ad in place, like steps 16 and 17 (make a dated backup first). --input reads
another object, --out writes elsewhere, --dry-run only reports. Deterministic (no seed). Env kpmp_repurposing.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import config as C

GLOM_PARQUET = Path(f"{KB_ANALYSIS}/in_glom_for_new_samples/in_glom_new_samples.parquet")


def as_bool(s: pd.Series) -> np.ndarray:
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "1.0"]).to_numpy()


def complete_in_glom(obs: pd.DataFrame, flags: pd.DataFrame) -> tuple[np.ndarray, pd.DataFrame]:
    """Completed `in_glom` (bool, in obs order) and a per-sample summary.

    obs   : needs columns sample, cell_id, in_glom
    flags : in_glom_new_samples.parquet (sample, cell_id, in_glom)
    """
    fkey = flags["sample"].astype(str) + "::" + flags["cell_id"].astype(str)
    if fkey.duplicated().any():
        raise SystemExit(f"{int(fkey.duplicated().sum())} duplicated (sample, cell_id) in {GLOM_PARQUET.name}")
    new = pd.Series(as_bool(flags["in_glom"]), index=fkey.to_numpy())
    key = obs["sample"].astype(str) + "::" + obs["cell_id"].astype(str)
    val = key.map(new)
    have = val.notna().to_numpy()
    if int(have.sum()) != len(flags):
        raise SystemExit(f"{len(flags) - int(have.sum())} of {len(flags)} parquet cells are not in the object")
    before = as_bool(obs["in_glom"])
    after = before.copy()
    after[have] = val[have].astype(bool).to_numpy()
    summ = (pd.DataFrame({"sample": obs["sample"].astype(str).to_numpy(), "before": before, "after": after,
                          "from_parquet": have})
              .groupby("sample").agg(cells=("after", "size"), from_parquet=("from_parquet", "sum"),
                                     true_before=("before", "sum"), true_after=("after", "sum")))
    return after, summ


def main() -> int:
    ap = argparse.ArgumentParser(description="Complete in_glom for the Region_1 samples (step 17b).")
    ap.add_argument("--input", default=str(C.XENIUM_LABELED), help="object to complete (default: %(default)s)")
    ap.add_argument("--out", default=None, help="write the completed object here instead of updating --input")
    ap.add_argument("--dry-run", action="store_true", help="report only, write nothing")
    args = ap.parse_args()

    a = sc.read_h5ad(args.input)
    flags = pd.read_parquet(GLOM_PARQUET)
    after, summ = complete_in_glom(a.obs, flags)
    print(f"[17b] {args.input}\n{summ.to_string()}")
    print(f"[17b] in_glom True {int(as_bool(a.obs['in_glom']).sum())} -> {int(after.sum())} "
          f"({int((after != as_bool(a.obs['in_glom'])).sum())} cells changed)", flush=True)
    if args.dry_run:
        return 0
    a.obs["in_glom"] = after
    out = args.out or args.input
    a.write_h5ad(out)
    print(f"[17b] wrote completed in_glom -> {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

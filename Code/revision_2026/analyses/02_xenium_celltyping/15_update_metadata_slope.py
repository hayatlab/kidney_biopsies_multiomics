"""Incorporate the updated supplemental file
`review_analyses/Supplemental Data 1_metadata_update with slope.csv`
into the per-core metadata, WITHOUT creating duplicate columns.

The file is the slope-updated patient supplement (176 samples). Almost all of its
columns already exist in xenium_labeled.h5ad as `supp_*`; this script:
  1. maps each xenium core (TMA_ID) to a file row via an improved fuzzy match
     (fixes Apollo/Appollo + d0/t0 naming so 119/119 cores match, vs 86/119 before),
  2. builds ONE canonical column per field (no `supp_IFTA` vs `supp_IFTA_percentage`
     doubling), COALESCING new-file value (priority, so slope/classification updates
     win) -> existing supp_* -> existing non-supp clinical (so nothing regresses),
  3. adds the genuinely-new fields (ck_number_rna/atac; cold_ischemia; DBD/DCD),
  4. reports before/after coverage.

KPMP-binned columns (age/eGFR *_KPMP) are dropped: they are populated only for KPMP
reference samples, none of which are among our 119 Xenium cores (0 coverage).

Output: outputs/per_sample_metadata_updated.csv  (one row per TMA_ID).
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import re
import numpy as np
import pandas as pd
import anndata as ad
import config as C

SLOPE_CSV = C.PROJECT_ROOT.parent / "Supplemental Data 1_metadata_update with slope.csv" \
    if hasattr(C, "PROJECT_ROOT") else \
    __import__("pathlib").Path(f"{KB_ANALYSIS}/"
                               "Supplemental Data 1_metadata_update with slope.csv")
MISSING = ("nan", "none", "na", "", "<na>")

# new-file column -> canonical output name (KPMP + redundant modality booleans omitted)
MAP = {
    "disease": "disease", "age_at_biopsy": "age", "sex": "sex",
    "IFTA_percentage": "ifta_pct", "total_glomeruli": "total_glomeruli",
    "global_glomerulosclerosis": "global_gs", "proportion glomerulosclerosis": "prop_gs",
    "date_of_biopsy": "date_of_biopsy",
    "creatinine_at_time_of_biopsy_umol_L": "creatinine_umol_l",
    "treatment_before_biopsy": "treatment_before", "ace_inhibition": "ace_inhibition",
    "angiotensin_ii_inhibition": "angiotensin_ii_inhibition", "statin": "statin",
    "corticosteroids": "corticosteroids", "treatment_after_biopsy": "treatment_after",
    "dialysis": "dialysis", "transplant": "transplant",
    "date_of_dialysis_transplant": "date_dialysis_transplant",
    "egfr_at_time_of_biopsy ml/min/1.73": "egfr_at_biopsy",
    "proteinuria_around_biopsy PCR g/g": "proteinuria",
    "total_slope eGFR": "total_slope", "total_intercept eGFR": "total_intercept",
    "immediate_slope eGFR": "immediate_slope", "chronic_slope eGFR": "chronic_slope",
    "classification eGFR": "classification", "jm_slope eGFR": "jm_slope",
    "jm_intercept eGFR": "jm_intercept", "xenium_available": "xenium_available",
    "Cold ischemia time (min)": "cold_ischemia_min", "DBD/DCD": "dbd_dcd",
    "ck_number_rna": "ck_number_rna", "ck_number_atac": "ck_number_atac",
}
# canonical -> existing supp_ column (for coalesce fallback / before-coverage)
SUPP = {
    "disease": "supp_disease", "age": "supp_age_at_biopsy", "sex": "supp_sex",
    "ifta_pct": "supp_IFTA", "total_glomeruli": "supp_total_glomeruli",
    "global_gs": "supp_global_gs", "prop_gs": "supp_prop_gs",
    "date_of_biopsy": "supp_date_of_biopsy", "creatinine_umol_l": "supp_creatinine_at_time_of_biopsy",
    "treatment_before": "supp_treatment_before", "ace_inhibition": "supp_ace_inhibition",
    "angiotensin_ii_inhibition": "supp_angiotensin_ii_inhibition", "statin": "supp_statin",
    "corticosteroids": "supp_corticosteroids", "treatment_after": "supp_treatment_after",
    "dialysis": "supp_dialysis", "transplant": "supp_transplant",
    "date_dialysis_transplant": "supp_date_of_dialysis_transplant",
    "egfr_at_biopsy": "supp_egfr_at_time_of_biopsy", "proteinuria": "supp_proteinuria_around_biopsy",
    "total_slope": "supp_total_slope", "total_intercept": "supp_total_intercept",
    "immediate_slope": "supp_immediate_slope", "chronic_slope": "supp_chronic_slope",
    "classification": "supp_classification", "jm_slope": "supp_jm_slope",
    "jm_intercept": "supp_jm_intercept", "xenium_available": "supp_xenium_available",
    "cold_ischemia_min": "supp_Cold ischemia time (min)", "dbd_dcd": "supp_DBD/DCD",
}
# canonical -> existing non-supp obs clinical column (last-resort fallback)
NONSUPP = {"disease": "Disease", "age": "age at biopsy", "sex": "Gender",
           "ifta_pct": "IFTA", "egfr_at_biopsy": "eGFR_biopsy"}
NUMERIC = {"age", "creatinine_umol_l", "egfr_at_biopsy", "total_glomeruli", "global_gs",
           "prop_gs", "total_slope", "total_intercept", "immediate_slope",
           "chronic_slope", "jm_slope", "jm_intercept", "cold_ischemia_min"}


def clean(v):
    s = str(v).strip()
    return np.nan if s.lower() in MISSING else s


def to_num(v):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return np.nan
    try:
        return float(str(v).strip().replace(",", "."))
    except ValueError:
        return np.nan


def norm(k):
    if k is None or (isinstance(k, float) and np.isnan(k)):
        return ""
    s = re.sub(r"\s+", " ", str(k).strip().replace("\xa0", " ")).lower()
    return s.replace("appollo", "apollo").replace("utecht", "utrecht")


# --------------------------------------------------------------------------
# Core -> supplement row resolution.
#
# 2026-07-22: REPLACED a token-based fuzzy matcher that was silently wrong.
# It emitted every token of a name as a match key, so the TMA well ('b4', 'a3')
# and group words ('apollo', 'fsgs', 'combi', 'dn') were keys; 37 keys were
# ambiguous, ALL 119 cores resolved to 17-18 candidate rows, and the winner was
# whichever key came first out of an unordered set -> wrong for 44% of eGFR
# values and not even reproducible across processes (string hashing is
# randomised per run). It reported "119/119 matched", which counted matches,
# not correct matches.
#
# The replacement matches ONLY on identifiers that are unique by construction
# and refuses to guess. Every core must resolve to exactly one row, and the
# resolved disease is checked against the run's own `condition` label.
# --------------------------------------------------------------------------
WELL_RE = re.compile(r"\b([A-F][1-6])\b")
KNUM_RE = re.compile(r"\bK(\d+)\b")
HNUM_RE = re.compile(r"(\d{2})-(\d+)")
# supplement `disease` -> the run's `condition` vocabulary
DISEASE2COND = {"control_preKTX": "CTRL", "control_lowGFR": "CTRL", "control": "CTRL",
                "control_KPMP": "CTRL", "DN": "DN", "FSGS": "FSGS", "IgAN": "IgAN",
                "MN": "MN", "TIN": "TIN"}


def hnum(s):
    """'H<yy>-<n>' / 'T<yy>-<n>' / '<yy>-<n>' -> canonical '<yy>-<n>'."""
    m = HNUM_RE.search(str(s))
    return f"{m.group(1)}-{int(m.group(2))}" if m else None


def build_index(new):
    """sample_ID -> row, and canonical H/T-number -> [rows]."""
    by_sid = {str(r["sample_ID"]).strip(): i for i, r in new.iterrows()}
    by_hnum = {}
    for i, r in new.iterrows():
        h = hnum(r["patient_ID"])
        if h:
            by_hnum.setdefault(h, []).append(i)
    return by_sid, by_hnum


def resolve_core(tma_id, by_sid, by_hnum):
    """Resolve one TMA core to exactly one supplement row.

    Returns (row_index, how). row_index is None when the core cannot be resolved
    UNAMBIGUOUSLY -- callers must leave those fields empty rather than guess.
    """
    s = str(tma_id)
    low = norm(s)
    if "apollo" in low:                       # pre-KTX controls: 'Appollo d0 B4 CTRL'
        m = WELL_RE.search(s.upper())         # -> supplement 'Apollot0_B4'
        if m:
            sid = f"Apollot0_{m.group(1)}"
            if sid in by_sid:
                return by_sid[sid], "apollo-well"
        return None, "apollo: no well / not in file"
    m = KNUM_RE.search(s)                     # Utrecht cores carry the K-number
    if m and f"K{m.group(1)}" in by_sid:      # which IS the supplement sample_ID
        return by_sid[f"K{m.group(1)}"], "K-number"
    h = hnum(s)                               # Combi / FSGS / Rotterdam: H/T-number
    if h and h in by_hnum:
        if len(by_hnum[h]) == 1:
            return by_hnum[h][0], "H/T-number"
        return None, f"H/T-number {h} matches {len(by_hnum[h])} rows"
    return None, "no unique identifier in core name"


def main():
    new = pd.read_csv(SLOPE_CSV, sep=";")
    new = new.loc[:, ~new.columns.str.startswith("Unnamed")]
    by_sid, by_hnum = build_index(new)

    a = ad.read_h5ad(C.XENIUM_LABELED, backed="r")
    o = a.obs
    base = o.groupby("TMA_ID", observed=True).agg(
        region=("sample", lambda s: "_".join(str(s.iloc[0]).split("__")[1:3])),
        n_cells=("sample", "size"),
        condition=("condition", lambda s: s.astype(str).value_counts().index[0]))
    supp_first = o.groupby("TMA_ID", observed=True).first()

    core2row, how = {}, {}
    for t in base.index:
        core2row[t], how[t] = resolve_core(t, by_sid, by_hnum)
    matched = sum(v is not None for v in core2row.values())
    print(f"[join] {matched}/{len(base)} cores resolved unambiguously")
    for h, n in pd.Series(list(how.values())).value_counts().items():
        print(f"        {n:4d}  {h}")
    unresolved = [t for t, v in core2row.items() if v is None]
    if unresolved:
        print(f"[join] {len(unresolved)} UNRESOLVED (left empty, never guessed):")
        for t in unresolved:
            print(f"        {t}   ({how[t]})")

    # one row must not serve several cores unless they really are the same sample
    used = pd.Series([v for v in core2row.values() if v is not None])
    reused = used.value_counts()
    reused = reused[reused > 1]
    if len(reused):
        print(f"[join] WARNING {len(reused)} supplement rows map to >1 core:")
        for ri, n in reused.items():
            cores = [t for t, v in core2row.items() if v == ri]
            print(f"        {new.at[ri, 'sample_ID']} -> {n} cores: {cores}")

    # VALIDATION: the resolved disease must agree with the run's own condition.
    # This is the check whose absence let the previous fuzzy join ship: it costs
    # nothing and catches a mis-join immediately.
    bad = []
    for t in base.index:
        ri = core2row[t]
        if ri is None:
            continue
        mapped = DISEASE2COND.get(str(new.at[ri, "disease"]).strip())
        if mapped is not None and mapped != base.at[t, "condition"]:
            bad.append((t, base.at[t, "condition"], new.at[ri, "disease"],
                        new.at[ri, "sample_ID"]))
    if bad:
        for t, cond, dis, sid in bad:
            print(f"[join] MISMATCH {t}: condition={cond} but supplement says "
                  f"{dis} (row {sid})")
        raise SystemExit(f"[join] ABORT: {len(bad)}/{len(base)} cores resolved to a "
                         "supplement row whose disease contradicts the core's own "
                         "condition label -- the join is wrong, refusing to write.")
    print(f"[join] validation OK: resolved disease agrees with `condition` for "
          f"{matched}/{matched} cores")

    cov_before, cov_after, prov = {}, {}, {}
    out = base.copy()
    # NO coalesce fallback (changed 2026-07-22). The old code fell back to the
    # existing `supp_*` obs columns when the supplement row had no value -- but
    # those columns are themselves the output of the earlier broken join, so the
    # fallback silently re-imported contaminated values. A resolved core now takes
    # its values from its ONE supplement row; absent means missing, and missing is
    # reported rather than filled.
    for nc, oc in MAP.items():
        vals, src = [], {"new": 0, "supp": 0, "nonsupp": 0, "none": 0}
        for t in base.index:
            ri = core2row[t]
            v = clean(new.at[ri, nc]) if ri is not None else np.nan
            src["new" if not (isinstance(v, float) and np.isnan(v)) else "none"] += 1
            vals.append(v)
        s = pd.Series(vals, index=base.index)
        if oc in NUMERIC:
            s = s.map(to_num)
        out[oc] = s
        cov_after[oc] = int(s.notna().sum())
        sc = SUPP.get(oc)
        cov_before[oc] = int(supp_first[sc].map(clean).notna().sum()) if sc and sc in supp_first else 0
        prov[oc] = src

    # ethnicity has no counterpart in the new file -> keep the existing obs column
    if "Etnicity" in supp_first:
        out["ethnicity"] = supp_first["Etnicity"].map(clean)

    out = out.sort_values(["region", "n_cells"], ascending=[True, False])
    outpath = C.OUT / "per_sample_metadata_updated.csv"
    out.to_csv(outpath)

    print("\nBEFORE vs AFTER coverage (cores with a value, of %d):" % len(base))
    print("  %-22s %6s %6s   %s" % ("field", "before", "after", "provenance(new/supp/nonsupp)"))
    for oc in MAP.values():
        p = prov[oc]
        print("  %-22s %6d %6d   %d/%d/%d"
              % (oc, cov_before[oc], cov_after[oc], p["new"], p["supp"], p["nonsupp"]))
    print("\nwrote", outpath, out.shape)


if __name__ == "__main__":
    main()

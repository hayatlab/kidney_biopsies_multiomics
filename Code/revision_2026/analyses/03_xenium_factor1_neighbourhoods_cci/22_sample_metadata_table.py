#!/usr/bin/env python
"""Per-sample metadata table for every sample in the snRNA object, with Milo inclusion.

Leon (2026-09-08): "give me the per sample_id with all the meta data including what was
considered fast progressing and so on like how it was classified. Do that for all available
samples in the single nuc data and make an extra column whether it was included in milo".

One row per `sample` in adata_20260717_with_clinical.h5ad — all 63, including the ones the
Milo contrast drops — carrying every obs column that is CONSTANT within a sample (46 of the
77; the other 31 are per-cell measurements like nCount_RNA and cannot be tabulated per
sample). Sample-level constancy is verified, not assumed.

Milo columns:
  in_milo_contrast        classification is Fast Progressor or stable
  in_milo_PRIMARY         and disease is not TIN (the PRIMARY run excludes TIN)
  milo_arm                Fast_Progressor (test) / stable (reference) / excluded
  excluded_because        why a sample is not in the PRIMARY run
  n_milo_pos, pct_milo_pos   Milo-positive cells contributed by this sample
  pct_of_all_milo_pos     this sample's share of every Milo+ cell in the run

How the classification relates to the data (measured here, not asserted): `jm_slope` is the
variable that tracks it. The script recomputes the jm_slope range of each class every run and
writes it to the JSON, so the statement stays true if the metadata changes.

CORRECTED MODE (Leon, 2026-09-10): `--sample-metadata outputs/metadata/sample_metadata_corrected.csv`
overlays the colleague's corrected per-sample table on the object's columns (classification,
the filled-in clinical data for K68/K95/K98, ...), keeps the colleague's
`classification_basis` text, and recomputes only the Milo bookkeeping. `in_milo_PRIMARY`
follows the standing decision (TIN excluded, Leon 2026-09-03 and again 2026-09-10), not the
TRUE the corrected sheet carries for TIN. Pair it with `--tag` of the corrected Milo run and
`--out-suffix _corrected` so the object-derived table (the reference 24 compares against) is
not overwritten.

Human GRCh38, HGNC symbols. Seed 42.
"""
from kbpaths import KB_ANALYSIS, KB_DATA  # repository path config (config/kbpaths.py)
import argparse, json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

HERE = Path(f"{KB_ANALYSIS}/vivien_takeover")
SNRNA = f"{KB_DATA}/charlotte_new_data/adata_20260717_with_clinical.h5ad"
OUT = HERE / "figures" / "milo_donor_composition"
SEED, ALPHA = 42, 0.05
GROUPS = ["Fast Progressor", "stable"]
LEAD = ["sample", "patient_ID", "disease", "classification", "milo_arm",
        "in_milo_contrast", "in_milo_PRIMARY", "excluded_because",
        "n_cells", "n_milo_pos", "pct_milo_pos", "pct_of_all_milo_pos",
        "jm_slope", "classification_basis", "egfr_at_time_of_biopsy"]
RECOMPUTED = ["milo_arm", "in_milo_contrast", "in_milo_PRIMARY", "n_in_milo", "n_milo_pos",
              "pct_milo_pos", "pct_of_all_milo_pos"]


def overlay_corrected(tab, path):
    """Replace the object-derived columns with the corrected table's, sample by sample."""
    corr = pd.read_csv(path)
    corr.index = corr["sample"].astype(str)
    missing = sorted(set(tab.index) - set(corr.index))
    if missing:
        raise SystemExit(f"[meta] object samples absent from {path}: {missing}")
    bad = [s for s in tab.index if int(corr.at[s, "n_cells"]) != int(tab.at[s, "n_cells"])]
    if bad:
        raise SystemExit(f"[meta] n_cells disagree with the object for {bad}")
    for c in corr.columns:
        if c not in RECOMPUTED + ["sample", "n_cells"]:
            tab[c] = corr[c].reindex(tab.index)
    print(f"[meta] overlaid {path}", flush=True)
    return tab


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-metadata", default=None,
                    help="corrected per-sample CSV overlaid on the object's columns")
    ap.add_argument("--tag", default="_scvi_noTIN", help="Milo run whose Milo+ cells are counted")
    ap.add_argument("--out-suffix", default="", help="appended to the output file names")
    args = ap.parse_args()
    corrected = bool(args.sample_metadata)
    np.random.seed(SEED)
    OUT.mkdir(parents=True, exist_ok=True)

    ad = sc.read_h5ad(SNRNA, backed="r")
    obs = ad.obs.copy()
    key = obs["sample"].astype(str)
    g = obs.groupby(key, observed=True)
    print(f"[load] {ad.shape}; {key.nunique()} samples, {len(obs.columns)} obs columns", flush=True)

    # keep only columns that are genuinely constant within a sample
    const = []
    for c in obs.columns:
        try:
            if g[c].apply(lambda v: v.astype(str).nunique()).max() <= 1:
                const.append(c)
        except Exception:
            continue
    dropped = [c for c in obs.columns if c not in const]
    print(f"[cols] {len(const)} sample-level columns kept; {len(dropped)} per-cell columns "
          f"dropped (e.g. {', '.join(dropped[:4])})", flush=True)

    tab = g.first()[const].copy()
    tab["sample"] = tab.index.astype(str)
    tab["n_cells"] = g.size()
    for c in const:
        if tab[c].dtype.name == "category":
            tab[c] = tab[c].astype(str)
    if corrected:
        tab = overlay_corrected(tab, args.sample_metadata)

    cl = tab["classification"].astype(str)
    tab["in_milo_contrast"] = cl.isin(GROUPS)
    tab["in_milo_PRIMARY"] = tab["in_milo_contrast"] & (tab["disease"].astype(str) != "TIN")
    tab["milo_arm"] = np.where(cl == "Fast Progressor", "Fast_Progressor (test)",
                        np.where(cl == "stable", "stable (reference)", "excluded"))
    given = (tab["excluded_because"].fillna("").astype(str)
             if corrected and "excluded_because" in tab else pd.Series("", index=tab.index))
    not_arm = np.where(given != "", given, "classification='" + cl + "' is not a contrast arm")
    tab["excluded_because"] = np.where(
        tab["in_milo_PRIMARY"], "",
        np.where(~tab["in_milo_contrast"], not_arm, "disease=TIN, excluded from the PRIMARY run"))

    # --- Milo+ contribution, from the chosen run -----------------------------------
    cells = pd.read_parquet(HERE / "outputs" / "milo" / f"milo_cell_scores{args.tag}.parquet")
    cells["signif"] = cells["milo_adj_log10pval"] > -np.log10(ALPHA)
    smp = cells["sample"].astype(str)
    pos = cells.groupby(smp, observed=True)["signif"].agg(n_in_milo="size", n_milo_pos="sum")
    tab = tab.drop(columns=[c for c in ["n_in_milo", "n_milo_pos"] if c in tab]).join(pos)
    tab["n_in_milo"] = tab["n_in_milo"].fillna(0).astype(int)
    tab["n_milo_pos"] = tab["n_milo_pos"].fillna(0).astype(int)
    tot = int(tab["n_milo_pos"].sum())
    tab["pct_milo_pos"] = np.where(tab["n_in_milo"] > 0,
                                   tab["n_milo_pos"] / tab["n_in_milo"].replace(0, np.nan) * 100,
                                   np.nan).round(4)
    tab["pct_of_all_milo_pos"] = (tab["n_milo_pos"] / tot * 100).round(3)
    arm_ok = tab.loc[tab["n_in_milo"] > 0, "in_milo_PRIMARY"].all() if "noTIN" in args.tag else True
    if not arm_ok:
        raise SystemExit(f"cells from samples outside in_milo_PRIMARY appear in {args.tag}")

    # --- what the classification actually tracks ---------------------------------
    jm = pd.to_numeric(tab["jm_slope"], errors="coerce")
    rng = {c: {"n": int(jm[cl == c].notna().sum()),
               "min": float(jm[cl == c].min()), "max": float(jm[cl == c].max()),
               "median": float(jm[cl == c].median())}
           for c in cl.unique() if jm[cl == c].notna().any()}
    fp_max = rng.get("Fast Progressor", {}).get("max", np.nan)
    st_min = rng.get("stable", {}).get("min", np.nan)
    boundary = (jm <= fp_max) & (cl == "stable")
    if not corrected:        # the corrected sheet carries its own, rule-based basis text
        tab["classification_basis"] = np.where(
            cl.isin(GROUPS),
            "jm_slope " + jm.round(3).astype(str) + " (Fast Progressor <= "
            + f"{fp_max:.3f}" + ", stable >= " + f"{st_min:.3f}" + ")",
            np.where(cl == "regenerate", "regenerate: all TIN, highest positive jm_slope",
                     "no progression class (controls / unlabelled)"))
        tab.loc[boundary, "classification_basis"] += "  [BOUNDARY CASE]"
    print(f"[class] jm_slope by class: " + "; ".join(
        f"{k} n={v['n']} [{v['min']:.3f}, {v['max']:.3f}]" for k, v in rng.items()), flush=True)
    print(f"[class] {int(boundary.sum())} stable sample(s) sit below the Fast Progressor "
          f"maximum ({fp_max:.3f}): {list(tab.index[boundary])}", flush=True)

    cols = [c for c in LEAD if c in tab.columns]
    cols += [c for c in tab.columns if c not in cols]
    tab = tab[cols].rename_axis(None).sort_values(
        ["in_milo_PRIMARY", "classification", "disease", "sample"],
        ascending=[False, True, True, True])

    stem = f"sample_metadata_milo{args.out_suffix}"
    tab.to_csv(OUT / f"{stem}.csv", index=False)
    try:
        tab.to_excel(OUT / f"{stem}.xlsx", index=False)
        print("[save] .xlsx written too", flush=True)
    except Exception as e:
        print(f"[save] no .xlsx ({type(e).__name__}: {e}); CSV only", flush=True)

    print(f"\n[counts] classification x disease (samples):\n"
          f"{pd.crosstab(tab['classification'].fillna('nan'), tab['disease'].fillna('nan'), margins=True)}",
          flush=True)
    print(f"\n[counts] Milo arm:\n"
          f"{tab.groupby('milo_arm', observed=True).agg(n_samples=('sample','size'), n_cells=('n_cells','sum'), milo_pos=('n_milo_pos','sum'))}",
          flush=True)
    print(f"\n[top contributors of Milo+ cells]\n"
          f"{tab.nlargest(8, 'n_milo_pos')[['sample','patient_ID','disease','classification','n_milo_pos','pct_milo_pos','pct_of_all_milo_pos']].to_string(index=False)}",
          flush=True)

    (OUT / f"22_sample_metadata_meta{args.out_suffix}.json").write_text(json.dumps({
        "source": SNRNA, "sample_metadata": args.sample_metadata, "milo_tag": args.tag,
        "n_samples": int(len(tab)), "n_columns": int(tab.shape[1]),
        "sample_level_columns": const, "per_cell_columns_dropped": dropped,
        "jm_slope_by_classification": rng,
        "boundary_samples": [str(i) for i in boundary.index[boundary]],
        "n_in_milo_PRIMARY": int(tab["in_milo_PRIMARY"].sum()),
        "total_milo_pos": tot, "alpha": ALPHA, "seed": SEED}, indent=2, default=str))
    print(f"\n[done] {OUT}/{stem}.csv  ({tab.shape[0]} rows x {tab.shape[1]} cols)", flush=True)


if __name__ == "__main__":
    main()

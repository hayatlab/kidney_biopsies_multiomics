#!/usr/bin/env python
"""Per-core Xenium quality metrics + biopsy age. See PLAN.md.
Human GRCh38, HGNC symbols (468-gene hKidney panel). No randomness. Login node, ~2 min."""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import sys, json, re
from pathlib import Path
import numpy as np, pandas as pd, h5py, scipy.sparse as sp
from scipy.spatial import ConvexHull
from scipy.special import gammaln
from anndata.io import read_elem

D = Path(f"{KB_ANALYSIS}/xenium_biopsy_age_qc")
XEN = f"{KB_ANALYSIS}/celltype_concordance/outputs/xenium_labeled.h5ad"
META = f"{KB_ANALYSIS}/celltype_concordance/outputs/per_sample_metadata_updated.csv"   # condition/disease (identifier-resolved join)
META_NEW = f"{KB_ANALYSIS}/xenium_biopsy_age_qc/Sample_clinical_Data_Metadata.xlsx"  # Leon, 2026-10-15: authoritative dates/clinical
FACT = f"{KB_ANALYSIS}/vivien_takeover/outputs/xenium_factor_scores.parquet"
RUN_DATE = pd.Timestamp("2025-04-16")      # experiment.xenium run_start_time, all 5 regions
# 2026-10-16: collaborator states that all control (Apollo pre-transplant) cores were biopsied in 2000, contradicting the
# 2021-2024 dates in Sample_clinical_Data_Metadata.xlsx. Applied as an explicit override (year only -> 1 July 2000),
# flagged per core in `biopsy_date_overridden`. Set to None to use the metadata file as is.
CTRL_YEAR_OVERRIDE = 2000
BLOCK = {"0041734_Region_1": "Combi", "0041734_Region_3": "Apollo_t0", "0042348_Region_1": "FSGS",
         "0042348_Region_2": "Utrecht2", "0042348_Region_3": "Utrecht1"}   # Sample overview xlsx, Tabelle1
RAREFY_DEPTH = 100
SMOKE = "--smoke" in sys.argv

with h5py.File(XEN, "r") as f:
    obs = read_elem(f["obs"])
    X = read_elem(f["layers/counts"]).tocsr()
    xy = read_elem(f["obsm/spatial"])
    var = read_elem(f["var"]).index
assert X.shape == (len(obs), 468), X.shape
obs["region"] = obs["sample"].astype(str).map(lambda s: "_".join(s.split("__")[1:3]))
obs["tma_block"] = obs["region"].map(BLOCK)
obs["TMA_ID"] = obs["TMA_ID"].astype(str)
if SMOKE:
    keep = obs.TMA_ID.isin(obs.TMA_ID.unique()[:6]); obs, X, xy = obs[keep], X[keep.values], xy[keep.values]

# ---- per-cell quantities --------------------------------------------------------------
Xd = X.astype(np.int64)
tot = np.asarray(Xd.sum(1)).ravel()
ngenes = np.diff(Xd.indptr)
obs["mat_total"] = tot; obs["n_genes"] = ngenes
# rarefaction: expected distinct genes at depth d (exact hypergeometric; Hurlbert 1971)
d = RAREFY_DEPTH
N_rep = np.repeat(tot, ngenes)
n = Xd.data
with np.errstate(invalid="ignore"):
    logp = gammaln(N_rep - n + 1) - gammaln(N_rep - n - d + 1) - gammaln(N_rep + 1) + gammaln(N_rep - d + 1)
term = np.where(N_rep - n >= d, 1 - np.exp(logp), 1.0)
term[N_rep < d] = np.nan
rare = np.full(len(tot), np.nan)
ok = tot >= d
rare_sum = np.add.reduceat(np.nan_to_num(term), Xd.indptr[:-1]) if len(n) else np.zeros(len(tot))
rare_sum[ngenes == 0] = 0
rare[ok] = rare_sum[ok]
obs["genes_at_100"] = rare
print("cells with >=100 transcripts (complexity denominators):", ok.sum(), "of", len(tot))
assert np.nanmax(rare) <= d + 1e-6
# sanity: cells-table transcript_counts vs matrix total
diff = (obs.transcript_counts.astype(float) - tot)
print("transcript_counts - matrix total: median", np.median(diff), "max|.|", np.abs(diff).max())

fs = pd.read_parquet(FACT)
fcols = [c for c in fs.columns if re.match(r"Factor\d_ulm$", c)]
fs = fs.loc[obs.index, fcols + ["view"]] if SMOKE else fs.reindex(obs.index)[fcols + ["view"]]
assert fs[fcols].notna().all().all(), "factor scores missing for some cells"
for c in fcols: obs[c] = fs[c].values
obs["view"] = fs["view"].values
obs["lvl2"] = obs["level2_refined"].astype(str)

# ---- per-core aggregation --------------------------------------------------------------
rows = []
for tid, g in obs.groupby("TMA_ID", sort=False):
    pts = xy[obs.index.get_indexer(g.index)]
    hull_mm2 = ConvexHull(pts).volume / 1e6 if len(g) >= 3 else np.nan
    tx = g.transcript_counts.astype(float).values
    r = dict(TMA_ID=tid, region=g.region.iloc[0], tma_block=g.tma_block.iloc[0], condition=str(g.condition.iloc[0]),
             n_cells=len(g), hull_area_mm2=hull_mm2, cell_density_per_mm2=len(g) / hull_mm2,
             cell_area_frac=g.cell_area.sum() / 1e6 / hull_mm2,
             median_tx_per_cell=np.median(tx), median_genes_per_cell=np.median(g.n_genes),
             tx_per_100um2=100 * tx.sum() / g.cell_area.sum(),
             genes_at_100=np.nanmedian(g.genes_at_100), n_cells_ge100=int((tx >= 100).sum()),
             frac_cells_lt20=np.mean(tx < 20), frac_cells_lt50=np.mean(tx < 50),
             negprobe_per_10k_tx=1e4 * g.control_probe_counts.sum() / (tx.sum() + g.control_probe_counts.sum()),
             unassigned_per_10k_tx=1e4 * g.unassigned_codeword_counts.sum() / (tx.sum() + g.unassigned_codeword_counts.sum()),
             negprobe_count=int(g.control_probe_counts.sum()), unassigned_count=int(g.unassigned_codeword_counts.sum()),
             median_cell_area=np.median(g.cell_area), median_nucleus_area=np.nanmedian(g.nucleus_area), frac_nucleus_na=np.mean(g.nucleus_area.isna()),
             median_conf_level2=np.median(g.conf_level2), frac_conf_lt05=np.mean(g.conf_level2 < 0.5))
    for ct in ["PT", "TAL"]:
        s = g[g.lvl2 == ct]
        r[f"n_{ct}"] = len(s)
        r[f"{ct}_median_tx_per_cell"] = np.median(s.transcript_counts) if len(s) >= 20 else np.nan
        r[f"{ct}_median_genes_per_cell"] = np.median(s.n_genes) if len(s) >= 20 else np.nan
    for c in fcols:
        r[f"{c}_mean_all"] = g[c].mean()
        pv = g[g.view == "PT"]
        r[f"{c}_mean_PTview"] = pv[c].mean() if len(pv) >= 20 else np.nan
    for lv in ["PT", "aPT", "dPT", "TAL", "Immune", "FIB", "POD", "PEC"]:
        r[f"frac_{lv}"] = np.mean(g.lvl2 == lv)
    rows.append(r)
core = pd.DataFrame(rows)

# ---- metadata: biopsy date / age ----------------------------------------------------
m_old = pd.read_csv(META)[["TMA_ID", "condition", "disease", "date_of_biopsy"]].rename(columns={"condition": "condition_meta", "date_of_biopsy": "date_old"})
mn = pd.read_excel(META_NEW, sheet_name="Samples")
mn = mn[mn.xenium].rename(columns={"xenium_TMA_IDs": "TMA_ID", "age_at_biopsy": "patient_age", "sex_clinical": "sex",
                                   "egfr_at_biopsy_final": "egfr", "IFTA": "ifta_raw"})
mn = mn[["TMA_ID", "date_of_biopsy", "patient_age", "sex", "egfr", "ifta_raw"]]
m = m_old.merge(mn, on="TMA_ID", how="left", validate="1:1")
assert m.date_of_biopsy.notna().all(), "new metadata lacks a biopsy date for some core"
m["ifta"] = pd.to_numeric(m.ifta_raw.astype(str).str.replace("[<>%]", "", regex=True), errors="coerce")
m["biopsy_year_differs_from_old"] = m.date_of_biopsy.astype(str).str[:4] != m.date_old.astype(str).str[:4]
print("biopsy year differs from the 2026-07 supplement for", int(m.biopsy_year_differs_from_old.sum()), "cores")
m = m.drop(columns=["date_old", "ifta_raw"])
m["biopsy_date_overridden"] = False
if CTRL_YEAR_OVERRIDE is not None:
    is_ctrl = m.condition_meta.eq("CTRL")
    m.loc[is_ctrl, "date_of_biopsy"] = str(CTRL_YEAR_OVERRIDE); m.loc[is_ctrl, "biopsy_date_overridden"] = True
    print(f"override: {int(is_ctrl.sum())} control cores set to biopsy year {CTRL_YEAR_OVERRIDE} (metadata file says 2021-2024)")
ds = m.date_of_biopsy.astype(str)
m["biopsy_date_is_year_only"] = ds.str.len() == 4
m["biopsy_year"] = ds.str.extract(r"(\d{4})")[0].astype(int)
m["biopsy_date"] = pd.to_datetime(np.where(m.biopsy_date_is_year_only, ds + "-07-01", ds), errors="raise")
m["biopsy_age_years"] = (RUN_DATE - m.biopsy_date).dt.days / 365.25
core = core.merge(m, on="TMA_ID", how="left", validate="1:1")
assert core.condition_meta.eq(core.condition).all()
assert core.biopsy_age_years.notna().all(), "cores without biopsy date"
core = core.drop(columns="condition_meta")

# ---- pseudonymise ---------------------------------------------------------------------
core = core.sort_values(["region", "TMA_ID"]).reset_index(drop=True)
core.insert(0, "core_id", [f"core_{i:03d}" for i in range(1, len(core) + 1)])
priv = D / "outputs/private"; priv.mkdir(parents=True, exist_ok=True)
core[["core_id", "TMA_ID", "region"]].to_csv(priv / ("core_id_map_smoke.csv" if SMOKE else "core_id_map.csv"), index=False)
pub = core.drop(columns=["TMA_ID", "biopsy_date", "date_of_biopsy"])   # no dates in the shareable table
out = D / "outputs" / ("core_qc_metrics_smoke.csv" if SMOKE else "core_qc_metrics.csv")
pub.to_csv(out, index=False)
meta = dict(n_cells=int(len(obs)), n_cores=int(len(core)), run_date=str(RUN_DATE.date()), rarefy_depth=d,
            cells_ge100=int(ok.sum()), factor_cols=fcols, year_only_dates=int(core.biopsy_date_is_year_only.sum()), ctrl_year_override=CTRL_YEAR_OVERRIDE, n_overridden=int(core.biopsy_date_overridden.sum()),
            missing_pct={c: round(100 * core[c].isna().mean(), 2) for c in pub.columns if core[c].isna().any()},
            versions={"numpy": np.__version__, "pandas": pd.__version__, "scipy": __import__("scipy").__version__},
            smoke=SMOKE, species="human GRCh38", gene_ids="HGNC symbols")
json.dump(meta, open(D / "outputs" / ("01_meta_smoke.json" if SMOKE else "01_meta.json"), "w"), indent=1)
print(json.dumps(meta, indent=1)); print(pub.describe().T[["min", "50%", "max"]].round(3).to_string())

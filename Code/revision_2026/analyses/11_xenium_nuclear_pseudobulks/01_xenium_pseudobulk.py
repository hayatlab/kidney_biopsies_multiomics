#!/usr/bin/env python
"""Cell type x sample pseudobulks of the Xenium data with decoupler's pseudobulk function, for the MCFA team's ambient-RNA
figure (reviewer comment 6): they want to test whether gene-expression changes removed by their strict feature mask are
genuinely shared across cell types on an imaging platform.

Input : celltype_concordance/outputs/xenium_labeled.h5ad (258,375 cells x 468 genes, layers['counts'] = raw transcript counts,
        HGNC symbols, human). Labels level2_refined (manuscript labels) and the MCFA view (level1 pooling used by the factor
        model: PT+aPT+dPT -> PT, FIB+MES+VSM-P -> Stromal, CNT+PC -> CNT_PC, IC-A+IC-B -> IC).
Sample: TMA core = one biopsy = one donor (119 cores, no donor has two cores). Sample id = snRNA-seq K number (ck_number_rna)
        where the core has one, otherwise a pseudonym xen_core_###; the key to the TMA IDs (pathology numbers) stays in
        outputs/private/ and is not shared.
Method: dc.pp.pseudobulk(mode='sum', min_cells=0, min_counts=0) -> unfiltered, so the recipient can filter with the
        psbulk_cells / psbulk_counts columns; a second, decoupler-default-filtered version (min_cells 10, min_counts 1000) is
        written alongside. eGFR per core from the corrected sources (workbook > per-sample metadata > object), as in
        mcfa_xenium_celltype_coordination_20261026.
Outputs (outputs/): pseudobulk_<grouping>_unfiltered.h5ad / .csv (counts, pseudobulks x genes), *_obs.csv, README.md.
Identifier scan at the end: no H-number, TMA_ID, patient_id or date may appear in any shared file.
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import re, json, os, sys, inspect, zipfile
import numpy as np, pandas as pd, scanpy as sc, anndata as ad, decoupler as dc
B = f"{KB_ANALYSIS}/mcfa_xenium_pseudobulk_20261026"
OBJ = f"{KB_ANALYSIS}/celltype_concordance/outputs/xenium_labeled.h5ad"
AUTH = f"{KB_ANALYSIS}/vivien_takeover/outputs/authoritative_core_egfr.csv"
META = f"{KB_ANALYSIS}/celltype_concordance/outputs/per_sample_metadata_updated.csv"
L2L1 = json.load(open(f"{KB_ANALYSIS}/vivien_takeover/outputs/01_factor_scoring_meta.json"))
print("decoupler", dc.__version__, "scanpy", sc.__version__, "anndata", ad.__version__)
print("pseudobulk signature:", inspect.signature(dc.pp.pseudobulk))

a = sc.read_h5ad(OBJ)
print(a)
a.X = a.layers["counts"].copy()
assert a.X.max() > 10 and np.allclose(a.X[:1000].toarray() % 1, 0), "counts layer is not integer counts"
obs = a.obs
obs["TMA_norm"] = obs["TMA_ID"].astype(str).str.replace(r"\s+", " ", regex=True).str.strip()
# ---- sample ids -------------------------------------------------------------------------------------------------------
def first_ck(s):
    v = [x for x in s.astype(str) if x not in ("nan", "None", "", "NaN")]
    return v[0] if v else ""
cores = obs.groupby("TMA_norm", observed=True).agg(ck=("ck_number_rna", first_ck), run=("sample", "first"), region=("region", "first"),
                                                   cond=("condition", "first")).reset_index()
# fill K numbers missing in the object from the per-sample metadata table (keyed by TMA id)
_m = pd.read_csv(META); _m["k"] = _m["TMA_ID"].astype(str).str.replace(r"\s+", " ", regex=True).str.strip()
_mk = _m.dropna(subset=["ck_number_rna"]).drop_duplicates("k").set_index("k")["ck_number_rna"].astype(str)
filled = cores.TMA_norm.map(_mk).fillna("")
conflict = cores[(cores.ck != "") & (filled != "") & (cores.ck != filled)]
# the object copies KB199 / KB135 onto two cores of the newer run whose diagnoses (DN, TIN) match the metadata's K numbers
# (CK404, CK436); the per-sample metadata table (updated 2026-07-22) is taken as authoritative where the two disagree
for i in conflict.index: print(f"  K-number conflict: object {cores.loc[i, 'ck']} vs metadata {filled[i]} (condition {cores.loc[i, 'cond']}) -> metadata")
n_before = (cores.ck != "").sum(); cores["ck"] = np.where(filled != "", filled, cores.ck)
print(f"K numbers: {n_before} from the object, {(cores.ck != '').sum() - n_before} added from per_sample_metadata")
pseud = {t: f"xen_core_{i+1:03d}" for i, t in enumerate(sorted(cores.TMA_norm))}
cores["core_pseudonym"] = cores.TMA_norm.map(pseud)
# a K number can label two cores (two cores of one biopsy on different slides, or a metadata error when the diagnoses differ):
# keep cores separate, make the sample id unique with a core suffix and carry the K number in donor_id
cores["donor_id"] = [c if c else pseud[t] for c, t in zip(cores.ck, cores.TMA_norm)]
dupk = cores[cores.ck != ""].groupby("ck").size(); dupk = dupk[dupk > 1].index.tolist()
cores["sample_id"] = cores.donor_id
for k in dupk:
    idx = cores.index[cores.ck == k]
    for j, i in enumerate(idx): cores.loc[i, "sample_id"] = f"{k}_core{j+1}"
    print(f"  K number {k} labels {len(idx)} cores (runs {cores.loc[idx, 'run'].str.slice(0, 20).tolist()}, conditions {cores.loc[idx, 'cond'].tolist()}) -> sample ids {cores.loc[idx, 'sample_id'].tolist()}")
cores[["TMA_norm", "sample_id", "donor_id", "core_pseudonym", "ck", "run", "region"]].rename(columns={"TMA_norm": "TMA_ID"}).to_csv(f"{B}/outputs/private/sample_key_PRIVATE.csv", index=False)
print(f"cores {len(cores)}; with snRNA K number {(cores.ck != '').sum()}; pseudonymised {(cores.ck == '').sum()}; K numbers on two cores: {dupk}")
obs = obs.merge(cores[["TMA_norm", "sample_id", "donor_id", "core_pseudonym"]], on="TMA_norm", how="left").set_index(a.obs.index)
# ---- eGFR per core (corrected) -----------------------------------------------------------------------------------------
auth = pd.read_csv(AUTH); auth["k"] = auth["core"].astype(str).str.replace(r"\s+", " ", regex=True).str.strip()
meta = pd.read_csv(META); meta["k"] = meta["TMA_ID"].astype(str).str.replace(r"\s+", " ", regex=True).str.strip()
eg = cores[["TMA_norm"]].merge(auth[["k", "eGFR"]].rename(columns={"eGFR": "e_auth"}), left_on="TMA_norm", right_on="k", how="left").drop(columns="k")
eg = eg.merge(meta[["k", "egfr_at_biopsy", "ifta_pct", "classification", "jm_slope"]], left_on="TMA_norm", right_on="k", how="left").drop(columns="k")
eobj = pd.to_numeric(obs.groupby("TMA_norm", observed=True)["EGFR"].first(), errors="coerce")
eg["e_obj"] = eg.TMA_norm.map(eobj)
eg["eGFR_at_biopsy"] = eg.e_auth.fillna(eg.egfr_at_biopsy).fillna(eg.e_obj)
eg["eGFR_source"] = np.select([eg.e_auth.notna(), eg.egfr_at_biopsy.notna(), eg.e_obj.notna()], ["workbook", "per_sample_metadata", "object"], "missing")
print("eGFR source:", eg.eGFR_source.value_counts().to_dict())
egm = eg.set_index("TMA_norm")
for c in ["eGFR_at_biopsy", "eGFR_source", "ifta_pct", "classification", "jm_slope"]: obs[c] = obs.TMA_norm.map(egm[c])
# ---- groupings -----------------------------------------------------------------------------------------------------------
l2l1 = L2L1.get("l2l1") or L2L1.get("level2_to_level1") or {k: v for k, v in L2L1.items() if isinstance(v, str)}
obs["level2_refined"] = obs["level2_refined"].astype(str)
obs["mcfa_view"] = obs["level2_refined"].map(lambda x: l2l1.get(x, x))
print("views:", obs.mcfa_view.value_counts().to_dict())
keep = ["sample_id", "donor_id", "core_pseudonym", "sample", "region", "condition", "Disease", "eGFR_at_biopsy", "eGFR_source", "ifta_pct", "ifta_group",
        "egfr_group", "classification", "jm_slope", "segmentation_method", "level2_refined", "mcfa_view"]
a.obs = obs[keep].rename(columns={"sample": "xenium_run"})
a.obs["xenium_run"] = a.obs["xenium_run"].astype(str)
for c in a.obs.columns:
    if a.obs[c].dtype.name == "category": a.obs[c] = a.obs[c].astype(str)
# ---- pseudobulk ----------------------------------------------------------------------------------------------------------
readme = []
for grp in ["mcfa_view", "level2_refined"]:
    pb_all = dc.pp.pseudobulk(a, sample_col="sample_id", groups_col=grp, layer="counts", mode="sum")
    for tag, mc, mcnt in [("unfiltered", 0, 0), ("min10cells_1000counts", 10, 1000)]:
        pb = pb_all[(pb_all.obs.psbulk_cells >= mc) & (pb_all.obs.psbulk_counts >= mcnt)].copy()
        pb.obs = pb.obs.rename(columns={grp: "cell_type"})
        pb.obs["grouping"] = grp
        print(f"{grp} {tag}: {pb.shape}, psbulk_cells median {pb.obs.psbulk_cells.median():.0f}, psbulk_counts median {pb.obs.psbulk_counts.median():.0f}")
        X = pd.DataFrame(pb.X.toarray() if hasattr(pb.X, "toarray") else pb.X, index=pb.obs_names, columns=pb.var_names).astype(int)
        pb.write_h5ad(f"{B}/outputs/pseudobulk_{grp}_{tag}.h5ad")
        X.to_csv(f"{B}/outputs/pseudobulk_{grp}_{tag}_counts.csv"); pb.obs.to_csv(f"{B}/outputs/pseudobulk_{grp}_{tag}_obs.csv")
        readme.append(f"- `pseudobulk_{grp}_{tag}`: {pb.shape[0]} pseudobulks x {pb.shape[1]} genes (min_cells {mc}, min_counts {mcnt}); samples {pb.obs.sample_id.nunique()}, cell types {pb.obs.cell_type.nunique()}")
# per-sample cell-type counts (all cells) for convenience
ct = pd.crosstab(a.obs["sample_id"], a.obs["level2_refined"]); ct.to_csv(f"{B}/outputs/cells_per_sample_level2_refined.csv")
pd.Series(a.var_names, name="gene").to_csv(f"{B}/outputs/xenium_panel_genes.csv", index=False)
with open(f"{B}/outputs/README.md", "w") as f:
    f.write(f"""# Xenium cell type x sample pseudobulks (decoupler {dc.__version__}, `dc.pp.pseudobulk`, mode = sum of raw transcript counts)

Built {pd.Timestamp.today():%A, %d %B %Y} from the kidney-biopsy Xenium object (258,375 segmented cells, 468-gene panel, HGNC symbols,
5 Xenium runs, 119 TMA cores; one core = one biopsy = one donor). Counts are raw transcript counts per cell summed per
cell type and sample; no normalisation.

**Sample id.** `sample_id` is the snRNA-seq K number of the same biopsy where one exists ({(cores.ck != '').sum()} cores), otherwise a
pseudonym `xen_core_###` ({(cores.ck == '').sum()} cores). K numbers come from the per-sample metadata table (which corrects two cores that the
Xenium object mislabels) and `donor_id` equals `sample_id` ({"no K number labels two cores" if not dupk else "K numbers on two cores: " + ", ".join(dupk) + ", kept separate as <K>_core1/_core2"}).
`core_pseudonym` is the pseudonym for every core. No pathology numbers or dates.
Pseudobulks were computed with `dc.pp.pseudobulk` without filtering; `*_min10cells_1000counts` is the same object restricted to
pseudobulks with >= 10 cells and >= 1,000 transcripts (the former decoupler defaults).

**Groupings.** `mcfa_view` = the cell-type views of the multicellular factor model (PT = PT + aPT + dPT; Stromal = FIB + MES + VSM-P;
CNT_PC = CNT + PC; IC = IC-A + IC-B; EC, Immune, TAL, DCT, POD, PEC); `level2_refined` = the manuscript's refined labels.

**Files** (each as .h5ad with counts in X, and as *_counts.csv (pseudobulks x genes) + *_obs.csv):
{chr(10).join(readme)}

**obs columns.** `sample_id`, `donor_id`, `core_pseudonym`, `cell_type`, `grouping`, `psbulk_cells` (cells summed), `psbulk_counts` (transcripts summed),
`xenium_run` (one of 5 runs; batch), `region`, `condition` (CTRL/DN/FSGS/IgAN/MN/TIN), `Disease`, `eGFR_at_biopsy` + `eGFR_source`
(workbook / per-sample metadata / object; 2 cores missing), `ifta_pct`, `ifta_group`, `egfr_group`, `classification` and `jm_slope`
(progression class and joint-model eGFR slope where available), `segmentation_method`.
`cells_per_sample_level2_refined.csv` = cells per sample and refined label; `xenium_panel_genes.csv` = the 468 panel genes.

**Caveat.** Cell labels are the scANVI label transfer from the snRNA-seq reference trained with the two AKI samples present
(the AKI samples were later excluded from the snRNA-seq analyses; the Xenium reference was not retrained). Segmentation
spill-over between neighbouring cells is the imaging analogue of ambient RNA; the per-cell nearest-neighbour distances are
available on request.
""")
# ---- identifier scan -------------------------------------------------------------------------------------------------------
bad = re.compile(r"H\d{2}-\d{3,}|Rotterdam|patient_id|\d{4}-\d{2}-\d{2}")
hits = {}
for fn in sorted(os.listdir(f"{B}/outputs")):
    p = f"{B}/outputs/{fn}"
    if os.path.isdir(p) or fn.endswith(".h5ad"): continue
    txt = open(p, errors="replace").read(); m = bad.findall(txt)
    if m: hits[fn] = sorted(set(m))[:5]
for fn in sorted(os.listdir(f"{B}/outputs")):
    if fn.endswith(".h5ad"):
        o = ad.read_h5ad(f"{B}/outputs/{fn}").obs
        m = [c for c in o.columns if o[c].astype(str).str.contains(bad).any()]
        if m: hits[fn] = m
print("identifier scan:", "CLEAN" if not hits else hits)
assert not hits, "identifiers found in shared outputs"
with zipfile.ZipFile(f"{B}/outputs/xenium_pseudobulks_for_mcfa_20261026.zip", "w", zipfile.ZIP_DEFLATED) as z:
    for fn in sorted(os.listdir(f"{B}/outputs")):
        if fn.endswith((".h5ad", ".csv", ".md")) and not fn.startswith("."): z.write(f"{B}/outputs/{fn}", fn)
print("zip written")

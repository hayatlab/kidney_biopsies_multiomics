#!/usr/bin/env python
"""Step 04 — nuclear-only cell type x sample pseudobulks of the Xenium data for the MCFA team (ambient-RNA figure, reviewer 6).

Transcript rule: from each run's transcripts.parquet keep gene transcripts (is_gene) with qv >= 20 (the 10x default; with this filter the
per-cell totals reproduce the object's counts layer exactly, checked) assigned to a cell of the analysed object, with overlaps_nucleus == 1.
Summed per cell -> AnnData (cells x 468 genes) with the same obs as the 26 Oct delivery -> decoupler dc.pp.pseudobulk(mode='sum',
no filtering) by (sample_id, mcfa_view) and (sample_id, level2_refined). Sample ids are the SAME as in the 26 Oct zip (snRNA K/CK/KB
number where one exists, else xen_core_###) so the files join; `zenodo_sample_id` adds the names of the public release.
Each pseudobulk carries psbulk_cells, psbulk_counts (= nuclear transcripts) and total_counts_same_cells (all Q>=20 transcripts of the same
cells, i.e. the previous pseudobulk value) so the nuclear fraction can be computed.
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import sys, os, json, zipfile, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import pyarrow.parquet as pq, scipy.sparse as sp, anndata as ad, scanpy as sc, decoupler as dc
PREV = f"{KB_ANALYSIS}/mcfa_xenium_pseudobulk_20261026/outputs"
D = os.environ.get("KB_NUCLEAR_PSEUDOBULK_OUT", f"{OUT}/for_mcfa_nuclear_20261027"); os.makedirs(D, exist_ok=True)  # set the variable to re-run without overwriting the delivered files
print("decoupler", dc.__version__, "scanpy", sc.__version__)
a = sc.read_h5ad(OBJ); o = object_obs(); genes = a.var_names
prev_obs = pd.read_csv(f"{PREV}/pseudobulk_mcfa_view_unfiltered_obs.csv", index_col=0)
l2l1 = json.load(open(f"{KB_ANALYSIS}/vivien_takeover/outputs/01_factor_scoring_meta.json"))
l2l1 = l2l1.get("l2l1") or l2l1.get("level2_to_level1") or {k: v for k, v in l2l1.items() if isinstance(v, str)}
rows = []; stats = {}
for r in RUNS:
    rs = run_short(r); m = (o.run == r).values; cid = o.cell_id.astype(str)[m].values
    t = pq.read_table(f"{RAW}/output-{r}/transcripts.parquet", columns=["cell_id", "feature_name", "qv", "is_gene", "overlaps_nucleus"]).to_pandas()
    t = t[t.is_gene & (t.qv >= QV_MIN) & t.cell_id.isin(cid)]
    tot = t.groupby("cell_id").size(); nuc = t[t.overlaps_nucleus == 1]
    ci = pd.Categorical(nuc.cell_id, categories=cid); gi = pd.Categorical(nuc.feature_name, categories=genes)
    assert not pd.isna(gi.codes).any() and (gi.codes >= 0).all(), "gene outside panel"
    X = sp.coo_matrix((np.ones(len(nuc), dtype=np.int32), (ci.codes, gi.codes)), shape=(len(cid), len(genes))).tocsr()
    # check: total (nuclear + non-nuclear) reproduces the object's counts layer for this run
    tot_all = pd.Categorical(t.cell_id, categories=cid); gt = pd.Categorical(t.feature_name, categories=genes)
    XT = sp.coo_matrix((np.ones(len(t), dtype=np.int32), (tot_all.codes, gt.codes)), shape=X.shape).tocsr()
    assert (XT != a[m].layers["counts"]).nnz == 0, f"{rs}: transcripts at qv>={QV_MIN} != counts layer"
    rows.append((X, XT, np.where(m)[0]))
    stats[rs] = dict(cells=int(len(cid)), transcripts_q20=int(len(t)), nuclear=int(len(nuc)), nuclear_fraction=round(len(nuc) / len(t), 4))
    print(rs, stats[rs], flush=True)
idx = np.concatenate([r[2] for r in rows]); Xn = sp.vstack([r[0] for r in rows]).tocsr(); Xt = sp.vstack([r[1] for r in rows]).tocsr()
obs = o.iloc[idx].copy()
obs["sample_id"] = obs.prev_sample_id.astype(str); obs["zenodo_sample_id"] = o.iloc[idx].sample_id.astype(str).values
# the 26 Oct delivery split two-core K numbers into <K>_core1/_core2 — reuse exactly its sample_id via the core pseudonym
pk = pd.read_csv(PREV_KEY); pk["TMA_norm"] = pk.TMA_ID.map(norm); obs["sample_id"] = obs.TMA_norm.map(pk.set_index("TMA_norm").sample_id).astype(str)
assert obs.sample_id.notna().all() and set(obs.sample_id) == set(prev_obs.sample_id), "sample ids differ from the 26 Oct delivery"
obs["level2_refined"] = obs["level2_refined"].astype(str); obs["mcfa_view"] = obs["level2_refined"].map(lambda x: l2l1.get(x, x))
keepc = ["sample_id", "zenodo_sample_id", "xenium_core_id", "run", "condition", "level2_refined", "mcfa_view"]
obs = obs[keepc].copy(); obs["run"] = obs.run.astype(str)
for c in obs.columns:
    if obs[c].dtype.name == "category": obs[c] = obs[c].astype(str)
cells = ad.AnnData(X=Xn, obs=obs, var=pd.DataFrame(index=genes), layers={"nuclear": Xn, "total": Xt})
cells.obs["nuclear_counts"] = np.asarray(Xn.sum(1)).ravel(); cells.obs["total_counts"] = np.asarray(Xt.sum(1)).ravel()
print("cells", cells.shape, "nuclear fraction overall", round(cells.obs.nuclear_counts.sum() / cells.obs.total_counts.sum(), 4))
readme = []
for grp in ["mcfa_view", "level2_refined"]:
    pb = dc.pp.pseudobulk(cells, sample_col="sample_id", groups_col=grp, layer="nuclear", mode="sum")
    pbt = dc.pp.pseudobulk(cells, sample_col="sample_id", groups_col=grp, layer="total", mode="sum")
    assert (pb.obs_names == pbt.obs_names).all()
    pb.obs = pb.obs.rename(columns={grp: "cell_type"}); pb.obs["grouping"] = grp; pb.obs["transcript_set"] = "nuclear_only (overlaps_nucleus == 1)"
    pb.obs["total_counts_same_cells"] = pbt.obs.psbulk_counts.values; pb.obs["nuclear_fraction"] = pb.obs.psbulk_counts / pb.obs.total_counts_same_cells
    pb.layers["total_counts_same_cells"] = pbt.X
    # cross-check against the 26 Oct totals for the mcfa_view grouping
    if grp == "mcfa_view":
        pv = sc.read_h5ad(f"{PREV}/pseudobulk_mcfa_view_unfiltered.h5ad"); common = pb.obs_names.intersection(pv.obs_names)
        assert len(common) == pb.n_obs == pv.n_obs, (pb.n_obs, pv.n_obs, len(common))
        d = np.abs(np.asarray(pbt[common].X.todense() if sp.issparse(pbt.X) else pbt[common].X) - np.asarray(pv[common].X.todense() if sp.issparse(pv.X) else pv[common].X)).max()
        print("total-count pseudobulks vs 26 Oct delivery: max abs diff", d); assert d == 0
    X = pd.DataFrame(pb.X.toarray() if sp.issparse(pb.X) else pb.X, index=pb.obs_names, columns=pb.var_names).astype(int)
    pb.write_h5ad(f"{D}/pseudobulk_nuclear_{grp}_unfiltered.h5ad"); X.to_csv(f"{D}/pseudobulk_nuclear_{grp}_unfiltered_counts.csv"); pb.obs.to_csv(f"{D}/pseudobulk_nuclear_{grp}_unfiltered_obs.csv")
    print(grp, pb.shape, "median psbulk_cells", pb.obs.psbulk_cells.median(), "median nuclear counts", pb.obs.psbulk_counts.median())
    readme.append(f"- `pseudobulk_nuclear_{grp}_unfiltered`: {pb.shape[0]} pseudobulks x {pb.shape[1]} genes; samples {pb.obs.sample_id.nunique()}, cell types {pb.obs.cell_type.nunique()}")
json.dump(stats, open(f"{D}/per_run_transcript_stats.json", "w"), indent=1)
with open(f"{D}/README.md", "w") as f:
    f.write(f"""# Xenium NUCLEAR-only cell type x sample pseudobulks (decoupler {dc.__version__}, `dc.pp.pseudobulk`, mode = sum)

Built {pd.Timestamp.today():%A, %d %B %Y}. Same object, cells, labels, groupings and sample ids as the 26 Oct zip (`xenium_pseudobulks_for_mcfa_20261026.zip`),
but counting only transcripts that overlap the nucleus: from each run's `transcripts.parquet`, gene transcripts with Q-score >= 20 (the 10x default;
with this filter the per-cell totals reproduce the cell-feature matrix / the previous pseudobulks exactly — asserted in the script) assigned to an
analysed cell and with `overlaps_nucleus == 1`. Overall nuclear fraction of transcripts {cells.obs.nuclear_counts.sum() / cells.obs.total_counts.sum():.3f}
(per run in `per_run_transcript_stats.json`). No filtering on cells or counts.

**obs columns.** `sample_id` (identical to the 26 Oct files), `zenodo_sample_id` (sample name of the public release), `xenium_core_id`, `run`, `condition`,
`cell_type`, `grouping`, `psbulk_cells` (cells summed), `psbulk_counts` (= nuclear transcripts summed), `total_counts_same_cells` (all Q>=20 transcripts of
the same cells = the previous pseudobulk's psbulk_counts), `nuclear_fraction`. `layers['total_counts_same_cells']` in the h5ad holds the per-gene total
counts of the same cells (= the 26 Oct matrix) for a direct nuclear-vs-total comparison.

**Groupings.** `mcfa_view` = the factor model's views (PT = PT + aPT + dPT; Stromal = FIB + MES + VSM-P; CNT_PC = CNT + PC; IC = IC-A + IC-B; EC, Immune, TAL,
DCT, POD, PEC); `level2_refined` = the manuscript's refined labels.

**Files** (h5ad with nuclear counts in X; *_counts.csv pseudobulks x genes; *_obs.csv):
{chr(10).join(readme)}

**Caveats.** Nuclear overlap is 10x's per-transcript flag against the DAPI-based nucleus segmentation; it reduces but does not remove misassignment
(nuclei of neighbouring cells can touch, and nuclear-expanded cells have no boundary stain). Labels are the scANVI transfer from the snRNA-seq reference
trained with the two AKI samples present.
""")
bad = re.compile(r"H\d{2}-\d{3,}|\b\d{2}-\d{4,}\b|T\d{2}-\d{4,}|Rotterdam|patient_id|\d{4}-\d{2}-\d{2}|Combi |Appollo")
hits = {}
for fn in sorted(os.listdir(D)):
    p = f"{D}/{fn}"
    if fn.endswith(".h5ad"):
        ob = ad.read_h5ad(p).obs; mm = [c for c in ob.columns if ob[c].astype(str).str.contains(bad).any()]
    else: mm = sorted(set(bad.findall(open(p, errors="replace").read())))[:5]
    if mm: hits[fn] = mm
print("identifier scan:", "CLEAN" if not hits else hits); assert not hits
with zipfile.ZipFile(f"{D}/xenium_nuclear_pseudobulks_for_mcfa_20261027.zip", "w", zipfile.ZIP_DEFLATED) as z:
    for fn in sorted(os.listdir(D)):
        if fn.endswith((".h5ad", ".csv", ".md", ".json")): z.write(f"{D}/{fn}", fn)
print("DONE zip", os.path.getsize(f"{D}/xenium_nuclear_pseudobulks_for_mcfa_20261027.zip") / 1e6, "MB")

#!/usr/bin/env python
"""Generalized Leiden-resolution sweep for a single level2 cell type (for manual annotation).

Mirrors pec_resubcluster_20260717/pec_leiden_sweep.py but:
  - parameterized by --celltype (a level2 label), subsetting the parent object directly;
  - crescent/ULM/custom-gene/provenance machinery REMOVED (PEC-specific);
  - a per-lineage CANONICAL marker panel added (dotplot per resolution) to aid annotation.

Embedding is rebuilt ON THE SUBSET (the incoming global X_harmony under-resolves within-type
structure), shared across all resolutions so partitions are directly comparable:
    log1p -> HVG(2000, batch=patient_ID) -> scale(10) -> PCA(30, seed42)
          -> harmony(patient_ID) -> neighbors(k=15) -> UMAP   [computed ONCE]
          -> leiden at each resolution (igraph, n_iterations=2, seed42)   [only this varies]

Per resolution -> <outroot>/<celltype>/res_<r>/:
  umap_res<r>.{pdf,png}              UMAP coloured by cluster
  markers_res<r>.csv                full Wilcoxon table (top 100/cluster)
  markers_dotplot_res<r>.{pdf,png}  top-5 Wilcoxon markers per cluster
  canonical_panel_res<r>.{pdf,png}  curated per-lineage marker dotplot (grouped by putative subtype)
  cluster_summary_res<r>.csv        n + top markers per cluster
Human, HGNC symbols. Seed 42 throughout. Nothing excluded -- exclusion is the user's call.
"""
from kbpaths import KB_ANALYSIS, KB_DATA  # repository path config (config/kbpaths.py)
import os, argparse, warnings
import numpy as np, pandas as pd, scanpy as sc
import scanpy.external as sce
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=FutureWarning)
sc.settings.verbosity = 1; sc.settings.autoshow = False
np.random.seed(42); SEED = 42

PARENT_DEFAULT = f"{KB_DATA}/charlotte_new_data/adata_20260717_with_clinical.h5ad"
OUTROOT_DEFAULT = f"{KB_ANALYSIS}/subcluster_level2_20260728"
RESOLUTIONS = [0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5]
BATCH = "patient_ID"; N_HVG = 2000; N_PCS = 30; CELLTYPE_COL = "level2"

# ---- per-lineage canonical panels (grouped by putative subtype/state); filtered to present genes ----
PANELS = {
 "PT":   {"PT-S1":["SLC5A12","SLC5A2","LRP2"],"PT-S2":["SLC22A6","SLC22A8","SLC13A3"],
          "PT-S3":["SLC7A13","SLC22A7","MIOX"],"general":["CUBN","GATM","SLC34A1"],
          "injury":["HAVCR1","VCAM1","DCDC2"]},
 "aPT":  {"aPT/injury":["VCAM1","HAVCR1","PROM1","DCDC2"],"dediff":["SOX9","VIM","CD24"],
          "PT-id":["LRP2","CUBN","SLC34A1"],"inflam":["SPP1","ITGB6","CST3"]},
 "dPT":  {"degen/stress":["MT1G","MT1H","MT2A"],"injury":["HAVCR1","VCAM1"],
          "PT-id":["LRP2","CUBN","SLC34A1"],"senescence":["CDKN1A","CDKN2A"]},
 "TAL":  {"TAL-core":["SLC12A1","UMOD","CASR"],"cTAL":["CLDN16","CLDN14","KCNJ1"],
          "mTAL":["CLDN10","ENOX1"],"MD":["NOS1","PAPPA2"],"injury":["PROM1","ITGB6","VCAM1"]},
 "DCT":  {"DCT1":["SLC12A3","TRPM6","WNK1"],"DCT2":["SLC8A1","CALB1","TRPV5"],
          "core":["KLHL3","EGF"],"injury":["VCAM1","ITGB6"]},
 "CNT":  {"CNT":["SLC8A1","CALB1","TRPV5","SCNN1G"],"CNT-PC":["AQP2","GATA3"]},
 "PC":   {"PC-core":["AQP2","AQP3","GATA3","FXYD4"],"cortical":["SCNN1G","SCNN1B"],
          "IMCD":["PROM1","SLC14A2"],"transitional":["SLC4A1","SLC26A4"]},
 "IC-A": {"IC-A":["SLC4A1","ATP6V0A4","ATP6V1G3","AQP6"],"core":["SLC26A7","KIT"],"injury":["VCAM1"]},
 "IC-B": {"IC-B":["SLC26A4","SLC4A9","INSRR"],"core":["SLC35F3","TMEM61"]},
 "EC":   {"GEC-glom":["EHD3","PLVAP","SOST"],"AEA/DVR":["SERPINE2","SLC14A1"],
          "lymphatic":["PROX1","PDPN","CCL21"],"core":["PECAM1","EMCN","FLT1","KDR"]},
 "POD":  {"POD-core":["NPHS1","NPHS2","PODXL","PTPRO"],"TF":["WT1","TCF21"],
          "receptor":["PLA2R1","THSD7A"],"injury":["VIM","CD44"]},
 "PEC":  {"PEC-core":["CFH","CLDN1","ALDH1A2"],"activated":["VCAM1","CD44","ITGB1"],
          "crescent":["MMP7","PDLIM7"],"id":["CAV2","PAX8"]},
 "FIB":  {"FIB-core":["PDGFRA","PDGFRB","COL1A1","DCN","LUM"],"myofib":["ACTA2","POSTN","TAGLN"],
          "medullary":["COL1A2","FN1"],"cortical":["MEG3","GREM1"]},
 "VSM-P":{"VSMC":["ACTA2","MYH11","TAGLN","CNN1"],"pericyte":["PDGFRB","RGS5","NOTCH3"],
          "JG/renin":["REN","AKR1B1"],"core":["MYL9"]},
 "MES":  {"MES-core":["PDGFRB","GATA3","ITGA8","POSTN"],"contractile":["ACTA2","MYH11"],
          "JG":["REN"],"id":["NOTCH3","RGS5"]},
 "Immune":{"T":["CD3D","CD3E","CD8A","IL7R"],"NK":["NKG7","GNLY","KLRD1"],"B":["CD79A","MS4A1"],
          "Plasma":["MZB1","IGHG1"],"Mono/Macro":["CD14","LYZ","CD68","C1QA","C1QB","CD163"],
          "DC":["FCER1A","CLEC9A","LILRA4"],"Mast":["CPA3","TPSAB1"],"Neutrophil":["FCGR3B","CSF3R"]},
}

def save(fig, stem):
    fig.patch.set_facecolor("white")
    for ax in fig.get_axes(): ax.set_facecolor("white")
    fig.savefig(stem + ".pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(stem + ".png", dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--celltype", required=True)
    ap.add_argument("--parent", default=PARENT_DEFAULT)
    ap.add_argument("--outroot", default=OUTROOT_DEFAULT)
    args = ap.parse_args()
    CT = args.celltype
    ROOT = f"{args.outroot}/{CT.replace('/','_')}"; os.makedirs(ROOT, exist_ok=True)

    # -------- subset from parent --------
    print(f"[{CT}] reading parent backed: {args.parent}")
    par = sc.read_h5ad(args.parent, backed="r")
    mask = (par.obs[CELLTYPE_COL].astype(str).values == CT)
    print(f"[{CT}] {CELLTYPE_COL}=={CT}: {int(mask.sum())} cells")
    assert mask.sum() >= 100, f"too few cells for {CT}"
    adata = par[mask].to_memory(); del par
    print(f"[{CT}] loaded {adata.shape}  obsm={list(adata.obsm)}  X.max={float(adata.X.max()):.1f}")
    if "X_harmony" in adata.obsm: adata.obsm["X_harmony_global"] = adata.obsm.pop("X_harmony")
    adata.layers["libnorm"] = adata.X.copy()
    sc.pp.log1p(adata); adata.layers["lognorm"] = adata.X.copy()
    print(f"[{CT}] log1p done; X.max={float(adata.X.max()):.2f}; {BATCH}: {adata.obs[BATCH].nunique()} levels")

    # -------- subset-specific embedding (shared across resolutions) --------
    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat", batch_key=BATCH)
    except Exception as e:
        print(f"[{CT}] batch-aware HVG failed ({e}); falling back to non-batch HVG")
        sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat")
    n_hvg = int(adata.var.highly_variable.sum()); print(f"[{CT}] HVGs: {n_hvg}")
    emb = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(emb, max_value=10)
    sc.tl.pca(emb, n_comps=N_PCS, svd_solver="arpack", random_state=SEED)
    vr = emb.uns["pca"]["variance_ratio"]
    sce.pp.harmony_integrate(emb, key=BATCH, basis="X_pca", adjusted_basis="X_pca_harmony",
                             random_state=SEED)
    adata.obsm["X_pca"] = emb.obsm["X_pca"]; adata.obsm["X_pca_harmony"] = emb.obsm["X_pca_harmony"]; del emb
    fig, ax = plt.subplots(figsize=(5.2,3.6))
    ax.plot(np.arange(1,N_PCS+1), vr, "o-", ms=4); ax.set_xlabel("PC"); ax.set_ylabel("variance ratio")
    ax.set_title(f"{CT} subset PCA ({n_hvg} HVGs, {N_PCS} PCs, cum={vr.sum():.2f})")
    fig.tight_layout(); save(fig, f"{ROOT}/pca_variance_ratio")
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=15, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    print(f"[{CT}] embedding done (shared across resolutions)")

    # -------- canonical panel (filter present genes) --------
    panel = PANELS.get(CT, {})
    panel = {grp:[g for g in gs if g in adata.var_names] for grp,gs in panel.items()}
    panel = {grp:gs for grp,gs in panel.items() if gs}
    absent = [g for grp,gs in PANELS.get(CT,{}).items() for g in gs if g not in adata.var_names]
    if absent: print(f"[{CT}] canonical panel genes absent from object: {sorted(set(absent))}")

    # -------- sweep --------
    overview = []
    for r in RESOLUTIONS:
        key = f"leiden_{r}"; tag = f"res{r}"; out = f"{ROOT}/res_{r}"; os.makedirs(out, exist_ok=True)
        sc.tl.leiden(adata, resolution=r, key_added=key, flavor="igraph",
                     n_iterations=2, directed=False, random_state=SEED)
        cats = list(adata.obs[key].cat.categories)
        counts = adata.obs[key].value_counts().sort_index()
        print(f"\n{'='*66}\n[{CT}] leiden {r} -> {len(cats)} clusters: {counts.to_dict()}")
        overview.append({"resolution": r, "n_clusters": len(cats),
                         "min_cluster": int(counts.min()), "max_cluster": int(counts.max()),
                         "sizes": ";".join(f"{c}:{n}" for c,n in counts.items())})
        # UMAP by cluster
        fig, ax = plt.subplots(figsize=(6.4,5.6))
        sc.pl.umap(adata, color=key, ax=ax, show=False, frameon=False, legend_loc="on data",
                   legend_fontsize=9, legend_fontoutline=2, title=f"{CT} leiden {r} ({len(cats)} clusters)")
        fig.tight_layout(); save(fig, f"{out}/umap_{tag}")
        # Wilcoxon markers
        sc.tl.rank_genes_groups(adata, key, method="wilcoxon", n_genes=100, key_added=f"rgg_{key}")
        mk = sc.get.rank_genes_groups_df(adata, group=None, key=f"rgg_{key}")
        mk.to_csv(f"{out}/markers_{tag}.csv", index=False)
        top5 = {c: mk[mk.group==c].sort_values("scores", ascending=False).names.head(5).tolist() for c in cats}
        ndot = len({g for v in top5.values() for g in v})
        dp = sc.pl.dotplot(adata, var_names=top5, groupby=key, standard_scale="var", return_fig=True,
                           figsize=(max(10,0.30*ndot), 0.55*len(cats)+2.2))
        dp.add_totals().style(cmap="Reds"); dp.make_figure()
        dp.fig.suptitle(f"{CT} top-5 Wilcoxon markers/cluster - leiden {r}", fontsize=12, y=1.02)
        save(dp.fig, f"{out}/markers_dotplot_{tag}")
        # canonical panel dotplot
        if panel:
            ng = sum(len(v) for v in panel.values())
            dp = sc.pl.dotplot(adata, var_names=panel, groupby=key, standard_scale="var", return_fig=True,
                               figsize=(0.34*ng+3, 0.55*len(cats)+2.4))
            dp.add_totals().style(cmap="Reds"); dp.make_figure()
            dp.fig.suptitle(f"{CT} canonical marker panel - leiden {r}", fontsize=11, y=1.03)
            save(dp.fig, f"{out}/canonical_panel_{tag}")
        # cluster summary
        summ = pd.DataFrame([{"cluster":c, "n":int((adata.obs[key]==c).sum()),
                              "top_markers":",".join(top5[c])} for c in cats])
        summ.to_csv(f"{out}/cluster_summary_{tag}.csv", index=False)
        print(summ.to_string(index=False)); print(f"[saved] {out}")

    ov = pd.DataFrame(overview); ov.to_csv(f"{ROOT}/resolution_overview.csv", index=False)
    print(f"\n[{CT}] RESOLUTION OVERVIEW\n{ov.to_string(index=False)}")
    # sweep overview figure
    ncol=5; nrow=int(np.ceil(len(RESOLUTIONS)/ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(5.6*ncol, 5.2*nrow)); axs=axs.ravel()
    for ax, r in zip(axs, RESOLUTIONS):
        sc.pl.umap(adata, color=f"leiden_{r}", ax=ax, show=False, frameon=False, legend_loc="on data",
                   legend_fontsize=8, legend_fontoutline=2,
                   title=f"leiden {r} ({adata.obs[f'leiden_{r}'].nunique()} cl)")
    for ax in axs[len(RESOLUTIONS):]: ax.axis("off")
    fig.suptitle(f"{CT} leiden sweep - subset HVG/PCA/harmony({BATCH}) (n={adata.n_obs}, seed {SEED})", fontsize=15)
    fig.tight_layout(); save(fig, f"{ROOT}/resolution_sweep_overview")
    adata.write_h5ad(f"{ROOT}/{CT.replace('/','_')}_leiden_sweep.h5ad")
    print(f"\n[{CT}] done -> {ROOT}")

if __name__ == "__main__":
    main()

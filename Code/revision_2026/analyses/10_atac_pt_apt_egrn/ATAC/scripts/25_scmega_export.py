#!/usr/bin/env python3
"""Export scMEGA inputs (PT -> aPT only; dPT branch excluded) in 10x-HDF5 format for Seurat.

RNA : adata_20260717.h5ad, level2 PT + aPT cells of the GLUE co-embedding; genes in >= 10 cells.
      X is ambient-corrected (non-integral) -- stored as float32; Seurat normalises downstream.
ATAC: peak_mat2.h5ad, ATAC clusters pt_C3 + pt_C6 (recommended_use == 'use'); standard
      chromosomes; peaks accessible in >= 1% of these cells. Names 'chr-start-end' for Signac.
Also: 30-dim GLUE latent (the co-embedding scMEGA pairs cells on) and per-cell metadata incl.
      pseudotime_joint / pseudotime_via. SMOKE copy = 5% stratified by modality x label.
Human hg38 | HGNC symbols | seed 42
"""
from kbpaths import KB_ANALYSIS, KB_DATA  # repository path config (config/kbpaths.py)
import anndata as ad, h5py, numpy as np, pandas as pd, scipy.sparse as sp, os
SEED=42; rng=np.random.default_rng(SEED)
A=f"{KB_ANALYSIS}/aki_c5_exclusion_20261022/atac_rerun/ATAC"; O=f"{A}/scmega/inputs"
RNA=f"{KB_DATA}/charlotte_new_data/adata_20260717.h5ad"
PEAK=f"{KB_DATA}/New_ATAC_data/peak_mat2.h5ad"
LAT=[str(i) for i in range(30)]; STD={f"chr{i}" for i in list(range(1,23))+["X","Y"]}

H=pd.read_csv(f"{A}/tables/HANDOFF_joint_pseudotime.csv",low_memory=False)
V=pd.read_csv(f"{A}/tables/via_joint_pseudotime.csv",usecols=["cellid","pseudotime_via"])
H=H.merge(V,on="cellid",how="left")
keep=((H.modality=="RNA")&H.level2_rna.isin(["PT","aPT"]))|((H.modality=="ATAC")&(H.recommended_use=="use"))
M=H[keep].copy()
assert not M["sample"].astype(str).isin(["K95", "K96", "K98"]).any(), "AKI cells in the hand-off table"
M["label"]=np.where(M.modality=="RNA",M.level2_rna,M.atac_cluster)
print(M.groupby(["modality","label"]).size().to_string(),flush=True)
G=pd.read_csv(f"{A}/glue_v3/glue_for_pseudotime.csv",usecols=["cellid"]+LAT,low_memory=False).set_index("cellid")
assert M.cellid.isin(G.index).all()

def write10x(path,X_cells_by_feat,barcodes,features,ftype):
    X=sp.csr_matrix(X_cells_by_feat)            # CSR (cells x features) == CSC (features x cells)
    with h5py.File(path,"w") as f:
        g=f.create_group("matrix")
        g.create_dataset("barcodes",data=np.array(barcodes,dtype="S"))
        g.create_dataset("data",data=X.data,compression="gzip")
        g.create_dataset("indices",data=X.indices.astype(np.int64),compression="gzip")
        g.create_dataset("indptr",data=X.indptr.astype(np.int64))
        g.create_dataset("shape",data=np.array([X.shape[1],X.shape[0]],dtype=np.int64))
        ft=g.create_group("features")
        ft.create_dataset("id",data=np.array(features,dtype="S")); ft.create_dataset("name",data=np.array(features,dtype="S"))
        ft.create_dataset("feature_type",data=np.array([ftype]*len(features),dtype="S"))
        ft.create_dataset("genome",data=np.array(["hg38"]*len(features),dtype="S"))
        ft.create_dataset("_all_tag_keys",data=np.array(["genome"],dtype="S"))

def export(sub,tag):
    os.makedirs(f"{O}{tag}",exist_ok=True); o=f"{O}{tag}"
    r=sub[sub.modality=="RNA"].cellid.tolist(); a_=sub[sub.modality=="ATAC"].cellid.tolist()
    R=ad.read_h5ad(RNA,backed="r")[r].to_memory()
    X=sp.csr_matrix(R.X,dtype=np.float32); ng=np.asarray((X>0).sum(0)).ravel(); gk=ng>=10
    write10x(f"{o}/rna_counts.h5",X[:,gk],r,R.var_names[gk].tolist(),"Gene Expression")
    print(f"[{tag or 'full'}] RNA {X.shape[0]:,} cells x {int(gk.sum()):,}/{X.shape[1]:,} genes (>=10 cells)",flush=True)
    del R,X
    P=ad.read_h5ad(PEAK,backed="r"); pn=P.var_names.to_numpy()
    Pm=P[a_].to_memory(); Y=sp.csr_matrix(Pm.X).astype(np.int32)
    ch=pd.Series(pn).str.split("_").str[0].isin(STD).to_numpy()
    det=np.asarray((Y>0).mean(0)).ravel(); pk=ch&(det>=0.01)
    names=pd.Series(pn[pk]).str.replace("_","-").tolist()
    write10x(f"{o}/atac_peaks.h5",Y[:,pk],a_,names,"Peaks")
    print(f"[{tag or 'full'}] ATAC {Y.shape[0]:,} cells x {int(pk.sum()):,}/{len(pn):,} peaks "
          f"(standard chr: {int(ch.sum()):,}; >=1% of cells: {int(pk.sum()):,})",flush=True)
    del P,Pm,Y
    G.loc[sub.cellid].to_csv(f"{o}/glue_latent.csv")
    sub[["cellid","modality","label","sample","pseudotime_joint","pseudotime_via","pseudotime_rna_anchored",
         "Factor1_rna","joint_dc1","joint_dc2"]].to_csv(f"{o}/metadata.csv",index=False)
    print(f"[{tag or 'full'}] wrote {o}",flush=True)

sm=(M.groupby(["modality","label"],group_keys=False).apply(lambda d:d.sample(frac=0.05,random_state=SEED)))
export(sm,"_smoke")
export(M,"")
print("done",flush=True)

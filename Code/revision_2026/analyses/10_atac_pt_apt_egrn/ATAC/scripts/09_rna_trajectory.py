#!/usr/bin/env python3
"""PT -> dPT -> aPT trajectory on the RNA half of the GLUE co-embedding.

This is the native version of everything we did on ATAC gene scores:
  * real measured expression, so no gene-length artifact
  * Factor 1 scored RNA-on-RNA (the loadings are RNA-derived) -- no hybrid caveat
  * dPT is resolved, so the geometry decides whether it is an intermediate or a branch

Source RNA : charlotte_new_data/adata_20260717.h5ad (363,656 x 33,000) -- the object GLUE
             was built from (100% barcode overlap with the GLUE RNA rows).
             X is non-integral with variable library size (ambient-corrected counts);
             handled as counts (normalize_total 1e4 + log1p). Provenance of the correction
             not established here.
Axis       : DiffMap DC1 on the GLUE latent, oriented so + = toward aPT.
Species human hg38 | HGNC symbols | seed 42
"""
from kbpaths import KB_ANALYSIS, KB_DATA  # repository path config (config/kbpaths.py)
import scanpy as sc, anndata as ad, numpy as np, pandas as pd, os, phate
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.stats import rankdata, spearmanr, t as tdist, mannwhitneyu
from statsmodels.stats.multitest import multipletests

SEED=42; np.random.seed(SEED); sc.settings.n_jobs=16
A=f"{KB_ANALYSIS}/aki_c5_exclusion_20261022/atac_rerun/ATAC"
RNA=f"{KB_DATA}/charlotte_new_data/adata_20260717.h5ad"
LOAD=f"{KB_ANALYSIS}/loadings~hvg~5.csv"
VER="glue_v3"; K=30; MIN_DET=0.10; MIN_CELLS_DONOR=50
os.makedirs(f"{A}/figures",exist_ok=True); os.makedirs(f"{A}/tables",exist_ok=True)
print("scanpy",sc.__version__,"| anndata",ad.__version__,"| phate",phate.__version__,flush=True)

g=pd.read_csv(f"{A}/{VER}/glue_for_pseudotime.csv",low_memory=False)
g=g[~g.cellid.astype(str).str.split("_").str[0].isin(["K95", "K96", "K98"])]  # AKI exclusion
g=g[g.modality=="RNA"].set_index("cellid")
LAT=[str(i) for i in range(30)]
print(f"GLUE RNA rows: {len(g):,} | level2: {g.level2.value_counts().to_dict()}",flush=True)

a=ad.read_h5ad(RNA,backed="r")
missing=g.index.difference(a.obs_names)
assert len(missing)==0, f"{len(missing)} GLUE RNA cells absent from the RNA object"
a=a[g.index].to_memory()
print(f"loaded {a.shape[0]:,} x {a.shape[1]:,}",flush=True)
a.obsm["X_glue"]=g[LAT].to_numpy(dtype=np.float32)
a.obs["level2_glue"]=pd.Categorical(g["level2"].values)
agree=(a.obs["level2"].astype(str).values==g["level2"].astype(str).values).mean()
print(f"level2 agreement object vs GLUE file: {100*agree:.2f}%",flush=True)

# ── Factor 1, scored natively on RNA ──────────────────────────────────────────
L=pd.read_csv(LOAD); W=L[(L.Factor=="Factor1")&(L.view=="PT")].set_index("variable")["Loading"]
W=W[~W.index.duplicated()]
b=a.copy(); b.X=b.X.astype(np.float32)
sc.pp.normalize_total(b,target_sum=1e4); sc.pp.log1p(b)
gl=[x for x in W.index if x in b.var_names]
print(f"Factor1 PT loadings present in RNA: {len(gl):,}/{len(W):,} "
      f"(missing {len(W)-len(gl):,}, {100*(len(W)-len(gl))/len(W):.1f}%)",flush=True)
Mf=b[:,gl].X; Mf=np.asarray(Mf.todense() if hasattr(Mf,"todense") else Mf,dtype=np.float32)
a.obs["Factor1"]=Mf@W.loc[gl].to_numpy(dtype=np.float32); del Mf
f1=a.obs["Factor1"].values
print("Factor1 by level2:\n"+a.obs.groupby("level2_glue",observed=True)["Factor1"]
      .agg(['size','mean','std']).round(2).to_string(),flush=True)

# ── embeddings on the GLUE latent ─────────────────────────────────────────────
sc.pp.neighbors(a,use_rep="X_glue",n_neighbors=K,random_state=SEED)
sc.tl.umap(a,random_state=SEED); sc.tl.diffmap(a,n_comps=15)
dc1=a.obsm["X_diffmap"][:,1].astype(float)
if spearmanr(dc1,f1).statistic<0: dc1=-dc1; a.obsm["X_diffmap"][:,1]*=-1
a.obs["DC1"]=dc1
a.obsm["X_phate"]=phate.PHATE(n_components=2,knn=K,random_state=SEED,n_jobs=16,
                              verbose=0).fit_transform(a.obsm["X_glue"])
lev=a.obs["level2_glue"].astype(str).values
isa=lev=="aPT"; isp=lev=="PT"
auc=mannwhitneyu(dc1[isa],dc1[isp],alternative="greater").statistic/(isa.sum()*isp.sum())
print(f"\nDC1: rho vs Factor1 {spearmanr(dc1,f1).statistic:+.3f} | AUC(aPT>PT) {auc:.3f}",flush=True)
print("DC1 mean by level2:\n"+pd.Series(dc1).groupby(lev).mean().round(4).to_string(),flush=True)
# is dPT intermediate or a branch?
print("\ndPT position test (DC1 percentile within the PT->aPT range):",flush=True)
lo,hi=np.median(dc1[isp]),np.median(dc1[isa])
md=np.median(dc1[lev=="dPT"])
print(f"   median DC1  PT {lo:.4f} | dPT {md:.4f} | aPT {hi:.4f}"
      f"  -> dPT sits at {100*(md-lo)/(hi-lo):.0f}% of the way",flush=True)

# ── genes vs DC1 ──────────────────────────────────────────────────────────────
det=np.asarray((b.X>0).mean(0)).ravel(); keep=det>=MIN_DET
genes=b.var_names[keep].to_numpy()
M=b[:,keep].X; M=np.asarray(M.todense() if hasattr(M,"todense") else M,dtype=np.float32); del b
n=M.shape[0]; print(f"\ngenes detected >= {MIN_DET:.0%}: {keep.sum():,}/{len(keep):,}",flush=True)
Rg=np.apply_along_axis(rankdata,0,M).astype(np.float32); ra=rankdata(dc1).astype(np.float32)
Rg-=Rg.mean(0); ra=ra-ra.mean()
rho=(Rg*ra[:,None]).sum(0)/(np.sqrt((Rg**2).sum(0))*np.sqrt((ra**2).sum())+1e-12); del Rg
rho=np.clip(rho,-0.999999,0.999999)
tt=rho*np.sqrt((n-2)/(1-rho**2)); pv=2*tdist.sf(np.abs(tt),df=n-2)
fdr=multipletests(pv,method="fdr_bh")[1]; z=np.arctanh(rho); se=1/np.sqrt(n-3)
res=pd.DataFrame({"gene":genes,"rho":rho,"ci95_lo":np.tanh(z-1.96*se),"ci95_hi":np.tanh(z+1.96*se),
                  "pval":pv,"fdr":fdr,"pct_detected":det[keep]}).sort_values("rho",ascending=False)
gi={x:i for i,x in enumerate(genes)}
top=res.reindex(res.rho.abs().sort_values(ascending=False).index).head(2000)
idx=np.array([gi[x] for x in top.gene]); don=a.obs["sample"].astype(str).values
ag=np.zeros(len(idx)); nd=0
for d in pd.unique(don):
    m=don==d
    if m.sum()<MIN_CELLS_DONOR: continue
    nd+=1
    Rd=np.apply_along_axis(rankdata,0,M[np.ix_(m,idx)]).astype(np.float32); pp=rankdata(dc1[m]).astype(np.float32)
    Rd-=Rd.mean(0); pp=pp-pp.mean()
    rd=(Rd*pp[:,None]).sum(0)/(np.sqrt((Rd**2).sum(0))*np.sqrt((pp**2).sum())+1e-12)
    ag+=np.sign(rd)==np.sign(top.rho.values)
res=res.merge(top.assign(donor_sign_agree=ag/nd)[["gene","donor_sign_agree"]],on="gene",how="left")
res.to_csv(f"{A}/tables/rna_gene_correlations.csv",index=False)
res[["gene","rho"]].to_csv(f"{A}/tables/rna_ranked_genes.rnk",sep="\t",index=False,header=False)
print(f"donors {nd} | median sign agreement {np.median(ag/nd):.2f} | FDR<0.05 {(res.fdr<0.05).sum():,}",flush=True)
print("\nTOP 15 UP toward aPT:\n"+res.head(15)[["gene","rho","fdr","donor_sign_agree"]].round(4).to_string(index=False),flush=True)
print("\nTOP 15 DOWN toward aPT:\n"+res.tail(15).iloc[::-1][["gene","rho","fdr","donor_sign_agree"]].round(4).to_string(index=False),flush=True)
try:
    d=pd.read_csv(f"{A}/tables/dc1_gene_correlations.csv")[["gene","rho"]].rename(columns={"rho":"rho_atac"})
    mm=res[["gene","rho"]].merge(d,on="gene")
    print(f"\nRNA vs ATAC gene-score ranking: Spearman {spearmanr(mm.rho,mm.rho_atac).statistic:.3f} on {len(mm):,} shared genes",flush=True)
except Exception as e: print("ATAC comparison skipped:",e,flush=True)

out=pd.DataFrame({"cellid":a.obs_names,"level2":lev,"sample":don,"Factor1":f1,"DC1":dc1,
                  "umap1":a.obsm["X_umap"][:,0],"umap2":a.obsm["X_umap"][:,1],
                  "phate1":a.obsm["X_phate"][:,0],"phate2":a.obsm["X_phate"][:,1],
                  "dc2":a.obsm["X_diffmap"][:,2]})
out.to_csv(f"{A}/tables/rna_trajectory_percell.csv",index=False)

CC={"PT":"#55A868","dPT":"#DD8452","aPT":"#937860"}
fig,axes=plt.subplots(2,3,figsize=(19,11.5))
o=np.random.default_rng(SEED).permutation(len(out))
n1=TwoSlopeNorm(vcenter=float(np.median(f1)),vmin=float(np.percentile(f1,1)),vmax=float(np.percentile(f1,99)))
for j,(x,y,nm) in enumerate([("umap1","umap2","UMAP"),("phate1","phate2","PHATE"),("DC1","dc2","DiffMap DC1/DC2")]):
    ax=axes[0,j]
    for cl in ["PT","dPT","aPT"]:
        m=out.level2==cl
        ax.scatter(out.loc[m,x],out.loc[m,y],s=.7,lw=0,alpha=.45,c=CC[cl],label=f"{cl} ({m.sum():,})",rasterized=True)
    ax.set_title(f"{nm} — level2",fontweight="bold")
    if j==0:
        lg=ax.legend(markerscale=14,fontsize=9)
        for h in lg.legend_handles: h.set_alpha(1)
    ax2=axes[1,j]
    s=ax2.scatter(out[x].values[o],out[y].values[o],s=.7,lw=0,c=f1[o],cmap="RdBu_r",norm=n1,rasterized=True)
    ax2.set_title(f"{nm} — Factor 1 (RNA-native)",fontweight="bold"); plt.colorbar(s,ax=ax2,fraction=.046)
    for aa in (ax,ax2): aa.set_xticks([]); aa.set_yticks([])
fig.suptitle(f"snRNA PT → dPT → aPT on the GLUE latent · {a.shape[0]:,} nuclei · "
             f"DC1 AUC(aPT>PT)={auc:.3f}",fontsize=15,fontweight="bold")
fig.tight_layout(rect=[0,0,1,.96])
fig.savefig(f"{A}/figures/rna_trajectory.png",dpi=140,bbox_inches="tight")
fig.savefig(f"{A}/figures/rna_trajectory.pdf",bbox_inches="tight")
print("\ndone",flush=True)

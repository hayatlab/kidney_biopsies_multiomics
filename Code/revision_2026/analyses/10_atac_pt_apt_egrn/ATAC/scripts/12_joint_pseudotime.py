#!/usr/bin/env python3
"""ONE pseudotime on the combined GLUE embedding (RNA + ATAC) -> hand-off table + split figure.

Design: pseudotime is computed once, on the combined embedding. Leon analyses the RNA half;
a colleague analyses the ATAC half. Both read the same axis from HANDOFF_joint_pseudotime.csv.

Cells      : glue_v3 = all RNA (PT/dPT/aPT) + all ATAC (pt_C1/2/3/4/6; pt_C5 removed upstream)
Graph      : kNN (K=30) on the 30-dim GLUE latent, all cells together
Pseudotime : rank-normalised DiffMap DC1 of the joint graph in [0,1], oriented so it rises with
             RNA Factor 1 (toward aPT). PRIMARY axis.
DPT        : scanpy DPT on the same graph, rooted at the RNA PT cell nearest the RNA-PT centroid.
             Reported alongside; on these data it compresses most cells into a narrow range.
RNA-anchored pseudotime (ATAC cells only): mean RNA-only pseudotime of each ATAC cell's 25
             nearest RNA neighbours in the GLUE latent. Use this for any "chromatin before
             expression" claim, because RNA alone defines it.
Human hg38 | HGNC symbols | seed 42.   SMOKE=1 -> 4% stratified subsample, outputs tagged _smoke.
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import scanpy as sc, anndata as ad, numpy as np, pandas as pd, os, json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.neighbors import NearestNeighbors
from scipy.stats import spearmanr, rankdata, mannwhitneyu

SEED=42; np.random.seed(SEED); sc.settings.n_jobs=16
A=f"{KB_ANALYSIS}/aki_c5_exclusion_20261022/atac_rerun/ATAC"
K=30; KNN_PROJ=25; N_DCS=10; LAT=[str(i) for i in range(30)]
EXCL=["pt_C1","pt_C2","pt_C4"]; TYPES=["PT","dPT","aPT"]
SMOKE=os.environ.get("SMOKE")=="1"; TAG="_smoke" if SMOKE else ""
print("scanpy",sc.__version__,"| anndata",ad.__version__,"| smoke",SMOKE,flush=True)

g=pd.read_csv(f"{A}/glue_v3/glue_for_pseudotime.csv",low_memory=False)
g=g[~g.cellid.astype(str).str.split(r"[_#]").str[0].isin(["K95", "K96", "K98"])].reset_index(drop=True)  # AKI exclusion
assert not g.cluster.astype(str).eq("pt_C5").any(), "pt_C5 present in the embedding"
assert g.cellid.is_unique, "cellids not unique across modalities"
g["label"]=np.where(g.modality=="RNA",g.level2.astype(str),g.cluster.astype(str))
if SMOKE:
    g=(g.groupby(["modality","label"],group_keys=False)
         .apply(lambda d:d.sample(frac=0.04,random_state=SEED)).reset_index(drop=True))
mod=g.modality.to_numpy(); isr=mod=="RNA"; isa=mod=="ATAC"
cid=g.cellid.astype(str).to_numpy()
level2=np.where(isr,g.level2.astype(str),"")
clus=np.where(isa,g.cluster.astype(str),"")
sample=np.where(isr,g.cellid.str.split("_").str[0],g.cellid.str.split("#").str[0])
X=g[LAT].to_numpy(np.float32)
print(f"cells {len(g):,}: RNA {isr.sum():,} | ATAC {isa.sum():,}",flush=True)

f1r=pd.read_csv(f"{A}/tables/rna_trajectory_percell.csv",usecols=["cellid","Factor1"]).set_index("cellid")["Factor1"]
f1a=pd.read_csv(f"{A}/tables/embeddings_glue_v3.csv",usecols=["cellid","Factor1"]).set_index("cellid")["Factor1"]
F1r=f1r.reindex(cid).to_numpy(float); F1a=f1a.reindex(cid).to_numpy(float)
assert np.isfinite(F1r[isr]).all(), "RNA Factor 1 missing for some RNA cells"
assert np.isfinite(F1a[isa]).all(), "ATAC hybrid Factor 1 missing for some ATAC cells"

# ── joint graph -> DiffMap -> pseudotime; DPT; UMAP ───────────────────────────
a=ad.AnnData(np.zeros((len(g),1),np.float32)); a.obs_names=cid; a.obsm["X_glue"]=X
sc.pp.neighbors(a,use_rep="X_glue",n_neighbors=K,random_state=SEED)
sc.tl.diffmap(a,n_comps=15); dm=a.obsm["X_diffmap"].copy()
sgn=np.sign(spearmanr(dm[isr,1],F1r[isr]).statistic); dc1=dm[:,1]*sgn; dc2=dm[:,2]
pt=(rankdata(dc1)-1)/(len(dc1)-1)
ipt=np.where(isr&(level2=="PT"))[0]
cen=X[ipt].mean(0); root=int(ipt[np.argmin(((X[ipt]-cen)**2).sum(1))])
a.uns["iroot"]=root; sc.tl.dpt(a,n_dcs=N_DCS)
dpt=a.obs["dpt_pseudotime"].to_numpy(float); fin=np.isfinite(dpt)
dptn=np.full(len(dpt),np.nan); dptn[fin]=(dpt[fin]-dpt[fin].min())/(dpt[fin].max()-dpt[fin].min())
sc.tl.umap(a,random_state=SEED); um=a.obsm["X_umap"].copy()
print("joint graph, DiffMap, DPT, UMAP done",flush=True)

# ── RNA-only axis (for the RNA-anchored projection of ATAC cells) ─────────────
r=ad.AnnData(np.zeros((int(isr.sum()),1),np.float32)); r.obsm["X_glue"]=X[isr]
sc.pp.neighbors(r,use_rep="X_glue",n_neighbors=K,random_state=SEED); sc.tl.diffmap(r,n_comps=15)
rd=r.obsm["X_diffmap"][:,1].astype(float)
if spearmanr(rd,F1r[isr]).statistic<0: rd=-rd
pt_ro=(rankdata(rd)-1)/(len(rd)-1)

nn=NearestNeighbors(n_neighbors=KNN_PROJ,n_jobs=16).fit(X[isr])
dist,ind=nn.kneighbors(X[isa])
pt_anch=np.full(len(g),np.nan); pt_anch[isa]=pt_ro[ind].mean(1)
nb=level2[isr][ind]
fr={t:np.full(len(g),np.nan) for t in TYPES}
for t in TYPES: fr[t][isa]=(nb==t).mean(1)
maj=np.full(len(g),np.nan,dtype=object)
maj[isa]=np.array(TYPES)[np.vstack([fr[t][isa] for t in TYPES]).argmax(0)]
d_ar=np.full(len(g),np.nan); d_ar[isa]=dist.mean(1)
rs=np.random.default_rng(SEED).choice(np.where(isr)[0],size=min(20000,int(isr.sum())),replace=False)
dr,_=nn.kneighbors(X[rs],n_neighbors=KNN_PROJ+1); d_rr=dr[:,1:].mean(1)

# ── summary statistics ────────────────────────────────────────────────────────
def auc(x,m1,m0):
    return float(mannwhitneyu(x[m1],x[m0],alternative="greater").statistic/(m1.sum()*m0.sum()))
S=dict(n_cells=int(len(g)),n_rna=int(isr.sum()),n_atac=int(isa.sum()),
       root_cell=str(cid[root]),root_pseudotime=float(pt[root]))
S["rho_jointDC1_vs_Factor1_rna"]=float(spearmanr(dc1[isr],F1r[isr]).statistic)
S["mean_pseudotime_rna"]={t:float(pt[isr&(level2==t)].mean()) for t in TYPES}
S["mean_pseudotime_atac"]={c:float(pt[isa&(clus==c)].mean()) for c in sorted(set(clus[isa]))}
S["auc_aPT_gt_PT_rna"]=auc(pt,isr&(level2=="aPT"),isr&(level2=="PT"))
S["auc_C6_gt_C3_atac_joint"]=auc(pt,isa&(clus=="pt_C6"),isa&(clus=="pt_C3"))
S["auc_C6_gt_C3_atac_rna_anchored"]=auc(pt_anch,isa&(clus=="pt_C6"),isa&(clus=="pt_C3"))
m_=S["mean_pseudotime_rna"]
S["dPT_position_pct_of_PT_to_aPT"]=100*(m_["dPT"]-m_["PT"])/(m_["aPT"]-m_["PT"])
S["dpt_nonfinite"]=int((~fin).sum()); S["dpt_frac_cells_below_0.2"]=float(np.mean(dptn[fin]<0.2))
S["rho_pseudotime_vs_dpt"]=float(spearmanr(pt[fin],dptn[fin]).statistic)
S["rho_dpt_vs_Factor1_rna"]=float(spearmanr(dptn[isr&fin],F1r[isr&fin]).statistic)
S["rho_joint_vs_rna_only_axis_on_rna_cells"]=float(spearmanr(pt[isr],pt_ro).statistic)
S["rho_joint_vs_rna_anchored_on_atac_cells"]=float(spearmanr(pt[isa],pt_anch[isa]).statistic)
S["median_knn_dist_atac_to_rna"]=float(np.median(d_ar[isa]))
S["median_knn_dist_rna_to_rna"]=float(np.median(d_rr))
S["atac_rna_neighbour_composition_pct"]={c:{t:round(100*float(np.nanmean(fr[t][isa&(clus==c)])),1)
                                           for t in TYPES} for c in sorted(set(clus[isa]))}
S["params"]=dict(input=f"{A}/glue_v3/glue_for_pseudotime.csv",K=K,n_dcs_dpt=N_DCS,
    knn_projection=KNN_PROJ,seed=SEED,
    orientation="sign of joint DC1 chosen so it correlates positively with RNA Factor 1",
    root_rule="RNA level2=PT cell nearest the RNA-PT centroid in the GLUE latent",
    pseudotime="(rank(joint DC1)-1)/(n-1) over all cells of both modalities",
    excluded_recommended=EXCL,scanpy=sc.__version__,anndata=ad.__version__,smoke=SMOKE)
for k,v in S.items():
    if k!="params": print(f"  {k}: {v}",flush=True)

# ── hand-off tables ───────────────────────────────────────────────────────────
use=np.where(isr,"use",np.where(np.isin(clus,EXCL),"exclude_recommended","use"))
nanif=lambda arr: pd.Series(arr).replace("",np.nan).to_numpy()
H=pd.DataFrame({"cellid":cid,"modality":mod,"sample":sample,
    "level2_rna":nanif(level2),"atac_cluster":nanif(clus),"recommended_use":use,
    "pseudotime_joint":pt,"dpt_joint":dptn,"joint_dc1":dc1,"joint_dc2":dc2,
    "umap1_joint":um[:,0],"umap2_joint":um[:,1],"pseudotime_rna_anchored":pt_anch,
    "rna_nbr_frac_PT":fr["PT"],"rna_nbr_frac_dPT":fr["dPT"],"rna_nbr_frac_aPT":fr["aPT"],
    "rna_nbr_majority":maj,"atac_to_rna_knn_dist":d_ar,
    "Factor1_rna":F1r,"Factor1_atac_genescore_hybrid":F1a})
H.to_csv(f"{A}/tables/HANDOFF_joint_pseudotime{TAG}.csv",index=False)
H[isa].to_csv(f"{A}/tables/HANDOFF_atac_cells_pseudotime{TAG}.csv",index=False)
with open(f"{A}/tables/HANDOFF_params{TAG}.json","w") as fh:
    json.dump(S,fh,indent=2,default=lambda o:o.item() if hasattr(o,"item") else str(o))
print(f"wrote HANDOFF tables ({len(H):,} rows; ATAC {int(isa.sum()):,})",flush=True)

# ── split figure ──────────────────────────────────────────────────────────────
CL={"PT":"#55A868","dPT":"#DD8452","aPT":"#937860"}
CA={"pt_C3":"#55A868","pt_C6":"#937860","pt_C1":"#4C72B0","pt_C2":"#8172B3","pt_C4":"#C44E52"}
ORD_A=["pt_C3","pt_C6","pt_C1","pt_C2","pt_C4"]
DM=np.c_[dc1,dc2]; rng=np.random.default_rng(SEED)
fig=plt.figure(figsize=(25,18)); G=fig.add_gridspec(3,4,hspace=.34,wspace=.36,top=.935)
axs=[[fig.add_subplot(G[i,j]) for j in range(4)] for i in range(3)]
def grey(ax,XY,m): ax.scatter(XY[m,0],XY[m,1],s=.25,lw=0,c="#dcdcdc",rasterized=True,zorder=1)
def lims(XY):
    p=.02*(XY.max(0)-XY.min(0)); return (XY[:,0].min()-p[0],XY[:,0].max()+p[0]),(XY[:,1].min()-p[1],XY[:,1].max()+p[1])
ulx,uly=lims(um); dlx,dly=lims(DM)
rows=[("RNA",isr,level2,CL,TYPES,.45,"level2"),("ATAC",isa,clus,CA,ORD_A,1.8,"ATAC cluster")]
for i,(nm,m,lab,pal,order,ps,labname) in enumerate(rows):
    o=rng.permutation(np.where(m)[0])
    ax=axs[i][0]; grey(ax,um,~m)
    for t in order:
        mm=m&(lab==t)
        if mm.any(): ax.scatter(um[mm,0],um[mm,1],s=ps,lw=0,alpha=.45 if t in EXCL else .6,
                                c=pal[t],rasterized=True,zorder=2)
    hs=[Line2D([],[],marker="o",ls="",ms=7,mfc=pal[t],mec="none",
               label=f"{t} ({int((m&(lab==t)).sum()):,})"+("  excl. recommended" if t in EXCL else ""))
        for t in order if (m&(lab==t)).any()]
    hs.append(Line2D([],[],marker="o",ls="",ms=7,mfc="#dcdcdc",mec="none",label="other modality"))
    ax.legend(handles=hs,fontsize=8.5,loc="best",framealpha=.9)
    ax.set_title(f"{nm} cells ({int(m.sum()):,}) — {labname}",fontweight="bold")
    for j,(vals,lbl,ttl) in enumerate([(pt,"joint pseudotime","joint pseudotime (primary)"),
                                       (dptn,"DPT (joint graph)","DPT rooted in RNA PT")],start=1):
        ax=axs[i][j]; grey(ax,um,~m); oo=o[np.isfinite(vals[o])]
        s_=ax.scatter(um[oo,0],um[oo,1],s=ps,lw=0,c=vals[oo],cmap="viridis",vmin=0,vmax=1,
                      rasterized=True,zorder=2)
        plt.colorbar(s_,ax=ax,fraction=.046,pad=.02,label=lbl)
        ax.set_title(f"{nm} — {ttl}",fontweight="bold")
    ax=axs[i][3]; grey(ax,DM,~m)
    s_=ax.scatter(DM[o,0],DM[o,1],s=ps,lw=0,c=pt[o],cmap="viridis",vmin=0,vmax=1,rasterized=True,zorder=2)
    plt.colorbar(s_,ax=ax,fraction=.046,pad=.02,label="joint pseudotime")
    ax.set_title(f"{nm} — joint pseudotime on DiffMap",fontweight="bold")
    ax.set_xlim(dlx); ax.set_ylim(dly); ax.set_xlabel("joint DC1"); ax.set_ylabel("joint DC2")
    for j in range(3):
        axs[i][j].set_xlim(ulx); axs[i][j].set_ylim(uly)
        axs[i][j].set_xticks([]); axs[i][j].set_yticks([])
        axs[i][j].set_xlabel("UMAP1 (joint)"); axs[i][j].set_ylabel("UMAP2 (joint)")
for j in (1,2): axs[0][j].scatter(um[root,0],um[root,1],marker="*",s=320,c="red",edgecolor="k",zorder=6)
axs[0][3].scatter(DM[root,0],DM[root,1],marker="*",s=320,c="red",edgecolor="k",zorder=6)
bins=np.linspace(0,1,51)
ax=axs[2][0]
for t in TYPES:
    mm=isr&(level2==t); ax.hist(pt[mm],bins=bins,density=True,alpha=.55,color=CL[t],
                                label=f"{t}  mean {pt[mm].mean():.3f}")
ax.legend(fontsize=9); ax.set_title("RNA — joint pseudotime by level2",fontweight="bold")
ax.set_xlabel("joint pseudotime"); ax.set_ylabel("density")
ax=axs[2][1]
for c in ORD_A:
    mm=isa&(clus==c)
    if not mm.any(): continue
    if c in EXCL: ax.hist(pt[mm],bins=bins,density=True,histtype="step",lw=1.7,color=CA[c],
                          label=f"{c}  mean {pt[mm].mean():.3f}  (excl.)")
    else: ax.hist(pt[mm],bins=bins,density=True,alpha=.55,color=CA[c],label=f"{c}  mean {pt[mm].mean():.3f}")
ax.legend(fontsize=8.5); ax.set_title("ATAC — joint pseudotime by cluster",fontweight="bold")
ax.set_xlabel("joint pseudotime"); ax.set_ylabel("density")
ax=axs[2][2]; ca=np.array([CA[c] for c in clus[isa]]); oo=rng.permutation(int(isa.sum()))
ax.scatter(pt[isa][oo],pt_anch[isa][oo],s=1.6,lw=0,alpha=.35,c=ca[oo],rasterized=True)
ax.plot([0,1],[0,1],"k--",lw=.8)
ax.set_title(f"ATAC — joint vs RNA-anchored pseudotime (rho {S['rho_joint_vs_rna_anchored_on_atac_cells']:.3f})",
             fontweight="bold")
ax.set_xlabel("joint pseudotime"); ax.set_ylabel("RNA-anchored pseudotime\n(mean of 25 nearest RNA cells)")
ax=axs[2][3]; ax.axis("off"); mr=S["mean_pseudotime_rna"]; ma=S["mean_pseudotime_atac"]
txt=("One pseudotime, both modalities\n\n"
     f"cells   {S['n_rna']:,} RNA + {S['n_atac']:,} ATAC\n"
     f"graph   kNN K={K} on 30-dim GLUE latent\n"
     f"axis    rank of joint DiffMap DC1\n"
     f"root    RNA PT centroid cell (red star)\n\n"
     f"RNA   PT {mr['PT']:.3f} > dPT {mr['dPT']:.3f} > aPT {mr['aPT']:.3f}\n"
     f"      AUC(aPT>PT) {S['auc_aPT_gt_PT_rna']:.3f}\n"
     f"ATAC  C3 {ma.get('pt_C3',float('nan')):.3f} > C6 {ma.get('pt_C6',float('nan')):.3f}\n"
     f"      AUC(C6>C3) {S['auc_C6_gt_C3_atac_joint']:.3f}\n\n"
     f"joint vs RNA-only axis   rho {S['rho_joint_vs_rna_only_axis_on_rna_cells']:.3f}\n"
     f"joint vs RNA-anchored    rho {S['rho_joint_vs_rna_anchored_on_atac_cells']:.3f}\n"
     f"DPT: {100*S['dpt_frac_cells_below_0.2']:.0f}% of cells < 0.2 (compressed)\n"
     f"kNN dist ATAC>RNA {S['median_knn_dist_atac_to_rna']:.3f}\n"
     f"         RNA>RNA  {S['median_knn_dist_rna_to_rna']:.3f}")
ax.text(0,1,txt,va="top",ha="left",fontsize=11.5,family="monospace")
fig.suptitle("One pseudotime on the combined GLUE embedding — RNA and ATAC cells shown separately "
             "(other modality in grey)",fontsize=16,fontweight="bold",y=.975)
fig.savefig(f"{A}/figures/TRAJECTORY_joint_split_rna_atac{TAG}.png",dpi=130,bbox_inches="tight")
fig.savefig(f"{A}/figures/TRAJECTORY_joint_split_rna_atac{TAG}.pdf",bbox_inches="tight")
print("wrote split figure\ndone",flush=True)

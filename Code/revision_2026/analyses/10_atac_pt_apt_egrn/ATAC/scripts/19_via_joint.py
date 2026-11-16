#!/usr/bin/env python3
"""VIA pseudotime on the combined GLUE embedding — does it agree with the DC1 ordering?

VIA (pyVIA) is root-anchored: it clusters the kNN graph, starts random walks at a root cluster and
measures Markov hitting time; terminal states are found automatically and each cell gets branch
probabilities towards them. Unlike DC1, its direction comes from the ROOT, not from Factor 1, so
the Factor 1 correlation below is not signed by construction.

Two runs on the same cells (glue_v3 RNA + ATAC, 30-dim GLUE latent, K=30), root = RNA level2 'PT':
  A  terminal states detected by VIA (the fair test; shows whether dPT is its own end state)
  B  terminal state enforced = aPT (a constrained PT -> aPT pseudotime)
Compared with pseudotime_joint (rank DC1) from tables/HANDOFF_joint_pseudotime.csv.
Human hg38 | seed 42 | SMOKE=1 -> 4% stratified subsample, outputs tagged _smoke
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import numpy as np, pandas as pd, os, json, time
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, mannwhitneyu
import pyVIA.core as via

SEED=42; np.random.seed(SEED)
A=f"{KB_ANALYSIS}/aki_c5_exclusion_20261022/atac_rerun/ATAC"; KNN=30; LAT=[str(i) for i in range(30)]
SMOKE=os.environ.get("SMOKE")=="1"; TAG="_smoke" if SMOKE else ""
g=pd.read_csv(f"{A}/glue_v3/glue_for_pseudotime.csv",low_memory=False)
g=g[~g.cellid.astype(str).str.split(r"[_#]").str[0].isin(["K95", "K96", "K98"])].reset_index(drop=True)  # AKI exclusion
g["label"]=np.where(g.modality=="RNA",g.level2.astype(str),g.cluster.astype(str))
if SMOKE:
    g=(g.groupby(["modality","label"],group_keys=False)
         .apply(lambda d:d.sample(frac=0.04,random_state=SEED)).reset_index(drop=True))
H=pd.read_csv(f"{A}/tables/HANDOFF_joint_pseudotime.csv",low_memory=False,
    usecols=["cellid","pseudotime_joint","dpt_joint","Factor1_rna","joint_dc2","umap1_joint","umap2_joint"]
    ).set_index("cellid").reindex(g.cellid)
assert H.pseudotime_joint.notna().all(), "cells missing from the hand-off table"
sub=pd.read_csv(f"{A}/tables/dPT_percell_scores.csv",usecols=["cellid","dPT_subset"]).set_index("cellid")["dPT_subset"]
X=g[LAT].to_numpy(np.float64); lab=g.label.tolist(); L=np.array(lab)
mod=g.modality.to_numpy(); isr=mod=="RNA"; isa=mod=="ATAC"
G2=L.astype(object); md=L=="dPT"; G2[md]="dPT_"+sub.reindex(g.cellid[md]).fillna("unscored").to_numpy()
dc=H.pseudotime_joint.to_numpy(); f1=H.Factor1_rna.to_numpy(); um=H[["umap1_joint","umap2_joint"]].to_numpy()
print(f"cells {len(g):,} (RNA {isr.sum():,} / ATAC {isa.sum():,}) | smoke {SMOKE}",flush=True)

def run(name,term):
    t=time.time()
    kw=dict(data=X,true_label=lab,knn=KNN,random_seed=SEED,num_threads=16,root_user=["PT"],
            do_compute_embedding=False)
    if term: kw["user_defined_terminal_group"]=term
    v=via.VIA(**kw); v.run_VIA()
    pt=np.asarray(v.single_cell_pt_markov,float)
    # pyVIA silently filters NaN entries out of single_cell_pt_markov (core.py ~l.2515); a shorter
    # list would shift every cell after the first NaN, so refuse rather than misalign
    assert len(pt)==len(g), f"VIA returned {len(pt)} pseudotimes for {len(g)} cells (NaNs dropped) - alignment lost"
    pt=(pt-np.nanmin(pt))/(np.nanmax(pt)-np.nanmin(pt))
    terms=[int(x) for x in getattr(v,"terminal_clusters",[])]; vl=np.asarray(v.labels)
    bp,bpn=None,None
    for at in ["single_cell_bp_rownormed","single_cell_bp"]:
        val=getattr(v,at,None)
        if val is not None:
            bp=np.asarray(val,float); bpn=at; break
    print(f"[{name}] {time.time()-t:.0f}s | VIA clusters {len(np.unique(vl))} | terminals {terms} | "
          f"branch probs {bpn} {None if bp is None else bp.shape}",flush=True)
    return pt,terms,vl,bp

def auc(x,m1,m0): return float(mannwhitneyu(x[m1],x[m0],alternative="greater").statistic/(m1.sum()*m0.sum()))
def compare(pt):
    S={"rho_vs_DC1_all":spearmanr(pt,dc).statistic,"rho_vs_DC1_rna":spearmanr(pt[isr],dc[isr]).statistic,
       "rho_vs_DC1_atac":spearmanr(pt[isa],dc[isa]).statistic,
       "rho_vs_Factor1_rna_unsigned_by_construction":spearmanr(pt[isr],f1[isr]).statistic,
       "auc_aPT_gt_PT":auc(pt,isr&(L=="aPT"),isr&(L=="PT")),"auc_C6_gt_C3":auc(pt,isa&(L=="pt_C6"),isa&(L=="pt_C3")),
       "frac_cells_below_0.2":float(np.mean(pt<0.2))}
    S["mean_pt"]={k:float(pt[G2==k].mean()) for k in ["PT","dPT_PT_lineage","dPT_distal_like","aPT",
                  "pt_C3","pt_C6","pt_C1","pt_C2","pt_C4"] if (G2==k).any()}
    m=pt[L=="dPT"].mean(); S["dPT_position_pct"]=100*(m-S["mean_pt"]["PT"])/(S["mean_pt"]["aPT"]-S["mean_pt"]["PT"])
    return {k:(float(v) if isinstance(v,(float,np.floating)) else v) for k,v in S.items()}
def terminals(vl,terms):
    rows=[]
    for t in terms:
        m=vl==t; vc=pd.Series(G2[m]).value_counts(normalize=True)
        rows.append({"via_cluster":t,"n":int(m.sum()),"majority":vc.index[0],"majority_pct":round(100*vc.iloc[0],1),
                     "composition":"; ".join(f"{k} {100*p:.0f}%" for k,p in vc.head(4).items()),
                     "mean_joint_dc2":float(H.joint_dc2.to_numpy()[m].mean())})
    return pd.DataFrame(rows)

out=pd.DataFrame({"cellid":g.cellid.values,"modality":mod,"label":L,"group":G2,"pseudotime_joint":dc})
R={}
for name,term in [("A_auto_terminals",None),("B_aPT_terminal",["aPT"])]:
    pt,terms,vl,bp=run(name,term)
    col="pseudotime_via" if name.startswith("A") else "pseudotime_via_aPTterminal"
    out[col]=pt; out[f"via_cluster_{name[0]}"]=vl
    T=terminals(vl,terms); S=compare(pt); S["terminals"]=T.to_dict("records")
    if bp is not None and bp.shape[1]==len(terms):
        B=pd.DataFrame(bp,columns=[f"bp_to_{t}" for t in terms]); B["group"]=G2
        S["branch_prob_by_group"]=B.groupby("group").mean().round(3).to_dict("index")
        if name.startswith("A"):
            for c in B.columns[:-1]: out[f"{c}_A"]=B[c].values
    R[name]=S
    print(f"\n=== {name} ===\n"+T.to_string(index=False),flush=True)
    for k,v in S.items():
        if k!="terminals": print(f"  {k}: {v}",flush=True)
out.to_csv(f"{A}/tables/via_joint_pseudotime{TAG}.csv",index=False)
json.dump(R,open(f"{A}/tables/via_joint_summary{TAG}.json","w"),indent=2,
          default=lambda o:o.item() if hasattr(o,"item") else str(o))

CL={"PT":"#55A868","dPT_PT_lineage":"#DD8452","dPT_distal_like":"#8172B3","aPT":"#937860"}
CA={"pt_C3":"#55A868","pt_C6":"#937860","pt_C1":"#4C72B0","pt_C2":"#8172B3","pt_C4":"#C44E52"}
fig,ax=plt.subplots(2,4,figsize=(26,11.5)); rng=np.random.default_rng(SEED)
pa=out.pseudotime_via.to_numpy(); pb=out.pseudotime_via_aPTterminal.to_numpy()
for j,(m,nm,ps) in enumerate([(isr,"RNA",.4),(isa,"ATAC",1.6)]):
    a_=ax[0,j]; a_.scatter(um[~m,0],um[~m,1],s=.2,c="#dcdcdc",lw=0,rasterized=True)
    o=rng.permutation(np.where(m)[0])
    s_=a_.scatter(um[o,0],um[o,1],s=ps,c=pa[o],cmap="viridis",vmin=0,vmax=1,lw=0,rasterized=True)
    plt.colorbar(s_,ax=a_,fraction=.046,label="VIA pseudotime (run A)")
    a_.set_title(f"{nm} cells — VIA pseudotime (auto terminals)",fontweight="bold"); a_.set_xticks([]); a_.set_yticks([])
for j,(pv,nm,rk) in enumerate([(pa,"A auto terminals","A_auto_terminals"),(pb,"B aPT terminal","B_aPT_terminal")],start=2):
    a_=ax[0,j]
    for k,c in CL.items():
        m=G2==k; o=rng.permutation(np.where(m)[0])[:20000]
        a_.scatter(dc[o],pv[o],s=.5,c=c,alpha=.35,lw=0,rasterized=True,label=k)
    a_.set_xlabel("rank-DC1 axis (pseudotime_joint)"); a_.set_ylabel(f"VIA pseudotime ({nm})")
    a_.set_title(f"VIA {nm} vs DC1 — RNA rho {R[rk]['rho_vs_DC1_rna']:.3f}",fontweight="bold")
    lg=a_.legend(markerscale=12,fontsize=8)
    for h in lg.legend_handles: h.set_alpha(1)
b_=np.linspace(0,1,51)
for k,c in CL.items():
    m=G2==k; ax[1,0].hist(pa[m],bins=b_,density=True,alpha=.5,color=c,label=f"{k} mean {pa[m].mean():.3f}")
ax[1,0].legend(fontsize=8); ax[1,0].set_title("RNA — VIA pseudotime (A) by group",fontweight="bold")
for k,c in CA.items():
    m=L==k
    if m.any(): ax[1,1].hist(pa[m],bins=b_,density=True,histtype="step" if k in("pt_C1","pt_C2","pt_C4") else "stepfilled",
                            alpha=.55,lw=1.6,color=c,label=f"{k} mean {pa[m].mean():.3f}")
ax[1,1].legend(fontsize=8); ax[1,1].set_title("ATAC — VIA pseudotime (A) by cluster",fontweight="bold")
nb=25; ed=np.linspace(0,1,nb+1); bid=np.clip(np.digitize(pa,ed)-1,0,nb-1); ctr=(ed[:-1]+ed[1:])/2
ok=isr&np.isfinite(f1); o=rng.permutation(np.where(ok)[0])[:40000]
ax[1,2].scatter(pa[o],f1[o],s=.4,alpha=.12,c="#34495e",lw=0,rasterized=True)
ax[1,2].plot(ctr,[np.nanmean(f1[ok&(bid==i)]) if (ok&(bid==i)).sum()>20 else np.nan for i in range(nb)],c="#c0392b",lw=3)
ax[1,2].set_title(f"Factor 1 along VIA pseudotime (A) — rho {R['A_auto_terminals']['rho_vs_Factor1_rna_unsigned_by_construction']:+.3f}\n(sign NOT imposed: direction comes from the PT root)",fontweight="bold")
ax[1,2].set_xlabel("VIA pseudotime (A)"); ax[1,2].set_ylabel("Factor 1 (RNA)")
bpg=R["A_auto_terminals"].get("branch_prob_by_group")
if bpg:
    Bd=pd.DataFrame(bpg).T.reindex([k for k in CL if k in bpg])
    tmap={f"bp_to_{d['via_cluster']}":f"{d['via_cluster']}: {d['majority']} {d['majority_pct']:.0f}%"
          for d in R["A_auto_terminals"]["terminals"]}
    Bd.rename(columns=tmap).plot(kind="bar",ax=ax[1,3],rot=12,width=.85)
    ax[1,3].set_ylabel("mean branch probability (row-normalised)")
    ax[1,3].set_title("Where each group is headed (run A terminals)",fontweight="bold")
    ax[1,3].legend(fontsize=7.5,title="terminal cluster: majority",title_fontsize=8,
                   loc="upper left",bbox_to_anchor=(1.01,1),borderaxespad=0)
else:
    ax[1,3].axis("off"); ax[1,3].text(0,.5,"branch probabilities not available",fontsize=12)
fig.suptitle("VIA pseudotime on the combined GLUE embedding vs the DC1 ordering",fontsize=16,fontweight="bold")
fig.tight_layout(rect=[0,0,1,.95])
fig.savefig(f"{A}/figures/VIA_vs_DC1{TAG}.png",dpi=130,bbox_inches="tight"); fig.savefig(f"{A}/figures/VIA_vs_DC1{TAG}.pdf",bbox_inches="tight")
print("\nwrote via_joint_pseudotime / via_joint_summary / VIA_vs_DC1\ndone",flush=True)

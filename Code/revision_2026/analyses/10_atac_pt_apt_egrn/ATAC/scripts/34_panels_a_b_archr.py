#!/usr/bin/env python3
"""Panels a and b (as scripts/26_panels_a_b.py) for the ArchR trajectories of 31_archr_trajectory.R. TRAJ=archr_dm | archr_glue.
Since 2026-09-14 the archr_dm set is the main figure (scmega/figures/); archr_glue is an alternative (scmega/figures_archr_glue/).
The ArchR pseudotime exists per PAIR of the scMEGA paired object, not per cell, so:
a  the embedding the ArchR curve was fitted in, with the pseudotime laid on it (Leon, 2026-09-14): archr_dm -> the joint
   diffusion map (DC1 / DC2; the fit used exactly these two, at pair means); archr_glue -> the joint UMAP of the GLUE
   latent (the fit used all 30 GLUE dimensions; the UMAP of the same kNN graph is only its 2-D view). Both cells of each
   pair are drawn at their own coordinates, coloured by the pair's ArchR value; cells of pairs ArchR left without a value
   in dark grey; all other cells light grey. Arrow = medians of the coloured cells' coordinates in 25 bins of the ArchR
   pseudotime, smoothed, pointing to the highest bin.
b  per-donor densities along the ArchR pseudotime of the RNA members of the pairs with a value (sample, disease and
   eGFR group of that RNA cell); donors with >= 50 such cells; mean +/- 95% t-interval across donors (as in 26). These
   are fewer cells than the DC1 panel b, which uses every RNA PT / aPT cell.
Out: scmega/{figures,figures_archr_glue}/panel_a_b_<TRAJ>.{png,pdf}; scmega/outputs/panel_b_donor_densities_<TRAJ>.csv
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import os, pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.patheffects as pe
from scipy.ndimage import gaussian_filter1d
from scipy.stats import t as tdist
TRAJ=os.environ.get("TRAJ",""); DIRS={"archr_dm":"figures","archr_glue":"figures_archr_glue"}
assert TRAJ in DIRS, "set TRAJ=archr_dm or TRAJ=archr_glue"
LAB={"archr_glue":"ArchR, GLUE latent","archr_dm":"ArchR, joint diffusion map"}[TRAJ]
A=f"{KB_ANALYSIS}/aki_c5_exclusion_20261022/atac_rerun/ATAC"; O=os.environ.get("FIGDIR") or f"{A}/scmega/{DIRS[TRAJ]}"; os.makedirs(O,exist_ok=True)
META=f"{KB_ANALYSIS}/vivien_takeover/figures/milo_donor_composition/sample_metadata_milo.csv"
H=pd.read_csv(f"{A}/tables/HANDOFF_joint_pseudotime.csv",low_memory=False)
Hi=H.set_index("cellid")
P=pd.read_csv(f"{A}/tables/archr_trajectories.csv")
P["t"]=P[f"traj_{TRAJ}"]/100                                        # ArchR scale (0, 100] -> (0, 1]
P["sample"]=Hi.loc[P.rna_cellid,"sample"].to_numpy()
has=P.t.notna().to_numpy()
print(f"{TRAJ}: {len(P):,} pairs, {int(has.sum()):,} with an ArchR value, {int((~has).sum()):,} without",flush=True)
cx,cy,ename,xl,yl={"archr_glue":("umap1_joint","umap2_joint","GLUE latent (joint UMAP view) + ArchR pseudotime","UMAP 1","UMAP 2"),
                   "archr_dm":("joint_dc1","joint_dc2","joint diffusion map + ArchR pseudotime","DC 1","DC 2")}[TRAJ]
Mb=pd.concat([P[["rna_cellid","t"]].set_axis(["cellid","t"],axis=1),P[["atac_cellid","t"]].set_axis(["cellid","t"],axis=1)],ignore_index=True)
Mb[cx]=Hi.loc[Mb.cellid,cx].to_numpy(); Mb[cy]=Hi.loc[Mb.cellid,cy].to_numpy()
R=H[H.modality=="RNA"]
fig,ax=plt.subplots(1,3,figsize=(21,6.2),gridspec_kw={"width_ratios":[1.05,1,1]})
# ── a ──
# ONLY_USED=1 (Leon, 2026-09-24): show the embedding with ONLY the cells the trajectory actually uses, i.e. the cells
# of pairs that carry an ArchR value; the other cells of the object are not drawn and the axes are refitted to them.
ONLY=os.environ.get("ONLY_USED")=="1"
a=ax[0]; bg=H[~H.cellid.isin(Mb.cellid)]
if not ONLY:
    a.scatter(bg[cx],bg[cy],s=.5,lw=0,c="#e2e2e2",rasterized=True,label=f"cells not in a pair ({len(bg):,})")
wv=Mb[Mb.t.notna()].sample(frac=1,random_state=42); nv=Mb[Mb.t.isna()]
sca=a.scatter(wv[cx],wv[cy],s=.8 if not ONLY else 1.6,lw=0,c=100*wv.t,cmap="plasma_r",vmin=0,vmax=100,rasterized=True)
if not ONLY:
    a.scatter(nv[cx],nv[cy],s=1.2,lw=0,c="#4a4a4a",rasterized=True,label=f"cells of pairs without an ArchR value ({len(nv):,})")
if ONLY:
    mx,my=wv[cx].to_numpy(),wv[cy].to_numpy(); px_,py_=(mx.max()-mx.min())*.04,(my.max()-my.min())*.04
    a.set_xlim(mx.min()-px_,mx.max()+px_); a.set_ylim(my.min()-py_,my.max()+py_)
bins=np.linspace(0,1,26); b=np.clip(np.digitize(wv.t,bins)-1,0,24)
px=gaussian_filter1d(wv.groupby(b)[cx].median().to_numpy(),1.2)
py=gaussian_filter1d(wv.groupby(b)[cy].median().to_numpy(),1.2)
a.plot(px[:-1],py[:-1],c="k",lw=2.6)
a.annotate("",xy=(px[-1],py[-1]),xytext=(px[-3],py[-3]),arrowprops=dict(arrowstyle="-|>",lw=2.6,color="k",mutation_scale=22))
def pos(t):
    m=(R.level2_rna==t)
    if cx=="joint_dc1":
        return {"PT":(np.percentile(R.joint_dc1[m],15),0.006),"aPT":(np.percentile(R.joint_dc1[m],90),0.006),
                "dPT":(R.joint_dc1[m].median(),R.joint_dc2[m].quantile(.97))}[t]
    # UMAP: the arrow runs through the cell-type medians, so put PT / aPT labels towards their outer edges
    return {"PT":(np.percentile(R[cx][m],8),R[cy][m].median()),"aPT":(np.percentile(R[cx][m],92),R[cy][m].median()),
            "dPT":(R[cx][m].median(),R[cy][m].median())}[t]
# in ONLY_USED mode the dPT branch is not drawn at all, and its label would sit far outside the new limits
labs=[("PT","PT"),("aPT","aPT")] if ONLY else [("PT","PT"),("aPT","aPT"),("dPT","dPT (branch,\nnot on trajectory)")]
for t,txt in labs:
    if ONLY:
        q=wv[wv.t<=wv.t.quantile(.1)] if t=="PT" else wv[wv.t>=wv.t.quantile(.9)]
        x,y=q[cx].median(),q[cy].median()
    else:
        x,y=pos(t)
    a.text(x,y,txt,fontsize=13,fontweight="bold",color="k",ha="center",path_effects=[pe.withStroke(linewidth=3,foreground="white")])
a.legend(loc="upper right",fontsize=9,frameon=False,markerscale=8)
fig.colorbar(sca,ax=a,fraction=.035,pad=.01).set_label("ArchR pseudotime (0–100)",fontsize=10)
a.set_xlabel(xl,fontsize=12); a.set_ylabel(yl,fontsize=12); a.set_xticks([]); a.set_yticks([])
a.set_title(f"a  {ename}",fontweight="bold",loc="left")
# ── b ──
Q=P[has]
R2=Q.merge(pd.read_csv(META)[["sample","disease","egfr_group"]],on="sample",how="left")
nod=R2.disease.isna(); print(f"panel b: {len(R2):,} RNA cells of pairs with an ArchR value; {int(nod.sum()):,} from "
      f"{R2.loc[nod,'sample'].nunique()} samples without a disease label excluded ({sorted(R2.loc[nod,'sample'].unique())})",flush=True)
NB=40; ed=np.linspace(0,1,NB+1); ctr=100*(ed[:-1]+ed[1:])/2
def donor_curves(df,grp):
    rows=[]
    for (g,s),d in df.groupby([grp,"sample"]):
        if len(d)<50: continue
        h=np.histogram(d.t,bins=ed)[0]/len(d)
        rows.append({"group":g,"sample":s,"n":len(d),**{f"b{i}":v for i,v in enumerate(gaussian_filter1d(h,1.2))}})
    return pd.DataFrame(rows)
GC={"control":"#2c3e50","MN":"#e74c3c","FSGS":"#3498db","IgA":"#9b59b6","DN":"#e67e22","TIN":"#16a085",
    ">90":"#2c3e50","60-90":"#3498db","30-60":"#e67e22","<30":"#e74c3c"}
tabs=[]
for k,(grp,title,order) in enumerate([("disease","b  disease group along the pseudotime",["control","MN","FSGS","IgA","DN","TIN"]),
                                      ("egfr_group","b′  eGFR group along the pseudotime",[">90","60-90","30-60","<30"])],start=1):
    D=donor_curves(R2.dropna(subset=[grp]),grp); D["grouping"]=grp; tabs.append(D)
    a=ax[k]; bc=[f"b{i}" for i in range(NB)]
    for g in order:
        d=D[D.group==g]
        if len(d)<2: continue
        Y=d[bc].to_numpy(); mu=Y.mean(0); se=Y.std(0,ddof=1)/np.sqrt(len(d)); hw=tdist.ppf(.975,len(d)-1)*se
        a.plot(ctr,mu,lw=2.3,c=GC[g],label=f"{g} ({len(d)} donors)"); a.fill_between(ctr,mu-hw,mu+hw,color=GC[g],alpha=.12)
    a.set_xlabel(f"Pseudotime ({LAB}; RNA cells of pairs)",fontsize=12)
    a.set_ylabel("Fraction of cells (per-donor mean ± 95% CI)",fontsize=11)
    a.set_title(title,fontweight="bold",loc="left"); a.legend(fontsize=9,frameon=False); a.set_xlim(0,100)
TD=os.environ.get("TABLEDIR") or f"{A}/scmega/outputs"  # output folder (release patch: the figure entry scripts in figures/ redirect it; default = original folder)
pd.concat(tabs).to_csv(f"{TD}/panel_b_donor_densities_{TRAJ}.csv",index=False)
fig.tight_layout()
fig.savefig(f"{O}/panel_a_b_{TRAJ}.png",dpi=160,bbox_inches="tight"); fig.savefig(f"{O}/panel_a_b_{TRAJ}.pdf",bbox_inches="tight")
print(f"wrote {O}/panel_a_b_{TRAJ}.png/.pdf",flush=True)

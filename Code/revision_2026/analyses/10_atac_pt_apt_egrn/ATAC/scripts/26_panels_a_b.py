#!/usr/bin/env python3
"""Panels a and b of the scOpen-Fig.5-style figure, for PT -> aPT (scMEGA scope; dPT excluded), DC1-rank set.
Since 2026-09-14 the main figure uses the ArchR diffusion-map trajectory (34_panels_a_b_archr.py); this set is an
alternative in scmega/figures_dc1rank/.

a  The embedding the pseudotime was computed on (Leon, 2026-09-14): the joint diffusion map (DC1 / DC2) of the kNN graph
   on the GLUE latent of all RNA + ATAC cells (script 12), with the pseudotime laid on the cells of the PT -> aPT set
   (RNA PT / aPT, ATAC C3 / C6); cells outside it (RNA dPT, ATAC C1 / C2 / C4) in grey. Arrow = medians of those cells'
   DC1 / DC2 in 25 pseudotime bins, smoothed, pointing to the highest bin.
b  Where each disease group's cells sit along the pseudotime (RNA PT + aPT). The paper pools cells;
   pooled curves are dominated by donors with many cells, so the primary version averages per-donor
   densities (donors with >= 50 cells) with a 95% t-interval across donors. eGFR group as a second
   grouping. Samples without a disease label are reported and excluded.
Pseudotime = pseudotime_joint (rank DC1, validated against VIA), re-ranked within the PT -> aPT set
and shown on a 0-100 scale.
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import pandas as pd, numpy as np, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.patheffects as pe
from scipy.ndimage import gaussian_filter1d
from scipy.stats import t as tdist
A=f"{KB_ANALYSIS}/aki_c5_exclusion_20261022/atac_rerun/ATAC"; O=f"{A}/scmega/figures_dc1rank"; os.makedirs(O,exist_ok=True)
META=f"{KB_ANALYSIS}/vivien_takeover/figures/milo_donor_composition/sample_metadata_milo.csv"
H=pd.read_csv(f"{A}/tables/HANDOFF_joint_pseudotime.csv",low_memory=False)
R=H[H.modality=="RNA"].copy()
sm=pd.read_csv(META)[["sample","disease","egfr_group"]]
R=R.merge(sm,on="sample",how="left")
# The joint pseudotime is ranked over ALL cells incl. dPT, which occupies ~0.60-0.70; plotting PT/aPT
# alone would show an empty stretch there. Re-rank within the PT -> aPT set (RNA PT/aPT + ATAC
# C3/C6): order unchanged, gap closed. scMEGA uses the same scale (tables/trajectory_pt_traj.csv).
T=H[((H.modality=="RNA")&H.level2_rna.isin(["PT","aPT"]))|((H.modality=="ATAC")&(H.recommended_use=="use"))].copy()
T["pt_traj"]=(T.pseudotime_joint.rank(method="average")-1)/(len(T)-1)
T[["cellid","modality","pt_traj"]].to_csv(f"{A}/tables/trajectory_pt_traj.csv",index=False)
fig,ax=plt.subplots(1,3,figsize=(21,6.2),gridspec_kw={"width_ratios":[1.05,1,1]})
# ── a ──
a=ax[0]; out=H[~H.cellid.isin(T.cellid)]
a.scatter(out.joint_dc1,out.joint_dc2,s=.5,lw=0,c="#cfcfcf",rasterized=True,
          label=f"not in PT → aPT set: RNA dPT, ATAC C1/C2/C4 ({len(out):,})")
o=T.sample(frac=1,random_state=42)
sca=a.scatter(o.joint_dc1,o.joint_dc2,s=.5,lw=0,c=100*o.pt_traj,cmap="plasma_r",vmin=0,vmax=100,rasterized=True)
bins=np.linspace(0,1,26); b=np.clip(np.digitize(T.pt_traj,bins)-1,0,24)
px=gaussian_filter1d(T.groupby(b).joint_dc1.median().to_numpy(),1.2)
py=gaussian_filter1d(T.groupby(b).joint_dc2.median().to_numpy(),1.2)
a.plot(px[:-1],py[:-1],c="k",lw=2.6)
a.annotate("",xy=(px[-1],py[-1]),xytext=(px[-3],py[-3]),arrowprops=dict(arrowstyle="-|>",lw=2.6,color="k",mutation_scale=22))
for t,(x,y) in {"PT":(np.percentile(R.joint_dc1[R.level2_rna=="PT"],15),0.006),
                "aPT":(np.percentile(R.joint_dc1[R.level2_rna=="aPT"],90),0.006),
                "dPT (branch,\nnot on trajectory)":(R.joint_dc1[R.level2_rna=="dPT"].median(),R.joint_dc2[R.level2_rna=="dPT"].quantile(.97))}.items():
    a.text(x,y,t,fontsize=13,fontweight="bold",color="k",ha="center",path_effects=[pe.withStroke(linewidth=3,foreground="white")])
a.legend(loc="upper right",fontsize=9,frameon=False,markerscale=8)
fig.colorbar(sca,ax=a,fraction=.035,pad=.01).set_label("pseudotime (0–100)",fontsize=10)
a.set_xlabel("DC 1",fontsize=12); a.set_ylabel("DC 2",fontsize=12); a.set_xticks([]); a.set_yticks([])
a.set_title("a  joint diffusion map (RNA + ATAC) + pseudotime",fontweight="bold",loc="left")
# ── b ──
R=R[R.level2_rna.isin(["PT","aPT"])].merge(T[["cellid","pt_traj"]],on="cellid",how="left")
assert R.pt_traj.notna().all()
nod=R.disease.isna(); print(f"panel b: {int(nod.sum()):,} RNA PT/aPT cells from "
      f"{R.loc[nod,'sample'].nunique()} samples without a disease label excluded "
      f"({sorted(R.loc[nod,'sample'].unique())})",flush=True)
NB=40; ed=np.linspace(0,1,NB+1); ctr=100*(ed[:-1]+ed[1:])/2
def donor_curves(df,grp):
    rows=[]
    for (g,s),d in df.groupby([grp,"sample"]):
        if len(d)<50: continue
        h=np.histogram(d.pt_traj,bins=ed)[0]/len(d)
        rows.append({"group":g,"sample":s,"n":len(d),**{f"b{i}":v for i,v in enumerate(gaussian_filter1d(h,1.2))}})
    return pd.DataFrame(rows)
GC={"control":"#2c3e50","MN":"#e74c3c","FSGS":"#3498db","IgA":"#9b59b6","DN":"#e67e22","TIN":"#16a085",
    ">90":"#2c3e50","60-90":"#3498db","30-60":"#e67e22","<30":"#e74c3c"}
tabs=[]
for k,(grp,title,order) in enumerate([("disease","b  disease group along the pseudotime",["control","MN","FSGS","IgA","DN","TIN"]),
                                      ("egfr_group","b′  eGFR group along the pseudotime",[">90","60-90","30-60","<30"])],start=1):
    D=donor_curves(R.dropna(subset=[grp]),grp); D["grouping"]=grp; tabs.append(D)
    a=ax[k]; bc=[f"b{i}" for i in range(NB)]
    for g in order:
        d=D[D.group==g]
        if len(d)<2: continue
        Y=d[bc].to_numpy(); mu=Y.mean(0); se=Y.std(0,ddof=1)/np.sqrt(len(d)); hw=tdist.ppf(.975,len(d)-1)*se
        a.plot(ctr,mu,lw=2.3,c=GC[g],label=f"{g} ({len(d)} donors)"); a.fill_between(ctr,mu-hw,mu+hw,color=GC[g],alpha=.12)
    a.set_xlabel("Pseudotime (re-ranked within PT → aPT)",fontsize=12); a.set_ylabel("Fraction of cells (per-donor mean ± 95% CI)",fontsize=11)
    a.set_title(title,fontweight="bold",loc="left"); a.legend(fontsize=9,frameon=False); a.set_xlim(0,100)
pd.concat(tabs).to_csv(f"{A}/scmega/outputs/panel_b_donor_densities.csv",index=False)
fig.tight_layout()
fig.savefig(f"{O}/panel_a_b.png",dpi=160,bbox_inches="tight"); fig.savefig(f"{O}/panel_a_b.pdf",bbox_inches="tight")
print("wrote scmega/figures_dc1rank/panel_a_b.png/.pdf",flush=True)

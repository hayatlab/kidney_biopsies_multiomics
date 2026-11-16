#!/usr/bin/env python3
"""Is dPT's distal-nephron signal real, ambient, doublets, or misannotation?

Trigger: along the joint pseudotime the top 'transient' genes all peak at ~0.65, where dPT sits,
and they are distal-nephron markers (UMOD, SLC12A1, PTGER3 = TAL; SLC12A3 = DCT; HSD11B2,
CDH16 = CNT/PC). If dPT carries distal transcripts, 'dPT is an intermediate between PT and aPT'
could partly reflect contamination rather than a PT state.

Per level2 (PT, dPT, aPT; TAL and DCT as positive references): % cells expressing distal and PT
markers, library size, genes per cell. Within dPT: does distal signal fall with library size
(ambient pattern) or rise with it (doublet pattern)? What fraction is distal-high AND PT-low
(misannotated distal cells) vs distal-high AND PT-high (mixed)?
X is already ambient-corrected (non-integral): distal signal here is RESIDUAL after correction.
Human hg38 | HGNC | seed 42
"""
from kbpaths import KB_ANALYSIS, KB_DATA  # repository path config (config/kbpaths.py)
import anndata as ad, scanpy as sc, numpy as np, pandas as pd
from scipy.stats import spearmanr
SEED=42; rng=np.random.default_rng(SEED)
A=f"{KB_ANALYSIS}/aki_c5_exclusion_20261022/atac_rerun/ATAC"
RNA=f"{KB_DATA}/charlotte_new_data/adata_20260717.h5ad"
DIST=["UMOD","SLC12A1","PTGER3","SLC12A3","HSD11B2","CDH16","MECOM"]
PTM=["CUBN","LRP2","SLC34A1","SLC22A6","ALDOB","SLC13A3"]
REF_N=15000
a=ad.read_h5ad(RNA,backed="r"); l2=a.obs["level2"].astype(str).to_numpy()
parts=[np.where(l2==t)[0] for t in ["PT","dPT","aPT"]]
for t in ["TAL","DCT"]:
    w=np.where(l2==t)[0]; parts.append(rng.choice(w,size=min(REF_N,len(w)),replace=False))
idx=np.sort(np.concatenate(parts)); a=a[idx].to_memory()
lv=a.obs["level2"].astype(str).to_numpy(); a.X=a.X.astype(np.float32)
tot=np.asarray(a.X.sum(1)).ravel(); ng=np.asarray((a.X>0).sum(1)).ravel()
sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
def g(x):
    v=a[:,x].X; return np.asarray(v.toarray() if hasattr(v,"toarray") else v).ravel()
D={x:g(x) for x in DIST if x in a.var_names}; P={x:g(x) for x in PTM if x in a.var_names}
print("distal markers:",list(D),"| PT markers:",list(P),flush=True)
ds=np.vstack(list(D.values())).mean(0); ps=np.vstack(list(P.values())).mean(0)
rows=[]
for t in ["PT","dPT","aPT","TAL","DCT"]:
    m=lv==t
    r={"level2":t,"n":int(m.sum()),"median_total_counts":float(np.median(tot[m])),
       "median_genes":float(np.median(ng[m])),"mean_distal_score":float(ds[m].mean()),
       "mean_PT_score":float(ps[m].mean())}
    for x,v in {**D,**P}.items(): r[f"pct_{x}"]=round(100*float((v[m]>0).mean()),1)
    rows.append(r)
T=pd.DataFrame(rows).set_index("level2"); T.to_csv(f"{A}/tables/dPT_distal_signal_check.csv")
print("\n"+T.T.round(3).to_string(),flush=True)
m=lv=="dPT"
print(f"\nwithin dPT: rho(distal score, total counts) {spearmanr(ds[m],tot[m]).statistic:+.3f}"
      f"   (negative = ambient-like; positive = doublet-like)",flush=True)
print(f"within dPT: rho(distal score, PT score)     {spearmanr(ds[m],ps[m]).statistic:+.3f}",flush=True)
dh=np.quantile(ds[np.isin(lv,["TAL","DCT"])],0.25); pl=np.quantile(ps[lv=="PT"],0.25)
print(f"\nthresholds: distal-high >= 25th pct of TAL/DCT distal score ({dh:.3f}); "
      f"PT-low < 25th pct of PT PT-score ({pl:.3f})",flush=True)
Q=[]
for t in ["PT","dPT","aPT","TAL","DCT"]:
    mm=lv==t; hi=ds[mm]>=dh; lo=ps[mm]<pl
    Q.append({"level2":t,"pct_distal_high":round(100*hi.mean(),1),
              "pct_distal_high_PT_low":round(100*(hi&lo).mean(),1),
              "pct_distal_high_PT_high":round(100*(hi&~lo).mean(),1)})
Q=pd.DataFrame(Q).set_index("level2"); Q.to_csv(f"{A}/tables/dPT_distal_high_fractions.csv")
print("\n"+Q.to_string(),flush=True)
pc=pd.DataFrame({"cellid":a.obs_names,"level2":lv,"distal_score":ds,"PT_score":ps,
                  "total_counts":tot,"n_genes":ng,"distal_high":ds>=dh,"PT_low":ps<pl})
pc["dPT_subset"]=np.where(lv!="dPT","",np.where(pc.distal_high&pc.PT_low,"distal_like","PT_lineage"))
pc[np.isin(lv,["PT","dPT","aPT"])].to_csv(f"{A}/tables/dPT_percell_scores.csv",index=False)
print("wrote tables/dPT_percell_scores.csv | dPT split:",pc.loc[lv=="dPT","dPT_subset"].value_counts().to_dict(),flush=True)
dec=pd.qcut(tot[m],10,labels=False)
print("\nwithin dPT, mean distal score by library-size decile (1 = smallest):",flush=True)
print(pd.Series(ds[m]).groupby(dec).mean().round(4).rename(lambda i:i+1).to_string(),flush=True)
print("done",flush=True)

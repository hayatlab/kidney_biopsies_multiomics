"""Step 17 - recover glomerular PEC cells that were mislabelled aPT.

Rationale
---------
scANVI level-1 transfer placed a set of parietal-epithelial (PEC) cells into the
PT lineage, so the level2 subclustering (step 09) called them aPT. These cells are
NOT proximal tubule: they carry no PT identity (MIOX/SLC22A6 ~ 0) and, for a subset,
sit physically in the glomerular niche with a positive PEC marker signature.

We re-embed the injury+parietal compartment (pred_level2 in {aPT,dPT,PEC}) exactly
like step 09's dPT recovery, then relabel to PEC ONLY the currently-aPT cells that
fall in a subcluster meeting THREE INDEPENDENT criteria - so the call rests on
spatial + marker evidence, not the (circular) prior label:
    (i)  spatially glomerular : in_glom fraction >= IN_GLOM_MIN   (aPT baseline ~0)
    (ii) PEC identity         : mean CFH        >= CFH_MIN         (aPT baseline ~0.7)
    (iii) prior-PEC concordant: prior-PEC frac  >= PRIORPEC_MIN

The tubulointerstitial CD24+ injured cells that carried an old "PEC" label but are
not glomerular (0% in_glom, CFH-low) are deliberately KEPT as aPT.

Updates obs['level2_refined'] + obs['subcluster_level2'] in outputs/xenium_labeled.h5ad
(a dated backup must be made BEFORE running). Writes a per-subcluster decision table
and a validation figure. seed 42.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import config as C

# ---- recovery criteria (subcluster-level; independent evidence axes) ----
COMPARTMENT = ["aPT", "dPT", "PEC"]      # pred_level2 pool to re-embed
IN_GLOM_MIN = 0.25                       # glomerular niche fraction (canonical PEC ~0.52, aPT ~0.00)
CFH_MIN     = 2.0                        # PEC marker (canonical PEC ~4.6, aPT ~0.7)
PRIORPEC_MIN = 0.50                      # concordance with prior PEC annotation
RES = 1.0

PRIOR_PEC = lambda predic, manual: (predic == "13_PEC") | (manual == "PEC/Glom")


def in_glom_bool(s: pd.Series) -> pd.Series:
    return s.astype(str).isin(["True", "1", "1.0"])


def main() -> int:
    a = sc.read_h5ad(C.XENIUM_LABELED)
    obs = a.obs
    assert "pred_level2" in obs and "level2_refined" in obs and "in_glom" in obs
    cats = list(obs["level2_refined"].cat.categories) if hasattr(obs["level2_refined"], "cat") \
        else sorted(obs["level2_refined"].astype(str).unique())

    l2_before = obs["level2_refined"].astype(str).copy()
    print("level2_refined BEFORE:")
    print(l2_before.value_counts().to_string())

    # ---- re-embed the injury+parietal compartment (same recipe as step 09) ----
    p2 = obs["pred_level2"].astype(str)
    sub = a[p2.isin(COMPARTMENT).to_numpy()].copy()
    print(f"\ncompartment (pred_level2 in {COMPARTMENT}): {sub.n_obs} cells")
    sub.X = sub.layers["counts"].copy().astype(np.float32)
    sc.pp.normalize_total(sub, target_sum=1e4); sc.pp.log1p(sub)
    sc.pp.pca(sub, n_comps=30, random_state=C.SEED)
    sc.pp.neighbors(sub, use_rep="X_pca", random_state=C.SEED)
    sc.tl.leiden(sub, resolution=RES, key_added="cl", flavor="igraph",
                 n_iterations=2, directed=False, random_state=C.SEED)

    predic = sub.obs["predic"].astype(str); manual = sub.obs["cell_type_manual"].astype(str)
    prior_pec = PRIOR_PEC(predic, manual).to_numpy()
    ig = in_glom_bool(sub.obs["in_glom"]).to_numpy()
    cur_apt = (sub.obs["level2_refined"].astype(str) == "aPT").to_numpy()
    cfh = np.asarray(sub[:, "CFH"].X.todense()).ravel() if hasattr(sub[:, "CFH"].X, "todense") \
        else np.asarray(sub[:, "CFH"].X).ravel()

    df = pd.DataFrame({"cl": sub.obs["cl"].astype(str).to_numpy(),
                       "prior_pec": prior_pec, "ig": ig, "cur_apt": cur_apt, "CFH": cfh},
                      index=sub.obs.index)
    g = df.groupby("cl", observed=True)
    stat = pd.DataFrame({"n": g.size(),
                         "in_glom_frac": g["ig"].mean().round(3),
                         "CFH_mean": g["CFH"].mean().round(2),
                         "prior_pec_frac": g["prior_pec"].mean().round(3),
                         "n_cur_aPT": g["cur_apt"].sum().astype(int)})
    stat["qualifies_PEC"] = ((stat["in_glom_frac"] >= IN_GLOM_MIN) &
                             (stat["CFH_mean"] >= CFH_MIN) &
                             (stat["prior_pec_frac"] >= PRIORPEC_MIN))
    stat = stat.sort_values("qualifies_PEC", ascending=False)
    print("\n=== per-subcluster decision (criteria: in_glom>=%.2f & CFH>=%.1f & priorPEC>=%.2f) ==="
          % (IN_GLOM_MIN, CFH_MIN, PRIORPEC_MIN))
    print(stat.to_string())
    C.ensure_dirs()
    stat.to_csv(C.MET / "pec_recovery_subcluster_decision.csv")

    qual = stat.index[stat["qualifies_PEC"]].tolist()
    recovered = df.index[df["cl"].isin(qual) & df["cur_apt"]]
    print(f"\nqualifying PEC subcluster(s): {qual}")
    print(f"currently-aPT cells recovered -> PEC: {len(recovered)}")

    # ---- write updated labels ----
    l2 = l2_before.copy()
    l2.loc[recovered] = "PEC"
    a.obs["level2_refined"] = pd.Categorical(l2, categories=cats)
    scl = a.obs["subcluster_level2"].astype(str).copy()
    scl.loc[recovered] = "aPTrec_PEC"
    a.obs["subcluster_level2"] = pd.Categorical(scl)

    print("\nlevel2_refined AFTER:")
    print(a.obs["level2_refined"].value_counts().to_string())
    n_nan = int(a.obs["level2_refined"].isna().sum())
    print(f"\n%% NaN in level2_refined after: {n_nan} ({100*n_nan/a.n_obs:.3f}%)")
    delta = (a.obs["level2_refined"].value_counts() - l2_before.value_counts()).dropna()
    print("delta:", delta[delta != 0].to_dict())

    # ---- validation figure: every subcluster placed by glomerular fraction vs CFH ----
    fig, ax = plt.subplots(figsize=(7.4, 5.8))
    sc_ = ax.scatter(stat["in_glom_frac"] * 100, stat["CFH_mean"], s=30 + stat["n"] / 6,
                     c=stat["prior_pec_frac"] * 100, cmap="viridis", vmin=0, vmax=100,
                     edgecolors="k", lw=0.6, zorder=3)
    for cl, r in stat.iterrows():
        ax.annotate(cl, (r["in_glom_frac"] * 100, r["CFH_mean"]), fontsize=8, ha="center",
                    va="center", fontweight="bold",
                    color="white" if r["prior_pec_frac"] < 0.55 else "black")
    ax.axvline(IN_GLOM_MIN * 100, ls="--", c="0.6", lw=1)
    ax.axhline(CFH_MIN, ls="--", c="0.6", lw=1)
    ax.text(60, ax.get_ylim()[1] * 0.93, "PEC recovery zone\n(glomerular + CFH+)",
            fontsize=9, color="#2ca02c", ha="center", fontweight="bold")
    ax.set_xlabel("% cells in glomerular niche (in_glom)")
    ax.set_ylabel("mean CFH (PEC marker, log-norm)")
    ax.set_title("PEC recovery from aPT: spatial location vs PEC identity\n"
                 "(size=n, colour=%% prior-PEC; recovered subcluster in the zone)")
    cb = plt.colorbar(sc_, ax=ax, fraction=0.045, pad=0.02); cb.set_label("% prior-PEC labelled")
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(C.FIG / "pec_recovery_decision.pdf", bbox_inches="tight")
    fig.savefig(C.FIG / "pec_recovery_decision.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    a.write_h5ad(C.XENIUM_LABELED)
    print(f"\nwrote updated level2_refined + subcluster_level2 -> {C.XENIUM_LABELED}")
    print(f"decision table -> {C.MET/'pec_recovery_subcluster_decision.csv'}")
    print(f"figure -> {C.FIG/'pec_recovery_decision.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

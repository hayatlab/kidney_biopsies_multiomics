#!/usr/bin/env python
"""Supervised pyVIA trajectories on the leiden_0.8 PEC subclustering (colleague request).

Four cluster configurations x three VIA input embeddings (umap / diffmap / phate):

  A_0to9_7cl          keep {0,2,4,8,9,11,12}   root  0            -> terminal 9
  B_0p2to9_7cl        keep {0,2,4,8,9,11,12}   root  0+2          -> terminal 9
  C_0to9_5cl          keep {0,2,4,9,11}        root  0            -> terminal 9
  D_0p2p4p11to9_5cl   keep {0,2,4,9,11}        root  0+2+4+11     -> terminal 9

Recipe is the one finalised in trajectory_new_embedding.py / pec_via_root0_finalplots.py:
subset-specific X_pca_harmony -> neighbors(15) -> diffmap(15); VIA with
edgepruning_clustering_resolution 0.15, cluster_graph_pruning 0.5, knn 15,
resolution_parameter 1.5, dataset 'group', seed 4, root_user + user_defined_terminal_group.
All three embeddings (UMAP, diffmap, PHATE) are computed ON THE SUBSET from the same
X_pca_harmony rep, so the embedding comparison is not confounded by the excluded clusters
(user decision 2026-07-23; the earlier leiden04 `umapbase` runs sliced the global UMAP).

Two implementation notes, both forced by the pyVIA source:
  * `root_user` does NOT support multiple simultaneous roots -- core.py:2251 loops the list
    and keeps the LAST match per graph component. Multi-root configs (B, D) are therefore
    implemented by merging the root clusters into one true_label group "ROOT", so VIA picks
    the least-connected VIA-cluster inside the union. Original per-cluster labels are kept
    in obs["cl"] for every plot.
  * root matching is SUBSTRING based (core.py:1744, `str(root_user) in str(majority_truth)`),
    so bare numeric labels would make root '2' match a cluster whose majority is '12' -- and
    12 is in the keep set. Labels are therefore prefixed: c0/c2/c4/c8/c9/c11/c12 (collision
    free for these sets). Terminal-group matching is exact (core.py:815) and needs no fix.

DPT is embedding independent (neighbors+diffmap on X_pca_harmony), so it is computed ONCE per
configuration and reported alongside all three VIA runs as a deterministic reference.

seed 42 (4 for VIA, as in the finalised run). VIA is not seed-deterministic.
Env: ViaEnv (needs LD_LIBRARY_PATH=$CONDA_PREFIX/lib).

Outputs: trajectories_leiden08/<config>/<embedding>/ plus <config>/_compare/.
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import os
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.sparse import issparse
from scipy.stats import spearmanr
import scanpy as sc
try: import pyVIA.core as via
except Exception: import via

np.random.seed(42)
DIR = f"{KB_ANALYSIS}/pec_resubcluster_20260717"
SRC = f"{DIR}/new_pec_annotation/pec_leiden_sweep.h5ad"
OUT = f"{DIR}/trajectories_leiden08"
SCORE = "crescent_lesion_score"
LABEL = "leiden_0.8"
SEED = 42
NJOBS = int(os.environ.get("OMP_NUM_THREADS", "8"))

KEEP7 = ["0", "2", "4", "8", "9", "11", "12"]
KEEP5 = ["0", "2", "4", "9", "11"]
KEEP8 = ["0", "2", "4", "8", "9", "10", "11", "13"]

CONFIGS = [
    dict(tag="A_0to9_7cl",        keep=KEEP7, root=["0"],                 terminal="9"),
    dict(tag="B_0p2to9_7cl",      keep=KEEP7, root=["0", "2"],            terminal="9"),
    dict(tag="C_0to9_5cl",        keep=KEEP5, root=["0"],                 terminal="9"),
    dict(tag="D_0p2p4p11to9_5cl", keep=KEEP5, root=["0", "2", "4", "11"], terminal="9"),
    # follow-up 2026-07-23: 8-cluster set, drops c12 vs A, adds c10 (n=38, endothelial
    # markers EMCN/LDB2 -- likely contaminant island) and c13 (n=359, retained_PEC)
    dict(tag="E_0to9_8cl",        keep=KEEP8, root=["0"],                 terminal="9"),
]
EMBEDDINGS = ["umap", "diffmap", "phate"]

# fixed cluster colours so the same cluster is the same colour in every config/embedding
CLCOL = {"c0": "#1f77b4", "c2": "#ff7f0e", "c4": "#2ca02c", "c8": "#d62728",
         "c9": "#9467bd", "c11": "#17becf", "c12": "#8c564b",
         "c10": "#e377c2", "c13": "#bcbd22"}

# gene-trend panels (filtered to genes present at runtime) -- same panels as the 0.4 runs
UP     = ["CD44", "LAMC2", "ITGB1", "PSD3"]            # crescent / activation program
DOWN   = ["SLC4A11", "CDH1", "NDRG1", "CFH"]           # PEC-identity / tubular genes
CUSTOM = ["CFH", "ITGB3", "SLC4A11", "PTPRQ", "ADAMTS19", "DDN"]  # colleague's panel


def _dense(x):
    return x.toarray() if issparse(x) else np.asarray(x)


def markers_of(label=LABEL):
    """cluster -> 'marker1/marker2' from the sweep cluster_summary (for legends)."""
    rstr = label.split("_", 1)[1]
    f = f"{DIR}/new_pec_annotation/res_{rstr}/cluster_summary_res{rstr}.csv"
    if not os.path.exists(f):
        return {}
    d = pd.read_csv(f)
    return {str(r["cluster"]): "/".join(str(r["top_markers"]).split(",")[:2]) for _, r in d.iterrows()}


ANN = markers_of()


def cl_label(c):
    """'0' -> 'c0  SLC8A1/NPHS1'"""
    return f"c{c}  {ANN.get(c, '')}".strip()


# ------------------------------------------------------------------ preparation
def prepare(cfg, P):
    """Subset, build the shared graph, and compute all three subset embeddings + DPT."""
    keep, root = cfg["keep"], cfg["root"]
    a = sc.read_h5ad(SRC)
    a = a[a.obs[LABEL].astype(str).isin(keep)].copy()
    a.obs["cl"] = ("c" + a.obs[LABEL].astype(str)).astype("category")

    # VIA true_label: prefixed cluster, with the root clusters merged when there are several
    is_root = a.obs[LABEL].astype(str).isin(root).to_numpy()
    if len(root) > 1:
        root_group = "ROOT"
        vl = a.obs["cl"].astype(str).to_numpy().copy()
        vl[is_root] = root_group
        a.obs["via_label"] = vl
        P(f"multi-root config: merged {['c'+r for r in root]} -> '{root_group}' "
          f"({int(is_root.sum())} cells) for VIA root finding")
    else:
        root_group = f"c{root[0]}"
        a.obs["via_label"] = a.obs["cl"].astype(str)
    term_group = f"c{cfg['terminal']}"

    P(f"cells={a.n_obs}; per-cluster n + crescent ULM mean:")
    for c in keep:
        m = (a.obs[LABEL].astype(str) == c).to_numpy()
        P(f"  c{c} ({ANN.get(c,'')}): n={int(m.sum())}  crescent={a.obs[SCORE][m].mean():.3f}")

    a.obsm["X_umap_global"] = a.obsm["X_umap"].copy()          # keep the sliced global UMAP

    # --- UMAP layout graph: n_neighbors=30 ---
    # At the default n_neighbors=15 the 7-cluster 2D UMAP layout throws off a spurious 57-cell
    # island (56x c4 + 1x c2, 1.6% of cells). VIA then builds a 2-component graph, the root c0 is
    # absent from the island so VIA sets a DUMMY root, and the terminal c9 lives in the other
    # component -> 0 terminal lineages, i.e. not the requested supervised trajectory.
    # The island is a 2D layout artifact, not structure: the same cells are a single component
    # under diffmap(15D), PHATE(2D) and X_pca_harmony(30D) itself. n_neighbors=30 (min_dist at
    # the 0.5 default) yields one component for BOTH the 5- and 7-cluster subsets. Applied
    # uniformly. VIA's own knn stays 15; only the UMAP layout graph differs.
    sc.pp.neighbors(a, use_rep="X_pca_harmony", n_neighbors=30, random_state=SEED)
    sc.tl.umap(a, random_state=SEED)
    a.obsm["X_umap_subset"] = a.obsm["X_umap"].copy()
    P("recomputed UMAP on the subset (layout graph n_neighbors=30; "
      "global sliced UMAP retained as X_umap_global)")

    # --- analysis graph (unchanged recipe): n_neighbors=15 -> diffmap -> DPT ---
    sc.pp.neighbors(a, use_rep="X_pca_harmony", n_neighbors=15, random_state=SEED)
    sc.tl.diffmap(a, n_comps=15)

    try:
        import phate
        ph = phate.PHATE(n_components=2, knn=15, random_state=SEED, n_jobs=NJOBS, verbose=0)
        a.obsm["X_phate"] = ph.fit_transform(a.obsm["X_pca_harmony"])
        P(f"PHATE ok (knn=15, n_jobs={NJOBS}, seed={SEED})")
    except Exception as e:
        P(f"!! PHATE failed ({type(e).__name__}: {e}) -- phate runs will be skipped")

    # connectivity check on every 2D embedding VIA will consume: a split graph silently
    # degrades the supervised run into "dummy root + no terminal" (see the UMAP note above)
    from scipy.sparse.csgraph import connected_components
    for nm, k in [("umap", "X_umap_subset"), ("phate", "X_phate")]:
        if k not in a.obsm: continue
        b = sc.AnnData(np.asarray(a.obsm[k], dtype=np.float32))
        sc.pp.neighbors(b, use_rep="X", n_neighbors=15, random_state=SEED)
        n, lab = connected_components(b.obsp["connectivities"], directed=False)
        sizes = pd.Series(lab).value_counts().sort_values(ascending=False).tolist()
        P(f"connectivity {nm} 2D @knn15: {n} component(s) sizes={sizes[:4]}"
          + ("" if n == 1 else "   *** SPLIT -> VIA root/terminal may be unreachable ***"))

    # ---- deterministic DPT, once per config (embedding independent) ----
    root_cells = np.where(is_root)[0]
    H = a.obsm["X_pca_harmony"]
    cent = H[root_cells].mean(0)
    iroot = int(root_cells[np.argmin(((H[root_cells] - cent) ** 2).sum(1))])
    a.uns["iroot"] = iroot
    sc.tl.dpt(a, n_dcs=15)
    dpt = a.obs["dpt_pseudotime"].to_numpy()
    dpt = np.clip(dpt, np.percentile(dpt, 0), np.percentile(dpt, 99.5))
    a.obs["pseudotime_dpt"] = (dpt - dpt.min()) / (dpt.max() - dpt.min() + 1e-9)
    rho_d, p_d = spearmanr(a.obs["pseudotime_dpt"], a.obs[SCORE])
    P(f"DPT (shared, iroot cell {iroot} in {'+'.join('c'+r for r in root)}): "
      f"Spearman(dpt, crescent) = {rho_d:+.3f} (p={p_d:.1e})")
    return a, root_group, term_group, rho_d


# ------------------------------------------------------------------ one VIA run
def run_embedding(a, cfg, emb, root_group, term_group, rho_dpt):
    keep, terminal = cfg["keep"], cfg["terminal"]
    od = f"{OUT}/{cfg['tag']}/{emb}"; os.makedirs(od, exist_ok=True)
    log = []
    def P(*m):
        s = " ".join(str(x) for x in m); print(s, flush=True); log.append(s)

    key = {"umap": "X_umap_subset", "diffmap": "X_diffmap", "phate": "X_phate"}[emb]
    if key not in a.obsm:
        P(f"!! {emb}: {key} missing -- skipped"); return None
    Vdata = a.obsm[key]
    P(f"\n----- {cfg['tag']} / {emb} -----")
    P(f"VIA base embedding: {emb} ({key}, shape {Vdata.shape}); "
      f"root_user=['{root_group}'] terminal=['{term_group}']")

    true_label = list(a.obs["via_label"].astype(str))
    # NB: pyVIA's preserve_disconnected=False is NOT a usable remedy for a split graph here --
    # connect_all_components() bridges via the MST, but the MST of a natively disconnected
    # adjacency is a forest with no cross-component edges, so it raises
    # "zero-size array to reduction operation minimum" (utils_via.py:195). Connectivity is
    # instead ensured upstream, in prepare(), by the UMAP layout graph. VIA params unchanged.
    v = via.VIA(data=Vdata, true_label=true_label, edgepruning_clustering_resolution=0.15,
                cluster_graph_pruning=0.5, knn=15, root_user=[root_group], resolution_parameter=1.5,
                dataset='group', random_seed=4, memory=5, do_compute_embedding=False,
                user_defined_terminal_group=[term_group])
    v.run_VIA()

    # ---- validity instrumentation: a supervised run is only meaningful if VIA actually
    # resolved the requested root AND at least one terminal lineage in ONE component ----
    def _aslist(x):     # v.root / terminal_clusters / comp_labels may be list OR ndarray
        return [] if x is None else np.asarray(x).ravel().tolist()
    n_term = len(_aslist(getattr(v, "terminal_clusters", None)))
    cc = _aslist(getattr(v, "connected_comp_labels", None))
    n_comp = len(set(cc)) if cc else 1
    labs, tl = np.asarray(v.labels), np.asarray(true_label)
    root_major = []
    for r in _aslist(getattr(v, "root", None)):
        mm = pd.Series(tl[labs == r]).mode()
        root_major.append(mm.iloc[0] if len(mm) else "?")
    root_ok = bool(root_major) and all(m == root_group for m in root_major)
    valid = bool(n_term > 0 and root_ok)
    P(f"VALIDITY: components={n_comp}  terminal_lineages={n_term}  "
      f"root_cluster_majority={root_major} (want all '{root_group}')  -> "
      f"{'OK' if valid else '*** INVALID: not the requested supervised trajectory ***'}")
    # for merged-root configs, root_major is always 'ROOT' by construction and hides WHERE inside
    # the union VIA anchored -- report the original leiden_0.8 composition of the root VIA-cluster
    cl_arr = a.obs["cl"].astype(str).to_numpy()
    root_orig = []
    for r in _aslist(getattr(v, "root", None)):
        vc = pd.Series(cl_arr[labs == r]).value_counts()
        root_orig.append("/".join(f"{k}:{int(n)}" for k, n in vc.head(3).items()))
    P(f"root VIA-cluster composition in original clusters: {root_orig}")

    pt = np.asarray(v.single_cell_pt_markov, dtype=float)
    pt = (pt - pt.min()) / (pt.max() - pt.min() + 1e-9)
    a.obs["pseudotime_via"] = pt
    dpt = a.obs["pseudotime_dpt"].to_numpy()

    rho_v, p_v = spearmanr(pt, a.obs[SCORE])
    rho_vd, _ = spearmanr(pt, dpt)
    P(f"VIA: Spearman(pt, crescent) = {rho_v:+.3f} (p={p_v:.1e})   |   "
      f"VIA-vs-DPT concordance = {rho_vd:+.3f}   |   DPT-vs-crescent = {rho_dpt:+.3f}")

    # per-cluster ordering + a cluster-level monotonicity check
    P("per-cluster median pseudotime (VIA / DPT) and mean crescent:")
    med_v, med_d, cres = [], [], []
    for c in keep:
        m = (a.obs[LABEL].astype(str) == c).to_numpy()
        med_v.append(float(np.median(pt[m]))); med_d.append(float(np.median(dpt[m])))
        cres.append(float(a.obs[SCORE][m].mean()))
        P(f"  c{c} ({ANN.get(c,'')}): VIA={med_v[-1]:.3f}  DPT={med_d[-1]:.3f}  crescent={cres[-1]:.2f}")
    rho_clu, _ = spearmanr(med_v, cres) if len(keep) > 2 else (np.nan, np.nan)
    rank_root = int(np.argsort(np.argsort(med_v))[keep.index(cfg["root"][0])])
    rank_term = int(np.argsort(np.argsort(med_v))[keep.index(terminal)])
    # strict-max check: argsort ranks break ties arbitrarily, so a saturated VIA pseudotime
    # (several clusters pinned at 1.0) can give the terminal top rank by luck
    ti = keep.index(terminal)
    others = [med_v[i] for i in range(len(keep)) if i != ti]
    term_strict = bool(med_v[ti] > max(others))
    P(f"cluster-level Spearman(median VIA pt, mean crescent) = {rho_clu:+.3f}")
    P(f"root c{cfg['root'][0]} pseudotime rank {rank_root}/{len(keep)-1} (want 0); "
      f"terminal c{terminal} rank {rank_term}/{len(keep)-1} (want {len(keep)-1}); "
      f"terminal median strictly highest = {term_strict}"
      + ("" if term_strict else "  [tie -> rank is not decisive]"))

    # ---------------- figures ----------------
    E = Vdata[:, :2] if emb != "diffmap" else a.obsm["X_diffmap"][:, 1:3]
    ax_names = {"umap": ("UMAP1", "UMAP2"), "phate": ("PHATE1", "PHATE2"),
                "diffmap": ("DC1", "DC2")}[emb]
    U = a.obsm["X_umap_subset"]
    cl = a.obs["cl"].astype(str).to_numpy()

    def scatter_panel(ax, XY, colour_by, cmap=None, title="", lab=("", "")):
        if colour_by is None:                       # categorical clusters
            for c in keep:
                m = cl == f"c{c}"
                ax.scatter(XY[m, 0], XY[m, 1], s=6, lw=0, color=CLCOL[f"c{c}"], label=cl_label(c))
            ax.legend(markerscale=2.5, fontsize=8, frameon=False)
            h = None
        else:
            h = ax.scatter(XY[:, 0], XY[:, 1], s=6, lw=0, c=colour_by, cmap=cmap)
        ax.set_title(title); ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel(lab[0]); ax.set_ylabel(lab[1])
        return h

    # (1) the VIA base embedding: clusters / VIA pt / DPT / crescent
    fig, ax = plt.subplots(1, 4, figsize=(27, 6.2)); fig.patch.set_facecolor("white")
    scatter_panel(ax[0], E, None, title=f"{LABEL} clusters", lab=ax_names)
    h1 = scatter_panel(ax[1], E, pt, "viridis", f"VIA pseudotime ({root_group}->{term_group})", ax_names)
    fig.colorbar(h1, ax=ax[1], shrink=.7)
    h2 = scatter_panel(ax[2], E, dpt, "viridis", "DPT pseudotime (deterministic ref)", ax_names)
    fig.colorbar(h2, ax=ax[2], shrink=.7)
    h3 = scatter_panel(ax[3], E, a.obs[SCORE], "magma", "crescent ULM score", ax_names)
    fig.colorbar(h3, ax=ax[3], shrink=.7)
    fig.suptitle(f"{cfg['tag']} — VIA on {emb.upper()}", fontsize=14)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{od}/{emb}_cluster_pt_crescent.{ext}", dpi=150 if ext == "png" else None,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # (2) the same pseudotime shown on the common subset UMAP (comparable across embeddings)
    fig, ax = plt.subplots(1, 3, figsize=(21, 6.2)); fig.patch.set_facecolor("white")
    scatter_panel(ax[0], U, None, title=f"{LABEL} clusters", lab=("UMAP1", "UMAP2"))
    h1 = scatter_panel(ax[1], U, pt, "viridis", f"VIA pseudotime (base={emb})", ("UMAP1", "UMAP2"))
    fig.colorbar(h1, ax=ax[1], shrink=.7)
    h3 = scatter_panel(ax[2], U, a.obs[SCORE], "magma", "crescent ULM score", ("UMAP1", "UMAP2"))
    fig.colorbar(h3, ax=ax[2], shrink=.7)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{od}/onUMAP_cluster_pt_crescent.{ext}", dpi=150 if ext == "png" else None,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # (3) VIA lineage curves on the base embedding
    try:
        v.embedding = E
        out = via.plot_trajectory_curves(via_object=v, embedding=E, draw_all_curves=False)
        figv = out[0] if isinstance(out, (list, tuple)) else plt.gcf()
        figv.patch.set_facecolor("white")
        for axx in figv.get_axes(): axx.set_facecolor("white")
        for ext in ("png", "pdf"):
            figv.savefig(f"{od}/via_lineage_curves.{ext}", dpi=150 if ext == "png" else None,
                         bbox_inches="tight", facecolor="white")
        plt.close(figv); P("saved via.plot_trajectory_curves")
    except Exception as e:
        P(f"!! plot_trajectory_curves failed: {type(e).__name__}: {e}")

    # (4) crescent score vs pseudotime
    fig, ax = plt.subplots(1, 2, figsize=(15, 5.6)); fig.patch.set_facecolor("white")
    for j, (ptv, name, rho) in enumerate([(pt, "VIA", rho_v), (dpt, "DPT", rho_dpt)]):
        order = np.argsort(ptv)
        ax[j].scatter(ptv, a.obs[SCORE], s=7, lw=0, alpha=0.35, color="0.5")
        w = 200; ys = a.obs[SCORE].to_numpy()[order]
        ax[j].plot(ptv[order], pd.Series(ys).rolling(w, center=True, min_periods=w // 4).mean(),
                   color="crimson", lw=2.2)
        ax[j].set_xlabel(f"{name} pseudotime"); ax[j].set_ylabel("crescent ULM score")
        ax[j].set_title(f"{name}: Spearman = {rho:+.3f}")
        for sp in ["top", "right"]: ax[j].spines[sp].set_visible(False)
    fig.suptitle(f"{cfg['tag']} / {emb}: crescent signature vs pseudotime", fontsize=13)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{od}/pseudotime_vs_crescent.{ext}", dpi=150 if ext == "png" else None,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # (5) per-cluster pseudotime boxplot, ordered by median VIA pt
    order_c = [keep[i] for i in np.argsort(med_v)]
    fig, ax = plt.subplots(figsize=(1.7 * len(keep) + 2.5, 4.8)); fig.patch.set_facecolor("white")
    bp = ax.boxplot([pt[(a.obs[LABEL].astype(str) == c).to_numpy()] for c in order_c],
                    tick_labels=[f"c{c}" for c in order_c], patch_artist=True, showfliers=False)
    for patch, c in zip(bp["boxes"], order_c):
        patch.set_facecolor(CLCOL[f"c{c}"]); patch.set_alpha(0.75)
    ax.set_ylabel("VIA pseudotime")
    ax.set_title(f"{cfg['tag']} / {emb}: pseudotime by cluster "
                 f"(root {root_group} -> terminal {term_group})")
    for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{od}/pseudotime_per_cluster_box.{ext}", dpi=150 if ext == "png" else None,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # (6) gene trends along VIA pseudotime (MAGIC-imputed)
    def gene_trend(genes, name):
        present = [g for g in genes if g in a.var_names]
        if not present:
            P(f"gene-trend {name}: none present"); return
        if n_term == 0:   # via.get_gene_expression does cm.get_cmap(cmap, 0) -> opaque ValueError
            P(f"gene-trend {name}: skipped, VIA resolved 0 terminal lineages"); return
        df = pd.DataFrame(_dense(a[:, present].X), columns=present, index=a.obs_names)
        try:
            dm = v.do_impute(df, magic_steps=3, gene_list=present)
            figg, _ = via.get_gene_expression(via_object=v, gene_exp=dm)
            figg.set_size_inches(14, 7); figg.patch.set_facecolor("white")
            for axx in figg.get_axes(): axx.set_facecolor("white")
            figg.tight_layout()
            for ext in ("png", "pdf"):
                figg.savefig(f"{od}/genetrend_{name}.{ext}", dpi=130 if ext == "png" else None,
                             bbox_inches="tight", facecolor="white")
            plt.close(figg); P(f"gene-trend {name}: {present}")
        except Exception as e:
            P(f"!! gene-trend {name} failed: {type(e).__name__}: {e}")
    gene_trend(UP, "UP"); gene_trend(DOWN, "DOWN"); gene_trend(CUSTOM, "CUSTOM")

    # tables
    pd.DataFrame({"cell": a.obs_names, "cluster": cl,
                  "pseudotime_via": pt, "pseudotime_dpt": dpt,
                  SCORE: a.obs[SCORE].to_numpy()}).to_csv(f"{od}/pseudotime.csv", index=False)
    pd.DataFrame({"cluster": [f"c{c}" for c in keep],
                  "markers": [ANN.get(c, "") for c in keep],
                  "n": [int((a.obs[LABEL].astype(str) == c).sum()) for c in keep],
                  "median_pt_via": np.round(med_v, 4), "median_pt_dpt": np.round(med_d, 4),
                  "mean_crescent": np.round(cres, 3)}
                 ).sort_values("median_pt_via").to_csv(f"{od}/per_cluster_pseudotime.csv", index=False)
    with open(f"{od}/summary.txt", "w") as fh:
        fh.write("\n".join(log) + "\n")

    return dict(config=cfg["tag"], embedding=emb, n_cells=a.n_obs,
                root="+".join("c" + r for r in cfg["root"]), terminal=term_group,
                rho_via_crescent=round(rho_v, 4), rho_dpt_crescent=round(rho_dpt, 4),
                rho_via_vs_dpt=round(rho_vd, 4), rho_cluster_level=round(float(rho_clu), 4),
                root_pt_rank=rank_root, terminal_pt_rank=rank_term,
                terminal_strict_max=term_strict, n_clusters=len(keep),
                via_valid=valid, n_components=n_comp, n_terminal_lineages=n_term,
                root_cluster_majority="|".join(root_major),
                root_cluster_origin="|".join(root_orig), pt_via=pt.copy())


# ------------------------------------------------------------------ per-config driver
def run_config(cfg):
    od = f"{OUT}/{cfg['tag']}"; os.makedirs(od, exist_ok=True)
    log = []
    def P(*m):
        s = " ".join(str(x) for x in m); print(s, flush=True); log.append(s)

    P(f"\n================ {cfg['tag']}: keep={cfg['keep']} "
      f"root={cfg['root']} -> terminal={cfg['terminal']} ================")
    a, root_group, term_group, rho_d = prepare(cfg, P)

    res = []
    for emb in EMBEDDINGS:
        r = run_embedding(a, cfg, emb, root_group, term_group, rho_d)
        if r is not None:
            a.obs[f"pseudotime_via_{emb}"] = r.pop("pt_via")
            res.append(r)

    # ---- cross-embedding comparison figure ----
    if res:
        cd = f"{od}/_compare"; os.makedirs(cd, exist_ok=True)
        embs = [r["embedding"] for r in res]
        U = a.obsm["X_umap_subset"]
        fig, ax = plt.subplots(2, len(embs), figsize=(7 * len(embs), 12.4), squeeze=False)
        fig.patch.set_facecolor("white")
        for j, e in enumerate(embs):
            XY = {"umap": a.obsm["X_umap_subset"], "phate": a.obsm.get("X_phate"),
                  "diffmap": a.obsm["X_diffmap"][:, 1:3]}[e]
            ptj = a.obs[f"pseudotime_via_{e}"].to_numpy()
            h = ax[0][j].scatter(XY[:, 0], XY[:, 1], s=6, lw=0, c=ptj, cmap="viridis")
            fig.colorbar(h, ax=ax[0][j], shrink=.7)
            ax[0][j].set_title(f"{e.upper()} space — VIA pseudotime\n"
                               f"rho(crescent)={res[j]['rho_via_crescent']:+.3f}  "
                               f"cluster-level={res[j]['rho_cluster_level']:+.3f}")
            h = ax[1][j].scatter(U[:, 0], U[:, 1], s=6, lw=0, c=ptj, cmap="viridis")
            fig.colorbar(h, ax=ax[1][j], shrink=.7)
            ax[1][j].set_title(f"same pseudotime on the common subset UMAP (base={e})")
            for i in (0, 1): ax[i][j].set_xticks([]); ax[i][j].set_yticks([])
        fig.suptitle(f"{cfg['tag']}: root {root_group} -> terminal {term_group} "
                     f"({a.n_obs} cells) — embedding comparison", fontsize=15)
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(f"{cd}/embedding_comparison.{ext}", dpi=140 if ext == "png" else None,
                        bbox_inches="tight", facecolor="white")
        plt.close(fig)

        # pairwise agreement between the three embeddings' pseudotimes
        M = pd.DataFrame(index=embs, columns=embs, dtype=float)
        for e1 in embs:
            for e2 in embs:
                M.loc[e1, e2] = spearmanr(a.obs[f"pseudotime_via_{e1}"],
                                          a.obs[f"pseudotime_via_{e2}"])[0]
        M.round(3).to_csv(f"{cd}/pseudotime_agreement_spearman.csv")
        P("\npairwise Spearman between embeddings' VIA pseudotimes:")
        P(M.round(3).to_string())

    a.write_h5ad(f"{od}/subset_all_pseudotimes.h5ad")
    with open(f"{od}/config_summary.txt", "w") as fh:
        fh.write("\n".join(log) + "\n")
    P(f"[done] {cfg['tag']} -> {od}")
    return res


def main(tags=None, embs=None):
    global EMBEDDINGS
    if embs: EMBEDDINGS = embs
    os.makedirs(OUT, exist_ok=True)
    cfgs = [c for c in CONFIGS if (not tags or c["tag"] in tags)]
    rows = [r for cfg in cfgs for r in run_config(cfg)]
    res = pd.DataFrame(rows)
    print("\n==================== SUMMARY ====================")
    print(res.to_string(index=False))
    master = f"{OUT}/trajectory_summary_leiden08.csv"
    if os.path.exists(master):
        old = pd.read_csv(master)
        keyold = old["config"].astype(str) + "|" + old["embedding"].astype(str)
        keynew = set(res["config"].astype(str) + "|" + res["embedding"].astype(str))
        res = pd.concat([old[~keyold.isin(keynew)], res], ignore_index=True)
    res.to_csv(master, index=False)
    print(f"\nwrote {master}")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--tags", default=None, help="comma-separated config tags (default: all 4)")
    ap.add_argument("--embeddings", default=None, help="comma-separated: umap,diffmap,phate")
    args = ap.parse_args()
    main(tags=args.tags.split(",") if args.tags else None,
         embs=args.embeddings.split(",") if args.embeddings else None)

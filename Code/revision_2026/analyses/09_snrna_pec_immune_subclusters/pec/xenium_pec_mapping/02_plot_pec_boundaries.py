"""Per-core + per-glomerulus Xenium cell-boundary plots of the mapped PEC leiden_0.4
subclusters, for visual plausibility review.

For every TMA core with >= MIN_PEC mapped PEC cells we render:
  1. per_core/<TMA>.{png,pdf}   full-core overview - where the PEC cells sit in the
                                 whole core (grey tissue context, podocytes pale blue).
  2. per_glomerulus/<TMA>.png   montage: the PEC cells are spatially clustered into
                                 glomeruli (DBSCAN on centroids); each glomerulus is
                                 shown zoomed so the leiden_0.4 cluster colours are
                                 obvious and you can judge whether they ring Bowman's
                                 capsule around the podocyte tuft.

Cells are drawn as their real Xenium segmentation polygons (matplotlib PolyCollection,
built straight from the boundary parquets - no geopandas). y-axis inverted to match
Xenium image orientation. PEC polygons are drawn on top of the grey context.

Inputs:
  ./pec_leiden04_predictions.csv                 (sample, cell_id -> pred_leiden04, conf)
  celltype_concordance/outputs/xenium_labeled.h5ad   (all cells: cell_id,sample,TMA_ID,level2_refined)
  $KB_XENIUM_RAW/output-<sample>/{cell,nucleus}_boundaries.parquet

Outputs (./boundary_plots_leiden04/): per_core/, per_glomerulus/, palette.csv,
  per_core_cluster_counts.csv, README.txt
"""
from __future__ import annotations
from kbpaths import KB_ANALYSIS, KB_XENIUM_RAW  # repository path config (config/kbpaths.py)
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.cluster import DBSCAN

from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
HERE = Path(f"{KB_ANALYSIS}/pec_resubcluster_20260717/xenium_pec_mapping")  # analysis workspace (outputs live here)
XENIUM_LABELED = Path(f"{KB_ANALYSIS}/celltype_concordance/"
                      "outputs/xenium_labeled.h5ad")
XENIUM_BASE = Path(f"{KB_XENIUM_RAW}")

SAMPLE_KEY, TMA_KEY, CT_KEY = "sample", "TMA_ID", "level2_refined"
PRED_COL = "pred_leiden04"        # overwritten per-run in main() based on --label
POD_COLOR = "#a6bddb"    # podocyte / glomerular tuft (anatomical anchor)
OTHER_COLOR = "#e8e8e8"  # generic tissue context
NUC_COLOR = "#9e9e9e"

# Cluster palette. First 13 are the hand-picked well-separated colours (stable across
# resolutions so leiden_0.1/0.4 keep their existing colours); the rest are drawn from
# tab20b/tab20c/tab20 with grey-ish and very light shades filtered out (they would clash
# with the grey tissue context and the pale-blue POD tint). Needed because leiden_0.8 has
# 19 clusters and leiden_1.5 has 28.
def _srgb_to_lab(rgb):
    rgb = np.asarray(rgb, float)
    lin = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    M = np.array([[0.4124, 0.3576, 0.1805], [0.2126, 0.7152, 0.0722], [0.0193, 0.1192, 0.9505]])
    xyz = lin @ M.T / np.array([0.95047, 1.0, 1.08883])
    d = 6 / 29
    f = np.where(xyz > d ** 3, np.cbrt(xyz), xyz / (3 * d * d) + 4 / 29)
    return np.stack([116 * f[..., 1] - 16, 500 * (f[..., 0] - f[..., 1]),
                     200 * (f[..., 1] - f[..., 2])], -1)


def _make_palette(n_needed=40):
    """Perceptually maximally-distinct categorical colours (greedy max-min CIELab).

    Built this way because the user reads cluster identity off the boundary plots by eye:
    a hand-picked tab-style list put two similar greens at indices 2 and 11, which are both
    strong PEC candidates at res 0.8. Colours near the grey tissue context / pale-blue POD
    tint / white are excluded so cluster fills never blend into the background.
    """
    import colorsys
    from matplotlib.colors import to_hex
    cand = np.array([colorsys.hsv_to_rgb(h, s, v)
                     for h in np.linspace(0, 1, 49)[:-1]
                     for s in (0.55, 0.80, 1.0)
                     for v in (0.45, 0.65, 0.85)])
    lab = _srgb_to_lab(cand)
    ctx = _srgb_to_lab(np.array([[0.910, 0.910, 0.910],   # OTHER_COLOR grey
                                 [0.651, 0.741, 0.859],   # POD pale blue
                                 [1.0, 1.0, 1.0],         # white bg
                                 [0.62, 0.62, 0.62]]))    # nucleus grey
    keep = np.min(np.linalg.norm(lab[:, None, :] - ctx[None, :, :], axis=-1), axis=1) > 30
    cand, lab = cand[keep], lab[keep]
    seed = int(np.argmin(np.linalg.norm(lab - _srgb_to_lab(np.array([0.12, 0.47, 0.71])), axis=1)))
    chosen = [seed]
    dmin = np.linalg.norm(lab - lab[seed], axis=1)
    while len(chosen) < min(n_needed, len(cand)):
        nxt = int(np.argmax(dmin))
        chosen.append(nxt)
        dmin = np.minimum(dmin, np.linalg.norm(lab - lab[nxt], axis=1))
    return [to_hex(cand[i]) for i in chosen]


CLUSTER_HEX = _make_palette(40)
PAL = {str(i): c for i, c in enumerate(CLUSTER_HEX)}

DBSCAN_EPS = 150.0       # microns; glomerular tuft ~150-250 um
DBSCAN_MIN = 3
GLOM_PAD = 70.0          # microns padding around a glomerulus crop
MAX_GLOM_PANELS = 12     # cap montage panels per core (biggest first)


def load_polys(parquet_path: Path):
    """Return (cell_ids[str], polys[list of (V,2)], centroids[(N,2)])."""
    df = pd.read_parquet(parquet_path, columns=["cell_id", "vertex_x", "vertex_y"])
    cid = df["cell_id"].to_numpy().astype(str)
    xs = df["vertex_x"].to_numpy(); ys = df["vertex_y"].to_numpy()
    change = np.r_[True, cid[1:] != cid[:-1]]
    starts = np.flatnonzero(change); ends = np.r_[starts[1:], len(cid)]
    polys = [np.column_stack((xs[s:e], ys[s:e])) for s, e in zip(starts, ends)]
    cents = np.array([p.mean(axis=0) for p in polys])
    return cid[starts], polys, cents


def draw(ax, polys, faces, edge, lw, alpha, zorder):
    if polys:
        ax.add_collection(PolyCollection(polys, facecolors=faces, edgecolors=edge,
                                         linewidths=lw, alpha=alpha, zorder=zorder))


def face_for(pred_val, ct_val):
    if pd.notna(pred_val):
        return PAL.get(str(int(pred_val)), OTHER_COLOR)
    return POD_COLOR if ct_val == "POD" else OTHER_COLOR


def render_box(ax, idx, polys, cents, pred, ct, nuc_polys, nuc_cents, xlim, ylim):
    """Draw all cells whose centroid is in [xlim,ylim]; PEC on top."""
    x0, x1 = xlim; y0, y1 = min(ylim), max(ylim)
    sel = np.where((cents[:, 0] >= x0) & (cents[:, 0] <= x1) &
                   (cents[:, 1] >= y0) & (cents[:, 1] <= y1))[0]
    bg_p, bg_f, fg_p, fg_f = [], [], [], []
    for i in sel:
        if pd.notna(pred[i]):
            fg_p.append(polys[i]); fg_f.append(PAL.get(str(int(pred[i])), OTHER_COLOR))
        else:
            bg_p.append(polys[i]); bg_f.append(POD_COLOR if ct[i] == "POD" else OTHER_COLOR)
    draw(ax, bg_p, bg_f, "white", 0.05, 0.9, 1)
    if nuc_polys is not None:
        nsel = np.where((nuc_cents[:, 0] >= x0) & (nuc_cents[:, 0] <= x1) &
                        (nuc_cents[:, 1] >= y0) & (nuc_cents[:, 1] <= y1))[0]
        draw(ax, [nuc_polys[i] for i in nsel], "none", NUC_COLOR, 0.12, 1.0, 2)
    draw(ax, fg_p, fg_f, "black", 0.35, 0.97, 3)
    ax.set_xlim(x0, x1); ax.set_ylim(y1, y0)          # inverted y
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def legend_handles(clusters_present):
    h = [Patch(facecolor=PAL.get(c, OTHER_COLOR), edgecolor="black", label=f"PEC cl {c}")
         for c in sorted(clusters_present, key=int)]
    h += [Patch(facecolor=POD_COLOR, edgecolor="none", label="POD (glom tuft)"),
          Patch(facecolor=OTHER_COLOR, edgecolor="none", label="other cells")]
    return h


def plot_core(tma, reg, out_core, out_glom, display=None, stem=None, overview=True):
    """reg: dict with cell_ids, polys, cents, pred(obj array), ct(str array),
    nuc_polys, nuc_cents, and a boolean mask 'in_core' aligned to region cells."""
    # display: text shown in the titles instead of the core label; stem: output file name (no extension)
    # instead of the core label; overview=False skips the full-core figure. Defaults = original behaviour;
    # the figure entry scripts pass a public core id because core labels carry pathology numbers.
    disp = tma if display is None else display
    fname = tma.replace(" ", "_").replace("/", "_") if stem is None else stem
    m = reg["in_core"]
    idx = np.where(m)[0]
    polys, cents, pred, ct = reg["polys"], reg["cents"], reg["pred"], reg["ct"]
    core_cents = cents[idx]
    pec_mask = np.array([pd.notna(pred[i]) for i in idx])
    n_pec = int(pec_mask.sum())
    clusters_present = {str(int(pred[i])) for i in idx[pec_mask]}

    # ---------- 1. full-core overview ----------
    if overview:
        x0, y0 = core_cents.min(0); x1, y1 = core_cents.max(0)
        fig, ax = plt.subplots(figsize=(11, 11))
        render_box(ax, idx, polys, cents, pred, ct, reg["nuc_polys"], reg["nuc_cents"],
                   (x0 - 30, x1 + 30), (y0 - 30, y1 + 30))
        ax.legend(handles=legend_handles(clusters_present), bbox_to_anchor=(1.01, 0.5),
                  loc="center left", fontsize=9, frameon=False, title="mapped label")
        ax.set_title(f"{disp}\n{n_pec} mapped PEC cells (full core)", fontsize=13)
        fig.tight_layout()
        base = out_core / fname
        fig.savefig(f"{base}.png", dpi=200, bbox_inches="tight")
        fig.savefig(f"{base}.pdf", bbox_inches="tight")
        plt.close(fig)

    # ---------- 2. per-glomerulus montage ----------
    pec_idx = idx[pec_mask]
    pec_cents = cents[pec_idx]
    pec_lab = np.array([int(pred[i]) for i in pec_idx])
    if len(pec_idx) >= DBSCAN_MIN:
        gl = DBSCAN(eps=DBSCAN_EPS, min_samples=DBSCAN_MIN).fit_predict(pec_cents)
    else:
        gl = np.zeros(len(pec_idx), int)
    glom_ids = [g for g in sorted(set(gl)) if g != -1]
    n_scattered = int((gl == -1).sum())
    # sort glomeruli by PEC count desc
    glom_ids = sorted(glom_ids, key=lambda g: -(gl == g).sum())
    shown = glom_ids[:MAX_GLOM_PANELS]
    n_extra = len(glom_ids) - len(shown)

    if shown:
        ncol = min(4, len(shown)); nrow = int(np.ceil(len(shown) / ncol))
        figm, axm = plt.subplots(nrow, ncol, figsize=(4.4 * ncol, 4.4 * nrow),
                                 squeeze=False)
        for k, g in enumerate(shown):
            ax = axm[k // ncol][k % ncol]
            gc = pec_cents[gl == g]
            bx0, by0 = gc.min(0) - GLOM_PAD; bx1, by1 = gc.max(0) + GLOM_PAD
            render_box(ax, idx, polys, cents, pred, ct, reg["nuc_polys"],
                       reg["nuc_cents"], (bx0, bx1), (by0, by1))
            lab = pec_lab[gl == g]
            comp = ", ".join(f"cl{c}:{int((lab == c).sum())}"
                             for c in sorted(set(lab)))
            ax.set_title(f"glom {k+1}  (n={len(gc)})\n{comp}", fontsize=9)
        for k in range(len(shown), nrow * ncol):
            axm[k // ncol][k % ncol].axis("off")
        sub = f"{disp}   -   {len(glom_ids)} glomeruli, {n_scattered} scattered PEC cells"
        if n_extra:
            sub += f"  (showing {len(shown)} largest of {len(glom_ids)})"
        figm.suptitle(sub, fontsize=12, y=1.005)
        figm.legend(handles=legend_handles(clusters_present), loc="lower center",
                    ncol=len(clusters_present) + 2, fontsize=9, frameon=False,
                    bbox_to_anchor=(0.5, -0.02))
        figm.tight_layout()
        figm.savefig(out_glom / (fname + ".png"),
                     dpi=160, bbox_inches="tight")
        plt.close(figm)
    return n_pec, len(glom_ids), n_scattered


def main(label="leiden_0.4", min_pec=5):
    global PRED_COL
    tag = "leiden" + label.split("_", 1)[1].replace(".", "")
    PRED_COL = f"pred_{tag}"
    print(f"[cfg] label={label}  tag={tag}  pred_col={PRED_COL}", flush=True)
    out = HERE / f"boundary_plots_{tag}"
    out_core = out / "per_core"; out_glom = out / "per_glomerulus"
    out_core.mkdir(parents=True, exist_ok=True); out_glom.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"cluster": list(PAL), "hex": list(PAL.values())}).to_csv(
        out / "palette.csv", index=False)

    pred_df = pd.read_csv(HERE / f"pec_{tag}_predictions.csv")
    pred_df["cell_id"] = pred_df["cell_id"].astype(str)
    pred_df["sample"] = pred_df["sample"].astype(str)
    pred_map = pred_df.set_index(["sample", "cell_id"])[PRED_COL]
    ALL_CL = sorted(pred_df[PRED_COL].astype(int).unique())   # clusters that received cells
    print(f"[cfg] {len(pred_df)} query cells; clusters with cells: {ALL_CL}", flush=True)

    obs = sc.read_h5ad(XENIUM_LABELED, backed="r").obs[
        ["cell_id", SAMPLE_KEY, TMA_KEY, CT_KEY]].copy()
    for c in obs.columns:
        obs[c] = obs[c].astype(str)

    samples = sorted(pred_df["sample"].unique())
    print(f"[run] {len(samples)} regions, min_pec={min_pec}", flush=True)
    core_rows, tot_core = [], 0
    for s in samples:
        folder = XENIUM_BASE / f"output-{s}"
        cb, nb = folder / "cell_boundaries.parquet", folder / "nucleus_boundaries.parquet"
        if not cb.exists():
            print(f"[skip region] {s}: no cell_boundaries.parquet", flush=True); continue
        print(f"[load] {s}", flush=True)
        cell_ids, polys, cents = load_polys(cb)
        nuc_polys = nuc_cents = None
        if nb.exists():
            _, nuc_polys, nuc_cents = load_polys(nb)

        obs_s = obs[obs[SAMPLE_KEY] == s].set_index("cell_id")
        lookup = pd.DataFrame(index=cell_ids)
        lookup[TMA_KEY] = obs_s[TMA_KEY].reindex(cell_ids).to_numpy()
        lookup[CT_KEY] = obs_s[CT_KEY].reindex(cell_ids).to_numpy()
        pmap_s = pred_map.loc[s] if s in pred_map.index.get_level_values(0) else pd.Series(dtype=object)
        lookup[PRED_COL] = pmap_s.reindex(cell_ids).to_numpy()

        pred_arr = lookup[PRED_COL].to_numpy()
        ct_arr = lookup[CT_KEY].astype(str).to_numpy()
        tma_arr = lookup[TMA_KEY].astype(str).to_numpy()

        pec_here = lookup.dropna(subset=[PRED_COL])
        counts = pec_here.groupby(TMA_KEY).size().sort_values(ascending=False)
        for tma, n in counts.items():
            if tma == "nan" or n < min_pec:
                continue
            reg = dict(cell_ids=cell_ids, polys=polys, cents=cents, pred=pred_arr,
                       ct=ct_arr, nuc_polys=nuc_polys, nuc_cents=nuc_cents,
                       in_core=(tma_arr == tma))
            n_pec, n_glom, n_scat = plot_core(tma, reg, out_core, out_glom)
            vc = pec_here.loc[pec_here[TMA_KEY] == tma, PRED_COL].astype(int).value_counts()
            row = {"TMA_ID": tma, "sample": s, "n_pec": int(n_pec),
                   "n_glomeruli": n_glom, "n_scattered": n_scat}
            row.update({f"cl{c}": int(vc.get(c, 0)) for c in ALL_CL})
            core_rows.append(row); tot_core += 1
            print(f"   [core] {tma}: {n_pec} PEC, {n_glom} glomeruli, {n_scat} scattered",
                  flush=True)

    pd.DataFrame(core_rows).to_csv(out / "per_core_cluster_counts.csv", index=False)
    clusters_present = sorted(pred_df[PRED_COL].astype(int).unique())
    (out / "README.txt").write_text(
        f"Xenium cell-boundary plots of mapped PEC {label} subclusters\n"
        "==============================================================\n\n"
        f"Cores rendered (>= {min_pec} mapped PEC cells): {tot_core}\n\n"
        "FOLDERS\n"
        "  per_core/<core>.{png,pdf}   full-core overview: where PEC cells sit in the\n"
        "                              whole core. PDF is vector - zoom in freely.\n"
        "  per_glomerulus/<core>.png   montage: each glomerulus zoomed so cluster colours\n"
        "                              are obvious (biggest glomeruli first, up to 12).\n\n"
        "COLOUR KEY\n"
        f"  filled colours transferred PEC {label} subcluster (palette.csv)\n"
        "  pale blue      podocytes (glomerular tuft) - anatomical anchor\n"
        "  light grey     all other cells (tissue context)\n"
        "  grey outlines  cell nuclei\n\n"
        "HOW LABELS WERE PRODUCED\n"
        f"  snRNA PEC subset (7,004 cells) {label} transferred onto the {len(pred_df)}\n"
        "  Xenium PEC cells (level2_refined==PEC) via scArches surgery\n"
        "  (SCVI->SCANVI->load_query_data), 468 shared panel genes, seed 42.\n"
        f"  Clusters that received Xenium cells: {clusters_present}. Per-cell confidence is\n"
        f"  in ../pec_{tag}_predictions.csv (low-confidence cells are still shown).\n\n"
        "WHAT TO LOOK FOR\n"
        "  Genuine PEC subclusters should ring Bowman's capsule around the pale-blue\n"
        "  podocyte tuft. Clusters sitting in tubular interstitium away from glomeruli are\n"
        "  likely tubular-contamination subclusters. per_core_cluster_counts.csv gives the\n"
        "  exact per-core composition (incl. n_glomeruli and n_scattered PEC cells).\n\n"
        "IMPORTANT - reading the montage panels\n"
        "  'glom N' is just a spatial DBSCAN cluster of PEC-LABELLED cells, not a verified\n"
        "  glomerulus. A TRUE glomerulus panel shows PEC cells arcing around a pale-blue\n"
        "  podocyte tuft. A panel with cells scattered across tissue and NO pale-blue tuft\n"
        "  is TUBULAR CONTAMINATION (the thin-limb / ATL-DTL subclusters are ~0% in_glom;\n"
        "  see pred_cluster_in_glom_summary.png). Use the pale-blue tuft as the arbiter.\n\n"
        "SUMMARY FIGURE\n"
        "  pred_cluster_in_glom_summary.{pdf,png} - per cluster, the % of its Xenium cells\n"
        "  inside a glomerulus. High % = glomerular (real PEC); ~0% = tubular contamination.\n")
    print(f"\n[done] {tot_core} cores -> {out_core} + {out_glom}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", default="leiden_0.4")
    ap.add_argument("--min-pec", type=int, default=5)
    a = ap.parse_args()
    main(label=a.label, min_pec=a.min_pec)

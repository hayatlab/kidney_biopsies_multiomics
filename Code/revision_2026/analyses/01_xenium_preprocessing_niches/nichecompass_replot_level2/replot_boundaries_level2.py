"""Cell-boundary two-panel per-core plots (LEFT: NicheCompass niches, RIGHT:
level2_refined cell type) for any of the NicheCompass versions.

Generalised from review_analyses/4_plot_niche_cell_boundaries_5regions.py:
 - cell-type panel uses the refined `level2_refined` annotation (fixed palette),
 - niche panel keeps the object's existing niche colours (uns) so it matches the
   predic version's left panel,
 - works for the overall object and the glom subsets (obs subset -> only those
   cells survive the left-join, so glom cores plot just the glom cells).

Runs under `spatial_env` (geopandas / shapely). One PDF+PNG per TMA core into
<dest>/cell_boundary_plots_per_core/ .
"""
from __future__ import annotations
from kbpaths import KB_ANALYSIS, KB_RESULTS, KB_XENIUM_RAW  # repository path config (config/kbpaths.py)
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.colors import ListedColormap
from shapely.geometry import Polygon

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from level2_common import merge_level2, L2_COLORS, ORDER   # noqa: E402

NEW = Path(f"{KB_RESULTS}/xenium/leon/new_data")
# the bulky predic overall source object was relocated to archived/
ARCHIVED = Path(f"{KB_RESULTS}/xenium/leon/archived")
XENIUM_BASE = Path(f"{KB_XENIUM_RAW}")
FULL_OBS_PATH = Path(f"{KB_ANALYSIS}/celltype_concordance/"
                     "outputs/xenium_labeled.h5ad")  # all QC-passed cells + TMA_ID (grey bg)
CT_KEY = "level2_refined"
SAMPLE_KEY = "sample"
TMA_KEY = "TMA_ID"

VERSIONS = {
    "overall": dict(
        adata=ARCHIVED / "overall_with_glom_04/xenium_with_niches_overall_with_glom.h5ad",
        niche_key="latent_leiden_0.4",
        dest=NEW / "overall_with_glom_04_level2/cell_boundary_plots_per_core"),
    "manual_glom": dict(
        adata=ARCHIVED / "manual_glom_niche_based_on_04/xenium_with_niches_manual_glom_niche.h5ad",
        niche_key="latent_leiden_0.6",
        dest=NEW / "manual_glom_niche_based_on_04_level2/cell_boundary_plots_per_core"),
    "in_glom": dict(
        adata=NEW / "in_glom_niche_based_on_04/xenium_with_niches_in_glom_niche.h5ad",
        niche_key="latent_leiden_0.6",
        dest=NEW / "in_glom_niche_based_on_04_level2/cell_boundary_plots_per_core"),
    "manual_glom_10": dict(
        adata=NEW / "manual_glom_niche_based_on_04_res10_level2/xenium_with_niches.h5ad",
        niche_key="latent_leiden_1.0",
        dest=NEW / "manual_glom_niche_based_on_04_res10_level2/cell_boundary_plots_per_core"),
}


def load_boundary_polygons(parquet_path: Path) -> gpd.GeoDataFrame:
    df = pd.read_parquet(parquet_path, columns=["cell_id", "vertex_x", "vertex_y"])
    cell_ids = df["cell_id"].to_numpy()
    xs = df["vertex_x"].to_numpy(); ys = df["vertex_y"].to_numpy()
    change = np.r_[True, cell_ids[1:] != cell_ids[:-1]]
    starts = np.flatnonzero(change)
    ends = np.r_[starts[1:], len(cell_ids)]
    polys = [Polygon(zip(xs[s:e], ys[s:e])) for s, e in zip(starts, ends)]
    return gpd.GeoDataFrame({"cell_id": cell_ids[starts]}, geometry=polys)


def niche_colors_from(adata, niche_key: str) -> dict:
    cats = sorted(adata.obs[niche_key].astype(str).unique(),
                  key=lambda v: int(v) if v.isdigit() else v)
    ck = f"{niche_key}_colors"
    if ck in adata.uns and len(adata.uns[ck]) == len(cats):
        return dict(zip(cats, list(adata.uns[ck])))
    # fall back to a generated palette
    sc.pl.umap(adata, color=[niche_key], show=False); plt.close("all")
    return dict(zip(cats, list(adata.uns.get(ck, plt.cm.tab20.colors[:len(cats)]))))


def safe_filename(s: str) -> str:
    return s.replace(" ", "_").replace("/", "_")


def plot_core(tma, gdf_cells, gdf_nuc, niche_key, niche_cmap_d, ct_cmap_d, out_dir,
              grey_bg: bool = False, label: str | None = None, out_base=None):
    # label / out_base (figure entry scripts): public core name for the title and the output file (without
    # extension); default None = the internal core label (contains a pathology number) / <out_dir>/<label>
    # niche-labelled (glom) cells; the rest of the core (grey_bg) has NaN niche
    gdf_niche = gdf_cells.dropna(subset=[niche_key]).copy()
    if gdf_niche.empty:
        print(f"   [skip] {tma if label is None else label}: no cells with niche label", flush=True)
        return
    gdf_grey = gdf_cells[gdf_cells[niche_key].isna()] if grey_bg else None
    present_n = set(gdf_niche[niche_key].astype(str).unique())
    n_here = [c for c in niche_cmap_d if c in present_n]
    gdf_niche[niche_key] = pd.Categorical(gdf_niche[niche_key].astype(str),
                                          categories=n_here, ordered=False)
    niche_cmap = ListedColormap([niche_cmap_d[c] for c in n_here])

    gdf_ct = gdf_niche.dropna(subset=[CT_KEY]).copy()
    ct_present = set(gdf_ct[CT_KEY].astype(str).unique())
    ct_here = [c for c in ct_cmap_d if c in ct_present]
    gdf_ct[CT_KEY] = pd.Categorical(gdf_ct[CT_KEY].astype(str),
                                    categories=ct_here, ordered=False)
    ct_cmap = ListedColormap([ct_cmap_d[c] for c in ct_here])

    fig, axes = plt.subplots(1, 2, figsize=(28, 12))
    # LEFT: grey background (non-niche cells) then niche-coloured cells
    if grey_bg and gdf_grey is not None and not gdf_grey.empty:
        gdf_grey.plot(color="0.86", edgecolor="0.72", linewidth=0.1, ax=axes[0], zorder=1)
    gdf_niche.plot(column=niche_key, categorical=True, legend=True, cmap=niche_cmap,
                   linewidth=0.3, alpha=0.85, edgecolor="black", ax=axes[0], zorder=2,
                   legend_kwds={"bbox_to_anchor": (1, 0.5), "loc": "center left",
                                "fontsize": 10, "title": "Niche"})
    if not gdf_nuc.empty:
        gdf_nuc.boundary.plot(color="grey", linewidth=0.15, ax=axes[0], zorder=10)
    axes[0].set_ylim(axes[0].get_ylim()[::-1])
    axes[0].set_title(f"NicheCompass niches ({niche_key})", fontsize=16)

    # RIGHT: grey background then cell-type-coloured cells
    if grey_bg and gdf_grey is not None and not gdf_grey.empty:
        gdf_grey.plot(color="0.86", edgecolor="0.72", linewidth=0.1, ax=axes[1], zorder=1)
    if not gdf_ct.empty:
        gdf_ct.plot(column=CT_KEY, categorical=True, legend=True, cmap=ct_cmap,
                    linewidth=0.3, alpha=0.85, edgecolor="black", ax=axes[1], zorder=2,
                    legend_kwds={"bbox_to_anchor": (1, 0.5), "loc": "center left",
                                 "fontsize": 10, "title": "Cell type"})
    if not gdf_nuc.empty:
        gdf_nuc.boundary.plot(color="grey", linewidth=0.15, ax=axes[1], zorder=10)
    axes[1].set_ylim(axes[1].get_ylim()[::-1])
    axes[1].set_title(f"Cell type ({CT_KEY})", fontsize=16)

    for ax in axes:
        ax.tick_params(axis="both", which="both", bottom=False, top=False,
                       left=False, right=False, labelbottom=False, labeltop=False,
                       labelleft=False, labelright=False)
        for spine in ax.spines.values():
            spine.set_visible(False)
    plt.subplots_adjust(wspace=0.001)
    ngrey = 0 if gdf_grey is None else len(gdf_grey)
    plt.suptitle(f"{tma if label is None else label}   ({len(gdf_niche)} niche cells"
                 + (f" + {ngrey} other cells (grey)" if grey_bg else "") + ")", fontsize=15)
    plt.tight_layout()
    base = out_dir / safe_filename(str(tma)) if out_base is None else Path(out_base)
    fig.savefig(f"{base}.pdf", bbox_inches="tight")
    fig.savefig(f"{base}.png", bbox_inches="tight", dpi=140)
    plt.close(fig)
    print(f"   saved {base}.pdf  ({len(gdf_niche)} niche / {ngrey} grey)", flush=True)


def process_region(sample, adata, niche_key, niche_cmap_d, ct_cmap_d, out_dir,
                   max_cores: int = 0, full_obs=None, cores: dict | None = None):
    folder = XENIUM_BASE / f"output-{sample}"
    cb, nb = folder / "cell_boundaries.parquet", folder / "nucleus_boundaries.parquet"
    if not cb.exists() or not nb.exists():
        print(f"[skip] missing boundaries for {sample} under {folder}", flush=True)
        return
    print(f"[load] {sample}: boundary polygons", flush=True)
    gdf_cells = load_boundary_polygons(cb)
    gdf_nuc = load_boundary_polygons(nb)
    gdf_cells["cell_id"] = gdf_cells["cell_id"].astype(str)
    gdf_nuc["cell_id"] = gdf_nuc["cell_id"].astype(str)
    obj_sub = adata.obs.loc[adata.obs[SAMPLE_KEY].astype(str) == sample,
                            ["cell_id", TMA_KEY, niche_key, CT_KEY]].copy()
    obj_sub["cell_id"] = obj_sub["cell_id"].astype(str)
    grey_bg = full_obs is not None
    if grey_bg:
        # TMA_ID for ALL cells of the region (so non-glom cells are placed in a core),
        # niche + cell type only for the object's (glom) cells -> NaN elsewhere = grey
        full_sub = full_obs.loc[full_obs[SAMPLE_KEY].astype(str) == sample,
                                ["cell_id", TMA_KEY]].copy()
        full_sub["cell_id"] = full_sub["cell_id"].astype(str)
        gdf_cells = gdf_cells.merge(full_sub, on="cell_id", how="left")
        gdf_cells = gdf_cells.merge(obj_sub[["cell_id", niche_key, CT_KEY]],
                                    on="cell_id", how="left")
        gdf_nuc = gdf_nuc.merge(full_sub, on="cell_id", how="left")
    else:
        gdf_cells = gdf_cells.merge(obj_sub, on="cell_id", how="left")
        gdf_nuc = gdf_nuc.merge(obj_sub[["cell_id", TMA_KEY]], on="cell_id", how="left")
    # only plot cores that actually contain glom (object) cells
    tmas = sorted(obj_sub[TMA_KEY].dropna().astype(str).unique())
    if max_cores:
        tmas = tmas[:max_cores]
    print(f"   {len(tmas)} TMA cores in this region"
          + ("  [full-context grey background]" if grey_bg else ""), flush=True)
    if cores is not None:   # figure entry scripts: only the requested cores, named by their public id
        tmas = [t for t in tmas if t in cores]
    for tma in tmas:
        try:
            plot_core(tma, gdf_cells[gdf_cells[TMA_KEY].astype(str) == tma],
                      gdf_nuc[gdf_nuc[TMA_KEY].astype(str) == tma],
                      niche_key, niche_cmap_d, ct_cmap_d, out_dir, grey_bg=grey_bg,
                      **(cores[tma] if cores is not None else {}))
        except Exception as exc:
            print(f"   [error] {tma if cores is None else cores[tma]['label']}: {exc!r}", flush=True)
            import traceback; traceback.print_exc()


def run(cfg: dict, max_cores: int = 0, full_context: bool = False, cores: dict | None = None) -> None:
    """cores (figure entry scripts): {internal core label: {"label": public name, "out_base": output path
    without extension}} -> plot only these cores and load only their regions; None = every core."""
    adata_path, niche_key, out_dir = cfg["adata"], cfg["niche_key"], cfg["dest"]
    out_dir.mkdir(parents=True, exist_ok=True)
    full_obs = None
    if full_context:
        import anndata as ad
        print(f"[full-context] loading all-cell obs (TMA_ID) from {FULL_OBS_PATH}", flush=True)
        fo = ad.read_h5ad(str(FULL_OBS_PATH), backed="r").obs
        full_obs = fo[[c for c in ["cell_id", SAMPLE_KEY, TMA_KEY] if c in fo.columns]].copy()
        print(f"[full-context] {len(full_obs)} reference cells", flush=True)
    print(f"[load] adata {adata_path}", flush=True)
    adata = sc.read_h5ad(adata_path)
    if CT_KEY not in adata.obs.columns:
        adata, _ = merge_level2(adata, key=CT_KEY)
    niche_vals = adata.obs[niche_key].astype(str)
    niche_cats = sorted(niche_vals.unique(), key=lambda v: int(v) if v.isdigit() else v)
    adata.obs[niche_key] = pd.Categorical(niche_vals, categories=niche_cats, ordered=False)
    niche_cmap_d = niche_colors_from(adata, niche_key)
    ct_cmap_d = {c: L2_COLORS.get(c, "#999999")
                 for c in ORDER if c in set(adata.obs[CT_KEY].astype(str))}
    print(f"[palette] niches={len(niche_cmap_d)}, cell types={len(ct_cmap_d)}", flush=True)
    for sample in sorted(adata.obs[SAMPLE_KEY].astype(str).unique()):
        if cores is not None and not (set(adata.obs.loc[adata.obs[SAMPLE_KEY].astype(str) == sample,
                                                        TMA_KEY].astype(str)) & set(cores)):
            continue
        try:
            process_region(sample, adata, niche_key, niche_cmap_d, ct_cmap_d, out_dir,
                           max_cores=max_cores, full_obs=full_obs, cores=cores)
        except Exception as exc:
            print(f"[error] region {sample}: {exc!r}", flush=True)
            import traceback; traceback.print_exc()
    n = len(list(out_dir.glob("*.pdf")))
    print(f"[done] wrote {n} per-core PDFs to {out_dir}", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--version", choices=list(VERSIONS) + ["all"], default="overall")
    ap.add_argument("--adata", type=str, default="")
    ap.add_argument("--niche-key", type=str, default="")
    ap.add_argument("--dest", type=str, default="")
    ap.add_argument("--max-cores", type=int, default=0, help="0=all; N=first N cores/region (smoke)")
    ap.add_argument("--dest-root", type=str, default="",
                    help="if set, remap dest under this base (scratch for smoke)")
    ap.add_argument("--full-context", action="store_true",
                    help="draw every cell in the core; cells without a niche in this run are grey "
                         "(useful for the glom subsets)")
    args = ap.parse_args()
    if args.adata:                       # ad-hoc override (e.g. Part B 0.5/0.6 runs)
        run(dict(adata=Path(args.adata), niche_key=args.niche_key, dest=Path(args.dest)),
            max_cores=args.max_cores, full_context=args.full_context)
    else:
        for v in (list(VERSIONS) if args.version == "all" else [args.version]):
            cfg = dict(VERSIONS[v])
            if args.dest_root:
                cfg["dest"] = Path(args.dest_root) / cfg["dest"].parent.name / cfg["dest"].name
            run(cfg, max_cores=args.max_cores, full_context=args.full_context)
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python
"""Fig. 6g (right) — snRNA-seq immune subclusters (Leiden 0.6) transferred onto Xenium immune cells, full-core maps
of the niche-20 / crescent core K50 (IgAN) and of the TIN core K18.
Panel        : Main Figures.pdf, Figure 6, panel g (right half, "predicted from snRNA-seq using immune subclusters")
Reads        : $KB_ANALYSIS/immune_mapping_leiden06/prediction/xenium_immune_with_prediction.h5ad
                 (analyses/09_snrna_pec_immune_subclusters/immune_mapping/map_immune_leiden06_xenium.py, wrapper
                  run_map_immune06.sh; query = Xenium Immune/PL/B cells of the 2025 three-region NicheCompass object,
                  immune_mapping/nichecompass_filtered_tma.ipynb)
               $KB_DATA/processed/xenium/geodataframes/{cell,nucleus}_boundaries.parquet (external: cell polygons of
                 the 2025 three-region object, with TMA_ID and cell_type_manual)
Writes       : $KB_FIG_OUT/Fig6/Fig6g_immune_subclusters_K50_XEN093.png, $KB_FIG_OUT/Fig6/Fig6g_immune_subclusters_K18_XEN116.png
Plotting code: analyses/09_snrna_pec_immune_subclusters/immune_mapping/plot_immune06_cell_boundaries.py::plot_core_categorical
               (imported unchanged; title and output path are its parameters); boundary/prediction join and palette
               extracted verbatim from ::main (original: predictions/predictions_<core label>.png)
Illustrator  : glomerular region cut out with a circular mask; legend codes renamed to subcluster names
               (0 T, 1 TOX+ T, 3 B, ... 15 DC); title removed
Environment  : conda env nichecompass
Verified     : 2026-11-12, DIFFERENT vs $KB_ANALYSIS/immune_mapping_leiden06/cell_boundary_plots_immune_0_60/predictions/
               predictions_<core label>.png (K50 and K18) — only the axes title differs (it named the core by its
               internal label; here by its public id); with the title band masked the images are pixel-identical
"""
import sys
from pathlib import Path

import anndata as ad
import geopandas as gpd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from niches_snrna import ANALYSES_09, CORES, core_display, load_script  # noqa: E402
from kbpaths import fig_out  # noqa: E402
from cores import core_label  # noqa: E402

im = load_script(ANALYSES_09 / "immune_mapping" / "plot_immune06_cell_boundaries.py")


def main():
    # ---- extracted from plot_immune06_cell_boundaries.py::main ----
    gdf_cells = gpd.read_parquet(im.CELL_BOUNDARIES_PATH)
    gdf_nuc = gpd.read_parquet(im.NUCLEUS_BOUNDARIES_PATH)
    adata = ad.read_h5ad(im.ANNDATA_PATH)
    if im.COL_PRED in gdf_cells.columns:
        gdf_cells = gdf_cells.drop(columns=[im.COL_PRED])
    obs_df = adata.obs[[im.COL_PRED]].copy()
    obs_df["join_id"] = im.norm_id(adata.obs[im.FORCE_OBS_CELL_ID_COL])
    obs_df = obs_df.dropna(subset=["join_id"])
    obs_df = obs_df[~obs_df["join_id"].duplicated(keep="first")]
    gdf_cells["join_id"] = im.norm_id(gdf_cells[im.GDF_JOIN_COL])
    merged = gdf_cells.merge(obs_df[["join_id", im.COL_PRED]], on="join_id", how="left")
    merged = gpd.GeoDataFrame(merged, geometry="geometry", crs=gdf_cells.crs).drop(columns=["join_id"])
    n_pred = int(merged[im.COL_PRED].notna().sum())
    print(f"[Fig6g immune] matched {n_pred}/{len(merged)} boundary cells to a prediction ({adata.n_obs} immune cells)",
          flush=True)
    # palette over all predicted categories of the whole object (as in the original run)
    pred_cats = im._catsort(im.unify_to_str(merged[im.COL_PRED]).dropna().unique().tolist())
    pal_pred = im.make_palette(pred_cats)

    for xen in CORES:
        core = core_label(xen)
        gdf_core = merged.loc[merged["TMA_ID"] == core].copy()
        if len(gdf_core) == 0:
            raise SystemExit(f"{xen}: no cells in the boundary table")
        gdf_core_nuc = gdf_nuc.loc[gdf_nuc["TMA_ID"] == core].copy() if "TMA_ID" in gdf_nuc.columns else None
        display, stem = core_display(xen)
        out = fig_out("Fig6", f"Fig6g_immune_subclusters_{stem}.png")
        im.plot_core_categorical(gdf_core, gdf_core_nuc, im.COL_PRED, pal_pred, out,
                                 title=f"{display} — immune leiden_0.6 predictions", legend_title="leiden_0.6")
        print(f"[Fig6g immune] {display}: {len(gdf_core)} cells, "
              f"{int(gdf_core[im.COL_PRED].notna().sum())} with a prediction -> {out}", flush=True)


if __name__ == "__main__":
    main()

#
# Created on Fri Aug 22 2025
#
# Copyright (c) 2025 Vivien Goepp
#

# --------------------------------------------------------------------------------------------
# This script preprocesses and merges spatial boundary data for Xenium kidney biopsy samples.
#
# Datasets used/loaded:
# - Raw Xenium spatial data (cell/nucleus boundaries, circles, morphology images) from multiple sample folders.
# - Preprocessed annotated cell metadata: xenium_processed.h5ad
#
# Datasets generated/saved:
# - GeoDataFrames of cell boundaries, cell circles, and nucleus boundaries with merged cell metadata:
#     - Saved as Parquet files in: /data/Kramann_lab/kidney_biopsy/data/processed/xenium/geodataframes/
#         - cell_boundaries.parquet
#         - cell_circles.parquet
#         - nucleus_boundaries.parquet
#
# Main processing steps:
# - Loads Xenium spatial data for each sample and merges with annotated cell metadata.
# - Extracts and merges cell boundaries, cell circles, and nucleus boundaries as GeoDataFrames.
# - Filters out cells with missing TMA_ID or condition.
# - Adds centroid coordinates to boundaries.
# - Reorders cell type categories for consistency.
# - Defines unique sample_cell_id and sets as index.
# - Assigns a planar CRS to all GeoDataFrames.
# - Saves processed GeoDataFrames for downstream spatial analysis and visualization.
#
# Plots:
# - (Commented out) Optionally renders and saves morphology images for each sample:
#     - Plots morphology_focus image using sdata.pl.render_images, saves as PNG in:
#         /data/Kramann_lab/From_Vivien/kidney_biopsy/cell_boundaries/morphology_image.png
#
# Note:
# - No other plots are generated or saved by this script.
# - The script prepares spatial boundary data for further analysis and visualization.
# --------------------------------------------------------------------------------------------

from spatialdata_io import xenium
import scanpy as sc
import geopandas as gpd
import pandas as pd
import pyarrow.parquet as pq

# Set directories ---------------------------

data_dir = "/data/zd775156/kidney_biopsy/data/"
xenium_path = "/data/Kramann_lab/From_Sidrah/new_data/output-XETG00229__0042348__Region_3__20250416__104148/"
xenium_path_dict = {
    "XETG00229__0042348__Region_2__20250416__104148":
        "/data/Kramann_lab/From_Sidrah/new_data/output-XETG00229__0042348__Region_2__20250416__104148/",
    "XETG00229__0042348__Region_3__20250416__104148":
        "/data/Kramann_lab/From_Sidrah/new_data/output-XETG00229__0042348__Region_3__20250416__104148/",
    "XETG00229__0041734__Region_3__20250416__104148":
        "/data/Kramann_lab/From_Sidrah/new_data/output-XETG00229__0041734__Region_3__20250416__104148/",
}
res_dir = "/data/Kramann_lab/kidney_biopsy/data/processed/xenium/"
fig_dir = "/data/Kramann_lab/From_Vivien/kidney_biopsy/cell_boundaries/"

# Define samples --------------------------------------------------------------------------------
samples = list(xenium_path_dict.keys())
sample = samples[0]

# Populate empty gdf objects ------------------------------------------------------------------
gdf_cell_boundaries = {sample: None for sample in samples}
gdf_cell_circles = {sample: None for sample in samples}
gdf_nucleus_boundaries = {sample: None for sample in samples}

# Monkey-patch the pandas function to fix importing issue -------------------------


def safe_read_parquet(*args, **kwargs):
    return pq.read_table(*args, **kwargs).to_pandas(types_mapper=None)


pd.read_parquet = safe_read_parquet

for sample in samples:

    # Load data ---------------------------
    sdata = xenium(xenium_path_dict[sample], cells_as_circles=True)

    del sdata["morphology_focus"]

    # SpatialData object
    # ├── Images
    # │     └── 'morphology_focus': DataTree[cyx](5, 34112, 31315), (5, 17056, 15657), (5, 8528, 7828),
    # (5, 4264, 3914), (5, 2132, 1957)
    # ├── Labels
    # │     ├── 'cell_labels': DataTree[yx](34112, 31315), (17056, 15657), (8528, 7828), (4264, 3914), (2132, 1957)
    # │     └── 'nucleus_labels': DataTree[yx](34112, 31315), (17056, 15657), (8528, 7828), (4264, 3914), (2132, 1957)
    # ├── Points
    # │     └── 'transcripts': DataFrame with shape: ( < Delayed > , 13) (3D points)
    # ├── Shapes
    # │     ├── 'cell_boundaries': GeoDataFrame shape: (63860, 1)(2D shapes)
    # │     ├── 'cell_circles': GeoDataFrame shape: (63860, 2)(2D shapes)
    # │     └── 'nucleus_boundaries': GeoDataFrame shape: (63549, 1)(2D shapes)
    # └── Tables
    # └── 'table': AnnData(63860, 468)
    # with coordinate systems:
    #     ▸ 'global', with elements:
    #         morphology_focus(Images), cell_labels(Labels), nucleus_labels(Labels), transcripts(
    #             Points), cell_boundaries(Shapes), cell_circles(Shapes), nucleus_boundaries(Shapes)

    sdata['table'].obs["cell_id"].value_counts().sort_values(ascending=False)

    # Load adata ------------------------------------------------

    adata = sc.read_h5ad("/data/Kramann_lab/kidney_biopsy/data/processed/xenium/xenium_processed.h5ad")

    adata.obs["cell_id"].value_counts().sort_values(ascending=False)
    adata.obs["sample_cell_id"] = adata.obs["TMA_ID"].astype(str) + "_" + adata.obs["cell_id"].astype(str)
    # Replace all " " with "_" in values of adata.obs["sample_cell_id"]
    adata.obs["sample_cell_id"] = adata.obs["sample_cell_id"].str.replace(" ", "_")
    adata.obs["sample_cell_id"].value_counts().sort_values(ascending=False)

    adata.obs['cell_type_manual'].replace('PT/aPT', 'PT1', inplace=True)
    adata.obs['cell_type_manual'].replace('PT', 'PT2', inplace=True)

    # Assert that the TMA_cell_id are unique
    assert len(adata.obs["sample_cell_id"].unique()) == adata.n_obs

    # Add obs from adata to sdata --------------------------------
    sdata["table"].obs["sample"] = sample

    cols_to_add = [c for c in adata.obs.columns
                   if c not in sdata["table"].obs.columns or c in ["cell_id", "sample"]]

    sdata["table"].obs = sdata["table"].obs.merge(
        adata.obs[cols_to_add],
        left_on=["cell_id", "sample"],
        right_on=["cell_id", "sample"],
        how="left"
    )

    # del adata

    assert len(sdata["table"].obs["cell_id"].unique()) == sdata["table"].n_obs

    # Access cell_boundaries in Shapes ------------------------
    sdata["table"].obs.columns
    sdata["nucleus_boundaries"]

    # Get geopandas info --------------------------------------------------
    # Merge sdata["table"].obs to each of the spatial objects and save as separate GeoDataFrames

    # 1. Cell boundaries
    gdf_cell_boundaries[sample] = sdata["cell_boundaries"].merge(
        sdata["table"].obs,
        left_index=True,
        right_on="cell_id",
        how="left"
    )
    gdf_cell_boundaries[sample] = gpd.GeoDataFrame(gdf_cell_boundaries[sample])
    # gdf_cell_boundaries[sample].to_parquet(res_dir + f"{sample}_cell_boundaries.parquet")

    # 2. Cell circles
    gdf_cell_circles[sample] = sdata["cell_circles"].merge(
        sdata["table"].obs,
        left_index=True,
        right_on="cell_id",
        how="left"
    )
    gdf_cell_circles[sample] = gpd.GeoDataFrame(gdf_cell_circles[sample])
    # gdf_cell_circles[sample].to_parquet(res_dir + f"{sample}_cell_circles.parquet")

    # 3. Nucleus boundaries
    gdf_nucleus_boundaries[sample] = sdata["nucleus_boundaries"].merge(
        sdata["table"].obs,
        left_index=True,
        right_on="cell_id",
        how="left"
    )
    gdf_nucleus_boundaries[sample] = gpd.GeoDataFrame(gdf_nucleus_boundaries[sample])

del sdata

# Merge gdf objects from the 3 TMAs --------------------------------------
gdf_nucleus_boundaries = pd.concat(gdf_nucleus_boundaries.values(), ignore_index=True)
gdf_cell_boundaries = pd.concat(gdf_cell_boundaries.values(), ignore_index=True)
gdf_cell_circles = pd.concat(gdf_cell_circles.values(), ignore_index=True)

# Remove cores with almost no cells ----------------------------------------
cores_to_rm = ['Utrecht 2 A5 T03-04110 K29 MN', 'Utrecht1 B5 T18-01374 K25 TIN', 'Utrecht2 C3 T17-11176 K24 TIN']
gdf_cell_boundaries = gdf_cell_boundaries[~gdf_cell_boundaries["TMA_ID"].isin(cores_to_rm)]
gdf_nucleus_boundaries = gdf_nucleus_boundaries[~gdf_nucleus_boundaries["TMA_ID"].isin(cores_to_rm)]
gdf_cell_circles = gdf_cell_circles[~gdf_cell_circles["TMA_ID"].isin(cores_to_rm)]

# Remove values with NaN TMA_ID ----------------------------------------
gdf_cell_boundaries["TMA_ID"].isna().value_counts(dropna=False)
gdf_cell_boundaries["condition"].isna().value_counts(dropna=False)
gdf_cell_circles["TMA_ID"].isna().value_counts(dropna=False)
gdf_cell_circles["condition"].isna().value_counts(dropna=False)
gdf_nucleus_boundaries["TMA_ID"].isna().value_counts(dropna=False)
gdf_nucleus_boundaries["condition"].isna().value_counts(dropna=False)

# Check that the indices where TMA_ID are NaN are the same then where condition is NaN
assert gdf_cell_boundaries[gdf_cell_boundaries["TMA_ID"].isna()].index.equals(
    gdf_cell_boundaries[gdf_cell_boundaries["condition"].isna()].index
)
assert gdf_cell_circles[gdf_cell_circles["TMA_ID"].isna()].index.equals(
    gdf_cell_circles[gdf_cell_circles["condition"].isna()].index
)
assert gdf_nucleus_boundaries[gdf_nucleus_boundaries["TMA_ID"].isna()].index.equals(
    gdf_nucleus_boundaries[gdf_nucleus_boundaries["condition"].isna()].index
)

# Filter out NaN TMA_IDs
gdf_cell_boundaries = gdf_cell_boundaries[gdf_cell_boundaries["TMA_ID"].notna()]
gdf_nucleus_boundaries = gdf_nucleus_boundaries[gdf_nucleus_boundaries["TMA_ID"].notna()]
gdf_cell_circles = gdf_cell_circles[gdf_cell_circles["TMA_ID"].notna()]

gdf_cell_boundaries["TMA_ID"].isna().value_counts()
gdf_nucleus_boundaries["TMA_ID"].isna().value_counts()
gdf_cell_circles["TMA_ID"].isna().value_counts()

# Add centroids  --------------------
gdf_cell_boundaries['centroid'] = gdf_cell_boundaries.geometry.centroid
gdf_nucleus_boundaries['centroid'] = gdf_nucleus_boundaries.geometry.centroid


# Reorder cell_type_manual categories ----------------------------------
cell_types = sorted(
    gdf_cell_boundaries["cell_type_manual"].unique(), key=str.lower)

gdf_cell_boundaries["cell_type_manual"] = gdf_cell_boundaries["cell_type_manual"].cat.reorder_categories(
    cell_types, ordered=True)
gdf_nucleus_boundaries["cell_type_manual"] = gdf_nucleus_boundaries["cell_type_manual"].cat.reorder_categories(
    cell_types, ordered=True)
gdf_cell_circles["cell_type_manual"] = gdf_cell_circles["cell_type_manual"].cat.reorder_categories(
    cell_types, ordered=True)

# Define sample_cell_id ----------------------
gdf_cell_boundaries["sample_cell_id"] = (gdf_cell_boundaries["TMA_ID"].astype(str) +
                                         '_' + gdf_cell_boundaries["cell_id"].astype(str))
gdf_nucleus_boundaries["sample_cell_id"] = (gdf_nucleus_boundaries["TMA_ID"].astype(str) +
                                            '_' + gdf_nucleus_boundaries["cell_id"].astype(str))
gdf_cell_circles["sample_cell_id"] = (gdf_cell_circles["TMA_ID"].astype(str) +
                                      '_' + gdf_cell_circles["cell_id"].astype(str))
# Set as index
gdf_cell_boundaries.index = gdf_cell_boundaries["sample_cell_id"]
gdf_nucleus_boundaries.index = gdf_nucleus_boundaries["sample_cell_id"]
gdf_cell_circles.index = gdf_cell_circles["sample_cell_id"]


# Set any arbitrary planar CRS to avoid warnings --------
gdf_cell_boundaries = gdf_cell_boundaries.set_crs("EPSG:3857")
gdf_nucleus_boundaries = gdf_nucleus_boundaries.set_crs("EPSG:3857")
gdf_cell_circles = gdf_cell_circles.set_crs("EPSG:3857")


# Save gdf objects ------------------------------------------------------
gdf_nucleus_boundaries.to_parquet(res_dir + "geodataframes/nucleus_boundaries.parquet")
gdf_cell_boundaries.to_parquet(res_dir + "geodataframes/cell_boundaries.parquet")
gdf_cell_circles.to_parquet(res_dir + "geodataframes/cell_circles.parquet")


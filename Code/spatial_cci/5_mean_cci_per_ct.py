#
# Created on Tue Nov 04 2025
#
# Copyright (c) 2025 Vivien Goepp
#

# --------------------------------------------------------------------------------------------
# This script computes and visualizes mean ligand–receptor (LR) interaction scores per
# cell-type pair and tissue core from Xenium kidney biopsy spatial transcriptomics data.
# It focuses on the APP–CD74 interaction across proximal tubule subtypes and immune cells.
#
# Datasets used/loaded:
# - Xenium processed AnnData object: xenium/xenium_processed.h5ad
# - Spatial cell and nucleus boundaries:
#       geodataframes/{cell_boundaries, nucleus_boundaries, cell_circles}.parquet
# - Precomputed CCI interaction scores (from all tissue cores):
#       local_cci_results/all_cores/interaction_scores_*.parquet
#
# Datasets/objects generated and saved:
# - mean_cci_score_df.csv: average CCI score per (core, cell-type pair, LR pair)
#
# Figures generated:
# mean_per_ct/mean_cci_per_ct_violin_APP_CD74_pvals.png
#   - Violin and strip plot of mean APP–CD74 CCI scores per sender cell type (aPT, PT1, PT2)
#   - Points colored by disease group
#   - Statistical annotations from Mann–Whitney tests between senders
#   - Title: “CCI score APP-CD74 with Immune as target (average per core)”
#
# The analysis pipeline:
# 1. Load spatial and transcriptomic data; harmonize cell-type labels (merge aPT/PEC → aPT).
# 2. Load precomputed cell–cell interaction scores and subset to relevant cores.
# 3. Compute mean interaction scores per (tissue core, cell-type pair, LR pair).
# 4. Annotate each core with its disease group using adata.obs["disease_id"].
# 5. Filter to cell-type pairs of interest:
#       - aPT→Immune, PT1→Immune, PT2→Immune
# 6. Restrict to the ligand–receptor pair “APP–CD74”.
# 7. Generate violin + strip plots stratified by sender cell type and disease group.
# 8. Perform Mann–Whitney tests between sender types and annotate statistical significance.
# 9. Save aggregated interaction table and publication-quality figure.
# --------------------------------------------------------------------------------------------


from statannotations.Annotator import Annotator
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import os
import geopandas as gpd
from matplotlib.colors import ListedColormap
import warnings
import sys
warnings.filterwarnings(
    "ignore", message="nopython is set for njit and is ignored")

# Load helper functions ----------------------------
sys.path.append(os.path.abspath("../utils/"))
from running_on_slurm import running_on_slurm
from spatial_cci_helpers import *


# Set directories ---------------------------
data_dir = "/data/zd775156/kidney_biopsy/data/"
xenium_path = "/data/Kramann_lab/From_Sidrah/new_data/output-XETG00229__0042348__Region_3__20250416__104148/"
res_dir = "/data/Kramann_lab/kidney_biopsy/data/processed/xenium/"
fig_dir = "/data/Kramann_lab/From_Vivien/kidney_biopsy/spatial_cci/"


# -------------------------------------------------------------------------------
# Load and prepare data ---------------------------------------------------------
# -------------------------------------------------------------------------------

# Load adata ---------------------------------------------------
adata = sc.read_h5ad(res_dir + "xenium/xenium_processed.h5ad")

# Replace aPT/PEC with aPT in adata.obs["cell_type_manual"]
adata.obs["cell_type_manual"] = adata.obs["cell_type_manual"].cat.rename_categories(
    {"aPT/PEC": "aPT"}
)

adata.obs["sample_cell_id"] = adata.obs["TMA_ID"].astype(str) + "_" + adata.obs["cell_id"].astype(str)
adata.obs.index = adata.obs["sample_cell_id"]

# Load xenium with boundaries ---------------------------------------------------
gdf_cell_boundaries = gpd.read_parquet(
    res_dir+"geodataframes/cell_boundaries.parquet")
gdf_nucleus_boundaries = gpd.read_parquet(
    res_dir+"geodataframes/nucleus_boundaries.parquet")
gdf_cell_circles = gpd.read_parquet(
    res_dir+"geodataframes/cell_circles.parquet")


# Replace aPT/PEC with aPT in gdf_cell_boundaries using cat.rename_categories
gdf_cell_boundaries["cell_type_manual"] = gdf_cell_boundaries["cell_type_manual"].cat.rename_categories(
    {"aPT/PEC": "aPT"}
)
gdf_nucleus_boundaries["cell_type_manual"] = gdf_nucleus_boundaries["cell_type_manual"].cat.rename_categories(
    {"aPT/PEC": "aPT"}
)

# Unname index
gdf_cell_boundaries.index.name = None
gdf_nucleus_boundaries.index.name = None

# Get f2 mapping and other info from gdf_cell_boundaries ----------------
cores = gdf_cell_boundaries["TMA_ID"].unique().tolist()
cell_types = gdf_cell_boundaries["cell_type_manual"].cat.categories.tolist()
factor2_map = gdf_cell_boundaries.set_index("sample_cell_id")["Factor2_ulm"]

# Create custom palette for cell_types ---------------------------------
colors_tab20 = plt.get_cmap("tab20")(range(20))
colors_husl = sns.color_palette('husl', len(cell_types) - 20)
my_pal = list(colors_tab20) + colors_husl
my_cmap = ListedColormap(my_pal, name="my_custom")


# If it does not exist yet, create and save lr_pairs.txt --------------

# Load CCI results ------
if running_on_slurm():
    end = 198361
else:
    end = 25000

print("Running on SLURM:", running_on_slurm())
print("End:", end)

cci_score_df = load_continuous_chunks(
    "interaction_scores",
    folder=fig_dir + "local_cci_results/all_cores/",
    start=0, end=end
)

# convert source and target as strings in cci_score_df
# Replace aPT/PEC with aPT in cci_score_df
cci_score_df["source"] = cci_score_df["source"].replace("aPT/PEC", "aPT")
cci_score_df["target"] = cci_score_df["target"].replace("aPT/PEC", "aPT")

cci_score_df["lr_pair"] = cci_score_df["gene_a"] + "_" + cci_score_df["gene_b"]
cci_score_df["source"] = cci_score_df["source"].astype(str)
cci_score_df["target"] = cci_score_df["target"].astype(str)

cci_score_df["ct_pair"] = cci_score_df["source"] + "_" + cci_score_df["target"]

# Define core by splitting root at "_" and taking first part
cci_score_df["core"] = cci_score_df["root"].str.split("_").str[0]

# Compute average cci_score_df per core and ct_pair ------------------------
ave_cci_score_df = cci_score_df.groupby(["core", "ct_pair", "lr_pair"])["score"].mean().reset_index()

# Add disease_id to cci_score_df by mapping from adata.obs["disease_id"]
adata_disease_map = adata.obs[["TMA_ID", "disease_id"]].drop_duplicates().set_index("TMA_ID")["disease_id"]
ave_cci_score_df["disease_id"] = ave_cci_score_df["core"].map(adata_disease_map)

# Save ave_cci_score_df
ave_cci_score_df.to_csv(os.path.join(fig_dir, "local_cci_results", "mean_cci_score_df.csv"), index=False)

# Define ct_pairs_of_interest ----------------------------------------
ct_pairs_of_interest = ["aPT_Immune", "PT1_Immune", "PT2_Immune"]
sources_of_interest = [ct.split("_")[0] for ct in ct_pairs_of_interest]

# Subset to ct_pairs_of_interest
ave_cci_score_df = ave_cci_score_df[ave_cci_score_df["ct_pair"].isin(ct_pairs_of_interest)]

# Subset to lr_pair = "APP_CD74"
ave_cci_score_df = ave_cci_score_df[ave_cci_score_df["lr_pair"] == "APP_CD74"]

ave_cci_score_df["disease_id"] = ave_cci_score_df["disease_id"].astype(str)

ave_cci_score_df["source"] = ave_cci_score_df["ct_pair"].str.split("_").str[0]


# Plot mean CCI score per ct_pair ------------------------------------------

sns.set(style="whitegrid")
plt.figure(figsize=(6, 4))
order = sources_of_interest
# draw violins with no inner elements; we'll remove the fill below to keep only outlines
ax = sns.violinplot(
    x="source",
    y="score",
    data=ave_cci_score_df,
    order=order,
    inner=None,
    linewidth=1.2,
    color="white",
)
# map color to disease_id using hue and a discrete palette
disease_ids = ave_cci_score_df["disease_id"].unique()
palette = dict(zip(disease_ids, sns.color_palette("tab10", n_colors=len(disease_ids))))

sns.stripplot(
    x="source",
    y="score",
    data=ave_cci_score_df,
    order=order,
    hue="disease_id",
    palette=palette,
    size=3,
    jitter=True,
    alpha=0.8,
    dodge=False
)

# place legend outside the plot
ax.legend(title="",
          bbox_to_anchor=(1.05, 0.75),
          borderaxespad=0.)

ax.set_xlabel("Sender")
ax.set_ylabel("")
ax.set_title("CCI score APP-CD74 with Immune as target (average per core)")
plt.tight_layout()

# Add pvalues
pairs = [("aPT", "PT1"), ("aPT", "PT2")]
annotator = Annotator(ax, pairs, data=ave_cci_score_df,
                      x="source", y="score", order=order)
annotator.configure(test="Mann-Whitney", text_format="star",
                    verbose=False)
annotator.apply_and_annotate()

os.makedirs(os.path.join(fig_dir, "mean_per_ct"), exist_ok=True)
if running_on_slurm():
    plt.savefig(os.path.join(fig_dir, "mean_per_ct",
                "mean_cci_per_ct_violin_APP_CD74_pvals.png"), dpi=300)
plt.show()

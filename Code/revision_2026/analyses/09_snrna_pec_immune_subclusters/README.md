# 09 — snRNA-seq PEC and immune subclusters: aPEC trajectory, aPEC differential expression, transfer to Xenium

## Purpose
Re-clusters the parietal epithelial cells (PECs, 7,004 nuclei) and the immune cells (19,486 nuclei) of the snRNA-seq
atlas on subset-specific embeddings, scores PEC subclusters with the crescent signature of glomerular niche 20
(folder `08_glomerular_niche20_crescent`), orders the PEC subclusters 0, 2, 4, 9 and 11 along a VIA/PHATE trajectory
towards the crescent-like aPEC subcluster 9, tests subcluster 9 against 0, 2, 4 and 11 (MAST, then GSEA), and transfers
the PEC (Leiden 0.8) and immune (Leiden 0.6) subcluster labels onto Xenium cells with scANVI query mapping.

## Analysis workspaces
| sub-folder here | workspace (outputs live there) |
|---|---|
| `pec/`, `pec/de_9vs0241/`, `pec/xenium_pec_mapping/` | `$KB_ANALYSIS/pec_resubcluster_20260717/` (`new_pec_annotation/`, `trajectories_leiden08/C_0to9_5cl/`, `de_9vs0241/de_results/9_vs_0_2_4_11/`, `xenium_pec_mapping/`) |
| `subcluster_sweep/` | `$KB_ANALYSIS/subcluster_level2_20260728/<celltype>/` (immune: `Immune/`) |
| `immune_mapping/` | `$KB_ANALYSIS/immune_mapping_leiden06/` (`prediction/`, `cell_boundary_plots_immune_0_60/`); the 2025 notebook wrote `$KB_RESULTS/xenium/leon/xenium_with_niches.h5ad` |

## Run order (source `config/env.sh` first; `.sh` files are SLURM wrappers, submit with `sbatch`)
| step | script | conda env | wrapper | output |
|---|---|---|---|---|
| P1 | `pec/build_pec_subset.py` | base | none | `pec_new_20260717_raw.h5ad` (level2 == PEC, counts, provenance of the earlier annotation) |
| P2 | `pec/pec_leiden_sweep.py` | base (`$CONDA_ROOT/bin/python`) | none (log `pec_leiden_sweep.log`) | `new_pec_annotation/pec_leiden_sweep.h5ad`, `res_<r>/` UMAPs, marker dotplots and tables; crescent ULM score |
| P3 | `pec/trajectory_leiden08.py --tags C_0to9_5cl` | ViaEnv | `pec/run_trajectory_leiden08.sh` | `trajectories_leiden08/C_0to9_5cl/` (subset UMAP/diffmap/PHATE, DPT, VIA per embedding; `subset_all_pseudotimes.h5ad`); published run: job 328231 |
| P4 | `pec/genetrend_niche_markers.py` | ViaEnv | `pec/run_genetrend_niche.sh` | `.../phate/genetrend_UP_top5.pdf`, `genetrend_DOWN_bottom5.pdf` (rebuilds VIA, see caveats) |
| P5 | `pec/heatmap_C5_phate_zbin.py` | ViaEnv (original env not recorded; reproduces in ViaEnv) | none | `.../phate/pseudotime_heatmap_zbinmean_HVGs.{pdf,png}` |
| P6 | `pec/prep_de_9vs0241.py` | base (not recorded; scanpy only) | none | `de_9vs0241/pec_de_input_leiden08_C5.h5ad` (2,543 nuclei, counts) |
| P7 | `pec/de_9vs0241/de_analysis_9vs0241.R` | de_env | `pec/de_9vs0241/run_de_9vs0241.sh` | `de_results/9_vs_0_2_4_11/DE/markers_MAST_9_vs_0_2_4_11.csv` (+ unfiltered volcano and GSEA) |
| P8 | `pec/de_9vs0241/filter_and_volcano_pecfiltered.py` | nichecompass | none | `..._pecfiltered.csv` (105 canonical non-PEC lineage markers removed), volcanos |
| P9 | `pec/de_9vs0241/gsea_pecfiltered.R` | de_env | `pec/de_9vs0241/run_gsea_pecfiltered.sh` | `GSEA_pecfiltered/{Hallmark,Reactome}/{gsea_results.csv,dotplot.pdf}` (published run: job 328394) |
| P10 | `pec/xenium_pec_mapping/01_map_leiden04_to_xenium.py --label leiden_0.8` | kpmp_repurposing | `pec/xenium_pec_mapping/run_map.sh --label leiden_0.8` | `pec_leiden08_predictions.csv`, `xenium_pec_leiden08.h5ad`, models |
| P11 | `pec/xenium_pec_mapping/02_plot_pec_boundaries.py --label leiden_0.8` | base (original run; kpmp_repurposing gives identical plots) | none (log `plot_leiden_0.8.log`) | `boundary_plots_leiden08/per_core/`, `per_glomerulus/` (one file per core) |
| I1 | `subcluster_sweep/subcluster_sweep.py --celltype Immune` | pathfire | `subcluster_sweep/run_sweep_array.sh` (array task 6 = Immune) | `Immune/Immune_leiden_sweep.h5ad`, `Immune/res_<r>/` |
| I2 | `immune_mapping/nichecompass_filtered_tma.ipynb` (2025) | nichecompass | none | `$KB_RESULTS/xenium/leon/xenium_with_niches.h5ad` (three-region NicheCompass object, `cell_type_manual`) |
| I3 | `immune_mapping/map_immune_leiden06_xenium.py` | kpmp_repurposing | `immune_mapping/run_map_immune06.sh` | `prediction/xenium_immune_with_prediction.h5ad` (+ csv) |
| I4 | `immune_mapping/plot_immune06_cell_boundaries.py` | nichecompass | none | `cell_boundary_plots_immune_0_60/{predictions,cell_type_manual,side_by_side...}/` (one file per core) |

P2 also needs step 1 of folder 08 (niche-20 MAST table); P10/P11 and I3/I4 need the Xenium objects of folders 01/02.
The wrappers `run_trajectory_leiden08.sh`, `run_genetrend_niche.sh`, `run_map.sh` and `run_sweep_array.sh` originally
called the workspace copy of their script by a relative path; in this repository they call the shipped copy (patch
rules in `build/patch_rules/fig_niches.py`). `filter_and_volcano_pecfiltered.py` likewise reads the `PANELS` literal
from the shipped `subcluster_sweep/subcluster_sweep.py` (identical literal).

## Inputs
| input | producer |
|---|---|
| `$KB_DATA/charlotte_new_data/adata_20260717_with_clinical.h5ad` (snRNA-seq atlas, `level2` labels, `patient_ID`) | `analyses/02_xenium_celltyping/18_add_clinical_to_snRNA_ref.py` (on the colleague's annotated object) |
| `$KB_DATA/scVI_kidney_biopsy_with_clinical.h5ad`, `$KB_PROJECT/trajectory_pec/pec_annotated.h5ad` | earlier atlas / earlier PEC annotation; read by P1 only to record provenance (not shipped) |
| `$KB_ANALYSIS/celltype_concordance/outputs/crescent_check/de/markers_manualglom10.csv` | `analyses/08_glomerular_niche20_crescent/de_crescent_pec_mast.R` (niche-20 PEC MAST; P2 uses its positive genes at padj < 0.05 as the crescent ULM signature) |
| `$KB_ANALYSIS/celltype_concordance/outputs/xenium_labeled.h5ad`, `.../shared_genes.csv` | `analyses/02_xenium_celltyping` (query cells with `level2_refined == "PEC"`, 2,577 cells; the 468 panel genes) |
| `$KB_XENIUM_RAW/output-<region>/{cell,nucleus}_boundaries.parquet` | external: Xenium Onboard Analysis output |
| `$KB_RESULTS/xenium/leon/xenium_with_niches.h5ad` | I2 (query: `cell_type_manual` in {Immune, PL/B, B}, 32,158 cells) |
| `$KB_DATA/processed/xenium/geodataframes/{cell,nucleus}_boundaries.parquet` | external: cell polygons of the 2025 three-region object (colleague's processing; carry `TMA_ID`, `cell_type_manual`) |

## Outputs and the panels that use them
| output | panel (entry script) |
|---|---|
| `xenium_pec_mapping/boundary_plots_leiden08/per_glomerulus/<K50>.png`, `<K18>.png` | Fig. 6g left (`figures/Fig6/Fig6g_pec_subclusters_xenium.py`) |
| `cell_boundary_plots_immune_0_60/predictions/predictions_<K50>.png`, `<K18>.png` | Fig. 6g right (`figures/Fig6/Fig6g_immune_subclusters_xenium.py`) |
| `trajectories_leiden08/C_0to9_5cl/phate/phate_cluster_pt_crescent.pdf` (panels 1 and 4) | Fig. 6h (`figures/Fig6/Fig6h_pec_phate_clusters_crescent.py`) |
| `.../phate/via_lineage_curves.pdf` (right panel) | Fig. 6i (`figures/Fig6/Fig6i_pec_via_pseudotime_curve.py`: pointer only, see caveats) |
| `.../phate/genetrend_UP_top5.pdf`, `genetrend_DOWN_bottom5.pdf` | Fig. 6j left (`figures/Fig6/Fig6j_genetrends.py`: pointer only) |
| `.../phate/pseudotime_heatmap_zbinmean_HVGs.pdf` | Fig. 6j right (`figures/Fig6/Fig6j_pseudotime_heatmap.py`) |
| `de_results/9_vs_0_2_4_11/DE/volcano_top10bottom10_9_vs_0_2_4_11_pecfiltered.pdf` | Fig. 6k (`figures/Fig6/Fig6k_volcano_aPEC_vs_PEC.py`) |
| `de_results/9_vs_0_2_4_11/GSEA_pecfiltered/Reactome/dotplot.pdf` | Fig. 6l (`figures/Fig6/Fig6l_gsea_reactome_dotplot.R`) |
| `new_pec_annotation/res_0.8/umap_res0.8.pdf`, `markers_dotplot_res0.8.pdf` | Supp. Fig. 18c-d (`figures/S18/S18c-d_pec_subclusters_umap_markers.py`) |
| `subcluster_level2_20260728/Immune/res_0.6/umap_res0.6.pdf`, `markers_dotplot_res0.6.pdf` | Supp. Fig. 18e-f (`figures/S18/S18e-f_immune_subclusters_umap_markers.py`) |

All entry scripts re-draw from the saved objects (layouts, labels, pseudotimes, Wilcoxon results, DE and GSEA tables);
none recomputes Leiden, UMAP, PHATE, VIA, MAST or GSEA. Except for 6i and 6j (trends), every regenerated panel matched
its original file pixel for pixel; the Fig. 6g plots differ only in the title, which now names the core by its public
id instead of its internal label (see `build/verification/niches.csv`).

## Caveats (read before reusing)
- **Leiden labels are the reference, not the code.** Subcluster numbers come from the saved objects
  (`pec_leiden_sweep.h5ad` `leiden_0.8`, `Immune_leiden_sweep.h5ad` `leiden_0.6`); a re-run with other package versions
  can renumber them, and every downstream step (trajectory subset {0, 2, 4, 9, 11}, the 9 vs 0/2/4/11 contrast, the
  Xenium transfer) refers to these numbers.
- **VIA is not seed-deterministic.** Repeated runs of the same configuration differ (job logs 328227/328230/328231;
  for C_0to9_5cl on PHATE, Spearman(pseudotime, crescent score) = 0.114 / 0.142 / 0.131). Fig. 6h/6i use job 328231
  (stored in `subset_all_pseudotimes.h5ad`). The VIA object was never saved, so the lineage curve of Fig. 6i cannot be
  redrawn, and P4 (Fig. 6j trends) rebuilt VIA: its pseudotime correlates rho = 0.985 with the stored one.
- **The PHATE-based VIA pseudotime agrees poorly with DPT** (Spearman +0.004 for C_0to9_5cl/phate; UMAP- and
  diffmap-based VIA pseudotimes of the same cells correlate 0.20 and 0.03 with it, `_compare/`), and its root cluster c0
  is not the earliest cluster (median pseudotime rank 1 of 4; c2 is earlier). The terminal cluster c9 is the latest.
- **Patched pyVIA.** ViaEnv carries pyVIA 0.2.4 with a local fix in `pyVIA/utils_via.py` (zip -> tuple(zip), needed
  for scipy >= 1.13; original kept as `utils_via.py.ORIG_backup`). ViaEnv also needs `LD_LIBRARY_PATH=$CONDA_PREFIX/lib`.
- **Fig. 6k-l**: MAST is a cell-level test (286 vs 2,257 nuclei) with `patient_ID` as a fixed latent covariate, not a
  pseudobulk or mixed model; P values are optimistic. 105 canonical non-PEC lineage markers were removed before the
  volcano and the GSEA (protected sets in P8). The model has no cellular-detection-rate covariate although cluster 9
  nuclei are 2.5x deeper (nCount 4,840 vs 1,946 as recorded; 1.85x genes per nucleus), so the ranking by adjusted P is
  depth-confounded (workspace STATUS.md, 24 Jul 2026); rank by fold change, as Fig. 6k does. Fig. 6l is a GSEA (fgsea
  multilevel, genes ranked by avg_log2FC, BH), not an over-representation test; 27 of its 93 Reactome terms, including
  all ten plotted, sit at the fgsea P floor (1e-10), so the colour scale of the panel carries no information.
- **Supp. Fig. 18e** is Leiden 0.6 (18 clusters); the caption says 0.8. Subcluster names of Fig. 6g (immune) and
  Supp. Fig. 18e were added in Illustrator; in the immune palette (tab20) subcluster 15 (DC) is light grey, close to the
  grey of unmapped cells.
- **Different Xenium cell sets in Fig. 6g.** The PEC transfer uses the current Xenium object (`xenium_labeled.h5ad`,
  `level2_refined == "PEC"`, five regions, Xenium-bundle polygons); the immune transfer uses the 2025 three-region
  NicheCompass object (`cell_type_manual` Immune/PL/B/B) and its own polygon tables. The immune reference was fed to
  scVI/scANVI as log1p values while the query is raw counts (stated in `map_immune_leiden06_xenium.py`).
- **Immune join by `cell_id` only.** `plot_immune06_cell_boundaries.py` joins predictions to polygons by `cell_id`
  without the region; 12 polygon rows share a `cell_id`, so 2 predictions are also painted onto a second cell (cores
  XEN006, XEN013, XEN091, XEN098; not the cores of Fig. 6g).
- **AKI samples.** The PEC and immune snRNA-seq objects still contain nuclei of K95 and K98 (AKI donor).
- **Identifiers.** Per-core output files of P11 and I4 are named after, and titled with, the internal core label
  (`TMA_ID`, which contains a pathology number); the objects carry `patient_ID`. Do not share those files; the entry
  scripts select cores with `config/cores.py::core_label("XEN###")` and name outputs by public id and K-number.

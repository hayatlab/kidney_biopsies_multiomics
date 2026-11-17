# 03 — Xenium: Factor 1 per cell, Voronoi neighbourhoods, local cell-cell communication

## Purpose
Projects the MCFA factors (the factor of interest is Factor 1) of the snRNA-seq model onto every Xenium cell
(decoupler ULM on per-view loadings), builds cell-boundary polygons, distance-bounded Voronoi graphs and the 2-hop
neighbourhood composition of every cell (by cell type and by NicheCompass niche), runs a local CellPhoneDB v5 analysis
inside each 2-hop neighbourhood, scores NF-kB activity (PROGENy) per cell, rebuilds a per-core eGFR table, and draws the
Xenium panels of Fig. 3, Extended Data Fig. 4 and 6 and Supplementary Fig. 12. It is a re-implementation, on the final
5-run Xenium object (258,375 cells, 119 cores, 16 `level2_refined` cell types), of the original paper code
(`Code/plot_egfr_f2_on_spatial`, `Code/nh_analysis`, `Code/spatial_cci`, `Code/cell_boundaries`), with Factor 1
replacing Factor 2. Human GRCh38, HGNC symbols, seed 42.

Function and file names still carry the panel letters of an earlier draft: `fig3d` = Fig. 3c, `fig3g` = Fig. 3f,
`fig3h` = Fig. 3g, `fig3i` = Fig. 3h, `fig3j` = Fig. 3i, `ed3b` = ED 4c, `ed3f` = ED 4f, `suppl7` = Supp. Fig. 12.

## Analysis workspace
`$KB_ANALYSIS/vivien_takeover/` (scripts were in `scripts/`; SLURM logs in `logs/`). Outputs in `outputs/`
(tables, parquet files, `geodataframes/`, `neighborhoods/`, `local_cci_results/`, `nets/`, `metadata/`) and `figures/`
(`xenium/`, `xenium/cci_lr_panels/`, `cci/`, `milo_donor_composition/`).

## Run order (source `config/env.sh` first; the wrappers are SLURM scripts that `cd` into the workspace)
| step | script | conda env | wrapper | main output |
|---|---|---|---|---|
| 1 | `01_score_factors_xenium.py` | nichecompass | `01_score_factors_xenium.sh` (200 GB) | `outputs/xenium_factor_scores.parquet` (Factor 1-5, ULM and weighted sum, per cell; obs columns) |
| 2 | `10_build_authoritative_egfr.py` | agent_env (openpyxl) | none | `outputs/authoritative_core_egfr.csv` (per-core eGFR from the source workbooks) |
| 3 | `02_build_geodataframes.py` | spatialdata | `02_build_geodataframes.sh` | `outputs/geodataframes/{cell_boundaries,nucleus_boundaries,cell_circles}.parquet` |
| 4 | `03_voronoi_neighborhoods.py` (`--threshold 10`) | spatialdata | `03_voronoi_neighborhoods.sh` | `outputs/neighborhoods/voronoi_graphs_10.pkl`, `nh_2hop_compo_10.csv`, `nh_2hop_niche_compo_10.csv` |
| 5 | `13_crop_justification.py` | not recorded | none | `figures/xenium/crop_justification_factor1.pdf`, `13_crop_justification_meta.json`: the evidence for the ribbon crop [-7.25, 5.25] that 05 and 09 hard-code |
| 6 | `05_xenium_figures.py --panels ...` | spatialdata | `05_xenium_figures.sh` | `figures/xenium/fig3d_*`, `fig3g_*`, `fig3h_*`, `ed3b_*`, `ed3f_*`, `suppl7_*` |
| 7 | `08_run_local_cci.py --start --end` | cci_env | `08_run_local_cci_array.sh` (array 0-15, 16 chunks of 16,149 roots, 180 GB each) | `outputs/local_cci_results/interaction_{scores,pvalues}_<start>_<end>.parquet` |
| 8 | `09_cci_postprocess.py` (full run ~70 min; `--reuse` re-plots from the saved table) | cci_env | `09_cci_postprocess.sh` (400 GB) | `local_cci_results/cci_score_df.parquet`, `cci_summary_all_pairs.parquet`; `figures/cci/fig3i_<sender>_to_<target>.pdf` |
| 9 | `34_cci_panel.py --pairs <18 pairs> --outsub cci_lr_panels` | spatialdata | `34b_cci_top_pairs.sh` (300 GB) | `figures/xenium/cci_lr_panels/fig3j_<pair>_apt_pt.{pdf,png,svg}`, `outputs/fig3j_<pair>_apt_pt_mean_cci_per_core.csv`, `fig3j_pairs_summary.csv` |
| 10 | `26_score_nfkb_xenium.py --fetch-nets` (login node, needs internet), then `26_score_nfkb_xenium.py` | nichecompass | `26_score_nfkb_xenium.sh` (scoring only) | `outputs/nets/` (PROGENy top 500, CollecTRI, retrieval date), `outputs/xenium_nfkb_scores.parquet` |
| 11 | `13_crop_justification.py --score NFkB_progeny --scores-parquet outputs/xenium_nfkb_scores.parquet --stem nfkb_progeny --bins-per-1-99 46` (and the same for `NFkB_collectri`) | not recorded | none | `figures/xenium/13_crop_justification_meta_nfkb_{progeny,collectri}.json` |
| 12 | `27_ed3f_nfkb.py` | spatialdata | none (SLURM job 329585) | `figures/xenium/fig3h_celltype_vs_nfkb*`, `ed3f_niche_vs_nfkb*`, `nfkb_score_context.pdf` |
| m1 | `22_sample_metadata_table.py` | not recorded | none | `figures/milo_donor_composition/sample_metadata_milo.csv` (per-sample table of the snRNA object) |
| m2 | `24_apply_corrected_metadata.py` | not recorded | none | `outputs/metadata/sample_metadata_corrected.csv` (+ `metadata_corrections_log.csv`) |

Steps 1-12 need the final `xenium_labeled.h5ad` (folder 02, after step 17b). Steps m1-m2 build the corrected per-sample
snRNA metadata that folders 04 and 05 read. Runs that drew the published panels (the `--panels` arguments are not
logged; each log lists the panels it drew): `05_xenium_figures.sh` jobs 329333/329334 (2026-09-02; ED 4c), 329425
(2026-09-07; Fig. 3c, 3g, ED 4f), 329429 (Fig. 3f) and 329443 (Supp. Fig. 12); `09_cci_postprocess.sh --reuse`
(Fig. 3h; job 329451); `34b_cci_top_pairs.sh` (Fig. 3i; job 330299); `python 27_ed3f_nfkb.py` (ED 6a; job 329585).

## Inputs
| input | producer |
|---|---|
| `$KB_ANALYSIS/celltype_concordance/outputs/xenium_labeled.h5ad` | `analyses/02_xenium_celltyping` (steps 02 -> 17b) |
| `$KB_ANALYSIS/loadings~hvg~5.csv` (5 factors x 10 views, 2026-08-04) | MCFA of the snRNA-seq pseudobulk; upstream, not in this repository |
| `$KB_DATA/charlotte_new_data/adata_20260717_with_clinical.h5ad` (level2 -> level1 map; m1) | `analyses/02_xenium_celltyping/18_add_clinical_to_snRNA_ref.py` |
| `$KB_XENIUM_RAW/output-<run>/{cell_boundaries,nucleus_boundaries,cells}.parquet` (5 runs) | external: Xenium Onboard Analysis output |
| `$KB_DATA/metadata/shared_for_xenium_analysis.xlsx`, `shared_for_xenium_analysis_with_kpmp_egfr.xlsx` | external: clinical metadata workbooks |
| `outputs/cpdb_db/cellphonedb.zip` (CellPhoneDB v5.0.0) | external: downloaded by hand, not scripted |
| PROGENy (human, top 500 per pathway) and CollecTRI via OmniPath | external: fetched at run time by `26 --fetch-nets` |
| `$KB_ANALYSIS/apt_immune_lr/outputs/apt_immune_pairs_panel.csv`, `rela_downstream_ligands.csv` | `apt_immune_lr/` (not shipped); used only for the cross-pair overview of 34, not for Fig. 3i |
| `outputs/milo/milo_cell_scores<tag>.parquet` (m1) | Milo analysis of the same workspace (`04_run_milo.py`, not shipped: not used by a panel here) |
| `sample_metadata_milo_korrigiert - sample_metadata_milo_korrigiert.csv` (m2) | external: the m1 table corrected by hand by a colleague (2026-09-10) |
| private core lists (`KB_EXCLUDED_CORES`, tag `dropped_low_cells`; step 3) | private file, see `config/cores.py` |

## Outputs and the panels that use them
| output | panel (entry script) |
|---|---|
| `figures/xenium/fig3d_factor1_vs_egfr.pdf` (`05::fig3d`) | Fig. 3c (`figures/Fig3/Fig3c_factor1_per_core_vs_egfr.py`) |
| `figures/xenium/fig3g_voronoi_2hop.pdf` (`05::fig3g`) | Fig. 3f (`figures/Fig3/Fig3f_voronoi_2hop_neighbourhood.py`) |
| `figures/xenium/fig3h_celltype_vs_factor1.pdf` (`05::ribbon`) | Fig. 3g (`figures/Fig3/Fig3g_celltype_composition_vs_factor1.py`) |
| `figures/cci/fig3i_aPT_to_Immune.pdf` (`09::plot`) | Fig. 3h (`figures/Fig3/Fig3h_cci_apt_to_immune_vs_factor1.py`) |
| `figures/xenium/cci_lr_panels/fig3j_app_cd74_apt_pt.pdf` (`34::draw_pair`) | Fig. 3i (`figures/Fig3/Fig3i_app_cd74_apt_vs_pt.py`) |
| `figures/xenium/ed3b_2hop_celltype_heatmap.pdf` (`05::ed3b`) | ED 4c (`figures/ED4/ED4c_2hop_celltype_heatmap.py`) |
| `figures/xenium/ed3f_niche_vs_factor1.pdf` (`05::ribbon`) | ED 4f (`figures/ED4/ED4f_niche_composition_vs_factor1.py`) |
| `figures/xenium/fig3h_celltype_vs_nfkb.pdf` (`27`, `05::ribbon`) | ED 6a (`figures/ED6/ED6a_celltype_composition_vs_nfkb.py`) |
| `figures/xenium/suppl7_<core>.pdf` (`05::suppl7`) | Supp. Fig. 12a-b (`figures/S12/S12a-b_core_maps.py`) |
| `outputs/xenium_factor_scores.parquet` | also folder 06 (Supp. Fig. 10d) |
| `outputs/authoritative_core_egfr.csv` via `05::corrected_core_egfr` | also folder 05 (Supp. Fig. 13; Xenium part) |
| `outputs/metadata/sample_metadata_corrected.csv` | folders 04 (Supp. Fig. 2, ED 2e-g) and 05 (Supp. Fig. 13b) |
| `outputs/geodataframes/`, `outputs/neighborhoods/`, `outputs/local_cci_results/` | inputs of steps 4-12 |

## Caveats (read before reusing)
- **eGFR of the Xenium object is corrupted.** Its `EGFR` column broadcasts one group value over many donors (22 CTRL
  cores read 67.0, 15 FSGS cores 58.0). Every figure uses `05::corrected_core_egfr()` = the workbook values from step 2
  for 85 cores and the object value for the 32 cores of the two newer runs. A second reconciliation (identifier join)
  disagrees for 4 cores; Fig. 3c shows r = -0.405 (p = 6.0e-06, 117 cores), the identifier-join eGFR gives r = -0.328.
  Two cores (K121, K110) still carry the broadcast value 58.0.
- **Cores are not donors.** Fig. 3c and 3i are per core ("per sample" in the figure); several punches can come from one
  donor. The donor-level check of Fig. 3c (111 donor keys, r = -0.400, p = 1.4e-05) is only printed by `fig3d`; donor
  keys are approximate (H.number known for 31 donors only, the rest grouped by condition and eGFR).
- **Crops are a rule, not the data range.** The ribbons and the CCI LOESS use [-7.25, 5.25] on Factor 1 (largest run of
  0.25-unit bins with >= 500 neighbourhoods; 99.07% of 258,375 neighbourhoods; step 5) and the same rule on NF-kB
  ([-1.52, 1.43]; step 11).
- **Fig. 3f was drawn by an earlier revision of `fig3g`** (2026-09-07 00:37). The shipped default draws the Factor 1
  colour bar at 55% height (`legend_kwds` shrink/pad/aspect, added afterwards); the entry script passes
  `factor_cbar_kwds={}` to draw it as published. Neighbourhood choice and data are unchanged (same core, root cell and
  19-cell composition). `fig3g` draws the neighbourhood edges in set order, which follows Python's string-hash seed,
  so the anti-aliasing where edges cross changes from run to run (~0.005% of pixels; set `PYTHONHASHSEED` to fix it).
  All other panels of this folder are redrawn pixel-identically by their entry scripts (core names in the titles of
  Fig. 3f and Supp. Fig. 12 aside).
- **Core names.** Internal core labels (`TMA_ID`) contain pathology numbers; the code selects cores through
  `config/cores.py::core_label("XEN###")`. Two defaults still name cores by their internal label and therefore do not
  resolve on the de-identified object: the control core of `05::SUPPL7_CORES` ("Appollo d0 B5 CTRL" = XEN008) and
  `34::DEFAULT_CORES` ("Utrecht2 F5" = XEN119, "Utrecht2 F2" = XEN116). The figure entry scripts pass public ids instead.
  `05::fig3d` writes `fig3d_per_core.csv` and `34::main` writes `<stem>_mean_cci_per_core.csv`, both indexed by the
  internal label; the entry scripts do not write them.
- **Not scripted:** the CellPhoneDB v5.0.0 database download; PROGENy/CollecTRI come from OmniPath at run time (date in
  `outputs/nets/nets_meta.json`), so a later fetch can differ.
- **Fig. 3i / 29:** `29_fig3j_app_cd74.py` (earlier violin version with dPT) is superseded by `34_cci_panel.py` and is not
  shipped. Fig. 3i tests aPT vs PT as senders with a paired Wilcoxon over the 116 cores with both senders (median
  paired difference 0.655, bootstrap 95% CI 0.131-1.016, p = 0.0022).
- Supp. Fig. 12b "high eGFR" is 60 ml/min/1.73 m2 (control core XEN008); ED 4d's "FSGS high eGFR" core (K82) has 63.

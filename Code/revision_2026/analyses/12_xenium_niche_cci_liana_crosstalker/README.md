# 12 — Cell-cell interactions within Xenium niches (LIANA + CrossTalkeR)

## Purpose
Ligand-receptor interactions between the cell types of one spatial niche of the Xenium data, for two niches:

- **glomerular niche 20** (manual-glomerulus NicheCompass model, Leiden 1.0; the crescent niche): Fig. 6d (network) and
  Fig. 6e (ligand-receptor pairs with PEC as sender);
- **niche 1** of the whole-core NicheCompass model (Leiden 0.4): Extended Data Fig. 4k (alluvial plot, aPT as sender).

Analysis by Sidrah Maryam (run 22 Jul 2026; scripts `crosstalker.R`, `cross_down.R` in her run folders, the notebook
`Fig_6e.ipynb`). The files here are her scripts with the run as a parameter and the paths taken from `config/`
(20 Nov 2026); the computation is unchanged. Step 3 and the checks in the figure scripts were added on 20 Nov 2026.

Method: the raw Xenium counts of all cells of the niche are log-normalised (Seurat `LogNormalize`, scale factor 10,000);
`liana::liana_wrap(method = "cellphonedb")` scores every ligand-receptor pair for every ordered pair of cell types
(`level2_refined`); interactions with permutation P ≤ 0.05 are kept; the LR score is `lr.mean`, the mean of ligand
expression in the sender and receptor expression in the receiver. `CrossTalkeR::generate_report` (single condition,
`threshold = 0`) builds the cell-cell graph: edge weight = sum of the LR scores of a cell-type pair, edge `inter` =
number of interactions of the pair rescaled to 0.1–1.1. Human, HGNC symbols, 468-gene Xenium panel. In the objects the
label `dPT` is the population called aPT2 in the paper.

## Analysis workspaces
`$KB_RESULTS/xenium/sidrah/whole_data_onlyglom_20/` (niche 20) and `$KB_RESULTS/xenium/sidrah/overall_glom_04_level2/`
(niche 1): input object, `CCI_CrossTalkeR/` with the LIANA result, the CrossTalkeR objects and the report.
`$KB_ANALYSIS/sidrah_scripts_20261120/outputs/cci/` (`$KB_CCI_EXPORT`): the CSV tables of step 3.

## Run order (source `config/env.sh` first)
| step | script | environment | output |
|---|---|---|---|
| 1 | `Rscript 01_liana_cellphonedb.R <niche20\|niche1>` | R with Seurat, SingleCellExperiment, liana | `<run>/CCI_CrossTalkeR/<object>_liana.rds` |
| 2 | `Rscript 02_crosstalker_report.R <niche20\|niche1>` | R with liana, CrossTalkeR | `<run>/CCI_CrossTalkeR/LR_data.Rds`, `LR_data_final.Rds`, `Single_overall_glom04_<niche>.html` |
| 3 | `Rscript 03_export_tables.R` | conda env de_env (igraph, SeuratObject; CrossTalkeR not needed) | `$KB_CCI_EXPORT/<run>_lr_table.csv`, `_graph_edges.csv`, `_cells_per_type.csv` |

`cci_runs.R` holds the two runs (folder, object, labels) and is sourced by all scripts. `plot_top_lr_per_receiver.py`
is the plotting function of Fig. 6e. Set `KB_CCI_OUT` to send the output of steps 1–2 of a re-run to another folder,
and `KB_CCI_EXPORT` for step 3.

Published runs: steps 1–2 on 22 Jul 2026 (interactive R sessions; package versions were not recorded). Results:

| run | cells | cell types in the result | interactions tested / with P ≤ 0.05 | edges |
|---|---|---|---|---|
| niche20 | 241 (PEC 134, Immune 56, FIB 28, aPT 7, EC 7, POD 5, MES 4) | 6 (the 4 MES cells are below LIANA's minimum of 5 cells) | 5,832 / 1,863 | 36 |
| niche1 | 35,773 (Immune 10,899, aPT 6,793, FIB 6,207, PT 3,239, TAL 3,223, EC 2,450, PC 1,346, aPT2 743, …) | 16 | 32,424 / 20,590 | 256 |

## Inputs
| input | producer |
|---|---|
| `whole_data_onlyglom_20/xenium_with_niches_glom10.rds` (Seurat; 241 cells × 468 genes; counts, `level2_refined`) | by its content the cells with `latent_leiden_1.0 == 20` of the manual-glomerulus NicheCompass object (`analyses/01_xenium_preprocessing_niches`), converted to Seurat by Sidrah Maryam; the conversion script is not in this repository |
| `overall_glom_04_level2/xenium_niches_overall_with_glom_niche1.rds` (Seurat; 35,773 cells × 468 genes) | by its content the cells with `latent_leiden_0.4 == 1` of the whole-core NicheCompass object, converted the same way |

## Outputs and the panels that use them
| output | panel | entry script |
|---|---|---|
| niche 20 cell-cell graph (`LR_data_final.Rds`) | Fig. 6d | `figures/Fig6/Fig6d_cci_network_niche20.R` |
| niche 20 interaction table (`niche20_lr_table.csv`) | Fig. 6e | `figures/Fig6/Fig6e_pec_sender_top_lr_dotplot.py` |
| niche 1 interaction table (`LR_data_final.Rds`) | Extended Data Fig. 4k | `figures/ED4/ED4k_apt_sender_alluvial_niche1.R` |

## Re-run check (20 Nov 2026)
- Step 3 and the three figure scripts were run from this folder. Fig. 6e was re-drawn and shows the printed rows and
  dots (4 pairs, 16 dots, LR score 4.03–5.55). For Fig. 6d and Extended Data Fig. 4k the plotted tables were
  reproduced and checked against the printed panels (Fig. 6d: maximum weight 484.29 for PEC → PEC, `inter` 0.1–1.1;
  Extended Data Fig. 4k: the 16 flows give the stratum sizes of the published PDF), but the plots were not re-drawn.
- **Steps 1 and 2 were not re-run**: liana and CrossTalkeR are not installed in the environments available for this
  check. That the published results were computed on all cells of each niche was confirmed from the data: the ligand
  expression stored by LIANA equals the mean log-normalised expression over all cells of the sender (APP in the 6,793
  aPT cells of niche 1: 4.893133; APP in the 134 PECs of niche 20: 4.849409).

## Caveats (read before reusing)
- **Niche 20 is small and dominated by one disease**: 241 cells, 204 of them from IgAN cores; aPT (7 cells), EC (7) and
  POD (5) are at or just above the minimum cell number. Scores involving these cell types rest on a handful of cells.
- LIANA treats cells as the unit; there is no donor or core level in the P values, and no test compares niches or
  conditions. The panels show descriptive scores.
- The LR score is an expression mean, so it favours highly expressed genes (APP, COL4A2, SPARC) and is bounded by the
  468-gene panel: a pair can only be scored if ligand and receptor are both on the panel.
- **Extended Data Fig. 4k shows both ends of the ranking**: `plot_sankey` with a row limit draws the highest- and the
  lowest-scoring interactions. Of the 1,753 significant aPT-sender interactions the panel shows the 8 highest (all APP)
  and the 8 lowest (EFNB2, BMP2, ADAM15).
- The copy of `crosstalker.R` in the niche 1 run folder refers to an undefined variable (`cond`) and cannot run as
  saved; step 1 uses the niche 20 version for both runs (see the check above).
- The published runs were unseeded; step 1 now sets a seed. Permutation P values of a re-run can differ slightly,
  and with them the set of interactions near P = 0.05.
- The internal core labels in the input objects (`TMA_ID`, `H.number`, `patient_id`) are not used by this analysis and
  must not be exported.

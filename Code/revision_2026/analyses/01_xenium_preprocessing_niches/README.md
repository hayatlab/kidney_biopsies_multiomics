# 01 — Xenium preprocessing and NicheCompass niches

## Purpose
Builds the 5-region Xenium object, trains the NicheCompass spatial niche models and draws every niche figure:

- **Region 1 additions and label transfer** — the two Region 1 slides added during review get core identities and
  whitelist metadata (`process_region1_new.py`) and the `predic` labels of the 3-region reference by joint PCA +
  Harmony + KNN (`phaseB_label_transfer_merge.py`).
- **Whole-core model** — NicheCompass on all 5 regions (`nichecompass_5regions.py`, run `single_sample/22042026_221229`,
  the parent run that carries `in_glom`), retrained with the `in_glom` annotation carried along
  (`nichecompass_overall_with_glom.py`, run `overall_with_glom/14052026_211003`; Leiden 0.4 = 10 niches, the
  niches of Fig. 3 / Extended Data Fig. 4).
- **Glomerular model** — NicheCompass on the cells inside manually drawn glomerulus polygons
  (`annotate_in_glom_new_samples.py` + `nichecompass_manual_glom_niche.py`, run `manual_glom_niche/14052026_211003`),
  clustered at Leiden 1.0 (21 niches; niche 20 = the crescent niche of Fig. 6).
- **Downstream battery** — `nichecompass_5regions_downstream.py` (composition, niche x cell-type heatmap, GP and
  Hallmark/Reactome activity, condition / eGFR / IFTA enrichment, ROI flags, within-niche cell-type markers).
- **Level-2 re-plots** (`nichecompass_replot_level2/`) — after the refined Xenium cell typing (`level2_refined`,
  analyses/02_xenium_celltyping) existed, the cell-type-dependent figures were re-drawn without retraining or
  re-clustering, the clinical-metadata figures were regenerated after the 2026-07-22 metadata join fix, and the
  Hallmark "dedup" heatmaps were recomputed on log-normalised counts.

## Analysis workspace
Scripts: `$KB_ANALYSIS/` (top level) and `$KB_ANALYSIS/nichecompass_replot_level2/` (sub-folder of the same name here).
Model artifacts: `$KB_PROJECT/niche_compass_results/artifacts/{single_sample,overall_with_glom,manual_glom_niche}/<timestamp>/`.
Results: `$KB_RESULTS/xenium/leon/new_data/` (current) and `$KB_RESULTS/xenium/leon/archived/` (cell-type `predic` era).

## Run order
Source `config/env.sh` first. The wrappers are SLURM scripts. Submit them from a folder that contains `logs/`,
because their `--output` paths are relative. Each one activates its conda env itself.

| step | script (wrapper) | conda env | what it does |
|---|---|---|---|
| 1 | `process_region1_new.py` (no wrapper saved) | not recorded (needs scanpy + shapely) | core identity + whitelist metadata for the two Region 1 h5ads -> `$KB_ANALYSIS/annotated_00*_Region_1.h5ad` |
| 2 | `phaseB_label_transfer_merge.py` (`phaseB.sh`) | omicverse | `predic` transfer (joint PCA, Harmony by sample, KNN k=15) and 5-region merge -> `$KB_XENIUM_RAW/h5ad_files_only_transcript/marker_genes/clustered_annotated_final_file_5regions.h5ad` |
| 3 | `nichecompass_5regions.py` (`nichecompass_5regions.sh`) | nichecompass | 5-region model, Leiden 0.4 -> `artifacts/single_sample/<ts>/` (the run used is `22042026_221229`) |
| 4 | `nichecompass_5regions_downstream.py` (`nichecompass_5regions_downstream.sh`) | nichecompass | downstream battery on step 3; adds the ROI flags (incl. `in_glom`) of the three original regions -> `<run>/xenium_with_niches.h5ad` |
| 5 | `annotate_in_glom_new_samples.py` (run inside the step 6/7 wrappers) | nichecompass | `in_glom` for the two Region 1 slides from the QuPath glomerulus polygons -> `$KB_ANALYSIS/in_glom_for_new_samples/in_glom_new_samples.parquet` |
| 6 | `nichecompass_overall_with_glom.py` (`nichecompass_overall_with_glom.sh`), then `overall_with_glom_finalize.py` + `glom_postprocess.py` + `copy_results_to_kramann.py` (`overall_glom_finalize_and_copy.sh`) | nichecompass | whole-core model with `in_glom` (training job 321201 hit the time limit after training; the finalize script redoes neighbours/UMAP/Leiden 0.4 and writes the adata) -> copied to `$KB_RESULTS/xenium/leon/new_data/overall_with_glom/`; since the 2026-07-03 reorganisation the cell-type `predic` folders are in `archived/` (`overall_with_glom_04`, `manual_glom_niche_based_on_04`) |
| 7 | `nichecompass_manual_glom_niche.py` (`nichecompass_manual_glom_niche.sh`; `copy_results_to_kramann.sh` copies the run) | nichecompass | glomerular model on `in_glom == True` cells (20,081 cells, 76 cores) |
| — | analyses/02_xenium_celltyping | — | refined cell typing `level2_refined` -> `$KB_ANALYSIS/celltype_concordance/outputs/xenium_labeled.h5ad` (it reads the step 6 object) |
| 8 | `nichecompass_replot_level2/replot_partA.py --version overall` and `replot_boundaries_level2.py --version overall` (`partE_pecrec_replot.sh`, job 327950; `partE_pecrec_boundaries.sh`, job 327951) | nichecompass / spatial_env | level2 re-plots of the whole-core niches -> `$KB_RESULTS/xenium/leon/new_data/overall_with_glom_04_level2/` |
| 9 | `nichecompass_replot_level2/bake_level2_into_niche_objects.py` (no wrapper) | nichecompass | writes `level2_refined`, `subcluster_level2` and the completed `in_glom` into the niche objects (backup `*.PRE_LEVEL2_20260709.h5ad`) |
| 10 | `nichecompass_5regions_downstream.py --leiden-res 1.0 --cell-type-key level2_refined` on the step 7 model (`nichecompass_replot_level2/partH_manualglom10_downstream.sh`, job 327992) | nichecompass | Leiden 1.0 re-clustering of the glomerular latent (no retraining) + downstream battery -> `$KB_RESULTS/xenium/leon/new_data/manual_glom_niche_based_on_04_res10_level2/` |
| 11 | `nichecompass_replot_level2/replot_boundaries_level2.py --version manual_glom_10 --full-context` (`partI_manualglom10_boundaries.sh`, job 327993) | spatial_env | per-core cell-boundary plots of the glomerular niches |
| 12 | `nichecompass_replot_level2/regenerate_hallmark_dedup.py` and `regenerate_metadata_figures.py`, then `sync_staged_figures.py` (`partJ_egfr_hallmark_fix.sh`, job 328172; `partK_joinfix_regen.sh`, job 328177) | nichecompass | Hallmark dedup heatmaps (log-normalised ULM) and the clinical-metadata figures after the identifier-join fix, staged in `$KB_ANALYSIS/figures_regen_20260722/` and copied into the results folders |

`level2_common.py` (merge of `level2_refined` + fixed palette) and `metadata_common.py` (per-core clinical table ->
obs) are helper modules for the level2 scripts.

## Inputs
| input | producer |
|---|---|
| `$KB_XENIUM_RAW/h5ad_files_only_transcript/output-*.h5ad` (per-region Xenium h5ads) and `marker_genes/clustered_annotated_final_file.h5ad` (3-region annotated reference, `predic`) | colleague (external) |
| `$KB_ANALYSIS/coordinates *.csv` (core polygons of the Region 1 slides), `$KB_ANALYSIS/samples_to_use.csv` (core whitelist) | manual (external) |
| `$KB_ANALYSIS/regions_of_interest/*_roi.csv` (ROI and `in_glom` flags of the three original regions) | colleague notebook `only_glom_roi/new_way_of_tying.ipynb` (external) |
| `$KB_ANALYSIS/in_glom_for_new_samples/*.csv` (QuPath glomerulus polygons of the two Region 1 slides) | manual annotation (external) |
| `$KB_PROJECT/niche_compass_results/data/gene_programs/` (OmniPath / NicheNet / MEBOCOST gene programs) | downloaded by NicheCompass on first use (external) |
| `$KB_PROJECT/niche_compass_results/Supplemental.csv` | legacy clinical table, used only as a fallback (external) |
| `$KB_ANALYSIS/celltype_concordance/outputs/xenium_labeled.h5ad` (`level2_refined`, `in_glom`) | analyses/02_xenium_celltyping |
| `$KB_ANALYSIS/celltype_concordance/outputs/per_sample_metadata_updated.csv` (per-core clinical table, identifier join) | analyses/02_xenium_celltyping/15_update_metadata_slope.py |
| `$KB_XENIUM_RAW/output-<region>/{cell,nucleus}_boundaries.parquet` | Xenium Onboard Analysis (external) |
| MSigDB Hallmark / Reactome gene sets | fetched at run time from OmniPath by decoupler (`dc.op.hallmark`) (external, not versioned) |
| private core lists (`KB_EXCLUDED_CORES`, read through `config/cores.py`) | private, not in the repository |

## Outputs (used downstream)
- `$KB_RESULTS/xenium/leon/new_data/manual_glom_niche_based_on_04_res10_level2/xenium_with_niches.h5ad`
  (20,081 glomerular cells; `latent_leiden_1.0`, `level2_refined`, clinical columns), `figures/`,
  `niche_marker_results/` (within-niche Wilcoxon markers per cell type), `cell_boundary_plots_per_core/`.
  Used by analyses/08_glomerular_niche20_crescent (niche 20).
- `$KB_RESULTS/xenium/leon/new_data/overall_with_glom_04_level2/xenium_with_niches_overall_with_glom.h5ad`
  (258,375 cells, 119 cores; `latent_leiden_0.4`, baked `level2_refined`), `cell_level2_and_niche.csv`,
  `figures/`, `niche_marker_results/`, `cell_boundary_plots_per_core/`. The path
  `archived/overall_with_glom_04/xenium_with_niches_overall_with_glom.h5ad` is a symlink to this object.
  The whole-core niches reach analyses/03 (Fig. 3, ED Fig. 4c,f) through `xenium_labeled.h5ad`.

## Figure panels (entry scripts in figures/)
| panel | entry script | code here |
|---|---|---|
| Fig. 6a | `figures/Fig6/Fig6a_niche_vs_celltype_heatmap.py` | `nichecompass_5regions_downstream.py::heatmap_niche_vs_celltype` (step 10) |
| Fig. 6b | `figures/Fig6/Fig6b_niche_maps_K50_K18.py` | `replot_boundaries_level2.py` (step 11) |
| Fig. 6c | `figures/Fig6/Fig6c_niche20_marker_dotplot.py` | `nichecompass_5regions_downstream.py::per_niche_celltype_markers` (step 10) |
| ED Fig. 4d-e | `figures/ED4/ED4d-e_celltype_niche_maps.py` | `replot_boundaries_level2.py --version overall` (step 8) |
| ED Fig. 4g | `figures/ED4/ED4g_niche_vs_celltype_heatmap.py` | `replot_partA.py::heatmap_noannot` (step 8) |
| ED Fig. 4h | `figures/ED4/ED4h_egfr_group_enrichment.py` | `regenerate_metadata_figures.py` -> `run_categorical_enrichment` (step 12) |
| ED Fig. 4i | `figures/ED4/ED4i_condition_stacked_bars.py` | `regenerate_metadata_figures.py` -> `plot_stacked_proportions` (step 12) |
| ED Fig. 4j | `figures/ED4/ED4j_hallmark_niches.py` | `regenerate_hallmark_dedup.py` -> `glom_postprocess.hallmark_pathway_dotplot` (step 12) |
| Supp. Fig. 19a | `figures/S19/S19a_hallmark_glom_niches.py` | same, config `manual_glom10` (step 12) |

The figure scripts read only the saved objects and tables above. They never retrain or re-cluster. Changes made
for the release (all backwards compatible, defaults = original behaviour) are listed in
`build/patch_rules/fig_niches.py`. Among them: per-panel output names, and selecting single cores by their public
id instead of the internal core label.

## Known caveats
- **Niche numbers are run-specific.** Leiden numbers clusters by size, and re-running the downstream step can
  renumber the niches. The crescent niche is `latent_leiden_1.0 == 20` *in the saved object*. The per-cell labels
  in the result h5ads are the reference, so do not regenerate them to reproduce figures. For the same reason
  `partK_joinfix_regen.sh` regenerated only the metadata figures.
- **Niche 20 is essentially one biopsy.** 203 of its 241 cells come from one core of K50 (IgAN). The rest come from
  K117 (19), Apollo t0 E1 (11), K110 (5), L14 (2) and L32 (1). By cell type: PEC 134, Immune 56, FIB 28, EC 7,
  aPT 7, POD 5, MES 4 (checked on the saved object, 2026-11-12). Statements about niche 20 are therefore
  statements about one crescentic glomerulus region.
- **`in_glom` completion was not scripted.** The `in_glom` flag in `xenium_labeled.h5ad` was completed outside
  any saved script (15,618 -> 20,081 cells according to the build trace; 20,081 `in_glom` cells in the current
  `xenium_labeled.h5ad`, matching the glomerular model). `bake_level2_into_niche_objects.py`
  copies it into the niche objects. The glomerular model itself was trained on the `in_glom` of step 4 plus the
  step 5 parquet.
- **Core exclusions.** The internal core labels contain pathology numbers. The NicheCompass exclusion lists
  (`TO_EXCLUDE`: QC-failed cores, unassigned cells, one Apollo control) are therefore read from a private list
  through `config/cores.py`. Without that file only `Unassigned` and the Apollo d0 A2 control are excluded. The
  other excluded cores are not part of the de-identified Zenodo object.
- **Two core labels differ from the private key by whitespace.** In both niche objects, the `TMA_ID` strings of
  XEN092 (K12) and XEN080 (K51) differ from the key labels in whitespace only. `cores.core_label()` matches exactly,
  so it cannot select these two cores without whitespace normalisation.
- **Clinical metadata.** Before 2026-07-22 the per-core clinical join was a token matcher. It assigned the wrong
  eGFR to 52 of 119 cores. Missing eGFR values also fell into the `>90` bin. Both bugs were fixed on 2026-07-22
  (`metadata_common.py`, `derive_ifta_egfr_groups`), and the figures were regenerated by job 328177. Any eGFR or
  IFTA figure older than that is stale. `attach_supplemental()` is kept only as a loudly-warning fallback.
- **Enrichment tests count cells.** The condition / eGFR / IFTA enrichments are one-sided Fisher tests on cells.
  Cells of a core are not independent, so the FDRs (ED Fig. 4h) are anti-conservative.
- **Hallmark sets are fetched at run time** from OmniPath (`dc.op.hallmark`) and are not versioned. A re-run on
  2026-11-12 got the same 7,318 pairs / 50 pathways and pixel-identical figures.
- **Staging.** The Kramann_lab quota was full in July 2026. The regeneration scripts therefore wrote to
  `$KB_ANALYSIS/figures_regen_20260722/` and `sync_staged_figures.py` copied the files with integrity checks.
- **Labels edited only in Illustrator:** `dPT` was renamed `aPT2`. The ED Fig. 4j colour bar reads "Z-scaled"
  but the values are min-max scaled. In Supp. Fig. 19a, niche 16 is labelled "10".

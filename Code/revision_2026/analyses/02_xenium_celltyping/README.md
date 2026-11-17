# 02 — Xenium cell typing: label transfer from snRNA-seq, level2 refinement, final Xenium object

## Purpose
Annotates the 258,375 Xenium cells (5 runs, 119 cores, 468-gene panel) with the 16 snRNA-seq level2 cell types and
writes the object every other Xenium analysis reads, `outputs/xenium_labeled.h5ad`: scVI/scANVI label transfer from the
snRNA-seq reference (scArches query surgery), a Leiden sweep, per-lineage subclustering with a consensus of the
transferred labels (`level2_refined`, all 16 types), per-core clinical metadata, recovery of glomerular PECs that the
transfer had called aPT, and completion of the manual-glomerulus flag `in_glom`. Also draws Supplementary Fig. 9a-d and
adds clinical metadata to the snRNA-seq reference. Human GRCh38, HGNC symbols, seed 42. Methods text:
`$KB_ANALYSIS/celltype_concordance/METHODS_celltyping.md`.

## Analysis workspace
`$KB_ANALYSIS/celltype_concordance/` (scripts at its root; `outputs/`, `outputs/metrics/`, `outputs/models/`,
`figures/`). Paths and constants: `config.py`; concordance metrics: `metrics.py`.

## Run order (source `config/env.sh` first; conda env kpmp_repurposing unless stated)
| step | script | wrapper | what it does / writes |
|---|---|---|---|
| 0 | `00_check_env.py` | none | environment + input check; appends versions to the workspace `RUN_LOG.md` |
| 1 | `01_prep_reference.py` | none (login node) | snRNA reference subset to the 468 panel genes, ambient-corrected counts rounded -> `outputs/humanized_ref.h5ad`, `metrics/rounding_report.csv` |
| 2 | `02_train_transfer.py` | `02_train_transfer.sh` (32 CPUs, 200 GB, no GPU) | scVI base + scANVI level1/level2 + query surgery -> `outputs/xenium_labeled.h5ad` with `pred_level1/2`, `conf_level1/2`; models in `outputs/models/` |
| 3 | `03_leiden_sweep.py` | `03_leiden_sweep.sh` | neighbours on `X_pca_harmony`, Leiden 0.2-2.0 -> `leiden_harmony_*` (in place) |
| 4 | `09_level2_subcluster.py` | none | per-lineage re-embedding + Leiden, consensus `pred_level2` per subcluster, hierarchical aPT/dPT step -> `level2_refined`, `subcluster_level2` (in place); `metrics/level2_subcluster_assignments.csv`, `level2_refined_distribution.csv` |
| 5 | `10_heatmap_pred_vs_refined.py`, `12_umap_level2_refined.py`, `13_dotplot_bubble.py` | none | Supp. Fig. 9a, 9b, 9c-d (drawn 2026-07-01, i.e. before steps 6-9) |
| 6 | `15_update_metadata_slope.py` | none | per-core clinical table -> `outputs/per_sample_metadata_updated.csv` (re-run 2026-07-22 with a join fix; the earlier table is kept as `per_sample_metadata_updated.PRE_JOINFIX_20260722.csv`) |
| 7 | `16_apply_metadata_to_h5ad.py` | none (back up the object first) | writes the per-core table into the object's obs (`supp_*`, `EGFR`, `egfr_group`, ...; in place, gzip). Imports `metadata_common.py` and `nichecompass_5regions_downstream.py` from `analyses/01_xenium_preprocessing_niches/` |
| 8 | `17_recover_pec_from_apt.py` | none (back up the object first) | re-embeds the aPT/dPT/PEC compartment, relabels 182 aPT cells PEC (subcluster with in_glom >= 0.25, CFH >= 2.0, prior PEC >= 0.5) -> `level2_refined`, `subcluster_level2` (in place); `metrics/pec_recovery_subcluster_decision.csv` |
| 9 | `17b_complete_in_glom.py` | none (back up the object first) | **reconstruction** (see below): merges the Region_1 glomerulus flags -> `in_glom` (in place; `--out`, `--input`, `--dry-run`) |
| 10 | `18_add_clinical_to_snRNA_ref.py` | none | donor-level clinical fields onto the snRNA reference -> `$KB_DATA/charlotte_new_data/adata_20260717_with_clinical.h5ad` (+ per-donor CSV, note) |

Backups kept in `outputs/`: `xenium_labeled.PRE_METADATA_UPDATE_20260702.h5ad` (the object after step 4, 2026-07-01, as
used for Supp. Fig. 9) and `xenium_labeled.PRE_PEC_RECOVERY_20260709.h5ad` (after step 7). Steps 6-9 were run on
2026-07-02 (6-7) and 2026-07-09 (8-9); the final object is dated 2026-07-09 15:58. Not shipped: steps 04-08, 11, 14 and 19
of the workspace (concordance metrics, report, label exports, a later lookup), which no panel uses.

### Step 17b (reconstruction)
The `in_glom` completion of 2026-07-09 was done by hand; `17b_complete_in_glom.py` re-implements it: the flags of the two
Region_1 samples from `in_glom_new_samples.parquet` are merged on (sample, cell_id); all other cells keep their value.
Check (2026-11-12): applied to the in_glom of the archived pre-step-8 object (step 8 does not change in_glom; the object
between steps 8 and 9 was not kept), it equals the final object's `in_glom` for all 258,375 cells (0 disagreements;
True 15,618 -> 20,081 = +2,587 and +1,876 in the two Region_1 samples). The final object differs from the pre-step-8
copy only in `in_glom` (4,463 cells), `level2_refined` and `subcluster_level2` (182 cells each).

## Inputs
| input | producer |
|---|---|
| `$KB_DATA/charlotte_new_data/adata_20260717.h5ad` (snRNA-seq reference, 363,656 nuclei, level1/level2) | upstream snRNA-seq processing; not in this repository |
| Xenium NicheCompass object `config.XENIUM` = `$KB_RESULTS/xenium/leon/new_data/overall_with_glom/xenium_with_niches_overall_with_glom.h5ad` | `analyses/01_xenium_preprocessing_niches` (`nichecompass_overall_with_glom.py` -> `nichecompass_5regions_downstream.py` -> `copy_results_to_kramann.py`). **That path no longer exists**: the object as read by step 2 is kept as `$KB_RESULTS/xenium/leon/archived/overall_with_glom_04/xenium_with_niches_overall_with_glom.PRE_LEVEL2_20260709.h5ad`; the current file in `new_data/overall_with_glom_04_level2/` has `level2_refined` etc. baked in |
| `$KB_ANALYSIS/Supplemental Data 1_metadata_update with slope.csv` (step 6) | external: clinical supplement (slope update) |
| `$KB_ANALYSIS/in_glom_for_new_samples/in_glom_new_samples.parquet` (step 9) | `analyses/01_xenium_preprocessing_niches/annotate_in_glom_new_samples.py` |
| `$KB_DATA/scVI_kidney_biopsy_with_clinical.h5ad` (step 10) | upstream snRNA-seq object with clinical metadata; not in this repository |

## Outputs and the panels that use them
| output | used by |
|---|---|
| `outputs/xenium_labeled.h5ad` (final) | folders 03 (Fig. 3c, 3f-i; ED 4c, 4f; ED 6a; Supp. Fig. 12), 05 (Supp. Fig. 13, Xenium column), 06 (Supp. Fig. 10d), 01 (`bake_level2_into_niche_objects.py`) and 09 |
| `figures/heatmap_pred_level2_vs_refined.pdf` (step 5, `10`) | Supp. Fig. 9a (`figures/S9/S9a_heatmap_pred_vs_refined.py`) |
| `figures/umap_level2_refined_ondata.pdf` (`12`) | Supp. Fig. 9b (`figures/S9/S9b_umap_level2_refined.py`) |
| `figures/dotplot_top3_markers_level2_refined.pdf`, `bubble_level2_by_condition.pdf` (`13`) | Supp. Fig. 9c-d (`figures/S9/S9c-d_dotplot_bubble.py`) |
| `outputs/per_sample_metadata_updated.csv` | folder 06 (step 01), folder 01 (`nichecompass_replot_level2/metadata_common.py`) |
| `adata_20260717_with_clinical.h5ad` (step 10) | folders 03 (steps 1, m1), 04 and 05 |

## Caveats (read before reusing)
- **Supp. Fig. 9a-d show the labels before the PEC recovery.** They were drawn on 2026-07-01 (aPT 15,072, PEC 2,395);
  step 8 later moved 182 cells from aPT to PEC (14,890 / 2,577), which every other panel uses. Run on the final object
  the S9 entry scripts differ from the published panels (heatmap: the transferred-aPT row moves 0.0095 and the
  transferred-PEC row 0.018 from aPT to PEC; bubble plot: aPT and PEC counts); with
  `--xenium outputs/xenium_labeled.PRE_METADATA_UPDATE_20260702.h5ad` they redraw the published panels pixel for pixel.
- **Order of steps 8 and 9.** Step 8 used `in_glom` (one of its three criteria) before step 9 completed it, i.e. with all
  60,141 Region_1 cells False. Of the 182 recovered PECs, 46 lie in the Region_1 samples; 103 were flagged glomerular
  when step 8 ran and 135 are after step 9. The workspace PLAN.md (2026-07-09) records the recovery decision as
  unchanged with the completed flag; that re-evaluation was not repeated here.
- `in_glom` is the manual QuPath glomerulus annotation (point in polygon), not a NicheCompass niche.
- The per-core clinical columns in the object are those written by step 7 on 2026-07-02; step 6 was re-run on
  2026-07-22 (join fix) and step 7 was not re-run, so the object's `supp_*` columns can differ from
  `per_sample_metadata_updated.csv`. The object's `EGFR` column is corrupted (group values broadcast across donors);
  the Xenium figures use the per-core eGFR of `analyses/03_.../10_build_authoritative_egfr.py` instead.
- `dPT` (Supp. Fig. 9b, file names, code) is labelled `aPT2` in most published panels (Illustrator edit).
- Supp. Fig. 9b is the UMAP stored in the object (`X_umap`, NicheCompass latent space), not a UMAP of the scVI latent.
- Supp. Fig. 9c colour is the log fold change of the scanpy Wilcoxon ranking, scaled per gene, not scaled mean expression.
- Step 2 ran on CPU (no GPU) with `scvi.settings.seed = 42`; whether a re-run reproduces the models bit for bit on
  other hardware or thread counts was not tested.

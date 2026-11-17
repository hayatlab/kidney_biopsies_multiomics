# 06 — Xenium quality per core vs biopsy age

## Purpose
Reviewer question: does the age of the archived biopsy (years between biopsy and the Xenium run on 2025-04-16) degrade
Xenium data quality? Per TMA core (119 cores) the script computes cell counts and density, transcripts and genes per
cell, transcripts per 100 um2 of cell area, complexity (genes at 100 transcripts, by rarefaction), negative-probe and
unassigned-codeword rates, segmentation and label-transfer confidence and the mean Factor 1 score, and relates each to
biopsy age: Spearman rho with a bootstrap 95% CI (2,000 resamples, seed 42), OLS on log-scale metrics (effect in % per
10 years) unadjusted and adjusted for condition or TMA block (HC3 SE), BH across metrics, a Huber robust-regression
sensitivity check, and age-bin medians. The main analysis uses the disease cores only (n = 86); control cores are shown.

## Analysis workspace
`$KB_ANALYSIS/xenium_biopsy_age_qc/` (`scripts/`, `outputs/`, `outputs/private/`, `figures/`, `report/`).

## Run order (source `config/env.sh` first; conda env agent_env; login node, minutes)
| step | script | wrapper | output |
|---|---|---|---|
| 1 | `01_core_qc_metrics.py` | none | `outputs/core_qc_metrics.csv` (one row per core, public `core_id`), `outputs/private/core_id_map.csv` (core_id -> internal label; private), `outputs/01_meta.json` |
| 2 | `02_association.py` | none | `outputs/associations_biopsy_age.csv` (23 metrics x subsets all / disease), `age_bin_medians.csv`, `02_meta.json`; figures `F1`-`F4` |
| 3 | `04_simple_plots.py` | none | `figures/S1_simple_scatter{,_nocolour}.{pdf,png}`, `figures/S2_simple_boxplots.{pdf,png}` (numbers read from step 2, not recomputed) |

Not shipped: `03_build_report.py` (LaTeX report of the workspace). Last run: 2026-10-16 (steps 1 -> 3 -> 2; step 3
reads only the disease rows of step 2, which did not change in that re-run).

## Inputs
| input | producer |
|---|---|
| `$KB_ANALYSIS/celltype_concordance/outputs/xenium_labeled.h5ad` (obs incl. negative-probe / unassigned-codeword counts, `layers["counts"]`, `obsm["spatial"]`) | `analyses/02_xenium_celltyping` |
| `$KB_ANALYSIS/celltype_concordance/outputs/per_sample_metadata_updated.csv` (condition / disease per core) | `analyses/02_xenium_celltyping/15_update_metadata_slope.py` |
| `$KB_ANALYSIS/xenium_biopsy_age_qc/Sample_clinical_Data_Metadata.xlsx`, sheet "Samples" (biopsy dates, clinical data) | external: clinical metadata workbook (2026-10-15) |
| `$KB_ANALYSIS/vivien_takeover/outputs/xenium_factor_scores.parquet` | `analyses/03_xenium_factor1_neighbourhoods_cci/01_score_factors_xenium.py` |

## Outputs and the panels that use them
| output | panel (entry script) |
|---|---|
| `figures/S1_simple_scatter_nocolour.pdf` (all six subpanels) from `outputs/core_qc_metrics.csv` + `associations_biopsy_age.csv` | Supp. Fig. 10d (`figures/S10/S10d_qc_metrics_vs_biopsy_age.py`) |

## Caveats (read before reusing)
- **Biopsy year of the control cores is an assumption.** `01_core_qc_metrics.py` sets `CTRL_YEAR_OVERRIDE = 2000`
  (collaborator statement, 2026-10-16: all Apollo pre-transplant control cores were biopsied in 2000; set to 1 July
  2000), overriding the workbook dates (2021-2024); `biopsy_date_overridden` flags the 33 cores. The disease-core
  analysis (n = 86) does not depend on it; the all-core estimates do (their unadjusted correlations change sign).
- **The published Supp. Fig. 10d predates that override.** In the published panel the open circles (controls) sit at
  biopsy ages of 0-5 years (workbook dates); the file on disk and the shipped code place all 33 at 24.8 years. The
  disease cores, the fitted lines and bands and every printed statistic are identical in both versions.
- All cores were measured in one Xenium run, so biopsy age is not confounded with run date. The control cores form
  their own TMA block, so for them condition and block cannot be separated.
- `outputs/private/core_id_map.csv` contains internal core labels with pathology numbers; it is git-ignored in the
  workspace and must not be shared.

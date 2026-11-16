# 05 — NELL1 (and APP) expression per cell type and association with eGFR (snRNA-seq; Xenium for APP)

## Purpose
Answers "how is NELL1 expressed across kidney cell types, and does it track kidney function?": per cell type the mean
log1p(CP10k) and the % of cells expressing NELL1 and APP in snRNA-seq (16 level2 types, 363,656 nuclei) and APP in
Xenium (NELL1 is not on the 468-gene panel; checked against the panel, not assumed), and per cell type the association
of donor-level pseudobulk expression with eGFR at biopsy. The unit is the donor, never the cell: per donor x cell type,
pseudobulk = log1p(sum of gene counts / sum of library sizes x 1e4), kept when the donor has >= 20 cells of that type;
Spearman rho with 95% CI (Fisher z, Bonett-Wright SE), BH-FDR across the cell types within gene x modality x unit, OLS
slope per 10 ml/min/1.73 m2 with Breusch-Pagan and Shapiro-Wilk checks, and a partial Spearman adjusted for disease.
Human GRCh38, HGNC symbols, seed 42.

## Analysis workspace
`$KB_ANALYSIS/nell1_app_egfr/` (`scripts/`, `tables/`, `figures/`, `logs/`).

## Run order (source `config/env.sh` first)
| step | script | conda env | wrapper | output |
|---|---|---|---|---|
| 1 | `01_nell1_app_egfr.py` | nichecompass | `01_nell1_app_egfr.sh` (8 CPUs, 250 GB) | `tables/expression_by_celltype.csv`, `pseudobulk_per_unit.csv`, `egfr_association.csv`, `01_meta.json`; `figures/nell1_app_expression_by_celltype.pdf`, `nell1_app_egfr_forest.pdf`, `{NELL1,APP}_snRNA_donor_egfr_scatter.pdf`, `APP_Xenium_donor_egfr_scatter.pdf` |

Published run: SLURM job 329592, 2026-09-10 (`--smoke` = stratified subsample of 40 cells per sample x cell type, for
testing only). Not shipped: `02_specificity_context.py` and `03_celltype_dotplot.py` of the workspace (no panel uses
them).

## Inputs
| input | producer |
|---|---|
| `$KB_DATA/charlotte_new_data/adata_20260717_with_clinical.h5ad` (snRNA-seq, ambient-corrected non-integer counts in `.X`, `level2`) | `analyses/02_xenium_celltyping/18_add_clinical_to_snRNA_ref.py` |
| `$KB_ANALYSIS/vivien_takeover/outputs/metadata/sample_metadata_corrected.csv` (eGFR at biopsy, donor = `patient_ID`, disease per sample) | `analyses/03_xenium_factor1_neighbourhoods_cci/24_apply_corrected_metadata.py` |
| `$KB_ANALYSIS/celltype_concordance/outputs/xenium_labeled.h5ad` (`layers["counts"]`, `level2_refined`) | `analyses/02_xenium_celltyping` |
| per-core eGFR via `05_xenium_figures.py::corrected_core_egfr` (imported from the copy in `analyses/03_.../`) -> `authoritative_core_egfr.csv` | `analyses/03_xenium_factor1_neighbourhoods_cci/10_build_authoritative_egfr.py` |

## Outputs and the panels that use them
| output | panel (entry script) |
|---|---|
| `figures/nell1_app_expression_by_celltype.pdf`, left column "snRNA NELL1" (`::dotplot`) | Supp. Fig. 8a (`figures/S8/S8a_nell1_dotplot.py`) |
| `figures/NELL1_snRNA_donor_egfr_scatter.pdf`, subpanels PT, POD, TAL, DCT (`::scatter_grid`) | Supp. Fig. 8b (`figures/S8/S8b_nell1_vs_egfr_per_donor.py`) |
| `tables/expression_by_celltype.csv`, `pseudobulk_per_unit.csv`, `egfr_association.csv`, `01_meta.json` | read by the two entry scripts (they recompute the light aggregation from the objects only if the tables are missing, or with `--from-objects`) |

## Caveats (read before reusing)
- **The tables contain patient identifiers.** `pseudobulk_per_unit.csv` lists every unit by its id: the snRNA donor
  (`patient_ID`) and the Xenium donor key, which is taken from the pathology number in the internal core label or from
  `H.number`. Keep `tables/` private; the figures carry no identifiers.
- snRNA donors: 61 (the donors of K94/K97 and K95/K98 give 2 samples each, pooled). The Xenium donor key is an
  approximation for the 33 control cores, which carry neither a T- nor an H-number (grouped by condition and eGFR).
- The q values in Supp. Fig. 8b are BH over the 17 rows (16 cell types + all cells) of snRNA / donor / NELL1; n per
  subpanel differs because a donor enters a cell type only with >= 20 nuclei of it (POD: 54 donors).
- eGFR differs systematically between diseases; the partial Spearman adjusted for disease is in `egfr_association.csv`
  (and the forest plot), not in the published panel.
- `dPT` is labelled `aPT2` in the published panel (Illustrator edit).

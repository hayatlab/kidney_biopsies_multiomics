# 14 — Xenium quality metrics per cell, by core and run

## Purpose
Box plots of three per-cell quality metrics for every analysed core of the five Xenium runs (Supplementary Fig. 11a–c:
a total transcripts, b unique transcripts, c cell area). Analysis by Sidrah Maryam (notebook `Supp_fig_4a-c.ipynb`);
`qc_per_core.py` applies the definitions of that notebook to the analysed object. The association of per-core quality
with biopsy age (Supplementary Fig. 11d) is in `analyses/06_xenium_biopsy_age_qc`.

Method: `scanpy.pp.calculate_qc_metrics` on the raw counts: total transcripts = sum over the 468 panel genes, unique
transcripts = genes with at least one count; cell area = `cell_area` of the Xenium cell table (µm²). `seaborn.boxplot`
per core: box = quartiles, whiskers = most extreme cell within 1.5 × IQR, circles = cells beyond. The cell is the unit;
no test is shown.

## Analysis workspace
None: the figure script reads the Xenium object and writes only the panels.

## Run (source `config/env.sh` first)
| script | conda env | output |
|---|---|---|
| `python "figures/S11/S11a-c_xenium_qc_boxplots_per_core.py"` (imports `qc_per_core.py`; login node, about 3 min) | nichecompass | `$KB_FIG_OUT/S11/S11a_total_transcripts_per_core.pdf`, `S11b_unique_transcripts_per_core.pdf`, `S11c_cell_area_per_core.pdf`, `S11a-c_box_statistics.csv` |

## Input
| input | producer |
|---|---|
| `$KB_RESULTS/xenium/leon/new_data/overall_with_glom_04_level2/xenium_with_niches_overall_with_glom.h5ad` (258,375 cells × 468 genes; obs `sample` = run, `cell_area`, `condition`) | `analyses/01_xenium_preprocessing_niches` |

## Result (20 Nov 2026, re-run of the figure script)
258,375 cells in 119 cores: 15, 33, 16, 28 and 27 cores in the five runs (rows of the figure, top to bottom). The
quartiles, medians and maxima of all 357 boxes equal those of the printed panels.

## Caveats (read before reusing)
- **Box order and tick labels differ from the figure:** this script orders the boxes of a run by public core id and
  takes label and box from the same table row. To compare with the figure, match boxes by their statistics
  (`S11a-c_box_statistics.csv`), not by their position.
- The notebook as received reads three runs with spatialdata-io and writes one PDF per run with the internal core
  labels; the published five-run figure was not made by the notebook as received. The script here reproduces its boxes
  from the analysed object and names cores by public ids (sample id, condition, `xenium_core_id`).
- Total transcripts are recomputed from the gene counts; they differ from the `total_counts` column of the Xenium
  cell table by at most 3 per cell (control probes and codewords are not counted).
- Internal core labels contain pathology numbers: never print `TMA_ID` in a figure or table.

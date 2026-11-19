# 04 — snRNA-seq overview, QC and CLR cell-type abundance

## Purpose
- **Overview and QC figures (Supp. Fig. 2a–b).** `01b_umap_from_csv.py` joins the project's UMAP coordinates to the obs
  table of the snRNA-seq object (363,656 nuclei, 63 samples); `03_overview_figures_v2.py` draws UMAPs by cell type,
  disease, sex and sample, per-sample violins of transcripts, genes and mitochondrial fraction per nucleus, nuclei per
  sample and the level1 / level2 composition per sample, and writes the numbers behind every figure.
- **CLR cell-type abundance (Extended Data Fig. 2e–g).** `01_build_clr_inputs.R` turns the per-sample level2 (and level1)
  composition into centred log-ratios and writes the sample metadata in the shape the saezlab MCFA figure script expects;
  `02_ext2e_clr_compositions_boxplot.R` is that script (`ext2e_clr_compositions_boxplot.R` of saezlab/kidney_biopsies_MCFA)
  re-run on the new annotation, with five changes marked `## CHANGED` in its header, and a sixth of 19 Nov 2026 that
  reproduces the panels printed in the final figures (all 63 samples, eGFR group of Supplementary Data 1; see caveats).
- `labels.py` is the palette module of the project's data explorer (`kidney-biopsy-explorer/export/lib/labels.py`), which
  wrote the palette file the original 03 run read. It is kept for reference; the figure entry script uses the frozen copy
  `figures/_lib/snrna_palette_20260917.json` (see caveats).

## Analysis workspaces
- `$KB_ANALYSIS/snrna_overview_20260831` — 01b and 03 (`umap_and_obs.parquet`, `figures/{png,pdf,svg}/`,
  `figures/tables/`, `figures/run_meta_v2.json`).
- `$KB_ANALYSIS/Charlotte_help` — 01_build_clr_inputs.R and 02_ext2e (scripts were in `Charlotte_help/scripts/`;
  outputs in `outputs/`, figures in `figures/`).

## Run order (source `config/env.sh` first)
| step | script (wrapper) | conda env | run | what it does |
|---|---|---|---|---|
| 1 | `01b_umap_from_csv.py` (no wrapper) | pathfire (inferred: numpy 2.0.2 / pandas 2.2.3 / h5py 3.12.1 in `umap_run_meta.json` match that env) | 2 Sep 2026 | obs columns of the object + UMAP coordinates -> `umap_and_obs.parquet` |
| 2 | `03_overview_figures_v2.py` (`03_overview_figures_v2.sh`, 4 CPUs / 48G) | agent_env | job 329965, 17 Sep 2026 | the 11 overview figures (PNG 300 dpi, PDF and SVG with rasterised points) + `figures/tables/` |
| 3 | `01_build_clr_inputs.R` (`Rscript`, no wrapper) | scmega2 | 25 Sep 2026 | `outputs/clr~level{2,1}~new_annotation.csv`, `counts~level{2,1}~new_annotation.csv`, `patient_metadata_new_annotation.csv` |
| 4 | `02_ext2e_clr_compositions_boxplot.R [level2 \| level1]` (`Rscript`, no wrapper) | scmega2 | 25 Sep 2026 (level2 = published) | `figures/ext2e_clr_compositions_boxplot_<level>.{pdf,png}`, `outputs/ext2e_stats_<level>.csv` |

03's wrapper `cd`s into the workspace and writes its log to `logs/` there. Steps 3–4 need the tables of step 2.

## Inputs
| input | producer |
|---|---|
| `$KB_DATA/charlotte_new_data/adata_20260717_with_clinical.h5ad` (obs read with h5py: sample, patient_ID, sex_clinical, disease, Study, level1, level2, nCount_RNA, nFeature_RNA, pct_counts_mt) | `analyses/02_xenium_celltyping/18_add_clinical_to_snRNA_ref.py` (clinical columns added to the snRNA reference `adata_20260717.h5ad`, external) |
| `$KB_ANALYSIS/umap_embeddings.csv` (cellid, umap_1, umap_2; 363,763 rows incl. 107 nuclei not in the object) | the project's UMAP of the snRNA object, supplied by Leon on 2 Sep 2026 (producer not recorded; external) |
| `$KB_ANALYSIS/vivien_takeover/outputs/metadata/sample_metadata_corrected.csv` | `vivien_takeover/scripts/24_apply_corrected_metadata.py` (not shipped). It has a `patient_ID` column: 03 uses only `disease` and `sex_clinical`; 01_build_clr_inputs.R writes only sample, disease, classification, eGFR, IF/TA and proteinuria |
| `$KB_SUPP_DATA1` (default `$KB_ANALYSIS/Supplementary Data 1.csv`; `;`-separated; columns sample_id, disease, egfr_at_biopsy, egfr_group, IFTA_pct) | Supplementary Data 1 of the paper (clinical table, no patient identifiers); read by `02_ext2e_clr_compositions_boxplot.R` for the final ED 2e–g panels |
| colour palette (disease, level1, level2) | original run: `kidney-biopsy-explorer/public/data/meta/labels.json` (written by `labels.py`); release: `figures/_lib/snrna_palette_20260917.json` |

## Outputs and the panels that use them
| output | used by | entry script |
|---|---|---|
| `figures/pdf/umap_sex`, `umap_disease`, `umap_sample` | Supp. Fig. 2a | `figures/S2/S2a-b_snrna_overview_qc.py` |
| `figures/pdf/violin_nCount_RNA`, `violin_nFeature_RNA`, `violin_pct_counts_mt`, `bar_nuclei_per_sample`, `stacked_level2` | Supp. Fig. 2b | same |
| `figures/tables/composition_level{1,2}_per_sample.csv`, `nuclei_per_sample.csv` | `01_build_clr_inputs.R` | — |
| `Charlotte_help/outputs/clr~level2~new_annotation.csv`, `patient_metadata_new_annotation.csv` | `02_ext2e_clr_compositions_boxplot.R` | — |
| `Charlotte_help/figures/ext2e_clr_compositions_boxplot_level2.pdf`, `outputs/ext2e_stats_level2.csv` | Extended Data Fig. 2e–g | `figures/ED2/ED2e-g_clr_compositions_boxplot.R` |

`umap_and_obs.parquet` keeps the object's `patient_ID` column; no script here writes it to any output.
Verification (12 Nov 2026, `build/verification/atac_snrna.csv`): all 11 overview figures pixel-identical and all six
tables identical; ED 2e–g (run of 25 Sep) identical except for the jittered points (statistics table byte-identical).
ED 2e–g as printed in the final figures: statistics checked on 19 Nov 2026 (see caveats).

## Caveats (read before reusing)
- **Palette.** 03 read `labels.json` of the data explorer, which was rewritten on 2 Oct 2026, after the Supp. Fig. 2
  figures (17 Sep). The rewrite only dropped `niche_glom.crescent`; the disease / level1 / level2 colours are the same in
  both versions (explorer commits 3b3d469 of 14 Sep and e9d8e7a of 2 Oct). The frozen copy holds the 14 Sep entries, and 03
  run with it reproduces all 11 figures pixel-identically.
- **ED 2e–g jitter.** `geom_jitter` ran without a seed, so the horizontal positions of the points in the published panels
  cannot be reproduced; boxes, brackets and p-values are identical. The entry script seeds the jitter (42).
- **ED 2e–g statistics.** The brackets show unadjusted Welch t-test p-values (ggpubr); `ext2e_stats_<level>.csv` adds the
  mean difference with 95 % CI and BH within panel. Samples are treated as independent: K94 and K97 (both control) are two
  biopsies of one person, so the control group (n = 7) has 6 donors; the other donor with two samples (K95 / K98) is
  control_AKI.
- **ED 2e–g, final panels (figure PDFs of 18 Nov 2026).** The printed panels were redrawn after the run of 25 Sep; the
  code that drew them is not on this system. `02_ext2e_clr_compositions_boxplot.R` (CHANGED 6, 19 Nov 2026) reproduces
  all 30 printed p-values with two differences from the 25 Sep run: the two control_AKI samples are kept (own group
  without a test in e; both in the "< 30" eGFR group in f, where they are two biopsies of one donor; no IF/TA value, so g
  is unchanged), and the eGFR group is the column `egfr_group` of Supplementary Data 1 (35 of 63 eGFR values differ from
  the corrected sample metadata of September; 9 of the 61 non-AKI samples change group; K5 and K53 have eGFR 30.0 and
  group "<30" in the table). Group sizes in f: 25 / 14 / 9 / 15 (>= 90, 60-90, 30-60, < 30). Without the control_AKI
  samples the "< 30" comparisons would be p = 2.1e-08 (aPT), 0.0014 (PT), 0.00013 (Immune) instead of the printed
  5.7e-06, 0.00065, 0.00099. `EXT2E_VERSION=20260925` restores the samples and statistics of the earlier run.
- **Reference group.** Upstream used 11 external KPMP controls as the reference; they are not in this object, so the
  reference is the study's own 7 controls. IF/TA is stored as a percentage and divided by 100 for the script's 10 % / 50 %
  cut-offs. CLR without pseudocount: all 16 level2 types occur in all 63 samples.
- **Corrected metadata.** Disease and sex come from the corrected sample table; in the object 13,456 nuclei (K68, K95,
  K98) carried no disease / sex.
- **Labels.** The level2 population `dPT` is printed as `dPT` here; other figures call it aPT2.

# 08 — Glomerular niche 20 ("crescent" niche): PEC differential expression, MCFA factor 1 overlap, external validation

## Purpose
Glomerular niche 20 of the manual-glomerulus NicheCompass model (Leiden resolution 1.0, folder
`analyses/01_xenium_preprocessing_niches`) is the crescent-like niche (241 cells: 134 PEC, 56 immune, 28 FIB, ...). This
folder characterises its parietal epithelial cells (PECs): MAST differential expression of niche-20 PECs versus PECs of
all other glomerular niches, marker dotplot and volcano, the overlap of the resulting gene sets with the MCFA factor 1
loadings of the PT and PEC views, and a projection of the niche-20 signature onto an independent glomerular spatial
dataset (GSE294965, Sultana et al., PMID 41028563).

## Analysis workspace
The scripts ran from the root of `$KB_ANALYSIS` (= `review_analyses/`). Outputs:
- `$KB_ANALYSIS/celltype_concordance/outputs/crescent_check/de/` — `pec_manualglom10.h5ad`, `markers_manualglom10.csv`, the 6f dotplots
- `$KB_ANALYSIS/niche20_pec/` — `DE/niche20_vs_others_pecOnly/` (volcano, sorted DE copy), `GSEA/niche20_vs_others_pecOnly/`, `venn_mofa_factor1/`
- `$KB_ANALYSIS/validation_data/niche20_score_validation/option1/` — external validation

## Run order (source `config/env.sh` first; wrappers are SLURM scripts)
| step | script | conda env | wrapper | output |
|---|---|---|---|---|
| 0 | `00_build_pec_manualglom10.py --out <...>/pec_manualglom10.h5ad` | nichecompass | none | PEC-only h5ad (see below) |
| 1 | `Rscript de_crescent_pec_mast.R <pec_manualglom10.h5ad> <markers_manualglom10.csv>` | de_env | none saved | MAST table, all 468 panel genes |
| 2 | `plot_niche20_pec_top5_bottom5.py` | nichecompass | none | top-5/bottom-5 dotplots (v1 and final) |
| 3 | `niche20_volcano_gsea.R` | de_env | `niche20_volcano_gsea.sh` | volcano, GSEA (Hallmark/Reactome/GO BP), sorted DE copy |
| 4 | `venn_niche20_mofa_factor1.py` | decoupler | none | Venn diagrams, overlap statistics, membership tables |
| 5 | `export_venn_tables_xlsx.py` | spatial_env (openpyxl) | none | `niche20_vs_MOFA_Factor1_3way_top200.xlsx` |
| 6 | `niche20_validation_pseudotime.py` | spatial_env | `niche20_validation_pseudotime.sh` (300 GB, 8 CPUs) | polygon-level PCA / DPT / k-means / scores, figures, tables |

Step 0 is a hand-written reconstruction: the script that built `pec_manualglom10.h5ad` was not saved. Its output was
verified against the file used in the paper (same 1,711 cells in the same order, identical counts, genes, obs names, `grp`
and `donor` columns including category order), and step 1 run on the rebuilt file reproduced `markers_manualglom10.csv`
byte for byte (12 Nov 2026, verification job 340733). Step 1 was run without a third argument, so the donor column is
`donor`; with one donor level in the crescent group the script drops the donor term (log line
`latent.vars = NONE (crescent is single-donor)`), i.e. the test is `FindMarkers(test.use = "MAST", latent.vars = NULL,
logfc.threshold = 0, min.pct = 0)` on log-normalised counts (Seurat 5.2.0, MAST 1.32.0).

## Inputs
| input | producer |
|---|---|
| `$KB_RESULTS/xenium/leon/new_data/manual_glom_niche_based_on_04_res10_level2/xenium_with_niches.h5ad` | `analyses/01_xenium_preprocessing_niches`: `nichecompass_manual_glom_niche.py`, then `nichecompass_5regions_downstream.py` at Leiden 1.0 with `level2_refined` (`nichecompass_replot_level2/partH_manualglom10_downstream.sh`) |
| `pec_manualglom10.h5ad` | step 0 |
| `markers_manualglom10.csv` | step 1 |
| `niche20_pec/DE/niche20_vs_others_pecOnly/markers_MAST_niche20_vs_others_pecOnly.csv` | step 3 (sorted copy of `markers_manualglom10.csv`) |
| `$KB_ANALYSIS/loadings~hvg~5.csv` | MCFA (multicellular factor analysis) loadings of the snRNA-seq model; produced outside this folder (not traced here) |
| `$KB_ANALYSIS/validation_data/GSE294965_processed_data.h5ad` | external: processed data of GEO GSE294965 (Sultana et al.) |
| `$KB_ANALYSIS/niche9_pec/DE/niche9_vs_others_pecOnly/markers_MAST_niche9_vs_others_pecOnly.tsv` | earlier niche-9 DE (not shipped); read by step 6 only for the overlap table `signature_overlap_vs_niche9.tsv` |

## Outputs and the panels that use them
| output | panel (entry script) |
|---|---|
| `crescent_check/de/dotplot_niche20_vs_other_pec_top5_bottom5_ndrg1_varscale.pdf` (ad hoc variant of step 2) | Fig. 6f (`figures/Fig6/Fig6f_pec_niche20_vs_other_dotplot.py`) |
| `niche20_pec/DE/niche20_vs_others_pecOnly/volcano_niche20_vs_others_pecOnly.pdf` | Supp. Fig. 13b (`figures/S13/S13b_niche20_pec_volcano.R`) |
| `niche20_pec/venn_mofa_factor1/main_restricted_universe/top200/venn_*.pdf` | Supp. Fig. 16a-f (`figures/S16/S16a-f_venn_niche20_mofa_factor1.py`) |
| `niche20_vs_MOFA_Factor1_3way_top200.xlsx` | Supplementary Data 21 (per the figure trace) |
| `validation_data/niche20_score_validation/option1/*.pdf` | Supp. Fig. 17a-i (`figures/S17/S17a-i_niche20_validation_pseudotime.py`) |

The niche-20 PEC DE table itself (`markers_manualglom10.csv`) also feeds Fig. 6f, Supp. Figs 16-17, the
supplementary DE table, and the crescent ULM score of the snRNA-seq PEC subclusters (its positive genes at padj < 0.05;
`analyses/09_snrna_pec_immune_subclusters/pec/pec_leiden_sweep.py`, used in Fig. 6h).

## Caveats (read before reusing)
- **Niche 20 is dominated by one biopsy.** 203 of its 241 cells and 114 of its 134 PECs come from K50 (IgAN); the other
  crescent PECs come from K117 (DN, 16) and K110 (FSGS, 4). The signature is therefore largely that lesion's.
- **The MAST contrast is not adjusted for donor or sample.** `donor` is the `H.number` column, which is filled only for
  the two regions added during review; every PEC of the three original regions, including all 134 niche-20 PECs, has
  donor `'nan'`. The script therefore found one donor level in the crescent group and dropped the term. This happened
  because of missing metadata, not because the crescent group is single-donor (it spans three biopsies), so statements
  such as "adjusted for sample" or "crescent PECs derive from a single donor" do not hold.
- **Cell-level test.** MAST treats the 1,711 cells as independent observations; P values are therefore optimistic.
- **Niche labels.** A Leiden re-run can renumber niches; the per-cell labels stored in the input object
  (`latent_leiden_1.0`) are the reference, and `pec_manualglom10.h5ad` keeps only the derived `grp`.
- **Private file.** `pec_manualglom10.h5ad` carries pathology numbers in `obs['donor']`: never share it or print that
  column. The de-identified Zenodo object has no `H.number` (identifier columns removed); its donor-level column is
  `sample_id` (one core = one biopsy). Using it as `donor` would NOT reproduce the published table: the crescent group
  would have three levels and de_crescent_pec_mast.R would then include the term.
- **Fig. 6f** is an ad hoc combination of the two figures of step 2: the direction-consistent gene selection (NDRG1
  instead of AFDN) with per-gene min-max scaling (`standard_scale='var'`) and the v1 title. The entry script reproduces it.
- **Supp. Fig. 13b**: gene-label placement (ggrepel) is not seeded, so labels move between runs. The GSEA of step 3 runs on
  a 468-gene universe: Hallmark returned one term (EPITHELIAL_MESENCHYMAL_TRANSITION, NES 1.98, adjusted P 0.012),
  Reactome and GO BP none.
- **Supp. Fig. 16**: "top 200" means the 200 genes with the largest absolute factor 1 loading in a view, split by sign
  (PT 80 positive / 120 negative, PEC 69 / 131), then restricted to the Xenium panel and the view's genes; circles are not
  area-proportional.
- **Supp. Fig. 17** is a consistency check, not an independent validation of the niche-20 signature: on the
  480-gene GSE294965 panel it is ~92% identical to the earlier niche-9 signature (23 up / 22 down genes; see
  `signature_overlap_vs_niche9.tsv`). The 95% CI whiskers of panels h-i are unseeded bootstraps.

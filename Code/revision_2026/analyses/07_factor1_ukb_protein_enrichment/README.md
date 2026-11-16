# 07 — Factor 1 tubular signature in UK Biobank plasma proteins: enrichment (Fig. 4e table)

## Purpose
Reviewer response (8–9 Oct 2026): which pathways and protein classes characterise the circulating proteins that carry
the Factor 1 tubular-mean signature weight and predict incident CKD in UK Biobank (Olink Explore 3072, adjusted Cox
models computed by the UK Biobank colleague).

- `01_build_protein_table.R` merges the per-protein Cox results (2,923 proteins) with the Factor 1 tubular-mean signature
  weights (984 proteins carry a weight; matched by case-insensitive HGNC symbol), defines the quadrants, and adds Human
  Protein Atlas classes and gene-family flags.
- `04_composite_and_lowerleft.R` adds the composite rank (loading rank + Cox-z rank within the 504 proteins with weight
  > 0 and z > 0), the lower-left quadrant sets with over-representation analyses, and **fgsea on two rankings for five
  MSigDB collections**: Cox z of all 2,923 proteins, and the Factor 1 weight of the 984 matched proteins. The rows with
  ranking "Factor 1 weight, 984 matched proteins" and collections Hallmark and Reactome are the data of Fig. 4e.
- `05_figures_addendum.R` draws Leon's report figures (fig5–fig7; fig7 = fgsea tails per collection). They are not the
  published Fig. 4e, which has a different layout (see below).

## Analysis workspace
`$KB_ANALYSIS/factor1_reviewer_comment` (`scripts/`, `outputs/`, `figures/`, `report/`, `resources/`). Its `raw/` folder
holds the files received from the UK Biobank colleague; it is git-ignored and is never copied or shared.

## Run order (source `config/env.sh` first; no wrappers, each script ran as `Rscript scripts/<name>` on the login node)
| step | script | conda env | run | outputs |
|---|---|---|---|---|
| 1 | `01_build_protein_table.R` | de_env | 8 Oct 2026 | `outputs/protein_table.csv` |
| 2 | `04_composite_and_lowerleft.R` | de_env | 8 Oct 2026 | `outputs/fgsea_all_collections.csv` (Fig. 4e), `ora_lowerleft.csv`, `class_tests_lowerleft.csv`, `top40_composite.csv`, `lowerleft_proteins.csv`; `protein_table.csv` extended in place |
| 3 | `05_figures_addendum.R` | de_env | 8–9 Oct 2026 | `figures/fig5_top_composite.*`, `fig6_lowerleft.*`, `fig7_fgsea_{hallmark,reactome}.*` |

de_env: R 4.4.1, data.table 1.16.4, msigdbr 7.5.1, clusterProfiler 4.14.4, fgsea 1.32.2, ggplot2 3.5.1, ggrepel 0.9.6,
patchwork 1.3.0; conda-list md5 `cb0f13172b4a` (`outputs/versions.txt`, `outputs/env_hash.txt`). Seed 42. Human, HGNC
symbols. `02_enrichment.R`, `03_figures.R` and the Ihara comparison (06–09) of the same workspace feed no published
panel and are not shipped.

## Inputs
| input | producer |
|---|---|
| `raw/F1_Tubular_mean_results/Factor1/final_results_df.rds` (per-protein HR, 95 % CI, p and z from adjusted Cox models for incident CKD) | UK Biobank colleague (external; UK Biobank data, access by application only) |
| `raw/F1_Tubular_mean_results/Factor1/signature_bundle.rds` (Factor 1 tubular-mean signature weights per gene) | same |
| `raw/F1_Tubular_mean_results/Factor1/top_biomarker_candidates_z_score.csv` (the colleague's 504 proteins with weight > 0 and z > 0; consistency check and composite rank) | same |
| `resources/Predicted_secreted_proteins.tsv`, `Predicted_membrane_proteins.tsv`, `Predicted_intracellular_proteins.tsv`, `Plasma_proteins.tsv` | Human Protein Atlas downloads, 8 Oct 2026 (`resources/HPA_download_date.txt`) |
| MSigDB Hallmark, Reactome, KEGG, GO BP, GO CC | R package msigdbr 7.5.1 |

## Outputs and the panels that use them
| output | panel | entry script |
|---|---|---|
| `outputs/fgsea_all_collections.csv` (ranking "Factor 1 weight, 984 matched proteins"; Hallmark + Reactome) | Fig. 4e | `figures/Fig4/Fig4e_fgsea_hallmark_reactome_factor1_weight.R` |

No script on this system drew the published Fig. 4e. The entry script is a **re-drawing** from the table: per sign of
NES, the 10 terms with the smallest adjusted p, Hallmark and Reactome pooled, ordered by NES, dark bars for adjusted
p < 0.05. Checked against the published page on 12 Nov 2026: all 20 labels and their order (PDF text layer), the shading,
and the NES of every bar (measured from the bar ends, within 0.004) match (`build/verification/atac_snrna.csv`).

## Caveats (read before reusing)
- **raw/ is not public.** 01 and 04 cannot be re-run without the UK Biobank-derived files; the fgsea table itself holds
  only gene-set-level statistics.
- **Signature naming (known discrepancy).** The manuscript methods describe this signature as factor 2, with 1,976
  proteins common to both datasets and median-split quadrants; the bundle used here is labelled Factor1 and only 984 of
  the 2,923 Olink proteins carry a weight (Fig. 4e says 984). To be reconciled with the colleague before publication.
- **Adjusted p in Fig. 4e.** fgsea ran per collection, so the BH adjustment is within Hallmark and within Reactome, not
  over the pooled 20 terms. Several shaded-light terms sit just above 0.05 (IL6 JAK STAT3 0.054, EMT 0.058,
  Spermatogenesis 0.058); fgsea p-values come from its multilevel sampling (seed 42, fgsea 1.32.2), so another seed or
  version can move them across 0.05.
- **Re-running 04.** 04 extends `protein_table.csv` in place and reads it at the start; re-run 01 before re-running 04,
  otherwise the merge duplicates the rank columns.
- 1,939 of the 2,923 proteins have no Factor 1 weight and are not in the weight ranking.

# Figure panels and their code

Every panel of the final figures (main, Extended Data, Supplementary; figure export of 16 Nov 2026). `made_by`: **this repository** = code here (entry script in `figures/`, analysis steps in `analyses/`); MCFA = [saezlab/kidney_biopsies_MCFA]; atlas = snRNA/snATAC atlas construction; colleague / UK Biobank = other authors' code (not in this folder; where public, the note names the folder under `Code/`); wet lab / image = no analysis code; not found = no code for this panel on the analysis server. `re-run check`: result of re-running the entry script on the saved analysis outputs (12 Nov 2026) and comparing it with the file used in the published figure.


## Fig. 1

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | study overview schematic | image |  |  |  |
| b | cohort overview (condition, modalities, biopsy age, follow-up) | not found |  |  | not produced by code on this system (likely the MCFA repository); clinical workbook |
| c | follow-up per disease | MCFA |  |  | 1c_1e_followup_boxplot_and_eGFR_vs_slope.R |
| d | clinical endpoints per disease | not found |  |  |  |
| e | eGFR slope vs eGFR at biopsy | MCFA |  |  | 1c_1e_followup_boxplot_and_eGFR_vs_slope.R |
| f | snRNA-seq UMAP (363,656 nuclei) | atlas |  |  |  |
| g | snATAC-seq UMAP (63,386 nuclei) | atlas |  |  | ArchR workflow |

## Fig. 2

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | MCFA schematic | image |  |  |  |
| b | MDS of samples by condition / baseline eGFR | MCFA |  |  |  |
| c | MDS of samples by condition / baseline eGFR | MCFA |  |  |  |
| d | factor x clinical and cell-type abundance associations | MCFA |  |  | 2d_plot_metadata_associations.R |
| e | factor 1 vs baseline eGFR and aPT/PT ratio | MCFA |  |  |  |
| f | motif enrichment in PT-open vs aPT-open peaks | not found |  |  | ArchR/GLUE collaborator |
| g | eGFR GWAS enrichment (gchromVAR) and disrupted TF motifs (motifbreakR) | colleague |  |  |  |
| h | eGFR GWAS enrichment (gchromVAR) and disrupted TF motifs (motifbreakR) | colleague |  |  |  |

## Fig. 3

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | factor 1 loadings per view | not found |  |  | likely MCFA repository |
| b | Xenium design schematic | image |  |  |  |
| c | mean factor 1 score (ULM) per core vs eGFR | this repository | `figures/Fig3/Fig3c_factor1_per_core_vs_egfr.py` | pixel-identical to the published source file |  |
| d | Xenium transcripts + PhenoCycler protein images | image |  |  |  |
| e | RNF152 / NFKBIZ in PT vs aPT per core | colleague |  |  | Xenium team notebook |
| f | distance-bounded Voronoi graph and one 2-hop neighbourhood (core K16) | this repository | `figures/Fig3/Fig3f_voronoi_2hop_neighbourhood.py` | identical except the title, which now shows the public core id (plus anti-aliasing specks) |  |
| g | cell-type composition of 2-hop neighbourhoods vs factor 1 | this repository | `figures/Fig3/Fig3g_celltype_composition_vs_factor1.py` | pixel-identical to the published source file |  |
| h | cell-cell interactions aPT -> immune vs factor 1 | this repository | `figures/Fig3/Fig3h_cci_apt_to_immune_vs_factor1.py` | pixel-identical to the published source file |  |
| i | APP / CD74 maps (K82, K18) and APP-CD74 score with aPT vs PT as sender | this repository | `figures/Fig3/Fig3i_app_cd74_apt_vs_pt.py` | pixel-identical to the published source file |  |

## Fig. 4

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | UK Biobank: phenome associations, CKD incidence, proteins, dose-response | UK Biobank |  |  | first-submission scripts in Code/ukbb_processing of this GitHub repository |
| b | UK Biobank: phenome associations, CKD incidence, proteins, dose-response | UK Biobank |  |  | first-submission scripts in Code/ukbb_processing of this GitHub repository |
| c | UK Biobank: phenome associations, CKD incidence, proteins, dose-response | UK Biobank |  |  | first-submission scripts in Code/ukbb_processing of this GitHub repository |
| d | UK Biobank: phenome associations, CKD incidence, proteins, dose-response | UK Biobank |  |  | first-submission scripts in Code/ukbb_processing of this GitHub repository |
| e | Hallmark & Reactome fgsea ranked by factor 1 tubular weight (984 proteins) | this repository | `figures/Fig4/Fig4e_fgsea_hallmark_reactome_factor1_weight.R` | re-drawn from the fgsea table (no original rendering code); terms, order, NES and shading match | fgsea table by Leon; the published rendering was re-drawn here |

## Fig. 5

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | Hallmark / PROGENy / CollecTRI on factor loadings | MCFA |  |  |  |
| b | Hallmark / PROGENy / CollecTRI on factor loadings | MCFA |  |  |  |
| c | Hallmark / PROGENy / CollecTRI on factor loadings | MCFA |  |  |  |
| d | RelA KO cells, bulk RNA-seq, organoids, mice | wet lab |  |  | RelA KO bulk RNA-seq analysis (g-j) not in this repository |
| e | RelA KO cells, bulk RNA-seq, organoids, mice | wet lab |  |  | RelA KO bulk RNA-seq analysis (g-j) not in this repository |
| f | RelA KO cells, bulk RNA-seq, organoids, mice | wet lab |  |  | RelA KO bulk RNA-seq analysis (g-j) not in this repository |
| g | RelA KO cells, bulk RNA-seq, organoids, mice | wet lab |  |  | RelA KO bulk RNA-seq analysis (g-j) not in this repository |
| h | RelA KO cells, bulk RNA-seq, organoids, mice | wet lab |  |  | RelA KO bulk RNA-seq analysis (g-j) not in this repository |
| i | RelA KO cells, bulk RNA-seq, organoids, mice | wet lab |  |  | RelA KO bulk RNA-seq analysis (g-j) not in this repository |
| j | RelA KO cells, bulk RNA-seq, organoids, mice | wet lab |  |  | RelA KO bulk RNA-seq analysis (g-j) not in this repository |
| k | RelA KO cells, bulk RNA-seq, organoids, mice | wet lab |  |  | RelA KO bulk RNA-seq analysis (g-j) not in this repository |
| l | RelA KO cells, bulk RNA-seq, organoids, mice | wet lab |  |  | RelA KO bulk RNA-seq analysis (g-j) not in this repository |
| m | RelA KO cells, bulk RNA-seq, organoids, mice | wet lab |  |  | RelA KO bulk RNA-seq analysis (g-j) not in this repository |
| n | RelA KO cells, bulk RNA-seq, organoids, mice | wet lab |  |  | RelA KO bulk RNA-seq analysis (g-j) not in this repository |
| o | RelA KO cells, bulk RNA-seq, organoids, mice | wet lab |  |  | RelA KO bulk RNA-seq analysis (g-j) not in this repository |
| p | RelA KO cells, bulk RNA-seq, organoids, mice | wet lab |  |  | RelA KO bulk RNA-seq analysis (g-j) not in this repository |
| q | RelA KO cells, bulk RNA-seq, organoids, mice | wet lab |  |  | RelA KO bulk RNA-seq analysis (g-j) not in this repository |
| r | RelA KO cells, bulk RNA-seq, organoids, mice | wet lab |  |  | RelA KO bulk RNA-seq analysis (g-j) not in this repository |

## Fig. 6

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | cell-type composition of the glomerular niches | this repository | `figures/Fig6/Fig6a_niche_vs_celltype_heatmap.py` | pixel-identical to the published source file |  |
| b | glomerular niche maps of K50 (crescent) and K18 | this repository | `figures/Fig6/Fig6b_niche_maps_K50_K18.py` | identical except the title, which now shows the public core id | H&E and protein images are microscopy |
| c | Xenium marker genes in glomerular niche 20 | this repository | `figures/Fig6/Fig6c_niche20_marker_dotplot.py` | byte-identical |  |
| d | CCI in glomerular niche 20; PEC as sender | colleague |  |  | LIANA + CrossTalkeR; see Code/CCI_Niches in this GitHub repository |
| e | CCI in glomerular niche 20; PEC as sender | colleague |  |  | LIANA + CrossTalkeR; see Code/CCI_Niches in this GitHub repository |
| f | PEC marker genes, niche 20 vs other niches | this repository | `figures/Fig6/Fig6f_pec_niche20_vs_other_dotplot.py` | pixel-identical to the published source file |  |
| g | snRNA-seq PEC and immune subclusters mapped onto Xenium (K50, K18) | this repository | `figures/Fig6/Fig6g_immune_subclusters_xenium.py`<br>`figures/Fig6/Fig6g_pec_subclusters_xenium.py` | identical except the title, which now shows the public core id |  |
| h | PHATE of the PEC subclusters and crescent score | this repository | `figures/Fig6/Fig6h_pec_phate_clusters_crescent.py` | pixel-identical to the published source file |  |
| i | VIA pseudotime of the PEC trajectory | this repository | `figures/Fig6/Fig6i_pec_via_pseudotime_curve.py` | not re-drawable: the VIA object was never saved; the script points to the analysis step |  |
| j | marker genes along pseudotime; pseudotime heatmap | this repository | `figures/Fig6/Fig6j_genetrends.py`<br>`figures/Fig6/Fig6j_pseudotime_heatmap.py` | heatmap pixel-identical; gene trends not re-drawable (VIA object not saved) |  |
| k | aPEC (subcluster 9) vs PEC subclusters 0/2/4/11, volcano | this repository | `figures/Fig6/Fig6k_volcano_aPEC_vs_PEC.py` | pixel-identical to the published source file |  |
| l | Reactome GSEA, aPEC vs other PEC | this repository | `figures/Fig6/Fig6l_gsea_reactome_dotplot.R` | pixel-identical to the published source file |  |

## Extended Data Fig. 1

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | cohort overlap, age, endpoints, markers, ATAC motifs, composition, co-embedding | atlas |  |  | a-c clinical tables, d-j atlas (snRNA/ArchR/GLUE); no code on this system |
| b | cohort overlap, age, endpoints, markers, ATAC motifs, composition, co-embedding | atlas |  |  | a-c clinical tables, d-j atlas (snRNA/ArchR/GLUE); no code on this system |
| c | cohort overlap, age, endpoints, markers, ATAC motifs, composition, co-embedding | atlas |  |  | a-c clinical tables, d-j atlas (snRNA/ArchR/GLUE); no code on this system |
| d | cohort overlap, age, endpoints, markers, ATAC motifs, composition, co-embedding | atlas |  |  | a-c clinical tables, d-j atlas (snRNA/ArchR/GLUE); no code on this system |
| e | cohort overlap, age, endpoints, markers, ATAC motifs, composition, co-embedding | atlas |  |  | a-c clinical tables, d-j atlas (snRNA/ArchR/GLUE); no code on this system |
| f | cohort overlap, age, endpoints, markers, ATAC motifs, composition, co-embedding | atlas |  |  | a-c clinical tables, d-j atlas (snRNA/ArchR/GLUE); no code on this system |
| g | cohort overlap, age, endpoints, markers, ATAC motifs, composition, co-embedding | atlas |  |  | a-c clinical tables, d-j atlas (snRNA/ArchR/GLUE); no code on this system |
| h | cohort overlap, age, endpoints, markers, ATAC motifs, composition, co-embedding | atlas |  |  | a-c clinical tables, d-j atlas (snRNA/ArchR/GLUE); no code on this system |
| i | cohort overlap, age, endpoints, markers, ATAC motifs, composition, co-embedding | atlas |  |  | a-c clinical tables, d-j atlas (snRNA/ArchR/GLUE); no code on this system |
| j | cohort overlap, age, endpoints, markers, ATAC motifs, composition, co-embedding | atlas |  |  | a-c clinical tables, d-j atlas (snRNA/ArchR/GLUE); no code on this system |

## Extended Data Fig. 2

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | MCFA model selection, variance explained, factor scores, robustness | MCFA |  |  |  |
| b | MCFA model selection, variance explained, factor scores, robustness | MCFA |  |  |  |
| c | MCFA model selection, variance explained, factor scores, robustness | MCFA |  |  |  |
| d | MCFA model selection, variance explained, factor scores, robustness | MCFA |  |  |  |
| e | CLR cell-type abundance by condition / eGFR / IF-TA | this repository | `figures/ED2/ED2e-g_clr_compositions_boxplot.R` | statistics identical; jittered points differ (original jitter unseeded) | plotting code adapted verbatim from the MCFA repository |
| f | CLR cell-type abundance by condition / eGFR / IF-TA | this repository | `figures/ED2/ED2e-g_clr_compositions_boxplot.R` | as e | plotting code adapted verbatim from the MCFA repository |
| g | CLR cell-type abundance by condition / eGFR / IF-TA | this repository | `figures/ED2/ED2e-g_clr_compositions_boxplot.R` | as e | plotting code adapted verbatim from the MCFA repository |
| h | HNF4A/RELA footprints; ATAC tracks at RHEX and VCAM1 | atlas |  |  | ArchR; no code on this system |
| i | HNF4A/RELA footprints; ATAC tracks at RHEX and VCAM1 | atlas |  |  | ArchR; no code on this system |
| j | HNF4A/RELA footprints; ATAC tracks at RHEX and VCAM1 | atlas |  |  | ArchR; no code on this system |

## Extended Data Fig. 3

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | PT/aPT embeddings, DE/DA volcanos, peak annotation, motif enrichment, GSEA | not found |  |  | ArchR/GLUE collaborator — ask Leon |
| b | PT/aPT embeddings, DE/DA volcanos, peak annotation, motif enrichment, GSEA | not found |  |  | ArchR/GLUE collaborator — ask Leon |
| c | PT/aPT embeddings, DE/DA volcanos, peak annotation, motif enrichment, GSEA | not found |  |  | ArchR/GLUE collaborator — ask Leon |
| d | PT/aPT embeddings, DE/DA volcanos, peak annotation, motif enrichment, GSEA | not found |  |  | ArchR/GLUE collaborator — ask Leon |
| e | PT/aPT embeddings, DE/DA volcanos, peak annotation, motif enrichment, GSEA | not found |  |  | ArchR/GLUE collaborator — ask Leon |
| f | PT/aPT embeddings, DE/DA volcanos, peak annotation, motif enrichment, GSEA | not found |  |  | ArchR/GLUE collaborator — ask Leon |
| g | PT/aPT embeddings, DE/DA volcanos, peak annotation, motif enrichment, GSEA | not found |  |  | ArchR/GLUE collaborator — ask Leon |
| h | PT/aPT embeddings, DE/DA volcanos, peak annotation, motif enrichment, GSEA | not found |  |  | ArchR/GLUE collaborator — ask Leon |
| i | PT/aPT embeddings, DE/DA volcanos, peak annotation, motif enrichment, GSEA | not found |  |  | ArchR/GLUE collaborator — ask Leon |
| j | scMEGA PT->aPT trajectory (joint diffusion map, ArchR pseudotime) | this repository | `figures/ED3/ED3j-k_archr_trajectory_egfr_density.py` | pixel-identical to the published source file |  |
| k | eGFR groups along the pseudotime | this repository | `figures/ED3/ED3j-k_archr_trajectory_egfr_density.py` | pixel-identical to the published source file |  |
| l | eGRN modules I and II with JUN / HNF4A insets | this repository | `figures/ED3/ED3l_egrn_two_modules.py`<br>`figures/ED3/ED3l_insets_scMEGA_TF_JUN_HNF4A.R` | pixel-identical to the published source file (network and insets) |  |
| m | TF activity, TF expression and target expression of the modules along pseudotime | this repository | `figures/ED3/ED3m_egrn_module_timing.py` | top and middle pixel-identical; the printed bottom panel comes from the September run, the script draws the AKI/pt_C5-excluded rerun |  |

## Extended Data Fig. 4

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | post-Xenium H&E; RNA + protein images | image |  |  |  |
| b | post-Xenium H&E; RNA + protein images | image |  |  |  |
| c | cell-type composition of 2-hop neighbourhoods | this repository | `figures/ED4/ED4c_2hop_celltype_heatmap.py` | pixel-identical to the published source file |  |
| d | cell-type maps (K82 FSGS, K18 TIN) | this repository | `figures/ED4/ED4d-e_celltype_niche_maps.py` | identical except the title, which now shows the public core id |  |
| e | niche map (K18) | this repository | `figures/ED4/ED4d-e_celltype_niche_maps.py` | identical except the title, which now shows the public core id |  |
| f | niche composition of 2-hop neighbourhoods vs factor 1 | this repository | `figures/ED4/ED4f_niche_composition_vs_factor1.py` | pixel-identical to the published source file |  |
| g | cell-type composition of the niches | this repository | `figures/ED4/ED4g_niche_vs_celltype_heatmap.py` | pixel-identical to the published source file |  |
| h | niche enrichment by eGFR group | this repository | `figures/ED4/ED4h_egfr_group_enrichment.py` | pixel-identical to the published source file |  |
| i | condition composition of the niches | this repository | `figures/ED4/ED4i_condition_stacked_bars.py` | pixel-identical to the published source file |  |
| j | Hallmark pathway scores per niche | this repository | `figures/ED4/ED4j_hallmark_niches.py` | pixel-identical to the published source file |  |
| k | top ligands in niche 1 with aPT as sender | colleague |  |  | LIANA + CrossTalkeR; see Code/CCI_Niches in this GitHub repository |

## Extended Data Fig. 5

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | NELL1 IHC/IF, overexpression, bulk RNA-seq | wet lab |  |  |  |
| b | NELL1 IHC/IF, overexpression, bulk RNA-seq | wet lab |  |  |  |
| c | NELL1 IHC/IF, overexpression, bulk RNA-seq | wet lab |  |  |  |
| d | NELL1 IHC/IF, overexpression, bulk RNA-seq | wet lab |  |  |  |
| e | NELL1 IHC/IF, overexpression, bulk RNA-seq | wet lab |  |  |  |
| f | NELL1 IHC/IF, overexpression, bulk RNA-seq | wet lab |  |  |  |
| g | NELL1 IHC/IF, overexpression, bulk RNA-seq | wet lab |  |  |  |
| h | NELL1 IHC/IF, overexpression, bulk RNA-seq | wet lab |  |  |  |
| i | NELL1 IHC/IF, overexpression, bulk RNA-seq | wet lab |  |  |  |
| j | NELL1 IHC/IF, overexpression, bulk RNA-seq | wet lab |  |  |  |

## Extended Data Fig. 6

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | 2-hop composition vs NF-kB activity (PROGENy) | this repository | `figures/ED6/ED6a_celltype_composition_vs_nfkb.py` | pixel-identical to the published source file |  |
| b | RelA KO construct, ICE, Sanger, bulk RNA-seq PCA/volcano/GSEA | wet lab |  |  | bulk analysis not in this repository |
| c | RelA KO construct, ICE, Sanger, bulk RNA-seq PCA/volcano/GSEA | wet lab |  |  | bulk analysis not in this repository |
| d | RelA KO construct, ICE, Sanger, bulk RNA-seq PCA/volcano/GSEA | wet lab |  |  | bulk analysis not in this repository |
| e | RelA KO construct, ICE, Sanger, bulk RNA-seq PCA/volcano/GSEA | wet lab |  |  | bulk analysis not in this repository |
| f | RelA KO construct, ICE, Sanger, bulk RNA-seq PCA/volcano/GSEA | wet lab |  |  | bulk analysis not in this repository |
| g | RelA KO construct, ICE, Sanger, bulk RNA-seq PCA/volcano/GSEA | wet lab |  |  | bulk analysis not in this repository |

## Extended Data Fig. 7

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | APP IHC, ciPEC APP overexpression, co-culture bulk RNA-seq | wet lab |  |  |  |
| b | APP IHC, ciPEC APP overexpression, co-culture bulk RNA-seq | wet lab |  |  |  |
| c | APP IHC, ciPEC APP overexpression, co-culture bulk RNA-seq | wet lab |  |  |  |
| d | APP IHC, ciPEC APP overexpression, co-culture bulk RNA-seq | wet lab |  |  |  |
| e | APP IHC, ciPEC APP overexpression, co-culture bulk RNA-seq | wet lab |  |  |  |
| f | APP IHC, ciPEC APP overexpression, co-culture bulk RNA-seq | wet lab |  |  |  |
| g | APP IHC, ciPEC APP overexpression, co-culture bulk RNA-seq | wet lab |  |  |  |
| h | APP IHC, ciPEC APP overexpression, co-culture bulk RNA-seq | wet lab |  |  |  |
| i | APP IHC, ciPEC APP overexpression, co-culture bulk RNA-seq | wet lab |  |  |  |
| j | APP IHC, ciPEC APP overexpression, co-culture bulk RNA-seq | wet lab |  |  |  |
| k | APP IHC, ciPEC APP overexpression, co-culture bulk RNA-seq | wet lab |  |  |  |
| l | APP IHC, ciPEC APP overexpression, co-culture bulk RNA-seq | wet lab |  |  |  |

## Supplementary Fig. 1

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | eGFR trajectories and slopes | not found |  |  | an earlier version (54-sample cohort) exists in egfr_trajectories_fig1e; the final version was not made on this system |
| b | eGFR trajectories and slopes | not found |  |  | an earlier version (54-sample cohort) exists in egfr_trajectories_fig1e; the final version was not made on this system |
| c | eGFR trajectories and slopes | not found |  |  | an earlier version (54-sample cohort) exists in egfr_trajectories_fig1e; the final version was not made on this system |
| d | eGFR trajectories and slopes | not found |  |  | an earlier version (54-sample cohort) exists in egfr_trajectories_fig1e; the final version was not made on this system |
| e | eGFR trajectories and slopes | not found |  |  | an earlier version (54-sample cohort) exists in egfr_trajectories_fig1e; the final version was not made on this system |
| f | eGFR trajectories and slopes | not found |  |  | an earlier version (54-sample cohort) exists in egfr_trajectories_fig1e; the final version was not made on this system |

## Supplementary Fig. 2

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | snRNA-seq QC: UMAPs, per-sample metrics, composition | this repository | `figures/S2/S2a-b_snrna_overview_qc.py` | pixel-identical to the published source file |  |
| b | snRNA-seq QC: UMAPs, per-sample metrics, composition | this repository | `figures/S2/S2a-b_snrna_overview_qc.py` | pixel-identical to the published source file |  |
| c | QC metrics vs biopsy age | not found |  |  |  |

## Supplementary Fig. 3

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | aPT2 nuclei: marker dot plots, per-cell-type QC metrics, aPT2 fraction vs biopsy age / SoupX fraction / eGFR | not found |  |  | new in the export of 16 Nov 2026; no code on this system |
| b | aPT2 nuclei: marker dot plots, per-cell-type QC metrics, aPT2 fraction vs biopsy age / SoupX fraction / eGFR | not found |  |  | new in the export of 16 Nov 2026; no code on this system |
| c | aPT2 nuclei: marker dot plots, per-cell-type QC metrics, aPT2 fraction vs biopsy age / SoupX fraction / eGFR | not found |  |  | new in the export of 16 Nov 2026; no code on this system |
| d | aPT2 nuclei: marker dot plots, per-cell-type QC metrics, aPT2 fraction vs biopsy age / SoupX fraction / eGFR | not found |  |  | new in the export of 16 Nov 2026; no code on this system |

## Supplementary Fig. 4

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | MCFA factors 2-5: scores by condition, clinical / composition / QC associations, model-size similarity, K68, soup fraction, UMIs, loadings | not found |  |  | new in the export of 16 Nov 2026; no code on this system (MCFA analysis) |
| b | MCFA factors 2-5: scores by condition, clinical / composition / QC associations, model-size similarity, K68, soup fraction, UMIs, loadings | not found |  |  | new in the export of 16 Nov 2026; no code on this system (MCFA analysis) |
| c | MCFA factors 2-5: scores by condition, clinical / composition / QC associations, model-size similarity, K68, soup fraction, UMIs, loadings | not found |  |  | new in the export of 16 Nov 2026; no code on this system (MCFA analysis) |
| d | MCFA factors 2-5: scores by condition, clinical / composition / QC associations, model-size similarity, K68, soup fraction, UMIs, loadings | not found |  |  | new in the export of 16 Nov 2026; no code on this system (MCFA analysis) |
| e | MCFA factors 2-5: scores by condition, clinical / composition / QC associations, model-size similarity, K68, soup fraction, UMIs, loadings | not found |  |  | new in the export of 16 Nov 2026; no code on this system (MCFA analysis) |
| f | MCFA factors 2-5: scores by condition, clinical / composition / QC associations, model-size similarity, K68, soup fraction, UMIs, loadings | not found |  |  | new in the export of 16 Nov 2026; no code on this system (MCFA analysis) |
| g | MCFA factors 2-5: scores by condition, clinical / composition / QC associations, model-size similarity, K68, soup fraction, UMIs, loadings | not found |  |  | new in the export of 16 Nov 2026; no code on this system (MCFA analysis) |
| h | MCFA factors 2-5: scores by condition, clinical / composition / QC associations, model-size similarity, K68, soup fraction, UMIs, loadings | not found |  |  | new in the export of 16 Nov 2026; no code on this system (MCFA analysis) |

## Supplementary Fig. 5

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | factor 1 with vs without ligand-receptor views (scores, loadings, per-view r) | not found |  |  | new in the export of 16 Nov 2026; no code on this system (MCFA analysis) |
| b | factor 1 with vs without ligand-receptor views (scores, loadings, per-view r) | not found |  |  | new in the export of 16 Nov 2026; no code on this system (MCFA analysis) |
| c | factor 1 with vs without ligand-receptor views (scores, loadings, per-view r) | not found |  |  | new in the export of 16 Nov 2026; no code on this system (MCFA analysis) |
| d | factor 1 with vs without ligand-receptor views (scores, loadings, per-view r) | not found |  |  | new in the export of 16 Nov 2026; no code on this system (MCFA analysis) |

## Supplementary Fig. 6

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | snRNA-seq QC metrics vs biopsy age (same plot as Supplementary Fig. 2c) | not found |  |  | new in the export of 16 Nov 2026; no code on this system |
| b | factor 1 vs biopsy age; MCFA refit on biopsies < 10 years | not found |  |  | new in the export of 16 Nov 2026; no code on this system (MCFA analysis) |
| c | factor 1 vs biopsy age; MCFA refit on biopsies < 10 years | not found |  |  | new in the export of 16 Nov 2026; no code on this system (MCFA analysis) |

## Supplementary Fig. 7

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | SoupX contamination vs factor 1; MCFA refit after masking ambient-dominated genes | not found |  |  | new in the export of 16 Nov 2026; no code on this system (MCFA analysis) |
| b | SoupX contamination vs factor 1; MCFA refit after masking ambient-dominated genes | not found |  |  | new in the export of 16 Nov 2026; no code on this system (MCFA analysis) |
| c | SoupX contamination vs factor 1; MCFA refit after masking ambient-dominated genes | not found |  |  | new in the export of 16 Nov 2026; no code on this system (MCFA analysis) |
| d | SoupX contamination vs factor 1; MCFA refit after masking ambient-dominated genes | not found |  |  | new in the export of 16 Nov 2026; no code on this system (MCFA analysis) |
| e | SoupX contamination vs factor 1; MCFA refit after masking ambient-dominated genes | not found |  |  | new in the export of 16 Nov 2026; no code on this system (MCFA analysis) |
| f | SoupX contamination vs factor 1; MCFA refit after masking ambient-dominated genes | not found |  |  | new in the export of 16 Nov 2026; no code on this system (MCFA analysis) |

## Supplementary Fig. 8

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | TF activity and expression along the PT-aPT trajectory | this repository | `figures/S8/S8a_tf_activity_expression_heatmap.R` | pixel-identical to the published source file |  |
| b | factor 1 top-positive genes in the eGRN, by number of regulating TFs | this repository | `figures/S8/S8b_egrn_factor1_targets_overview.py` | pixel-identical to the published source file | panel c in the PDFs of 12 Nov |
| c | TFs upstream of the factor 1 genes in the eGRN | this repository | `figures/S8/S8c_egrn_factor1_tf_ranking.py` | pixel-identical to the published source file | panel b in the PDFs of 12 Nov |
| d | FOSL2 regulon | this repository | `figures/S8/S8d_egrn_tf_centric_FOSL2.py` | pixel-identical to the published source file |  |

## Supplementary Fig. 9

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | Xenium label transfer concordance, UMAP, markers, cells per condition | this repository | `figures/S9/S9a_heatmap_pred_vs_refined.py` | pixel-identical when run on the archived 1 Jul object (--xenium); the final object has 182 aPT cells moved to PEC |  |
| b | Xenium label transfer concordance, UMAP, markers, cells per condition | this repository | `figures/S9/S9b_umap_level2_refined.py` | pixel-identical when run on the archived 1 Jul object (--xenium); the final object has 182 aPT cells moved to PEC |  |
| c | Xenium label transfer concordance, UMAP, markers, cells per condition | this repository | `figures/S9/S9c-d_dotplot_bubble.py` | pixel-identical when run on the archived 1 Jul object (--xenium); the final object has 182 aPT cells moved to PEC |  |
| d | Xenium label transfer concordance, UMAP, markers, cells per condition | this repository | `figures/S9/S9c-d_dotplot_bubble.py` | pixel-identical when run on the archived 1 Jul object (--xenium); the final object has 182 aPT cells moved to PEC |  |

## Supplementary Fig. 10

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | Xenium QC per core | colleague |  |  | Xenium team notebook; see Code/Preprocessing_Xenium in this GitHub repository |
| b | Xenium QC per core | colleague |  |  | Xenium team notebook; see Code/Preprocessing_Xenium in this GitHub repository |
| c | Xenium QC per core | colleague |  |  | Xenium team notebook; see Code/Preprocessing_Xenium in this GitHub repository |
| d | Xenium per-core QC vs biopsy age | this repository | `figures/S10/S10d_qc_metrics_vs_biopsy_age.py` | pixel-identical to the file on disk; the published version places control cores at 0-5 y (code: 24.8 y) |  |

## Supplementary Fig. 11

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | Xenium Explorer images of factor 1 genes (no panel letters) | image |  |  |  |

## Supplementary Fig. 12

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | per-core cell types, factor 1 score, niches | this repository | `figures/S12/S12a-b_core_maps.py` | identical except the title, which now shows the public core id |  |
| b | per-core cell types, factor 1 score, niches | this repository | `figures/S12/S12a-b_core_maps.py` | identical except the title, which now shows the public core id |  |

## Supplementary Fig. 13

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | NELL1 per cell type and vs eGFR | this repository | `figures/S13/S13a_nell1_dotplot.py` | pixel-identical to the published source file |  |
| b | NELL1 per cell type and vs eGFR | this repository | `figures/S13/S13b_nell1_vs_egfr_per_donor.py` | pixel-identical to the published source file |  |
| c | NELL1 immunohistochemistry (Human Protein Atlas) | image |  |  |  |

## Supplementary Fig. 14

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | lentiviral NELL1 construct and sorting strategy (wet lab) | wet lab |  |  |  |
| b | lentiviral NELL1 construct and sorting strategy (wet lab) | wet lab |  |  |  |
| c | lentiviral NELL1 construct and sorting strategy (wet lab) | wet lab |  |  |  |
| d | lentiviral NELL1 construct and sorting strategy (wet lab) | wet lab |  |  |  |

## Supplementary Fig. 15

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | uncropped immunoblots (wet lab) | wet lab |  |  |  |
| b | uncropped immunoblots (wet lab) | wet lab |  |  |  |
| c | uncropped immunoblots (wet lab) | wet lab |  |  |  |

## Supplementary Fig. 16

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | RELA KO design, ICE, blots, top DE genes | wet lab |  |  | DESeq2 analysis of e not in this repository |
| b | RELA KO design, ICE, blots, top DE genes | wet lab |  |  | DESeq2 analysis of e not in this repository |
| c | RELA KO design, ICE, blots, top DE genes | wet lab |  |  | DESeq2 analysis of e not in this repository |
| d | RELA KO design, ICE, blots, top DE genes | wet lab |  |  | DESeq2 analysis of e not in this repository |
| e | RELA KO design, ICE, blots, top DE genes | wet lab |  |  | DESeq2 analysis of e not in this repository |

## Supplementary Fig. 17

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | Rela deletion in mouse kidneys, RNAscope (wet lab) | wet lab |  |  |  |
| b | Rela deletion in mouse kidneys, RNAscope (wet lab) | wet lab |  |  |  |
| c | Rela deletion in mouse kidneys, RNAscope (wet lab) | wet lab |  |  |  |

## Supplementary Fig. 18

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | Hallmark per glomerular niche; niche-20 PEC DE; snRNA PEC / immune subclusters | this repository | `figures/S18/S18a_hallmark_glom_niches.py` | pixel-identical to the published source file |  |
| b | Hallmark per glomerular niche; niche-20 PEC DE; snRNA PEC / immune subclusters | this repository | `figures/S18/S18b_niche20_pec_volcano.R` | points, axes and text identical; gene-label placement differs (ggrepel unseeded) |  |
| c | Hallmark per glomerular niche; niche-20 PEC DE; snRNA PEC / immune subclusters | this repository | `figures/S18/S18c-d_pec_subclusters_umap_markers.py` | pixel-identical to the published source file |  |
| d | Hallmark per glomerular niche; niche-20 PEC DE; snRNA PEC / immune subclusters | this repository | `figures/S18/S18c-d_pec_subclusters_umap_markers.py` | pixel-identical to the published source file |  |
| e | Hallmark per glomerular niche; niche-20 PEC DE; snRNA PEC / immune subclusters | this repository | `figures/S18/S18e-f_immune_subclusters_umap_markers.py` | pixel-identical to the published source file |  |
| f | Hallmark per glomerular niche; niche-20 PEC DE; snRNA PEC / immune subclusters | this repository | `figures/S18/S18e-f_immune_subclusters_umap_markers.py` | pixel-identical to the published source file |  |

## Supplementary Fig. 19

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | APP770-overexpressing ciPECs (wet lab) | wet lab |  |  |  |
| b | APP770-overexpressing ciPECs (wet lab) | wet lab |  |  |  |
| c | APP770-overexpressing ciPECs (wet lab) | wet lab |  |  |  |
| d | APP770-overexpressing ciPECs (wet lab) | wet lab |  |  |  |

## Supplementary Fig. 20

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | flow cytometry of ciPECs and iPSC-derived immune cells (wet lab) | wet lab |  |  |  |
| b | flow cytometry of ciPECs and iPSC-derived immune cells (wet lab) | wet lab |  |  |  |
| c | flow cytometry of ciPECs and iPSC-derived immune cells (wet lab) | wet lab |  |  |  |
| d | flow cytometry of ciPECs and iPSC-derived immune cells (wet lab) | wet lab |  |  |  |

## Supplementary Fig. 21

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | niche-20 genes vs factor 1 loadings (Venn) | this repository | `figures/S21/S21a-f_venn_niche20_mofa_factor1.py` | pixel-identical to the published source file |  |
| b | niche-20 genes vs factor 1 loadings (Venn) | this repository | `figures/S21/S21a-f_venn_niche20_mofa_factor1.py` | pixel-identical to the published source file |  |
| c | niche-20 genes vs factor 1 loadings (Venn) | this repository | `figures/S21/S21a-f_venn_niche20_mofa_factor1.py` | pixel-identical to the published source file |  |
| d | niche-20 genes vs factor 1 loadings (Venn) | this repository | `figures/S21/S21a-f_venn_niche20_mofa_factor1.py` | pixel-identical to the published source file |  |
| e | niche-20 genes vs factor 1 loadings (Venn) | this repository | `figures/S21/S21a-f_venn_niche20_mofa_factor1.py` | pixel-identical to the published source file |  |
| f | niche-20 genes vs factor 1 loadings (Venn) | this repository | `figures/S21/S21a-f_venn_niche20_mofa_factor1.py` | pixel-identical to the published source file |  |

## Supplementary Fig. 22

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | niche-20 score along the crescent trajectory (Sultana et al.) | this repository | `figures/S22/S22a-i_niche20_validation_pseudotime.py` | pixel-identical to the published source file |  |
| b | niche-20 score along the crescent trajectory (Sultana et al.) | this repository | `figures/S22/S22a-i_niche20_validation_pseudotime.py` | pixel-identical to the published source file |  |
| c | niche-20 score along the crescent trajectory (Sultana et al.) | this repository | `figures/S22/S22a-i_niche20_validation_pseudotime.py` | pixel-identical to the published source file |  |
| d | niche-20 score along the crescent trajectory (Sultana et al.) | this repository | `figures/S22/S22a-i_niche20_validation_pseudotime.py` | pixel-identical to the published source file |  |
| e | niche-20 score along the crescent trajectory (Sultana et al.) | this repository | `figures/S22/S22a-i_niche20_validation_pseudotime.py` | pixel-identical to the published source file |  |
| f | niche-20 score along the crescent trajectory (Sultana et al.) | this repository | `figures/S22/S22a-i_niche20_validation_pseudotime.py` | pixel-identical to the published source file |  |
| g | niche-20 score along the crescent trajectory (Sultana et al.) | this repository | `figures/S22/S22a-i_niche20_validation_pseudotime.py` | pixel-identical to the published source file |  |
| h | niche-20 score along the crescent trajectory (Sultana et al.) | this repository | `figures/S22/S22a-i_niche20_validation_pseudotime.py` | near-identical (unseeded bootstrap CI whiskers) |  |
| i | niche-20 score along the crescent trajectory (Sultana et al.) | this repository | `figures/S22/S22a-i_niche20_validation_pseudotime.py` | near-identical (unseeded bootstrap CI whiskers) |  |

## Supplementary Fig. 23

| panel | content | made by | entry script(s) | re-run check | note |
|---|---|---|---|---|---|
| a | nuclei sorting gates (wet lab; no panel letters) | wet lab |  |  |  |

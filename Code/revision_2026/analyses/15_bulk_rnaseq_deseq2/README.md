# 15 — Bulk RNA-seq of cultured proximal tubule cells (DESeq2, fgsea)

## Purpose
Differential expression and gene-set enrichment for three bulk RNA-seq experiments on cultured human proximal tubule
cells. Script of Sidrah Maryam (`Bulk.R`, received 20 Nov 2026), unchanged except that its two folders can be set
through environment variables.

| experiment | groups (3 replicates each) | contrasts | panels |
|---|---|---|---|
| 1 | DMSO or etoposide pretreatment × vehicle or soluble APPα | sAPP vs vehicle (in DMSO and in etoposide), etoposide vs DMSO | by its design the soluble APPα experiment; not a figure panel of the paper |
| 2 | `Ca_rela_5_1`, `Ca_rela_5_2`, `Ca_rela_3_3`, `CaNT` (by their names: RELA guide clones and a non-targeting control) | each clone vs `CaNT` (`Ca_rela_3_3` vs `CaNT` plus one further control sample) | not confirmed |
| 3 | `CaNT1_rela`, `CaNT2_rela`, `Ca_rela_5_1`, `Ca_rela_5_2`, `Ca_rela_3_3` | each clone vs each control | not confirmed |

**Which contrasts of experiments 2 and 3 are the RELA knockout panels of the paper (Fig. 5g–j, Extended Data Fig. 6e–g)
has not been confirmed with the authors of those panels.** `FIGURE_PANELS.md` therefore does not list this script as
the source of a panel.

Method: gene-level counts (`star_salmon/salmon.merged.gene_counts.tsv`, the output layout of nf-core/rnaseq), summed per gene symbol; genes
with more than 10 counts over all samples are kept (RELA always). DESeq2 with design `~ condition` over all groups of an
experiment (separate two-group model when a contrast uses a group with added samples), Wald test, Benjamini-Hochberg
(`alpha = 0.05`); volcano plots with |log2 fold change| > 0.5 and adjusted P < 0.05; PCA on variance-stabilised counts.
`fgseaMultilevel` on the Wald statistic for MSigDB v2026.1 Hallmark, GO BP and Reactome; g:Profiler GO BP for the 500
genes with the highest and the lowest statistic; heat map of the 50 genes with the lowest adjusted P. Seed 42. Human,
gene symbols.

## Run
| script | environment | output |
|---|---|---|
| `Rscript Bulk.R <1\|2\|3>` | R with DESeq2, fgsea, gprofiler2, EnhancedVolcano, gplots, pheatmap, tidyverse | `$KB_BULK_DIR/experiment<n>/DEG_results/`: `pca.experiment<n>.pdf` and, per contrast, `DESEQ2.<A>.vs.<B>.txt`, `volcano.*.pdf`, `fgsea.<collection>.*.txt`, `GSEA_<collection>.*.pdf`, `GO_BP.<up\|dn>.*.pdf`, `heatmap.top50.*.pdf`, `pvalue_hist.*.pdf` |

`KB_BULK_DIR` (default `/beegfs/data/KramannLab/RelA_Bulk_TS/Processed`) holds `star_salmon/` and receives the results;
`KB_BULK_GMT_DIR` (default `/beegfs/data/KramannLab/RelA_Bulk_TS/reference`) holds the three MSigDB `.gmt` files.

## Inputs
| input | producer |
|---|---|
| `$KB_BULK_DIR/star_salmon/salmon.merged.gene_counts.tsv` (samples TS390–TS428) | read alignment and quantification (file layout of nf-core/rnaseq); external, not in this repository |
| `h.all.v2026.1.Hs.symbols.gmt`, `c5.all.v2026.1.Hs.symbols.gmt`, `c2.all.v2026.1.Hs.symbols.gmt` | MSigDB download |

## Not checked
**This script was not run for this repository**: the count matrix is on another cluster and is not available where the
other analyses were checked. The file parses; nothing else was verified, and no output of it was compared with a
figure.

## Caveats (read before reusing)
- Three replicates per group. With `design = ~ condition` over all groups, the dispersion estimates of a contrast
  borrow from the other groups of the experiment.
- Control sample TS414 is used both as a fourth `CaNT` sample in experiment 2 and as `CaNT1_rela` in experiment 3.
- RELA is kept through the count filter and, if DESeq2's independent filtering removes it, re-inserted into the volcano
  plot from a run without independent filtering and Cook's cutoff; its adjusted P then comes from a different
  correction than that of the other genes.
- g:Profiler needs internet access and its results depend on the database release at run time.
- The paths of the published run are those of another cluster; set `KB_BULK_DIR` and `KB_BULK_GMT_DIR`.

# Code for calculating CCI for condition vs rest all diseases in different niches

#Load libraries
library(Seurat)
library(tidyverse)
library(dplyr)
library(zellkonverter)
require(biomaRt)
require(SingleCellExperiment)
require(liana)
require(purrr)
require(reticulate)
library(magrittr)
library(SeuratObject)
require(tibble)
library(patchwork)
library(rmarkdown)
library(tinytex)

#Read Seurat file 
S = readRDS('input_file.rds')
Idents(S) = S$latent_leiden_0.4
niches = unique(S$latent_leiden_0.4)
condition = c('CTRL', 'DN', 'TIN', 'MN', 'IgAN', 'FSGS')

for (ni in niches) {
  #for each niche
  niche_dir = file.path('path/to_each/niche/output/folder/creating/', paste0('niche_', ni)) 
  dir.create(niche_dir, showWarnings = FALSE)
  s_sub = subset(S, idents = ni)
  data = as.SingleCellExperiment(s_sub)
  for (cond in condition) {
    ctr <-assays(data)[["counts"]][,colData(data)$condition!=cond]
    metadata_ctr <-tibble::tibble(rownames(colData(data))[colData(data)$condition!=cond],
                                  as.character(colData(data)$cell_type_manual[colData(data)$condition!=cond]))
    #create seurat object for running liana for the rest condition
    control = CreateSeuratObject(counts = ctr)
    Idents(control) = metadata_ctr$"as.character(...)"
    #Run liana for CCI using cellphonedb
    control = NormalizeData(control, normalization.method = "LogNormalize", scale.factor = 10000)
    liana_test <-  liana_wrap(control, method = 'cellphonedb')
    saveRDS(liana_test, file = paste0('not', cond, '_', ni,'.rds'))
    #### for positive
    ctr <-assays(data)[["counts"]][,colData(data)$condition==cond]
    metadata_ctr <-tibble::tibble(rownames(colData(data))[colData(data)$condition==cond],
                                  as.character(colData(data)$cell_type_manual[colData(data)$condition==cond]))
    #create seurat object for running liana for conditions separately
    control = CreateSeuratObject(counts = ctr)
    Idents(control) = metadata_ctr$"as.character(...)"
    #Run liana for both conditions
    control = NormalizeData(control, normalization.method = "LogNormalize", scale.factor = 10000)
    liana_test <-  liana_wrap(control, method = 'cellphonedb')
    saveRDS(liana_test, file = paste0(cond,'_',ni, '.rds'))
  }
}

# Downstream with Crosstalker

library(CrossTalkeR)
library(rmarkdown)
library(tinytex)

liana_test = readRDS('Reference.rds') # the condition file which is taken as the baseline to compare the interactions with 
liana_test2 = readRDS('Condition.rds') # condition file for which the interaction has to be checked

gud <-liana_test%>%
  filter(pvalue<=0.05) %>%
  mutate(gene_A=gsub('_','*',ligand.complex)) %>%
  mutate(gene_B=gsub('_','*',receptor.complex)) %>%
  mutate(source=gsub("_","",source)) %>%
  mutate(target=gsub("_","",target)) %>%
  mutate(type_gene_A = rep('Ligand',length(.data$ligand))) %>%
  mutate(type_gene_B = rep('Receptor',length(.data$receptor))) 

pur <-liana_test2 %>% 
  filter(pvalue<=0.05) %>%
  mutate(gene_A=gsub('_','*',ligand.complex)) %>%
  mutate(gene_B=gsub('_','*',receptor.complex)) %>%
  mutate(source=gsub("_","",source)) %>%
  mutate(target=gsub("_","",target)) %>%
  mutate(type_gene_A = rep('Ligand',length(.data$ligand))) %>%
  mutate(type_gene_B = rep('Receptor',length(.data$receptor))) 

#create a list for input in generating CCI
l = list()
l[['Reference']] = gud
l[['Condition']] = pur

data_dis <- generate_report(lrpaths = l,
                            genes=NULL,
                            out_path='/output/folder/',
                            threshold=0,
                            out_file = 'ctr_vs_rest.html',
                            output_fmt = "html_document",
                            report = T,
                            sel_columns = c('source','target','gene_A','gene_B','type_gene_A','type_gene_B','lr.mean'))

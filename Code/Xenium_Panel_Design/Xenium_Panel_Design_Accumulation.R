library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

setwd("/Users/tbleckwehl/Documents/Kidney_Biopsy_Project/Xenium/")

#SpaPros - probeset.csv and manually added markers
#MOFA - Factor2_top50loadingsperview_transcriptome.csv
#LIANA+ - Factor2_top50loadingsperview_LIANA.csv and 
#DEG - top_10_unique_genes_by_celltype_disease.csv and Disease specific LIANA results from Sidrah "Disease"_LR.csv

# Spapros Panel
SpaPros_Panel <- read.csv("probeset.csv")

#Filter for lncRNA and Xist
SpaPros_Panel <- SpaPros_Panel %>%
  rename("Gene" = X) %>% 
  mutate(Method = "SpaPros", view = "Subclustering") %>%
  filter(selection=="True") %>%
  filter(Gene != "XIST") %>%
  filter(!str_starts(Gene, "LINC")) %>%
  select(Gene, Method, view)

# MOFA Factor 2
MOFA_F2 <- read.csv("Factor2_top50loadingsperview_transcriptome.csv")
MOFA_F2 <- MOFA_F2 %>% 
  filter(view == "CNT_PC" | view =="DCT" | view =="TAL" | view == "PEC_POD" | view == "PT") %>%
  rename("Gene" = variable) %>% 
  mutate(Method = "MOFA") %>%
  filter(Loading>0) %>%
  filter(!str_starts(Gene, "LINC")) %>%
  select(Gene, Method, view)

#Selected cell types
LIANA_CCI <- c("PEC_POD&EC", "PEC_POD&EC-GC", "PEC_POD&Fib_VSMC", "PEC_POD&Immune", "PEC_POD&aPT_aTAL_aDCT", "PEC_POD&PT",
               "CNT_PC&DCT", "CNT_PC&TAL",
               "Immune&EC", "Immune&EC-GC", "Immune&Fib_VSMC", "Immune&PEC_POD", "Immune&PT", "Immune&TAL",
               "PT&Immune", "PT&PEC_POD", "aPT_aTAL_aDCT&PEC_POD", "dPT_dTAL_dDCT&PEC_POD")

LIANA_LR <- read.csv("Factor2_top50loadingsperview_LIANA.csv")
LIANA_LR <- LIANA_LR %>% 
  filter(view %in% LIANA_CCI) %>%
  rename("Gene" = variable) %>% 
  mutate(Method = "LIANA") %>%
  filter(Loading>0) %>%
  separate_rows(Gene, sep = "\\^|_") %>%
  filter(!str_starts(Gene, "LINC")) %>%
  select(Gene, Method, view)

LIANA_Sidrah <- data.frame()
for(Disease in Diseases){
  print(Disease)
  Kidney_Disease <- read.csv(paste0("./LIANA_Sidrah/", Disease, "_LR.csv"))
  Kidney_Disease$view = Disease
  Kidney_Disease$Method = "LIANA_Sidrah"
  
  Kidney_Disease <- Kidney_Disease %>%
    pivot_longer(
      cols = c(ligand, receptor),
      names_to = "CT",
      values_to = "Gene"
    )
  
  Kidney_Disease <- Kidney_Disease %>% 
    arrange(desc(interaction_stat)) %>%
    filter(interaction_pvalue <0.01) %>%
    slice(1:200) %>% 
    select(Gene, Method, view)
  LIANA_Sidrah <- bind_rows(Kidney_Disease, LIANA_Sidrah)
}

Crosstalk <- data.frame()
for(Disease in Diseases){
  print(Disease)
  Kidney_Disease <- read.csv(paste0("./CrossTalkR/LR_Report_", Disease, ".csv"))
  Kidney_Disease$view = Disease
  Kidney_Disease$Method = "CrosstalkR"
  
  Kidney_Disease <- Kidney_Disease %>%
    pivot_longer(
      cols = c(gene_A, gene_B),
      names_to = "CT",
      values_to = "Gene"
    ) %>%
    mutate(Gene = sub("\\|.*", "", Gene))
  
  Kidney_Disease <- Kidney_Disease %>% 
    arrange(desc(LRScore)) %>% 
    slice(1:200) %>% 
    select(Gene, Method, view)
  LIANA_Sidrah <- bind_rows(Kidney_Disease, LIANA_Sidrah)
}

Marker_Turgay <- data.frame(Gene = c("SLC5A12", "SLC13A3", "CASR", "UMOD", "ESRRB", "SCNN1G", "SLC12A3", "KLHL3", "FLT1", "TEK", "PECAM1", "MEIS1", "PTPRC", "CLNK", "SLC26A7", "ATP6V1C2", "SLC44A5", "IGFBP5", "PLA2R1", "CFH", "SLC4A9", "PRKX", "IGKC", "CD38", "CDH16", "MYH1"))
Marker_Turgay$Method = "Marker_Turgay"
Marker_Turgay$view ="All"

celltype_disease <- read.csv("top_10_unique_genes_by_celltype_disease.csv")
celltype_disease <- celltype_disease %>% 
  rename("Gene" = gene) %>%
  mutate(Method = "celltype_disease", view = sub("_.*", "", disease)) %>%
  select(Gene, Method, view)

Panel <- rbind(SpaPros_Panel, MOFA_F2, LIANA_LR, LIANA_Sidrah, Crosstalk, Marker_Turgay, celltype_disease) #, DEG, DEG_2)
P <- Panel %>% distinct(Gene, .keep_all = TRUE) %>% filter(Gene != "XIST") %>% count(name = "Genes", Method, view)

color_palette <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(P$view)))
ggplot(P, aes(x=view, y=Genes, fill=view))+geom_bar(position="stack", stat="identity")+
  theme(panel.background = element_blank(), axis.text.x = element_text(size=12, colour = "black", angle = 90, vjust = 0.5, hjust = 0.8))+scale_y_continuous(expand = expansion(mult = c(0, .1)))+labs(y="Donor", x="")+facet_wrap(. ~ Method, scales = "free")+
  scale_fill_manual(values = color_palette)+theme(strip.text = element_text(size = 20))+theme(legend.position = "", legend.text=element_text(size=16), legend.title=element_text(size=18), axis.title=element_text(size=14))

Panel %>% distinct(Gene, .keep_all = TRUE) %>% count(name = "Genes", Method)
nrow(Panel %>% distinct(Gene, .keep_all = TRUE))
write.csv(Panel %>% distinct(Gene, .keep_all = TRUE), file = "Xenium_5k_Panel_Kidney_biospy.csv")

Panel <- read.csv("Xenium_480_Panel_Kidney_biospy.csv", row.names = 1)

#Genes excluded by the Xenium Designer
Excluded_genes <- c("ABCA9", "ABI3BP", "ACSM2A", "ACTR2", "ADAM9", "ADGRE2", "ADORA1", "AGRN", "ALDH1A2", "ALDOB", "ALOX5", "ANGPTL5", "ANO9", "APOE", "AQP1", "ARHGAP29", "ARHGAP4", "ASPG",
                    "BCL11B", "BRIP1", "CABLES2", "CACNA1E", "CADM1", "CADM3", "CALB1", "CALM1", "CALM2", "CALR", "CAPN8", "CCN1", "CCN2", "CDH6", "CEP135", "CHRD", "CLCNKB", "CNKSR2", "CNTN3", "CNTN4", "COL1A1", "COL4A6", "COL6A1", "COL8A1", "COMP", "COPA", "CUBN", "CYP21A2",
                    "DCDC2", "DCN", "DDX3Y", "DLGAP1", "DTHD1", "EEF1A1", "EFNA1", "EGFLAM", "EGR1", "EMB", "EMCN", "EPHA2", "ERBB3", "ERBB4", "ETV1", "FCRL5", "FCRL6", "FGF1", "FGF2", "FGFR1", "FKBP5", "FN1", "FOS", "FOSB",
                    "GAS7", "GNAI2", "GPC3", "GPC5", "GPX3", "GREM1", "GRID2", "GRN", "GSDMD", "HAVCR2", "HMGB1", "HSPA8", "ICOS", "IFI44L", "IGF1R", "IGHA1", "IGHG1", "IGHG3", "IGHG4", "IGHGP", "IGHM", "IGLC2", "IGLC3", "IL32", "ITGA4", "ITGA6", "ITGAX", "ITGB2-AS1", "ITGB8", "ITIH4", "JUN",
                    "KCNB2", "KCNIP4", "KCNMA1", "KCNQ3", "KDM5D", "KIF26B-AS1", "KLHL4", "KLK1", "LAMA1", "LAMA3", "LAMA4", "LDB2", "LEF1-AS1", "LEP", "LHCGR", "LPL", "LPP", "LRP2", "LRPAP1", "LRRC3B", "LTB", "LTBP2", "MAML2", "MCOLN2", "MEIS3", "MIR3681HG", "MIR4500HG", "MRC2", "MYLK", "MYO7A",
                    "NAMPT", "NLGN4Y", "NOTCH4", "NPDC1", "NPIPB13", "NPNT", "NRG1", "NRG3", "NRP2", "NTNG1", "ODF3B", "PARD3", "PARP10", "PAX8", "PDGFD", "PET117", "PGGHG", "PKM", "PLAU", "PLEKHS1", "PLG", "PNPLA3", "PODXL", "POTEG", "PPFIBP1", "PRKCB", "PRKY", "PSEN1", "PTH2R", "PTPRA", "PTPRG", "PTPRM", "PTPRQ",
                    "RAG2", "RASSF6", "RBFOX1", "RBFOX3", "RBPMS", "RHEX", "ROBO1", "ROR1", "ROR1-AS1", "RPS4Y1", "RPSA", "RXFP1", "SDC1", "SDC2", "SDC4", "SEMA3E", "SERPING1", "SFRP2", "SHH", "SLC12A1", "SLC17A3", "SLC22A10", "SLC5A12", "SLPI", "SPON2", "SPTAN1", "SULT1C3", "SVIP",
                    "TBC1D10C", "TBC1D3D", "TBC1D8-AS1", "TBL1Y", "TCIRG1", "TDO2", "TEX9", "TGFA", "TGM2", "TMC6", "TMC8", "TMEM132D", "TMEM150B", "TMEM150C", "TNC", "TNFRSF10B", "TNFRSF11B", "TNFRSF1A", "TNFSF10", "TOX2", "TPM1", "TRPM6", "TSPAN12", "TTC29", "TTTY10",
                    "VCAN", "VEGFA", "VEGFC", "VIM", "XAF1", "ZC4H2", "ZCWPW2", "ZFHX4", "ZMYND15", "ZNF534")
Panel$excluded <- ifelse(Panel$Gene %in% Excluded_genes, TRUE, FALSE)
P <- Panel %>% distinct(Gene, .keep_all = TRUE) %>% filter(Gene != "XIST") %>% count(name = "Genes", Method, excluded, view)
ggplot(P, aes(x=view, y=Genes, fill=excluded))+geom_bar(position="stack", stat="identity")+
  theme(panel.background = element_blank(), axis.text.x = element_text(size=12, colour = "black", angle = 90, vjust = 0.5, hjust = 0.8))+scale_y_continuous(expand = expansion(mult = c(0, .1)))+labs(y="Donor", x="")+facet_wrap(. ~ Method, scales = "free")+
  theme(strip.text = element_text(size = 20))+theme(legend.position = "right", legend.text=element_text(size=16), legend.title=element_text(size=18), axis.title=element_text(size=14))
write.csv(Panel %>% distinct(Gene, .keep_all = TRUE), file = "Xenium_Kidney_biospy_Xenium_Design_exclusion.csv")

#Final version of the Xenium Panel
Panel <- read.csv("Xenium_Kidney_biospy_Xenium_Design_V3.csv", row.names = 1)
CDKN1A <- data.frame(Gene = "CDKN1A", Method = "Marker_Turgay", view = "All", excluded = "FALSE")
Panel <- rbind(Panel, CDKN1A)

Xenium_list <- read.csv("./Final_Order/Xenium_Design_Files/GT3ZEN_hKidney_468g_gene_list.csv")
Final_Panel <- Panel %>% filter(Gene %in% Xenium_list$Gene)
P <- Final_Panel %>% distinct(Gene, .keep_all = TRUE) %>% filter(Gene != "XIST") %>% count(name = "Gene", Method, excluded, view)

ggplot(P, aes(x=view, y=Genes, fill=excluded))+geom_bar(position="stack", stat="identity")+
  theme(panel.background = element_blank(), axis.text.x = element_text(size=12, colour = "black", angle = 90, vjust = 0.5, hjust = 0.8))+scale_y_continuous(expand = expansion(mult = c(0, .1)))+labs(y="Donor", x="")+facet_wrap(. ~ Method, scales = "free")+
  theme(strip.text = element_text(size = 20))+theme(legend.position = "right", legend.text=element_text(size=16), legend.title=element_text(size=18), axis.title=element_text(size=14))


Final_Panel$Method[Final_Panel$Method == "LIANA_Sidrah"] <- "LIANA+"
Final_Panel$Method[Final_Panel$Method == "LIANA"] <- "LIANA+"
Final_Panel$Method[Final_Panel$Method == "celltype_disease"] <- "DEG"
Final_Panel$Method[Final_Panel$Method == "Marker_Turgay"] <- "SpaPros"
Final_Panel$Method[Final_Panel$Method == "MOFA"] <- "MFA"
P <- Final_Panel %>% distinct(Gene, .keep_all = TRUE) %>% filter(Gene != "XIST") %>% count(name = "Genes", Method, excluded, view)

ggplot(P, aes(x=Genes, y=Method, fill=Method))+geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = c("lightcoral", "indianred", "firebrick", "darkred")) +
  theme(panel.background = element_blank(), axis.text.x = element_text(size=14, colour = "black"))+labs(y="", x="# Gene with Xenium probes")+
  theme(axis.text.y = element_text(size = 14), axis.title=element_text(size=14, color="black"), legend.position = "")
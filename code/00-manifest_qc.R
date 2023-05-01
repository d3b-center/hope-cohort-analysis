suppressPackageStartupMessages({
  library(tidyverse)
  library(reshape2)
})

# Bo's manifest
bo_manifest = readxl::read_xlsx("data/Project Hope Data Clearup.xlsx")
bo_manifest = bo_manifest[,1:6]

# hope clinical file
proteomics_dat = read_tsv("data/hopeonly_clinical_table_011823.tsv")
bo_manifest$proteomics = bo_manifest$biospecimen_id_harvest %in% proteomics_dat$Sample_ID

# tumor-normal paired
# MAF
maf = readRDS("data/merged_files/snv_merged.rds")
bo_manifest$tumor_normal_consensus_maf <- bo_manifest$T_BS %in% maf$Kids_First_Biospecimen_ID

# CNV
cnv = readRDS("data/merged_files/cnv_merged.rds")
bo_manifest$tumor_normal_cnv <- bo_manifest$T_BS %in% cnv$Kids_First_Biospecimen_ID

# RNA
rna = readRDS("data/merged_files/gene-expression-rsem-tpm-collapsed.rds")
bo_manifest$rna_expression <- bo_manifest$RNA_BS %in% colnames(rna)

# Fusion
fusion = readRDS("data/merged_files/fusions_merged.rds")
bo_manifest$rna_fusion <- bo_manifest$RNA_BS %in% fusion$Kids_First_Biospecimen_ID

# MAF tumor only
maf = readRDS("data/merged_files/snv_merged_tumor_only.rds")
bo_manifest$tumor_only_mutect_maf <- bo_manifest$T_BS %in% maf$Kids_First_Biospecimen_ID

# CNV tumor only
cnv = readRDS("data/merged_files/cnv_merged_tumor_only.rds")
bo_manifest$tumor_only_cnv <- bo_manifest$T_BS %in% cnv$Kids_First_Biospecimen_ID

# CRAM (for MSI analysis)
msi_analysis <- data.frame()

# get cram file name
cram_manifest <- read_tsv("data/manifest/manifest_20230428_070806_cram.tsv")
cram_manifest <- cram_manifest %>% 
  dplyr::select(case_id, `Kids First Biospecimen ID`, sample_id, sample_type, name)
cram_manifest$sample_type[cram_manifest$`Kids First Biospecimen ID` == "BS_0TCRV9AC"] <- "Normal"
cram_manifest1 <- dcast(cram_manifest, sample_id ~ sample_type, value.var = "name")
colnames(cram_manifest1)[2:4] <- c("NonTumor_cram", "N_cram", "T_cram")

# get bs id
cram_manifest2 <- dcast(cram_manifest, sample_id ~ sample_type, value.var = "Kids First Biospecimen ID")
colnames(cram_manifest2)[2:4] <- c("NonTumor_BS", "N_BS", "T_BS")

msi_analysis <- cram_manifest1 %>%
  inner_join(cram_manifest2) %>%
  filter(sample_id %in% bo_manifest$`SDG ID QC`)

# add it to Bo's manifest
bo_manifest <- bo_manifest %>%
  left_join(msi_analysis, by = c("SDG ID QC" = "sample_id", "N_BS", "T_BS")) %>%
  dplyr::select(-c(NonTumor_cram, NonTumor_BS))
write_tsv(bo_manifest, file = "data/Project Hope Data Clearup_add.tsv")

# issues
# BS_0TCRV9AC sample is a Normal (in Bo's clean upfile) but in the Cavatica manifest, it says "Tumor"
# BS_7P378T0E cram is missing from Cavatica (Tumor-only) so no MAF/CNV/MSI
# 7316-3000 has two MSI outputs: 
# 25bb3e6c-40cd-420a-a080-8bf49f1157a0_tumor_only_msisensor_pro was run recently and has no C-id or BS-id associated with it but should be 
# C744765 and 7316-3000
# BS_SGS26NXP multiple runs

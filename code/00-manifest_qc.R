suppressPackageStartupMessages({
  library(tidyverse)
  library(reshape2)
})

# input directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")

# Bo's manifest
bo_manifest = readxl::read_xlsx("data/Project Hope Data Clearup.xlsx")
bo_manifest = bo_manifest[,1:6]

# add cram file name
lf <- list.files(file.path(data_dir, "manifest"), full.names = T)
lf <- lf[grep("cram", lf)]
cram_manifest <- read_tsv(lf)
cram_manifest <- cram_manifest %>% 
  dplyr::select(case_id, `Kids First Biospecimen ID`, sample_id, sample_type, name)
cram_manifest$sample_type[cram_manifest$`Kids First Biospecimen ID` == "BS_0TCRV9AC"] <- "Normal"
cram_manifest_tumor <- dcast(cram_manifest %>% filter(sample_type == "Tumor"), `Kids First Biospecimen ID` ~ sample_type, value.var = "name")
colnames(cram_manifest_tumor)[2] <- "T_cram"
cram_manifest_normal <- dcast(cram_manifest %>% filter(sample_type == "Normal"), `Kids First Biospecimen ID` ~ sample_type, value.var = "name")
colnames(cram_manifest_normal)[2] <- "N_cram"

# add it to Bo's manifest
bo_manifest <- bo_manifest %>%
  left_join(cram_manifest_tumor, by = c("T_BS" = "Kids First Biospecimen ID")) %>%
  left_join(cram_manifest_normal, by = c("N_BS" = "Kids First Biospecimen ID"))
bo_manifest$N_cram <- gsub("harmonized-data/aligned-reads/", "", bo_manifest$N_cram)
bo_manifest$T_cram <- gsub("harmonized-data/aligned-reads/", "", bo_manifest$T_cram)

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

# MSI
msi = readRDS("data/merged_files/msi_merged.rds")
bo_manifest$tumor_normal_msi <- bo_manifest$biospecimen_id_harvest %in% msi$sample_id

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

# MSI tumor only
msi = readRDS("data/merged_files/msi_merged_tumor_only.rds")
bo_manifest$tumor_only_msi <- bo_manifest$biospecimen_id_harvest %in% msi$sample_id

# save additions to Bo's manifest
write_tsv(bo_manifest, file = "data/Project Hope Data Clearup_add.tsv")

# tumor-only check
bo_manifest %>% 
  filter(!is.na(T_BS)) %>%
  dplyr::select(tumor_only_mutect_maf, tumor_only_cnv, tumor_only_msi) %>%
  unique()

# tumor-normal check
bo_manifest %>% 
  filter(!is.na(N_BS), !is.na(T_BS)) %>%
  dplyr::select(tumor_normal_consensus_maf, tumor_normal_cnv, tumor_normal_msi) %>%
  unique()

# BS_SNWRBH0J missing RNA expression and fusion
bo_manifest %>%
  filter(!is.na(RNA_BS)) %>%
  dplyr::select(rna_expression, rna_fusion) %>%
  unique()

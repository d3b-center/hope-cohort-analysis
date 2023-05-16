# Script to create master histology for HOPE cohort
suppressPackageStartupMessages({
  library(tidyverse)
})

# input directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")

# clinical info
annot <- readr::read_tsv(file.path(data_dir, "hopeonly_clinical_table_011823.tsv"))

# add proteomics cluster 
proteomics <- read_tsv(file.path(data_dir, "cluster_data_090722.tsv"))
proteomics <- proteomics %>%
  dplyr::select(id, rdt.cc) 
dev_clusters <- read_tsv(file.path(data_dir, "cluster_data101922.tsv"))
dev_clusters <- dev_clusters %>%
  dplyr::select(id, dtt.cc, age, rdt.name) %>%
  dplyr::rename("age_three_groups" = "age") # three age groups
cluster_meta <- proteomics %>%
  inner_join(dev_clusters)
annot <- annot %>%
  left_join(cluster_meta, by = c("Sample_ID" = "id")) 

# add two age groups
age_info <- readxl::read_xlsx(file.path(data_dir, "clini_m_030722-for_Komal.xlsx")) %>%
  dplyr::rename("age_two_groups" = "age.class")
annot <- annot %>%
  left_join(age_info %>% dplyr::select(age_two_groups, id), by = c("Sample_ID" = "id"))

# msi sensor (paired)
msi_paired <- readRDS(file.path(data_dir, "merged_files", "msi_merged.rds"))
msi_paired <- msi_paired %>%
  dplyr::rename("msi_paired" = "Percent") %>%
  dplyr::select(sample_id, msi_paired)

# msi sensor (tumor-only)
msi_tumor_only <- readRDS(file.path(data_dir, "merged_files", "msi_merged_tumor_only.rds"))
msi_tumor_only <- msi_tumor_only %>%
  dplyr::rename("msi_tumor_only" = "Percent") %>%
  dplyr::select(sample_id, msi_tumor_only)

# add both to annotation
annot <- annot %>%
  left_join(msi_paired, by = c("Sample_ID" = "sample_id")) %>%
  left_join(msi_tumor_only, by = c("Sample_ID" = "sample_id"))

# add ALT status
alt_status <- read_tsv(file.path(data_dir, "hope_cohort_alt_status.txt")) 
alt_status <- alt_status %>%
  filter(sample_id %in% annot$Sample_ID) %>%
  dplyr::select(sample_id, t_n_telomere_content, ALT_status) # (n = 68)
annot <- annot %>%
  left_join(alt_status, by = c("Sample_ID" = "sample_id"))

# add TMB values (paired)
tmb_output <- read_tsv('code/tmb-calculation/results/snv-mutation-tmb-coding.tsv')
tmb_output <- tmb_output %>%
  dplyr::rename("TMB_paired" = "tmb") %>%
  dplyr::select(Tumor_Sample_Barcode, TMB_paired)

# add sample_id to TMB
lf <- list.files(file.path(data_dir, "manifest"), full.names = T)
lf <- lf[grep("snv", lf)]
hist_df = read_tsv(lf)
hist_df = hist_df %>%
  filter(sample_id %in% annot$Sample_ID)
tmb_output = hist_df %>% 
  dplyr::select(`Kids First Biospecimen ID`, sample_id) %>%
  inner_join(tmb_output, by = c("Kids First Biospecimen ID" = "Tumor_Sample_Barcode"))
annot <- annot %>%
  left_join(tmb_output %>% dplyr::select(-c(`Kids First Biospecimen ID`)), by = c("Sample_ID" = "sample_id"))

# add TMB values (tumor-only)
tmb_output <- read_tsv('code/tmb-calculation/results/snv-tumor-only-mutation-tmb-coding.tsv')
tmb_output <- tmb_output %>%
  dplyr::rename("TMB_tumor_only" = "tmb") %>%
  dplyr::select(Tumor_Sample_Barcode, TMB_tumor_only)

# add sample_id to TMB
lf <- list.files(file.path(data_dir, "manifest"), full.names = T)
lf <- lf[grep("snv", lf)]
hist_df = read_tsv(lf)
hist_df = hist_df %>%
  filter(sample_id %in% annot$Sample_ID)
tmb_output = hist_df %>% 
  dplyr::select(`Kids First Biospecimen ID`, sample_id) %>%
  inner_join(tmb_output, by = c("Kids First Biospecimen ID" = "Tumor_Sample_Barcode"))
annot <- annot %>%
  left_join(tmb_output %>% dplyr::select(-c(`Kids First Biospecimen ID`)), by = c("Sample_ID" = "sample_id"))
write_tsv(annot, file = file.path(data_dir, "master_histology_hope_cohort.tsv"))

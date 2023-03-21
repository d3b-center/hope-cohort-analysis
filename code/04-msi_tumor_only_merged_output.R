# msisensor pro summary output file for Mateusz
suppressPackageStartupMessages({
  library(tidyverse)
})

# output directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")

# hope cohort manifest
hope_manifest <- read_tsv('data/hopeonly_clinical_table_011823.tsv')

# MSI output (n = 88)
fname <- 'results/msisensor-pro-tumor-only/hope_cohort_msi_sensor_output.tsv'
msi_output <- read_tsv(fname)
msi_output <- msi_output %>%
  filter(sample_id %in% hope_manifest$Sample_ID) %>%
  dplyr::rename("MSI_Percent" = "Percent")

# TMB output
tmb_output <- read_tsv('code/tmb-calculation/results/snv-tumor-only-mutation-tmb-coding.tsv')
tmb_output <- tmb_output %>%
  dplyr::rename("TMB" = "tmb") %>%
  dplyr::select(Tumor_Sample_Barcode, TMB)

# add sample_id to TMB
hist_df = read_tsv(file = "data/manifest_tumor_only/manifest_20230216_162009_snv.tsv")
hist_df = hist_df %>%
  filter(sample_id %in% hope_manifest$Sample_ID)
tmb_output = hist_df %>% 
  dplyr::select(`Kids First Biospecimen ID`, sample_id) %>%
  inner_join(tmb_output, by = c("Kids First Biospecimen ID" = "Tumor_Sample_Barcode")) %>%
  dplyr::arrange(desc(TMB)) %>%
  distinct(sample_id, .keep_all = T)
  
msi_output <- msi_output %>%
  left_join(tmb_output, by = "sample_id")

# add ALT status from Mateusz
pos_controls <- read_tsv('data/hope_cohort_alt_status.txt')
pos_controls <- pos_controls %>%
  dplyr::select(sample_id, ALT_status) 
msi_output <- msi_output %>%
  left_join(pos_controls, by = "sample_id")

# add proteomics cluster 
proteomics <- read_tsv('../d3b-miRNA-analysis/analyses/actmir-analysis/input/cluster_data_090722.tsv')
proteomics <- proteomics %>%
  dplyr::select(id, rdt.cc) 
dev_clusters <- read_tsv('data/cluster_data101922.tsv')
dev_clusters <- dev_clusters %>%
  dplyr::select(id, dtt.cc, age, rdt.name) %>%
  dplyr::rename("age_three_groups" = "age") # three age groups
cluster_meta <- proteomics %>%
  inner_join(dev_clusters)
msi_output <- msi_output %>%
  left_join(cluster_meta, by = c("sample_id" = "id")) 

# add two age groups
age_info <- readxl::read_xlsx(file.path(data_dir, "clini_m_030722-for_Komal.xlsx"))
msi_output <- msi_output %>%
  left_join(age_info %>% dplyr::select(age.class, id), by = c("sample_id" = "id")) %>%
  dplyr::rename("age_two_groups" = "age.class")
write_tsv(msi_output, file = "results/msisensor-pro-tumor-only/msi_output_merged_tumor_only.tsv")

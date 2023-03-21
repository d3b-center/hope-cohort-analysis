# script to summarise msi-sensor pro output
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
})

# output directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
input_dir <- file.path(data_dir, "msisensor_pro_tumor_only")
results_dir <- file.path(root_dir, "results", "msisensor-pro-tumor-only")
dir.create(results_dir, showWarnings = F, recursive = T)

# read msisensor pro files
lf <- list.files(path = input_dir, pattern = 'msisensor_pro', recursive = T, full.names = T)
df <- lapply(lf, read_tsv)
df <- do.call(rbind,df)
colnames(df)[3] <- "Percent"
df$name <- gsub('.*/', '', lf)

# read manifest
manifest <- read_tsv(file.path(data_dir, "manifest", "manifest_20230321_111302_msi.tsv"))
manifest <- manifest %>%
  dplyr::select(sample_id, name)

# add Gender info from updated clinical data from Mateusz (n = 88)
annot <- readr::read_tsv(file.path(data_dir, "hopeonly_clinical_table_011823.tsv"))
annot <- annot %>% 
  dplyr::select(Sample_ID, Gender) %>%
  dplyr::rename("gender" = "Gender")

# combine 
manifest <- manifest %>%
  inner_join(annot, by = c("sample_id" = "Sample_ID"))
manifest <- manifest %>% 
  inner_join(df, by = "name") %>%
  dplyr::arrange(desc(Percent)) %>%
  distinct(sample_id, .keep_all = T) %>%
  dplyr::select(sample_id, gender, Percent)

# modify
output_df <- manifest %>%
  mutate(Type = ifelse(Percent >= 3.5, "High", "Low")) %>%
  arrange(Percent)
fname <- file.path(results_dir, "hope_cohort_msi_sensor_output.tsv")
write_tsv(output_df, file = fname)

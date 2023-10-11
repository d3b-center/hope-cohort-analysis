# script to correlate ALT status to clinical variables like
# Protein clusters, Age (2 and 3 groups), Gender and TMB
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(ggforce)
})

# output directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
# input_dir <- file.path(root_dir, "analyses", "master-annotation", "results")
analyses_dir <- file.path(root_dir, "analyses", "alt-analysis")
output_dir <- file.path(analyses_dir, "results")
dir.create(output_dir, showWarnings = F, recursive = T)

# master histology file
hist_df <- read_tsv(file = file.path(root_dir, "data", "v1", "Hope-GBM-histologies.tsv"))
alt_status <- read_tsv(file.path("../alt-analysis/results/alt_status_aya_hgg.tsv"))
tmp <- hist_df %>%
  filter(!is.na(molecular_subtype)) %>%
  inner_join(alt_status) %>%
  dplyr::select(sample_id, experimental_strategy, ALT_status, t_n_telomere_content, HARMONY_age_class_derived) %>%
  group_by(sample_id) %>%
  mutate(experimental_strategy = toString(experimental_strategy)) %>%
  unique()
tmp <- tmp[grep("WGS", tmp$experimental_strategy),]

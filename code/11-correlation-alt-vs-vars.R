# script to correlation ALT status to clinical variables
# protein clusters, age, sex and MSI and TMB
suppressPackageStartupMessages({
  library(tidyverse)
})

# output directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")

# function to compute correlations
compute_corr <- function(){
  
}

# 1) tumor-normal-paired
# read combined file
alt_status <- read_tsv(file.path("results/msisensor-pro/msi_output_merged.tsv"))

# 2) tumor-only 
# read combined file
alt_status_tumor_only <- read_tsv(file.path("results/msisensor-pro-tumor-only/msi_output_merged_tumor_only.tsv"))



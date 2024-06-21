# Function: script to compare TMB values from overlapping samples between OpenPedCan and HOPE cohort

# load libraries
suppressPackageStartupMessages({
  library(tidyverse)
})

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "analyses", "tmb-calculation")
results_dir <- file.path(module_dir, "results")

# TMB from HOPE (T/N paired)
tmb_hope <- read_tsv(file.path(results_dir, "wgs_paired", "snv-mutation-tmb-coding.tsv"))

# TMB from OpenPedCan
tmb_openpedcan <- read_tsv('~/Projects/OpenPedCan-analysis/analyses/tmb-calculation/results/snv-mutation-tmb-coding.tsv')
tmb_openpedcan <- tmb_openpedcan %>%
  dplyr::select(Tumor_Sample_Barcode, tmb) %>%
  dplyr::rename("tmb_openpedcan" = "tmb")
tmb_hope <- tmb_hope %>%
  dplyr::select(Tumor_Sample_Barcode, tmb) %>%
  dplyr::rename("tmb_hope" = "tmb")
dat <- tmb_openpedcan %>%
  inner_join(tmb_hope)

# correlation between OpenPedCan and HOPE cohort TMB values
print(cor(dat$tmb_openpedcan, dat$tmb_hope))
# 1

# correlation between HOPE cohort T/N paired and tumor-only values
tmb_hope_tumor_only = read_tsv(file.path(results_dir, "wgs_tumor_only", "snv-mutation-tmb-coding.tsv"))
tmb_hope_tumor_only = tmb_hope_tumor_only %>%
  dplyr::select(Tumor_Sample_Barcode, tmb) %>%
  dplyr::rename("tmb_hope_tumor_only" = "tmb")
dat = dat %>%
  inner_join(tmb_hope_tumor_only)
print(cor(dat$tmb_hope, dat$tmb_hope_tumor_only))
# 0.9203342

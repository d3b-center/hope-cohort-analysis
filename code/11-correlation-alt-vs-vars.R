# script to correlation ALT status to clinical variables
# protein clusters, age, sex and MSI and TMB
suppressPackageStartupMessages({
  library(tidyverse)
})

# output directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")

# function to compute correlations
compute_corr <- function(x){

  # for ALT vs Age (two groups)
  age_two_group_pvalue <- chisq.test(x = factor(x$ALT_status), y = factor(x$age_two_groups))$p.value
  age_two_group_pvalue <- round(age_two_group_pvalue, digits = 3)
  # age_pvalue <- kruskal.test(formula = factor(alt_status$ALT_status) ~ factor(alt_status$age))$p.value
  
  # for ALT vs Age (three groups)
  age_three_group_pvalue <- chisq.test(x = factor(x$ALT_status), y = factor(x$age_three_groups))$p.value
  age_three_group_pvalue <- round(age_three_group_pvalue, digits = 3)
  # age_pvalue <- kruskal.test(formula = factor(alt_status$ALT_status) ~ factor(alt_status$age))$p.value
  
  # Sex
  gender_pvalue <- chisq.test(x = factor(x$ALT_status), y = factor(x$gender))$p.value
  gender_pvalue <- round(gender_pvalue, digits = 3)
  # gender_pvalue <- kruskal.test(formula = factor(alt_status$ALT_status) ~ factor(alt_status$gender))$p.value
  
  # Protein cluster
  rdt_cc_pvalue <- chisq.test(x = factor(x$ALT_status), y = factor(x$rdt.cc))$p.value
  rdt_cc_pvalue <- round(rdt_cc_pvalue, digits = 3)
  # rdt_cc_pvalue <- kruskal.test(formula = factor(alt_status$ALT_status) ~ factor(alt_status$rdt.cc))$p.value
  
  # Protein cluster name
  rdt_name_pvalue <- chisq.test(x = factor(x$ALT_status), y = factor(x$rdt.name))$p.value
  rdt_name_pvalue <- round(rdt_name_pvalue, digits = 3)
  # rdt_name_pvalue <- kruskal.test(formula = factor(alt_status$ALT_status) ~ factor(alt_status$rdt.name))$p.value
  
  # MSI
  msi_pvalue <- chisq.test(x = factor(x$ALT_status), y = x$MSI_Percent)$p.value
  msi_pvalue <- round(msi_pvalue, digits = 3)
  # msi_pvalue <- kruskal.test(formula = factor(alt_status$ALT_status) ~ alt_status$MSI_Percent)$p.value
  
  # TMB
  tmb_pvalue <- chisq.test(x = factor(x$ALT_status), y = x$TMB)$p.value
  tmb_pvalue <- round(tmb_pvalue, digits = 3)
  # tmb_pvalue <- kruskal.test(formula = factor(alt_status$ALT_status) ~ alt_status$TMB)$p.value

  output_df <- data.frame(age_two_group_pvalue, age_three_group_pvalue, 
                          gender_pvalue, 
                          rdt_cc_pvalue, rdt_name_pvalue,
                          msi_pvalue, tmb_pvalue)
  return(output_df)
}

# 1) tumor-normal-paired
# read combined file
alt_status <- read_tsv(file.path("results/msisensor-pro/msi_output_merged.tsv"))
df <- compute_corr(x = alt_status)
rownames(df) <- "Tumor-Normal"

# 2) tumor-only 
# read combined file
alt_status_tumor_only <- read_tsv(file.path("results/msisensor-pro-tumor-only/msi_output_merged_tumor_only.tsv"))
df_tumor_only <- compute_corr(x = alt_status_tumor_only)
rownames(df_tumor_only) <- "Tumor-only"

# write output
rbind(df, df_tumor_only) %>%
  rownames_to_column("Analysis_Type") %>%
  write_tsv(file = "results/alt_correlations.tsv")


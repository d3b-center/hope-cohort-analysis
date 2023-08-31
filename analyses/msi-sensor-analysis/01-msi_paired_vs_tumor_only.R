# script to compare paired and tumor-only msi-sensor output
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
})

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
input_dir <- file.path(root_dir, "data")
analyses_dir <- file.path(root_dir, "analyses", "msi-sensor-analysis")
output_dir <- file.path(analyses_dir, "results", "msisensor-pro-combined")
dir.create(output_dir, showWarnings = F, recursive = T)

# master histology
annot <- read_tsv(file = file.path(input_dir, "master_histology_hope_cohort.tsv"))

# combined barplot of common sample ids
annot <- annot %>%
  dplyr::select(Sample_ID, msi_paired, msi_tumor_only) %>%
  filter(!is.na(msi_tumor_only), !is.na(msi_paired))

# scatter plot
pdf(file = file.path(output_dir, "tumor_only_vs_paired_analysis.pdf"), width = 8, height = 6)
p <- ggplot(annot, aes(x = msi_paired, y = msi_tumor_only)) +
  geom_point(pch = 21, size = 4) +
  ggpubr::theme_pubr(base_size = 10, legend = "bottom") + 
  ggtitle(paste0("Paired vs Tumor-only analysis (n = ", nrow(annot) , ")")) +
  stat_cor(method = "pearson") + ylab("MSI (tumor-only)") + xlab("MSI (paired)")
print(p)
dev.off()

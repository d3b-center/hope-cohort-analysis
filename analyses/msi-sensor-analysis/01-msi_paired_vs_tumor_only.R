# Function: compare paired and tumor-only msi-sensor outputs

# load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
})

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analyses_dir <- file.path(root_dir, "analyses", "msi-sensor-analysis")
input_dir <- file.path(analyses_dir, "input")
output_dir <- file.path(analyses_dir, "results")
dir.create(output_dir, showWarnings = F, recursive = T)

# master histology
msi_paired <- file.path(output_dir, "msisensor-pro-paired", "Hope-msi-paired.tsv") %>% 
  read_tsv() %>%
  dplyr::rename("MSI_paired" = "Percent") %>%
  dplyr::select(sample_id, MSI_paired)
msi_tumor_only <- file.path(output_dir, "msisensor-pro-tumor-only", "Hope-msi-tumor_only.tsv") %>% 
  read_tsv() %>%
  dplyr::rename("MSI_tumor_only" = "Percent") %>%
  dplyr::select(sample_id, MSI_tumor_only)

# combined barplot of common sample ids
annot <- msi_paired %>%
  inner_join(msi_tumor_only, by = "sample_id")

# scatter plot
pdf(file = file.path(output_dir, "msisensor-pro-combined", "tumor_only_vs_paired_analysis.pdf"), width = 8, height = 6)
p <- ggplot(annot, aes(x = MSI_paired, y = MSI_tumor_only)) +
  geom_point(pch = 21, size = 4) +
  ggpubr::theme_pubr(base_size = 10, legend = "bottom") + 
  ggtitle(paste0("Paired vs Tumor-only analysis (n = ", nrow(annot) , ")")) +
  stat_cor(method = "pearson") + ylab("MSI (tumor-only)") + xlab("MSI (paired)")
print(p)
dev.off()

# script to summarise msi-sensor pro output
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
})

# Paired (n = 69)
paired_results <- read_tsv('results/msisensor-pro/hope_cohort_msi_sensor_output.tsv')
paired_results$analysis_type <- "paired"

# Tumor-only (n = 88)
tumor_only_results <- read_tsv('results/msisensor-pro-tumor-only/hope_cohort_msi_sensor_output.tsv')
tumor_only_results$analysis_type <- "tumor_only"
combined_df <- rbind(paired_results, tumor_only_results)

# combined barplot of common sample ids
unique_ids <- combined_df[which(duplicated(combined_df$sample_id)),"sample_id"] %>% 
  unique() %>%
  pull(sample_id)
combined_df <- combined_df %>% 
  filter(sample_id %in% unique_ids)
pdf(file = "results/msisensor-pro-combined/tumor_only_vs_paired_analysis.pdf", width = 12, height = 6)
ggplot(combined_df, 
       aes(x = sample_id, y = Percent, fill = analysis_type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  ggpubr::theme_pubr(base_size = 10, legend = "bottom") + 
  ggtitle(paste0("Comparison of samples processed using both Paired analysis and Tumor-only analysis (n = ", length(unique_ids), ")")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()

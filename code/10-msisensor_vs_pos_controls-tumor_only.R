# compare somatic MMR status to paper: https://www.biorxiv.org/content/10.1101/2022.08.05.502870v1.supplementary-material
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
})

# msisensor pro
fname <- 'results/msisensor-pro-tumor-only/hope_cohort_msi_sensor_output.tsv'
output_df <- read_tsv(fname)
output_df <- output_df %>%
  dplyr::select(-c(Type))

# positive controls
pos_controls <- readxl::read_xlsx('data/media-4.xlsx')
pos_controls <- pos_controls %>%
  dplyr::select(`SAMPLE ID`, MMR_SOMATIC, MMR_GERMLINE)
output_df <- output_df %>%
  inner_join(pos_controls, by = c("sample_id" = "SAMPLE ID"))
write_tsv(output_df, file = 'results/msisensor-pro-tumor-only/msisensorpro_vs_pos_controls.tsv')

# plot
output_df$label <- ifelse(output_df$MMR_GERMLINE != "none", output_df$sample_id, "")
ggplot(output_df, aes(x = "", Percent, color = MMR_SOMATIC, label = label)) + 
  geom_point(pch = 21, size = 4) + 
  geom_text(hjust = -0.2, vjust = 0.5) +
  theme_pubr(legend = "right") + 
  ylab("% MSI") + 
  xlab(paste0("Samples: n = ",nrow(output_df))) +
  ggtitle("% MSI vs MMR_Somatic")
ggsave(filename = "results/msisensor-pro-tumor-only/msisensorpro_vs_pos_controls.png", height = 6, width = 6)

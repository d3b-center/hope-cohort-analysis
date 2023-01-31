# compare somatic MMR status to paper: https://www.biorxiv.org/content/10.1101/2022.08.05.502870v1.supplementary-material
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
})

# msisensor pro
fname <- 'results/msisensor-pro-tumor-only/hope_cohort_msi_sensor_tumor_only_output.tsv'
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
ggplot(output_df, aes(x = "", Percent, color = MMR_SOMATIC)) + 
  geom_point(pch = 21, size = 4) + 
  theme_pubr(legend = "right") + ylab("% MSI") + xlab("") +
  ggtitle("% MSI vs MMR_Somatic") 
ggsave(filename = "results/msisensor-pro-tumor-only/msisensorpro_vs_pos_controls.png", height = 6, width = 6)

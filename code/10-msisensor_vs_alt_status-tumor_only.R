# compare somatic MMR status to paper: https://www.biorxiv.org/content/10.1101/2022.08.05.502870v1.supplementary-material
# issue that discusses this: https://github.com/d3b-center/bixu-tracker/issues/1534
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(ggpubr)
})

# msisensor pro
fname <- 'results/msisensor-pro-tumor-only/msi_output_merged_tumor_only.tsv'
output_df <- read_tsv(fname)

# boxplot
output_df <- output_df %>%
  filter(!is.na(ALT_status))
ggplot(output_df, aes(x = factor(ALT_status), y = MSI_Percent, color = ALT_status)) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.5, outlier.size = 1) +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", size = 4) +
  xlab("") + 
  ylab("% Microsatellite Instability") +
  ggtitle("% Microsatellite Instability vs. ALT Status") +
  geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
  theme(legend.position = "none") 
ggsave(filename = "results/msisensor-pro-tumor-only/msisensorpro_vs_alt_status.png", height = 6, width = 6)

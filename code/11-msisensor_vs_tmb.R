# msi sensor pro vs tmb scores 
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(ggpubr)
})

# msisensor pro
fname <- 'results/msisensor-pro/hope_cohort_msi_sensor_output.tsv'
output_df <- read_tsv(fname)
output_df <- output_df %>%
  dplyr::select(-c(Type))

# get TMB from OT
dat <- read_tsv('~/Projects/OpenPedCan-analysis/analyses/tmb-calculation/results/snv-mutation-tmb-coding.tsv')
output_df <- output_df %>%
  inner_join(dat, by = c("Kids First Biospecimen ID" = "Tumor_Sample_Barcode"))

# do correlation
ggplot(output_df, aes(Percent, tmb)) + 
  geom_point() + theme_pubr() + xlab("% MSI") + ylab("TMB") + ggtitle("% MSI vs TMB") +
  geom_smooth() + 
  stat_cor(method = "pearson", label.x = max(output_df$Percent)-2, label.y = max(output_df$tmb))
ggsave(filename = "results/msisensor-pro/msisensorpro_vs_tmb.png", height = 6, width = 6)

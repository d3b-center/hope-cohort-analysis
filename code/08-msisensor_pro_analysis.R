# script to summarise msi-sensor pro output
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
})

# read msisensor pro files
lf <- list.files(path = 'data/msisensor_pro/', pattern = 'msisensor_pro', recursive = T, full.names = T)
df <- lapply(lf, read_tsv)
df <- do.call(rbind,df)
colnames(df)[3] <- "Percent"
df$name <- gsub('.*/', '', lf)
manifest <- read_tsv('data/manifest_20221014_120538.tsv')
tumor_normal_pair <- read_tsv('data/Tumour_normal_participate.tsv')
output_df <- manifest %>%
  dplyr::select(name, sample_id, `Kids First Biospecimen ID`, sample_type) %>%
  mutate(name = gsub("_somatic", "", name)) %>%
  inner_join(tumor_normal_pair, by = c("Kids First Biospecimen ID" = "Tumour")) %>%
  inner_join(df, by = c("name")) %>%
  dplyr::select(-c(name, Normal, participate)) %>%
  mutate(Type = ifelse(Percent >= 3.5, "High", "Low")) %>%
  arrange(Percent)
write_tsv(output_df, file = 'results/hope_cohort_msi_sensor_output.tsv')

# read clinical data from proteomics
proteomics <- read_tsv('../d3b-miRNA-analysis/analyses/actmir-analysis/input/cluster_data_090722.tsv')
output_df <- output_df %>%
  inner_join(proteomics, by = c("sample_id" = "id")) %>%
  filter(!is.na(rdt.cc))

# plot boxplot of proteomics clusters vs clusters
output_df$rdt.cc <- as.character(output_df$rdt.cc)
output_df <- output_df %>%
  group_by(rdt.cc) %>%
  mutate(n = n()) %>%
  mutate(rdt.cc = paste0(rdt.cc, "\n(n = ",n,")")) 
p <- ggplot(output_df, aes(x = as.character(rdt.cc), y = Percent, color = as.character(rdt.cc))) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.7, outlier.size = 1) +
  ggpubr::theme_pubr() + ylab("") + 
  stat_compare_means(method = "wilcox.test", color = "red", 
                     comparisons = list(c("1\n(n = 20)", "2\n(n = 33)"), c("2\n(n = 33)", "3\n(n = 7)"), c("1\n(n = 20)", "3\n(n = 7)"))) +
  xlab("Cluster") + 
  ylab("% Microsatellite Instability") +
  geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
  theme(legend.position = "none") 
ggsave(filename = 'results/msi_sensor_vs_proteomics_clusters.pdf', width = 4, height = 4)

# plot histogram
ggplot(df, aes(Percent, y = "")) +
  geom_jitter(pch = 21, size = 4) +
  ggpubr::theme_pubr() + ylab("") + xlab("% Somatic Sites") + 
  ggtitle("Distribution of % Somatic Sites") +
  scale_x_continuous(breaks=seq(0, 6, 1)) +
  coord_flip()

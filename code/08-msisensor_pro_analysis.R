# script to summarise msi-sensor pro output
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
})

# read msisensor pro files
fname <- 'results/msisensor-pro/hope_cohort_msi_sensor_output.tsv'
if(!file.exists(fname)){
  lf <- list.files(path = 'data/msisensor_pro/', pattern = 'msisensor_pro', recursive = T, full.names = T)
  df <- lapply(lf, read_tsv)
  df <- do.call(rbind,df)
  colnames(df)[3] <- "Percent"
  df$name <- gsub('.*/', '', lf)
  manifest <- read_tsv('data/manifest_20221014_120538.tsv')
  tumor_normal_pair <- read_tsv('data/Tumour_normal_participate.tsv')
  output_df <- manifest %>%
    dplyr::select(name, sample_id, `Kids First Biospecimen ID`, gender, sample_type) %>%
    mutate(name = gsub("_somatic", "", name)) %>%
    inner_join(tumor_normal_pair, by = c("Kids First Biospecimen ID" = "Tumour")) %>%
    inner_join(df, by = c("name")) %>%
    dplyr::select(-c(name, Normal, participate)) %>%
    mutate(Type = ifelse(Percent >= 3.5, "High", "Low")) %>%
    arrange(Percent)
  write_tsv(output_df, file = fname)
} else {
  output_df <- read_tsv(fname)
}

output_df <- output_df %>%
  filter(Percent < 1)

# read clinical data from proteomics
proteomics <- read_tsv('../d3b-miRNA-analysis/analyses/actmir-analysis/input/cluster_data_090722.tsv')
proteomics <- proteomics %>%
  dplyr::select(id, rdt.cc) %>%
  dplyr::rename("proteomics_rdt_cc" = "rdt.cc")
dev_clusters <- read_tsv('data/cluster_data101922.tsv')
dev_clusters <- dev_clusters %>%
  dplyr::select(id, dtt.cc, age, rdt.name) 
cluster_meta <- proteomics %>%
  inner_join(dev_clusters)
output_df <- output_df %>%
  left_join(cluster_meta, by = c("sample_id" = "id")) 

# proteomics clusters
plot_data <- output_df %>%
  filter(!is.na(proteomics_rdt_cc)) %>%
  mutate(proteomics_rdt_cc = as.character(proteomics_rdt_cc)) %>%
  group_by(proteomics_rdt_cc) %>%
  mutate(n = n()) %>%
  mutate(proteomics_rdt_cc = paste0(proteomics_rdt_cc, "\n(n = ",n,")")) 
p <- ggplot(plot_data, aes(x = as.character(proteomics_rdt_cc), y = Percent, color = as.character(proteomics_rdt_cc))) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.5, outlier.size = 1) +
  ggpubr::theme_pubr(base_size = 8) + ylab("") + 
  stat_compare_means(color = "red", 
                     comparisons = combn(unique(plot_data$proteomics_rdt_cc), m = 2, simplify = F)) +
  stat_compare_means(color = "blue") +
  xlab("") + 
  ylab("% Microsatellite Instability") +
  ggtitle("% Microsatellite Instability vs. Proteomics Cluster") +
  # geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
  theme(legend.position = "none") 
ggsave(filename = 'results/msisensor-pro/msi_sensor_vs_proteomics_clusters.png', width = 6, height = 6)

# developmental clusters
plot_data <- output_df %>%
  filter(!is.na(dtt.cc)) %>%
  mutate(dtt.cc = as.character(dtt.cc)) %>%
  group_by(dtt.cc) %>%
  mutate(n = n()) %>%
  mutate(dtt.cc = paste0(dtt.cc, "\n(n = ",n,")")) 
q <- ggplot(plot_data, aes(x = as.character(dtt.cc), y = Percent, color = as.character(dtt.cc))) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.5, outlier.size = 1) +
  ggpubr::theme_pubr(base_size = 8) + ylab("") + 
  stat_compare_means(color = "red", 
                     comparisons = combn(unique(plot_data$dtt.cc), m = 2, simplify = F)) +
  stat_compare_means(color = "blue") +
  xlab("") + 
  ylab("% Microsatellite Instability") +
  ggtitle("% Microsatellite Instability vs. Developmental Cluster") +
  # geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
  theme(legend.position = "none") 
ggsave(filename = 'results/msisensor-pro/msi_sensor_vs_dev_clusters.png', width = 6, height = 6)

# age
plot_data <- output_df %>%
  filter(!is.na(age)) %>%
  mutate(age = as.character(age)) %>%
  group_by(age) %>%
  mutate(n = n()) %>%
  mutate(age = paste0(age, "\n(n = ",n,")")) 
r <- ggplot(plot_data, aes(x = as.character(age), y = Percent, color = as.character(age))) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.5, outlier.size = 1) +
  ggpubr::theme_pubr(base_size = 8) + ylab("") + 
  stat_compare_means(color = "red", 
                     comparisons = combn(unique(plot_data$age), m = 2, simplify = F)) +
  stat_compare_means(color = "blue") +
  xlab("") + 
  ylab("% Microsatellite Instability") +
  ggtitle("% Microsatellite Instability vs. Age") +
  # geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
  theme(legend.position = "none") 
ggsave(filename = 'results/msisensor-pro/msi_sensor_vs_age.png', width = 6, height = 6)

# developmental cluster name
plot_data <- output_df %>%
  filter(!is.na(rdt.name)) %>%
  mutate(rdt.name = as.character(rdt.name)) %>%
  group_by(rdt.name) %>%
  mutate(n = n()) %>%
  mutate(rdt.name = paste0(rdt.name, "\n(n = ",n,")")) 
s <- ggplot(plot_data, aes(x = as.character(rdt.name), y = Percent, color = as.character(rdt.name))) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.5, outlier.size = 1) +
  ggpubr::theme_pubr(base_size = 8) + ylab("") + 
  stat_compare_means(color = "red", 
                     comparisons = combn(unique(plot_data$rdt.name), m = 2, simplify = F)) +
  stat_compare_means(color = "blue") +
  xlab("") + 
  ylab("% Microsatellite Instability") +
  ggtitle("% Microsatellite Instability vs. Dev. Cluster Name") +
  # geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
  theme(legend.position = "none") 
ggsave(filename = 'results/msisensor-pro/msi_sensor_vs_dev_cluster_name.png', width = 6, height = 6)

# gender
plot_data <- output_df %>%
  filter(!is.na(gender),
         !gender %in% c("Not Reported")) %>%
  mutate(gender = as.character(gender)) %>%
  group_by(gender) %>%
  mutate(n = n()) %>%
  mutate(gender = paste0(gender, "\n(n = ",n,")")) 
p <- ggplot(plot_data, aes(x = as.character(gender), y = Percent, color = as.character(gender))) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.5, outlier.size = 1) +
  ggpubr::theme_pubr(base_size = 8) + ylab("") + 
  stat_compare_means(color = "red", 
                     comparisons = combn(unique(plot_data$gender), m = 2, simplify = F)) +
  stat_compare_means(color = "blue") +
  xlab("") + 
  ylab("% Microsatellite Instability") +
  ggtitle("% Microsatellite Instability vs. Gender") +
  geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
  theme(legend.position = "none") 
ggsave(filename = 'results/msisensor-pro/msi_sensor_vs_gender.png', width = 6, height = 6)

# # add cluster from miRNA clustering (shiny)
# mirna_clusters <- read.csv('data/PBTA_Sample_Clusters.csv', check.names = F)
# mirna_clusters$sample_id <- gsub(".*-7316", "7316", mirna_clusters$sample_id)
# mirna_clusters$sample_id <- gsub("_.*", "", mirna_clusters$sample_id)
# output_df <- output_df %>%
#   inner_join(mirna_clusters, by = "sample_id")
# 
# output_df <- output_df %>%
#   group_by(`5%_cluster`) %>%
#   mutate(n = n()) %>%
#   mutate(cluster_5 = paste0(`5%_cluster`, "\n(n = ",n,")")) 
# output_df <- output_df %>%
#   group_by(`10%_cluster`) %>%
#   mutate(n = n()) %>%
#   mutate(cluster_10 = paste0(`10%_cluster`, "\n(n = ",n,")")) 
# 
# p <- ggplot(output_df, aes(x = as.character(cluster_5), y = Percent, color = as.character(cluster_5))) +
#   stat_boxplot(geom ='errorbar', width = 0.2) +
#   geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.5, outlier.size = 1) +
#   ggpubr::theme_pubr() + ylab("") + 
#   stat_compare_means(method = "wilcox.test", color = "red", 
#                      comparisons = list(c("1\n(n = 2)", "2\n(n = 2)"), 
#                                         c("1\n(n = 2)", "3\n(n = 16)"),
#                                         c("1\n(n = 2)", "4\n(n = 4)"),
#                                         c("2\n(n = 2)", "3\n(n = 16)"), 
#                                         c("2\n(n = 2)", "4\n(n = 4)"),
#                                         c("3\n(n = 16)", "4\n(n = 4)"))) +
#   xlab("5% Clusters") + 
#   ylab("% Microsatellite Instability") +
#   geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
#   theme(legend.position = "none") 
# ggsave(filename = 'results/msi_sensor_vs_mirna_clusters_5.pdf', width = 4, height = 4)
# 
# q <- ggplot(output_df, aes(x = as.character(cluster_10), y = Percent, color = as.character(cluster_10))) +
#   stat_boxplot(geom ='errorbar', width = 0.2) +
#   geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.5, outlier.size = 1) +
#   ggpubr::theme_pubr() + ylab("") + 
#   stat_compare_means(method = "wilcox.test", color = "red", 
#                      comparisons = list(c("1\n(n = 1)", "2\n(n = 3)"), 
#                                         c("1\n(n = 1)", "3\n(n = 16)"),
#                                         c("1\n(n = 1)", "4\n(n = 4)"),
#                                         c("2\n(n = 3)", "3\n(n = 16)"), 
#                                         c("2\n(n = 3)", "4\n(n = 4)"),
#                                         c("3\n(n = 16)", "4\n(n = 4)"))) +
#   xlab("10% Clusters") + 
#   ylab("% Microsatellite Instability") +
#   geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
#   theme(legend.position = "none") 
# ggsave(filename = 'results/msi_sensor_vs_mirna_clusters_10.pdf', width = 4, height = 4)

# correlation
# broom::tidy(chisq.test(x = output_df$`5%_cluster`, y = output_df$Percent))
# broom::tidy(lm(output_df$Percent~0 + as.factor(output_df$`5%_cluster`)))


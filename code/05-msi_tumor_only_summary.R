# msi sensor pro vs tmb scores 
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(ggpubr)
  library(ggrepel)
})

# output directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
results_dir <- file.path(root_dir, "results" , "msisensor-pro-tumor-only")
dir.create(results_dir, recursive = T, showWarnings = F)

# MSI output (n = 88)
output_df <- read_tsv(file.path(results_dir, "msi_output_merged_tumor_only.tsv"))

# modify type
output_df <- output_df %>%
  dplyr::mutate(Type = ifelse(Type == "High", sample_id, ""))

# 1) MSI vs TMB
ggplot(output_df, aes(MSI_Percent, TMB)) + 
  geom_point() + 
  geom_text_repel(aes(label = Type), na.rm=TRUE, hjust=0, vjust=0, size = 3, color = "red") +
  theme_pubr() + 
  xlab("% MSI") + ylab("TMB") + ggtitle("% MSI vs TMB") +
  stat_cor(method = "pearson")
ggsave(filename = file.path(results_dir, "msisensorpro_vs_tmb.png"), height = 6, width = 8)

# 2) proteomics clusters
plot_data <- output_df %>%
  dplyr::rename("proteomics_rdt_cc" = "rdt.cc") %>%
  filter(!is.na(proteomics_rdt_cc)) %>%
  mutate(proteomics_rdt_cc = as.character(proteomics_rdt_cc)) %>%
  group_by(proteomics_rdt_cc) %>%
  mutate(n = n()) %>%
  mutate(proteomics_rdt_cc = paste0(proteomics_rdt_cc, "\n(n = ",n,")")) 
p <- ggplot(plot_data, aes(x = as.character(proteomics_rdt_cc), y = MSI_Percent, color = as.character(proteomics_rdt_cc))) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.5, outlier.size = 1) +
  geom_text_repel(aes(label = Type), na.rm=TRUE, hjust=0, vjust=0, size = 2, color = "black") +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", 
                     label = "p.format",
                     ref.group = ".all.",
                     size = 4) +
  xlab("") + 
  ylab("% Microsatellite Instability") +
  ggtitle("% Microsatellite Instability vs. Proteomics Cluster") +
  geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
  theme(legend.position = "none") 
ggsave(filename = file.path(results_dir, "msi_sensor_vs_proteomics_clusters.png"), width = 6, height = 6)

# 3) developmental clusters
plot_data <- output_df %>%
  filter(!is.na(dtt.cc)) %>%
  mutate(dtt.cc = as.character(dtt.cc)) %>%
  group_by(dtt.cc) %>%
  mutate(n = n()) %>%
  mutate(dtt.cc = paste0(dtt.cc, "\n(n = ",n,")")) 
q <- ggplot(plot_data, aes(x = as.character(dtt.cc), y = MSI_Percent, color = as.character(dtt.cc))) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.5, outlier.size = 1) +
  geom_text_repel(aes(label = Type), na.rm=TRUE, hjust=0, vjust=0, size = 2, color = "black") +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", 
                     label = "p.format",
                     ref.group = ".all.",
                     size = 4) +
  xlab("") + 
  ylab("% Microsatellite Instability") +
  ggtitle("% Microsatellite Instability vs. Developmental Cluster") +
  geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
  theme(legend.position = "none") 
ggsave(filename = file.path(results_dir, "msi_sensor_vs_dev_clusters.png"), width = 6, height = 6)

# 4) age (three groups)
plot_data <- output_df %>%
  filter(!is.na(age_three_groups)) %>%
  mutate(age_three_groups = as.character(age_three_groups)) %>%
  group_by(age_three_groups) %>%
  mutate(n = n()) %>%
  mutate(age_three_groups = paste0(age_three_groups, "\n(n = ",n,")")) 
plot_data$age_three_groups <- factor(plot_data$age_three_groups, levels = c("[0,15]\n(n = 58)", "(15,26]\n(n = 22)", "(26,40]\n(n = 8)"))
r <- ggplot(plot_data, aes(x = age_three_groups, y = MSI_Percent, color = as.character(age_three_groups))) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.5, outlier.size = 1) +
  geom_text_repel(aes(label = Type), na.rm=TRUE, hjust=0, vjust=0, size = 2, color = "black") +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", 
                     label = "p.format",
                     ref.group = ".all.",
                     size = 4) +
  xlab("") + 
  ylab("% Microsatellite Instability") +
  ggtitle("% Microsatellite Instability vs. Age") +
  geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
  theme(legend.position = "none") 
ggsave(filename = file.path(results_dir, "msi_sensor_vs_age_three_groups.png"), width = 6, height = 6)

# 5) age (two groups)
plot_data <- output_df %>%
  filter(!is.na(age_two_groups)) %>%
  mutate(age_two_groups = as.character(age_two_groups)) %>%
  group_by(age_two_groups) %>%
  mutate(n = n()) %>%
  mutate(age_two_groups = paste0(age_two_groups, "\n(n = ",n,")")) 
plot_data$age_two_groups <- factor(plot_data$age_two_groups, levels = c("[0,15]\n(n = 58)", "(15,40]\n(n = 30)"))
r <- ggplot(plot_data, aes(x = age_two_groups, y = MSI_Percent, color = as.character(age_two_groups))) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.5, outlier.size = 1) +
  geom_text_repel(aes(label = Type), na.rm=TRUE, hjust=0, vjust=0, size = 2, color = "black") +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", 
                     label = "p.format",
                     ref.group = ".all.",
                     size = 4) +
  xlab("") + 
  ylab("% Microsatellite Instability") +
  ggtitle("% Microsatellite Instability vs. Age") +
  geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
  theme(legend.position = "none") 
ggsave(filename = file.path(results_dir, "msi_sensor_vs_age_two_groups.png"), width = 6, height = 6)

# 6) developmental cluster name
plot_data <- output_df %>%
  filter(!is.na(rdt.name)) %>%
  mutate(rdt.name = as.character(rdt.name)) %>%
  group_by(rdt.name) %>%
  mutate(n = n()) %>%
  mutate(rdt.name = paste0(rdt.name, "\n(n = ",n,")")) 
s <- ggplot(plot_data, aes(x = as.character(rdt.name), y = MSI_Percent, color = as.character(rdt.name))) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.5, outlier.size = 1) +
  geom_text_repel(aes(label = Type), na.rm=TRUE, hjust=0, vjust=0, size = 2, color = "black") +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", 
                     label = "p.format",
                     ref.group = ".all.",
                     size = 4) +
  xlab("") + 
  ylab("% Microsatellite Instability") +
  ggtitle("% Microsatellite Instability vs. Dev. Cluster Name") +
  geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
  theme(legend.position = "none") 
ggsave(filename = file.path(results_dir, "msi_sensor_vs_dev_cluster_name.png"), width = 8, height = 6)

# 7) gender
plot_data <- output_df %>%
  filter(!is.na(gender),
         !gender %in% c("Not Reported")) %>%
  mutate(gender = as.character(gender)) %>%
  group_by(gender) %>%
  mutate(n = n()) %>%
  mutate(gender = paste0(gender, "\n(n = ",n,")")) 
p <- ggplot(plot_data, aes(x = as.character(gender), y = MSI_Percent, color = as.character(gender))) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.5, outlier.size = 1) +
  geom_text_repel(aes(label = Type), na.rm=TRUE, hjust=0, vjust=0, size = 2, color = "black") +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", size = 4) +
  xlab("") + 
  ylab("% Microsatellite Instability") +
  ggtitle("% Microsatellite Instability vs. Gender") +
  geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
  theme(legend.position = "none") 
ggsave(filename = file.path(results_dir, "msi_sensor_vs_gender.png"), width = 6, height = 6)

# 8) ALT 
output_df <- output_df %>%
  filter(!is.na(ALT_status))
p <- ggplot(output_df, aes(x = ALT_status, MSI_Percent, color = ALT_status)) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.5, outlier.size = 1) +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", size = 4) +
  xlab("") + 
  ylab("% Microsatellite Instability") +
  ggtitle("% Microsatellite Instability vs. ALT Status") +
  geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
  theme(legend.position = "none") 
ggsave(plot = p, filename = file.path(results_dir, "msisensorpro_vs_alt_status.png"), height = 6, width = 6)

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


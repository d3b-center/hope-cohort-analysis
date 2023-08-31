# comparison of msi sensor pro (paired) vs different variables 
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(ggpubr)
  library(ggrepel)
})

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
input_dir <- file.path(root_dir, "data")
analyses_dir <- file.path(root_dir, "analyses", "msi-sensor-analysis")
output_dir <- file.path(analyses_dir, "results", "msisensor-pro-paired")
dir.create(output_dir, showWarnings = F, recursive = T)

# master histology
annot <- read_tsv(file = file.path(input_dir, "master_histology_hope_cohort.tsv"))
output_df <- annot %>%
  filter(!is.na(msi_paired))

# modify type
output_df <- output_df %>%
  dplyr::mutate(Type = ifelse(msi_paired > 3.5, Sample_ID, ""))

# 1) MSI vs TMB
p <- ggplot(output_df, aes(msi_paired, TMB_paired)) + 
  geom_point() + 
  geom_text_repel(aes(label = Type), na.rm = TRUE, hjust = 0, vjust = 0, size = 3, color = "red") +
  theme_pubr() + 
  xlab("% MSI") + ylab("TMB") + ggtitle("% MSI vs TMB") +
  stat_cor(method = "pearson", color = "red")
ggsave(plot = p, filename = file.path(output_dir, "msi_vs_tmb.png"), height = 6, width = 8)

# 2) proteomics clusters
plot_data <- output_df %>%
  dplyr::rename("proteomics_rdt_cc" = "rdt.cc") %>%
  mutate(proteomics_rdt_cc = as.character(proteomics_rdt_cc)) %>%
  group_by(proteomics_rdt_cc) %>%
  mutate(n = n()) %>%
  mutate(proteomics_rdt_cc = paste0(proteomics_rdt_cc, "\n(n = ",n,")")) 
p <- ggplot(plot_data, aes(x = as.character(proteomics_rdt_cc), y = msi_paired, color = as.character(proteomics_rdt_cc))) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.4, outlier.size = 1) +
  geom_text_repel(aes(label = Type), na.rm = TRUE, hjust = 0, vjust = 0, size = 3, color = "black") +
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
q <- ggplot(plot_data %>% filter(!grepl("NA", proteomics_rdt_cc)), aes(x = as.character(proteomics_rdt_cc), y = msi_paired, color = as.character(proteomics_rdt_cc))) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.4, outlier.size = 1) +
  geom_text_repel(aes(label = Type), na.rm = TRUE, hjust = 0, vjust = 0, size = 3, color = "black") +
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
ggsave(plot = ggarrange(plotlist = list(p, q), ncol = 2), filename = file.path(output_dir, "msi_vs_proteomics_clusters.png"), width = 10, height = 6)

# 3) developmental clusters
plot_data <- output_df %>%
  mutate(dtt.cc = as.character(dtt.cc)) %>%
  group_by(dtt.cc) %>%
  mutate(n = n()) %>%
  mutate(dtt.cc = paste0(dtt.cc, "\n(n = ",n,")")) 
p <- ggplot(plot_data, aes(x = as.character(dtt.cc), y = msi_paired, color = as.character(dtt.cc))) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.4, outlier.size = 1) +
  geom_text_repel(aes(label = Type), na.rm = TRUE, hjust = 0, vjust = 0, size = 3, color = "black") +
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
q <- ggplot(plot_data %>% filter(!grepl("NA", dtt.cc)), aes(x = as.character(dtt.cc), y = msi_paired, color = as.character(dtt.cc))) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.4, outlier.size = 1) +
  geom_text_repel(aes(label = Type), na.rm = TRUE, hjust = 0, vjust = 0, size = 3, color = "black") +
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
ggsave(plot = ggarrange(plotlist = list(p, q), ncol = 2), filename = file.path(output_dir, "msi_vs_dev_clusters.png"), width = 10, height = 6)

# 4) age (three groups)
plot_data <- output_df %>%
  mutate(age_three_groups = as.character(age_three_groups)) %>%
  group_by(age_three_groups) %>%
  mutate(n = n()) %>%
  mutate(age_three_groups = paste0(age_three_groups, "\n(n = ",n,")")) 
plot_data$age_three_groups <- factor(plot_data$age_three_groups, levels = c("[0,15]\n(n = 45)", "(15,26]\n(n = 20)", "(26,40]\n(n = 8)"))
p <- ggplot(plot_data, aes(x = age_three_groups, y = msi_paired, color = as.character(age_three_groups))) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.4, outlier.size = 1) +
  geom_text_repel(aes(label = Type), na.rm = TRUE, hjust = 0, vjust = 0, size = 3, color = "black") +
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
ggsave(plot = p, filename = file.path(output_dir, "msi_vs_age_three_groups.png"), width = 6, height = 6)

# 5) age (two groups)
plot_data <- output_df %>%
  mutate(age_two_groups = as.character(age_two_groups)) %>%
  group_by(age_two_groups) %>%
  mutate(n = n()) %>%
  mutate(age_two_groups = paste0(age_two_groups, "\n(n = ",n,")")) 
plot_data$age_two_groups <- factor(plot_data$age_two_groups, levels = c("[0,15]\n(n = 45)", "(15,40]\n(n = 28)"))
p <- ggplot(plot_data, aes(x = age_two_groups, y = msi_paired, color = as.character(age_two_groups))) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.4, outlier.size = 1) +
  geom_text_repel(aes(label = Type), na.rm = TRUE, hjust = 0, vjust = 0, size = 3, color = "black") +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", size = 4) +
  xlab("") + 
  ylab("% Microsatellite Instability") +
  ggtitle("% Microsatellite Instability vs. Age") +
  geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
  theme(legend.position = "none") 
ggsave(plot = p, filename = file.path(output_dir, "msi_vs_age_two_groups.png"), width = 6, height = 6)

# 6) developmental cluster name
plot_data <- output_df %>%
  mutate(rdt.name = as.character(rdt.name)) %>%
  group_by(rdt.name) %>%
  mutate(n = n()) %>%
  mutate(rdt.name = paste0(rdt.name, "\n(n = ",n,")")) 
p <- ggplot(plot_data, aes(x = as.character(rdt.name), y = msi_paired, color = as.character(rdt.name))) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.4, outlier.size = 1) +
  geom_text_repel(aes(label = Type), na.rm = TRUE, hjust = 0, vjust = 0, size = 3, color = "black") +
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
q <- ggplot(plot_data %>% filter(!grepl("NA", rdt.name)), aes(x = as.character(rdt.name), y = msi_paired, color = as.character(rdt.name))) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.4, outlier.size = 1) +
  geom_text_repel(aes(label = Type), na.rm = TRUE, hjust = 0, vjust = 0, size = 3, color = "black") +
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
ggsave(plot = ggarrange(plotlist = list(p, q), ncol = 2), filename = file.path(output_dir, "msi_vs_dev_cluster_name.png"), width = 18, height = 6)

# 7) gender
plot_data <- output_df %>%
  mutate(Gender = as.character(Gender)) %>%
  group_by(Gender) %>%
  mutate(n = n()) %>%
  mutate(Gender = paste0(Gender, "\n(n = ",n,")")) 
p <- ggplot(plot_data, aes(x = as.character(Gender), y = msi_paired, color = as.character(Gender))) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.4, outlier.size = 1) +
  geom_text_repel(aes(label = Type), na.rm = TRUE, hjust = 0, vjust = 0, size = 3, color = "black") +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", size = 4) +
  xlab("") + 
  ylab("% Microsatellite Instability") +
  ggtitle("% Microsatellite Instability vs. Gender") +
  geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
  theme(legend.position = "none") 
ggsave(plot = p, filename = file.path(output_dir, "msi_vs_gender.png"), width = 6, height = 6)

# 8) ALT 
output_df = output_df %>% 
  group_by(ALT_status) %>% 
  mutate(n = n(), ALT_status = paste0(ALT_status, "\n(n = ", n, ")"))
p <- ggplot(output_df, aes(x = ALT_status, msi_paired, color = ALT_status)) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.4, outlier.size = 1) +
  geom_text_repel(aes(label = Type), na.rm = TRUE, hjust = 0, vjust = 0, size = 3, color = "black") +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", 
                     label = "p.format",
                     ref.group = ".all.",
                     size = 4) +
  xlab("") + 
  ylab("% Microsatellite Instability") +
  ggtitle("% Microsatellite Instability vs. ALT Status") +
  geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
  theme(legend.position = "none") 
ggsave(plot = p, filename = file.path(output_dir, "msi_vs_alt_status.png"), height = 6, width = 6)

# 9) MSI vs ALT telomere content
p <- ggplot(output_df, aes(x = msi_paired, y = t_n_telomere_content)) +
  geom_point() +
  geom_text_repel(aes(label = Type), na.rm = TRUE, hjust = 0, vjust = 0, size = 3, color = "red") +
  theme_pubr() + 
  xlab("% MSI") + ylab("TMB") + ggtitle("% MSI vs Telomere content") +
  stat_cor(method = "pearson", color = "red")
ggsave(plot = p, filename = file.path(output_dir, "msi_vs_telomere_content.png"))

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


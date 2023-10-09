# Function: comparison of msi sensor pro (tumor-only) vs clinical variables 

# load libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(ggpubr)
  library(ggrepel)
})

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data", "v1")
analyses_dir <- file.path(root_dir, "analyses", "msi-sensor-analysis")
input_dir <- file.path(analyses_dir, "input")
output_dir <- file.path(analyses_dir, "results", "msisensor-pro-tumor-only")
dir.create(output_dir, showWarnings = F, recursive = T)

# read histologies 
annot <- read_tsv(file = file.path(data_dir, "Hope-GBM-histologies.tsv"))
annot <- annot %>% 
  filter(!is.na(molecular_subtype))

# read MSI tumor-only output 
msi_tumor_only_output <- read_tsv(file.path(output_dir, "Hope-msi-tumor_only.tsv")) %>%
  dplyr::mutate(Type = ifelse(Percent > 3.5, sample_id, "")) %>%
  dplyr::rename("msi_tumor_only" = "Percent")

# read TMB tumor-only output 
tmb_tumor_only_output <- read_tsv("../tmb-calculation/results/wgs_tumor_only/snv-mutation-tmb-coding.tsv") %>%
  dplyr::rename("tmb_tumor_only" = "tmb")

# add proteomics cluster data
proteomic_data <- read_tsv(file.path(input_dir, "cluster_data_090722.tsv")) %>%
  dplyr::select(id, rdt.cc)

# add developmental cluster data
dev_data <- read_tsv(file.path(input_dir, "cluster_data_101922.tsv")) %>%
  dplyr::select(id, dtt.cc, rdt.name)

# combine all
output_df <- msi_tumor_only_output %>%
  inner_join(tmb_tumor_only_output, by = c("Kids_First_Biospecimen_ID" = "Tumor_Sample_Barcode")) %>%
  inner_join(proteomic_data, by = c("sample_id" = "id")) %>%
  inner_join(dev_data, by = c("sample_id" = "id")) %>%
  inner_join(annot)

# 1) MSI vs TMB_tumor_only
ggplot(output_df, aes(msi_tumor_only, tmb_tumor_only)) + 
  geom_point() + 
  geom_text_repel(aes(label = Type), na.rm = TRUE, hjust = 0, vjust = 0, size = 3, color = "red") +
  theme_pubr() + 
  xlab("% MSI") + ylab("TMB") + ggtitle("% MSI vs TMB (tumor-only)") +
  stat_cor(method = "pearson", color = "red")
ggsave(filename = file.path(output_dir, "msi_vs_tmb.png"), height = 6, width = 8)

# 2) proteomics clusters
plot_data <- output_df %>%
  dplyr::rename("proteomics_rdt_cc" = "rdt.cc") %>%
  mutate(proteomics_rdt_cc = as.character(proteomics_rdt_cc)) %>%
  group_by(proteomics_rdt_cc) %>%
  mutate(n = n()) %>%
  mutate(proteomics_rdt_cc = paste0(proteomics_rdt_cc, "\n(n = ",n,")")) 
p <- ggplot(plot_data, aes(x = as.character(proteomics_rdt_cc), y = msi_tumor_only, color = as.character(proteomics_rdt_cc))) +
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
  ggtitle("% Microsatellite Instability vs. Proteomics Cluster") +
  geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
  theme(legend.position = "none") 
q <- ggplot(plot_data %>% filter(!grepl("NA", proteomics_rdt_cc)), aes(x = as.character(proteomics_rdt_cc), y = msi_tumor_only, color = as.character(proteomics_rdt_cc))) +
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
p <- ggplot(plot_data, aes(x = as.character(dtt.cc), y = msi_tumor_only, color = as.character(dtt.cc))) +
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
q <- ggplot(plot_data %>% filter(!grepl("NA", dtt.cc)), aes(x = as.character(dtt.cc), y = msi_tumor_only, color = as.character(dtt.cc))) +
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
  mutate(HARMONY_age_class_derived = as.character(HARMONY_age_class_derived)) %>%
  group_by(HARMONY_age_class_derived) %>%
  mutate(n = n()) %>%
  mutate(HARMONY_age_class_derived = paste0(HARMONY_age_class_derived, "\n(n = ",n,")")) 
plot_data$HARMONY_age_class_derived <- factor(plot_data$HARMONY_age_class_derived, levels = c("[0,15]\n(n = 55)", "(15,26]\n(n = 23)", "(26,40]\n(n = 8)"))
p <- ggplot(plot_data, aes(x = HARMONY_age_class_derived, y = msi_tumor_only, color = as.character(HARMONY_age_class_derived))) +
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
  mutate(HARMONY_age_class_derived = as.character(HARMONY_age_class_derived)) %>%
  mutate(HARMONY_age_class_derived = ifelse(HARMONY_age_class_derived %in% c("(15,26]", "(26,40]"), "(15,40]", HARMONY_age_class_derived)) %>%
  group_by(HARMONY_age_class_derived) %>%
  mutate(n = n()) %>%
  mutate(HARMONY_age_class_derived = paste0(HARMONY_age_class_derived, "\n(n = ",n,")")) 
plot_data$HARMONY_age_class_derived <- factor(plot_data$HARMONY_age_class_derived, levels = c("[0,15]\n(n = 55)", "(15,40]\n(n = 31)"))
p <- ggplot(plot_data, aes(x = HARMONY_age_class_derived, y = msi_tumor_only, color = as.character(HARMONY_age_class_derived))) +
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
p <- ggplot(plot_data, aes(x = as.character(rdt.name), y = msi_tumor_only, color = as.character(rdt.name))) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.4, outlier.size = 1) +
  geom_text_repel(aes(label = Type), na.rm=TRUE, hjust = 0, vjust = 0, size = 3, color = "black") +
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
q <- ggplot(plot_data %>% filter(!grepl("NA", rdt.name)), aes(x = as.character(rdt.name), y = msi_tumor_only, color = as.character(rdt.name))) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.4, outlier.size = 1) +
  geom_text_repel(aes(label = Type), na.rm=TRUE, hjust = 0, vjust = 0, size = 3, color = "black") +
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

# 7) Gender
plot_data <- output_df %>%
  mutate(HARMONY_Gender = as.character(HARMONY_Gender)) %>%
  group_by(HARMONY_Gender) %>%
  mutate(n = n()) %>%
  mutate(HARMONY_Gender = paste0(HARMONY_Gender, "\n(n = ",n,")")) 
p <- ggplot(plot_data, aes(x = as.character(HARMONY_Gender), y = msi_tumor_only, color = as.character(HARMONY_Gender))) +
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

# # 8) ALT 
# plot_data <- output_df %>% 
#   group_by(ALT_status) %>% 
#   mutate(n = n(), ALT_status = paste0(ALT_status, "\n(n = ", n, ")"))
# p <- ggplot(plot_data, aes(x = ALT_status, msi_tumor_only, color = ALT_status)) +
#   stat_boxplot(geom ='errorbar', width = 0.2) +
#   geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.4, outlier.size = 1) +
#   geom_text_repel(aes(label = Type), na.rm = TRUE, hjust = 0, vjust = 0, size = 3, color = "black") +
#   ggpubr::theme_pubr(base_size = 10) + ylab("") + 
#   stat_compare_means(color = "red", 
#                      label = "p.format",
#                      ref.group = ".all.",
#                      size = 4) +
#   xlab("") + 
#   ylab("% Microsatellite Instability") +
#   ggtitle("% Microsatellite Instability vs. ALT Status") +
#   geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
#   theme(legend.position = "none") 
# q <- ggplot(plot_data %>% filter(!grepl("NA", ALT_status)), aes(x = ALT_status, msi_tumor_only, color = ALT_status)) +
#   stat_boxplot(geom ='errorbar', width = 0.2) +
#   geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.4, outlier.size = 1) +
#   geom_text_repel(aes(label = Type), na.rm = TRUE, hjust = 0, vjust = 0, size = 3, color = "black") +
#   ggpubr::theme_pubr(base_size = 10) + ylab("") + 
#   stat_compare_means(color = "red", size = 4) +
#   xlab("") + 
#   ylab("% Microsatellite Instability") +
#   ggtitle("% Microsatellite Instability vs. ALT Status") +
#   geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
#   theme(legend.position = "none") 
# ggsave(plot = ggarrange(plotlist = list(p, q), ncol = 2), filename = file.path(output_dir, "msi_vs_alt_status.png"), height = 6, width = 10)
# 
# # 9) MSI vs ALT telomere content
# p <- ggplot(output_df, aes(x = msi_paired, y = t_n_telomere_content)) +
#   geom_point() +
#   geom_text_repel(aes(label = Type), na.rm = TRUE, hjust = 0, vjust = 0, size = 3, color = "red") +
#   theme_pubr() + 
#   xlab("% MSI") + ylab("TMB") + ggtitle("% MSI vs Telomere content") +
#   stat_cor(method = "pearson", color = "red")
# ggsave(plot = p, filename = file.path(output_dir, "msi_vs_telomere_content.png"))


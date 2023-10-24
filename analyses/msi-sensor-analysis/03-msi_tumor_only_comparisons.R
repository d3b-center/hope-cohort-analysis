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
  filter(!is.na(HOPE_diagnosis))

# read MSI tumor-only output 
msi_tumor_only_output <- read_tsv(file.path(output_dir, "Hope-msi-tumor_only.tsv")) %>%
  dplyr::mutate(Type = ifelse(Percent > 3.5, sample_id, "")) %>%
  dplyr::rename("msi_tumor_only" = "Percent")

# read TMB tumor-only output 
tmb_tumor_only_output <- read_tsv("../tmb-calculation/results/wgs_tumor_only/snv-mutation-tmb-coding.tsv") %>%
  dplyr::rename("tmb_tumor_only" = "tmb")

# add proteomic/developmental cluster data
cluster_data <- read_tsv(file.path(input_dir, "cluster_data_101922.tsv")) %>%
  dplyr::select(id, dtt.cc, rdt.name)

# combine all
output_df <- msi_tumor_only_output %>%
  inner_join(tmb_tumor_only_output, by = c("Kids_First_Biospecimen_ID" = "Tumor_Sample_Barcode")) %>%
  inner_join(cluster_data, by = c("sample_id" = "id")) %>%
  inner_join(annot)

# 1) MSI vs TMB (tumor only)
p <- ggplot(output_df, aes(msi_tumor_only, tmb_tumor_only)) + 
  geom_point(shape = 21) + 
  geom_text_repel(aes(label = Type), na.rm = TRUE, hjust = 0, vjust = 0, size = 3, color = "red") +
  theme_pubr() + 
  xlab("% Microsatellite Instability") + ylab("TMB") + ggtitle("% Microsatellite Instability vs TMB (tumor-only)") +
  stat_cor(method = "pearson", color = "red")
ggsave(plot = p, filename = file.path(output_dir, "msi_vs_tmb.pdf"))

# 2) MSI vs Age (three groups)
plot_data <- output_df %>%
  mutate(HARMONY_age_class_derived = as.character(HARMONY_age_class_derived)) %>%
  group_by(HARMONY_age_class_derived) %>%
  mutate(n = n()) %>%
  mutate(HARMONY_age_class_derived = paste0(HARMONY_age_class_derived, "\n(n = ",n,")")) 
plot_data$HARMONY_age_class_derived <- factor(plot_data$HARMONY_age_class_derived, levels = c("[0,15]\n(n = 58)", "(15,26]\n(n = 23)", "(26,40]\n(n = 8)"))
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
  theme(legend.position = "none") +
  scale_color_manual(values = c("[0,15]\n(n = 58)" = "#C7E9C0",
                                "(15,26]\n(n = 23)" = "#74C476",
                                "(26,40]\n(n = 8)" = "#238B45"))
ggsave(plot = p, filename = file.path(output_dir, "msi_vs_age_three_groups.pdf"), width = 6, height = 6)

# 5) MSI vs Age (two groups)
plot_data <- output_df %>%
  mutate(HARMONY_age_class_derived = as.character(HARMONY_age_class_derived)) %>%
  mutate(HARMONY_age_class_derived = ifelse(HARMONY_age_class_derived %in% c("(15,26]", "(26,40]"), "(15,40]", HARMONY_age_class_derived)) %>%
  group_by(HARMONY_age_class_derived) %>%
  mutate(n = n()) %>%
  mutate(HARMONY_age_class_derived = paste0(HARMONY_age_class_derived, "\n(n = ",n,")")) 
plot_data$HARMONY_age_class_derived <- factor(plot_data$HARMONY_age_class_derived, levels = c("[0,15]\n(n = 58)", "(15,40]\n(n = 31)"))
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
  theme(legend.position = "none") +
  scale_color_manual(values = c("[0,15]\n(n = 58)" = "#C7E9C0",
                                "(15,40]\n(n = 31)" = "#238B45"))
ggsave(plot = p, filename = file.path(output_dir, "msi_vs_age_two_groups.pdf"), width = 6, height = 6)

# 4) MSI vs Developmental cluster name
plot_data <- output_df %>%
  mutate(rdt.name = as.character(rdt.name)) %>%
  group_by(rdt.name) %>%
  mutate(n = n()) %>%
  mutate(rdt.name = paste0(rdt.name, "\n(n = ",n,")")) 
plot_data$rdt.name <- factor(plot_data$rdt.name, levels = plot_data %>% arrange(dtt.cc) %>% pull(rdt.name) %>% unique())
pdf(file = file.path(output_dir, "msi_vs_dev_cluster_name.pdf"), width = 10, height = 6)
p <- ggplot(plot_data, aes(x = rdt.name, y = msi_tumor_only, color = rdt.name)) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.4, outlier.size = 1) +
  geom_text_repel(aes(label = Type), na.rm=TRUE, hjust = 0, vjust = 0, size = 3, color = "black") +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", size = 4) +
  xlab("") + 
  ylab("% Microsatellite Instability") +
  ggtitle("% Microsatellite Instability vs. Dev. Cluster Name") +
  geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
  theme(legend.position = "none") +
  scale_color_manual(values = c("Classical\n(n = 23)" = "#88BED8",
                                "Mesenchymal-IDHMutant\n(n = 40)" = "#89A544",
                                "Mesenchymal-IDHWT\n(n = 13)" = "#CE9D21",
                                "Pro-neural\n(n = 11)" = "#CE61A2",
                                "NA\n(n = 2)" = "gray"))
p
q <- ggplot(plot_data %>% filter(!grepl("NA", rdt.name)), aes(x = rdt.name, y = msi_tumor_only, color = rdt.name)) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.4, outlier.size = 1) +
  geom_text_repel(aes(label = Type), na.rm=TRUE, hjust = 0, vjust = 0, size = 3, color = "black") +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", size = 4) +
  xlab("") + 
  ylab("% Microsatellite Instability") +
  ggtitle("% Microsatellite Instability vs. Dev. Cluster Name") +
  geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
  theme(legend.position = "none") +
  scale_color_manual(values = c("Classical\n(n = 23)" = "#88BED8",
                                "Mesenchymal-IDHMutant\n(n = 40)" = "#89A544",
                                "Mesenchymal-IDHWT\n(n = 13)" = "#CE9D21",
                                "Pro-neural\n(n = 11)" = "#CE61A2"))
q
dev.off()

# 5) MSI vs Gender
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
  theme(legend.position = "none") +
  scale_color_manual(values = c("Male\n(n = 48)" = "#0707CF",
                                "Female\n(n = 41)" = "#CC0303"))
ggsave(plot = p, filename = file.path(output_dir, "msi_vs_gender.pdf"), width = 6, height = 6)

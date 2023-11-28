# Function: correlate ALT status to clinical variables like Protein clusters, Age (2 and 3 groups), Gender and tmb_paired

# load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(ggforce)
  library(ggrepel)
})

# output directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analyses_dir <- file.path(root_dir, "analyses", "alt-analysis")
input_dir <- file.path(analyses_dir, "input")
output_dir <- file.path(analyses_dir, "results")
dir.create(output_dir, showWarnings = F, recursive = T)

# read histologies 
annot <- read_tsv(file = file.path(data_dir, "Hope-GBM-histologies.tsv"))
annot <- annot %>% 
  filter(!is.na(HOPE_diagnosis))

# read MSI paired output 
msi_paired_output <- read_tsv(file.path("../msi-sensor-analysis/results/msisensor-pro-paired/Hope-msi-paired.tsv")) %>%
  dplyr::mutate(Type = ifelse(Percent > 3.5, sample_id, "")) %>%
  dplyr::rename("msi_paired" = "Percent")

# read tmb_paired paired output 
tmb_paired_paired_output <- read_tsv("../tmb-calculation/results/wgs_paired/snv-mutation-tmb-coding.tsv") %>%
  dplyr::rename("tmb_paired" = "tmb")

# read ALT status
alt_status_output <- read_tsv(file.path(output_dir, "alt_status_aya_hgg.tsv")) %>%
  dplyr::select(sample_id, t_n_telomere_content, ALT_status)

# add proteomic/developmental cluster data
cluster_data <- read_tsv(file.path(input_dir, "cluster_data_101922.tsv")) %>%
  dplyr::select(id, dtt.cc, rdt.name)

# combine all (n = 70)
output_df <- msi_paired_output %>%
  inner_join(tmb_paired_paired_output, by = c("Kids_First_Biospecimen_ID" = "Tumor_Sample_Barcode")) %>%
  inner_join(alt_status_output) %>%
  inner_join(cluster_data, by = c("sample_id" = "id")) %>%
  inner_join(annot)

# split Age into 2 and 3 groups
output_df <- output_df %>%
  mutate(HARMONY_age_class_derived = as.character(HARMONY_age_class_derived)) %>%
  mutate(age_three_groups = HARMONY_age_class_derived,
         age_two_groups = ifelse(HARMONY_age_class_derived %in% c("(15,26]", "(26,40]"), "(15,40]", HARMONY_age_class_derived))

# function to compute correlations
compute_corr <- function(x){

  # ALT status vs Age (two groups)
  alt_status_vs_age_two_group_chisq_pvalue <- chisq.test(x = factor(x$ALT_status), y = factor(x$age_two_groups))$p.value
  alt_status_vs_age_two_group_chisq_pvalue <- round(alt_status_vs_age_two_group_chisq_pvalue, digits = 3)
  
  # ALT status vs Age (three groups)
  alt_status_vs_age_three_group_chisq_pvalue <- chisq.test(x = factor(x$ALT_status), y = factor(x$age_three_groups))$p.value
  alt_status_vs_age_three_group_chisq_pvalue <- round(alt_status_vs_age_three_group_chisq_pvalue, digits = 3)
  
  # ALT status vs Sex
  alt_status_vs_gender_chisq_pvalue <- chisq.test(x = factor(x$ALT_status), y = factor(x$HARMONY_Gender))$p.value
  alt_status_vs_gender_chisq_pvalue <- round(alt_status_vs_gender_chisq_pvalue, digits = 3)
  
  # ALT status vs cluster name
  alt_status_vs_protein_cluster_name_chisq_pvalue <- chisq.test(x = factor(x$ALT_status), y = factor(x$rdt.name))$p.value
  alt_status_vs_protein_cluster_name_chisq_pvalue <- round(alt_status_vs_protein_cluster_name_chisq_pvalue, digits = 3)
  
  # summarize
  output_df <- data.frame(alt_status_vs_age_two_group_chisq_pvalue, 
                          alt_status_vs_age_three_group_chisq_pvalue, 
                          alt_status_vs_gender_chisq_pvalue, 
                          alt_status_vs_protein_cluster_name_chisq_pvalue)
  return(output_df)
}

# 1) ALT status chisq comparisons
df <- compute_corr(x = output_df)
df <- reshape2::melt(df, variable.name = "comparison", value.name = "chisq_pvalue")
write_tsv(x = df, file = file.path(output_dir, "alt_status_chisq_output.tsv"))

# 2) ALT status vs TMB (paired) boxplot
plot_data <- output_df %>% 
  group_by(ALT_status) %>% 
  mutate(n = n(), ALT_status = paste0(ALT_status, "\n(n = ", n, ")"))
p <- ggplot(plot_data, aes(x = ALT_status, y = tmb_paired, color = ALT_status)) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.4, outlier.size = 1) +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", size = 4) +
  xlab("") + 
  ylab("TMB") +
  ggtitle("ALT Status vs TMB") +
  theme(legend.position = "none") 
ggsave(plot = p, filename = file.path(output_dir, "alt_status_vs_tmb.pdf"), height = 6, width = 6)

# 3) ALT status vs MSI (paired)
p <- ggplot(plot_data, aes(x = ALT_status, y = msi_paired, color = ALT_status)) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.4, outlier.size = 1) +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", size = 4) +
  xlab("") + 
  ylab("% Microsatellite instability") +
  ggtitle("ALT Status vs % Microsatellite instability") +
  theme(legend.position = "none") 
ggsave(plot = p, filename = file.path(output_dir, "alt_status_vs_msi.pdf"), height = 6, width = 6)

# 4) ALT Telomere content vs ALT status
p <- ggplot(plot_data, aes(x = ALT_status, y = t_n_telomere_content, color = ALT_status)) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.4, outlier.size = 1) +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", size = 4) +
  xlab("") + 
  ylab("Telomere Content") +
  ggtitle("Telomere Content vs ALT Status") +
  theme(legend.position = "none") 
ggsave(plot = p, filename = file.path(output_dir, "telomere_content_vs_alt_status.pdf"), height = 6, width = 4)

# 5) ALT telomere content vs Age three groups
plot_df <- output_df %>%
  group_by(age_three_groups) %>%
  mutate(n = n(), age_three_groups = paste0(age_three_groups, "\n(n = ", n, ")"))
plot_df$age_three_groups <- factor(plot_df$age_three_groups,  levels = c("[0,15]\n(n = 43)", "(15,26]\n(n = 19)", "(26,40]\n(n = 8)"))
p <- ggplot(plot_df, aes(x = age_three_groups, y = t_n_telomere_content, color = age_three_groups)) +
  geom_boxplot(coef=0) +
  geom_sina() +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", size = 4, method = "kruskal.test") +
  xlab("") + 
  ylab("Telomere Content\n") +
  ggtitle("ALT telomere content vs. Age") + 
  theme(legend.position = "none") +
  scale_color_manual(values = c("[0,15]\n(n = 43)" = "#C7E9C0",
                                "(15,26]\n(n = 19)" = "#74C476",
                                "(26,40]\n(n = 8)" = "#238B45"))
ggsave(plot = p, filename = file.path(output_dir, "telomere_content_vs_age_three_groups.pdf"), height = 6, width = 6)

# 6) ALT telomere content vs Age two groups 
plot_df <- output_df %>%
  group_by(age_two_groups) %>%
  mutate(n = n(), age_two_groups = paste0(age_two_groups, "\n(n = ", n, ")"))
plot_df$age_two_groups <- factor(plot_df$age_two_groups, levels = c("[0,15]\n(n = 43)", "(15,40]\n(n = 27)"))
q <- ggplot(plot_df, aes(x = age_two_groups, y = t_n_telomere_content, color = age_two_groups)) +
  geom_boxplot(coef=0) +
  geom_sina() +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", size = 4, method = "wilcox.test") +
  xlab("") + 
  ylab("Telomere Content\n") +
  ggtitle("ALT telomere content vs. Age") + 
  theme(legend.position = "none") +
  scale_color_manual(values = c("[0,15]\n(n = 43)" = "#C7E9C0",
                                "(15,40]\n(n = 27)" = "#238B45"))
ggsave(plot = q, filename = file.path(output_dir, "telomere_content_vs_age_two_groups.pdf"), height = 6, width = 6)

# colored by tumor location
# q <- ggplot(plot_df, aes(x = age_two_groups, y = t_n_telomere_content)) +
#   geom_boxplot(coef = 0, outlier.shape = NA) +
#   geom_jitter(aes(color = HOPE_Tumor.Location.condensed)) +
#   ggpubr::theme_pubr(base_size = 10, legend = "right") + ylab("") + 
#   stat_compare_means(color = "red", size = 4, method = "wilcox.test") +
#   xlab("") + 
#   ylab("Telomere Content\n") +
#   ggtitle("ALT telomere content vs. Age") 
q <- ggplot(plot_df, aes(x = age_two_groups, y = t_n_telomere_content, color = age_two_groups)) +
  geom_boxplot(coef = 0, outlier.shape = NA) +
  geom_jitter(pch = 21, aes(fill = HOPE_Tumor.Location.condensed), size = 3) +
  ggpubr::theme_pubr(base_size = 10, legend = "right") + ylab("") + 
  stat_compare_means(color = "red", size = 4, method = "wilcox.test") +
  xlab("") + 
  ylab("Telomere Content\n") +
  ggtitle("ALT telomere content vs. Age") +
  scale_color_manual(values = c("[0,15]\n(n = 43)" = "#C7E9C0",
                                "(15,40]\n(n = 27)" = "#238B45")) +
  guides(color = "none")
ggsave(plot = q, filename = file.path(output_dir, "telomere_content_vs_age_two_groups_by_tumor_loc.pdf"), height = 6, width = 8)

# 7) ALT telomere content vs Gender
plot_data <- output_df %>%
  mutate(HARMONY_Gender = as.character(HARMONY_Gender)) %>%
  group_by(HARMONY_Gender) %>%
  mutate(n = n()) %>%
  mutate(HARMONY_Gender = paste0(HARMONY_Gender, "\n(n = ",n,")")) 
p <- ggplot(plot_data, aes(x = HARMONY_Gender, y = t_n_telomere_content, color = HARMONY_Gender)) +
  geom_boxplot(coef=0) +
  geom_sina() +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", size = 4, method = "wilcox.test") +
  xlab("") + 
  ylab("Telomere Content\n") +
  ggtitle("ALT telomere content vs. Gender") + 
  theme(legend.position = "none") +
  scale_color_manual(values = c("Male\n(n = 41)" = "#0707CF",
                                "Female\n(n = 29)" = "#CC0303"))
ggsave(plot = p, filename = file.path(output_dir, "telomere_content_vs_gender.pdf"), height = 6, width = 6)

# 8) ALT telomere content vs Developmental cluster name
plot_data <- output_df %>%
  mutate(rdt.name = as.character(rdt.name)) %>%
  group_by(rdt.name) %>%
  mutate(n = n()) %>%
  mutate(rdt.name = paste0(rdt.name, "\n(n = ",n,")")) 
plot_data$rdt.name <- factor(plot_data$rdt.name, levels = plot_data %>% arrange(dtt.cc) %>% pull(rdt.name) %>% unique())
p <- ggplot(plot_data, aes(x = rdt.name, y = t_n_telomere_content, color = rdt.name)) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.4, outlier.size = 1) +
  geom_text_repel(aes(label = Type), na.rm = TRUE, hjust = 0, vjust = 0, size = 3, color = "black") +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", size = 4) +
  xlab("") + 
  ylab("Telomere Content") +
  ggtitle("Telomere Content vs. Dev. Cluster") +
  geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
  theme(legend.position = "none") +
  scale_color_manual(values = c("Classical\n(n = 17)" = "#88BED8",
                                "Mesenchymal-IDHMutant\n(n = 35)" = "#89A544",
                                "Mesenchymal-IDHWT\n(n = 9)" = "#CE9D21",
                                "Pro-neural\n(n = 9)" = "#CE61A2"))
ggsave(plot = p, filename = file.path(output_dir, "telomere_content_vs_dev_cluster_name.pdf"), width = 10, height = 6)

# 8) ALT telomere content vs TMB (paired) scatter plot
p <- ggplot(output_df, aes(x = t_n_telomere_content, y = tmb_paired)) +
  geom_point(pch = 21) +
  theme_pubr() + 
  xlab("Telomere Content") + ylab("TMB") + ggtitle("Telomere content vs TMB") +
  stat_cor(method = "pearson", color = "red")
ggsave(plot = p, filename = file.path(output_dir, "telomere_content_vs_tmb.pdf"))

# 9) ALT telomere content vs MSI (paired)
p <- ggplot(output_df, aes(x = t_n_telomere_content, y = msi_paired)) +
  geom_point(pch = 21) +
  theme_pubr() + 
  xlab("Telomere Content") + ylab("% MSI") + ggtitle("Telomere content vs % MSI") +
  stat_cor(method = "pearson", color = "red")
ggsave(plot = p, filename = file.path(output_dir, "telomere_content_vs_msi.pdf"))

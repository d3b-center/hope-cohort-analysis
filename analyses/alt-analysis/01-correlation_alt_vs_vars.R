# Function: correlate ALT status to clinical variables like Protein clusters, Age (2 and 3 groups), Gender and tmb_paired

# load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(ggforce)
})

# output directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data", "v1")
analyses_dir <- file.path(root_dir, "analyses", "alt-analysis")
input_dir <- file.path(analyses_dir, "input")
output_dir <- file.path(analyses_dir, "results")
dir.create(output_dir, showWarnings = F, recursive = T)

# read histologies 
annot <- read_tsv(file = file.path(data_dir, "Hope-GBM-histologies.tsv"))
annot <- annot %>% 
  filter(!is.na(molecular_subtype))

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

# add proteomics cluster data
proteomic_data <- read_tsv(file.path(input_dir, "cluster_data_090722.tsv")) %>%
  dplyr::select(id, rdt.cc)

# add developmental cluster data
dev_data <- read_tsv(file.path(input_dir, "cluster_data_101922.tsv")) %>%
  dplyr::select(id, dtt.cc, rdt.name)

# combine all (n = 71)
output_df <- msi_paired_output %>%
  inner_join(tmb_paired_paired_output, by = c("Kids_First_Biospecimen_ID" = "Tumor_Sample_Barcode")) %>%
  inner_join(alt_status_output) %>%
  left_join(proteomic_data, by = c("sample_id" = "id")) %>%
  left_join(dev_data, by = c("sample_id" = "id")) %>%
  inner_join(annot)

# split Age into 2 and 3 groups
output_df <- output_df %>%
  mutate(HARMONY_age_class_derived = as.character(HARMONY_age_class_derived)) %>%
  mutate(age_three_groups = HARMONY_age_class_derived,
         age_two_groups = ifelse(HARMONY_age_class_derived %in% c("(15,26]", "(26,40]"), "(15,40]", HARMONY_age_class_derived))

# function to compute correlations
compute_corr <- function(x){

  # t_n_telomere_content vs MSI
  tel_content_vs_msi_paired_correlation = cor.test(x = x$t_n_telomere_content, y = x$msi_paired)$estimate
  tel_content_vs_msi_paired_pvalue = cor.test(x = x$t_n_telomere_content, y = x$msi_paired)$p.value
  
  # t_n_telomere_content vs tmb_paired
  tel_content_vs_tmb_paired_correlation = cor.test(x = x$t_n_telomere_content, y = x$tmb_paired)$estimate
  tel_content_vs_tmb_paired_pvalue = cor.test(x = x$t_n_telomere_content, y = x$tmb_paired)$p.value
  
  # ALT status vs Age (two groups)
  alt_status_vs_age_two_group_chisq_pvalue <- chisq.test(x = factor(x$ALT_status), y = factor(x$age_two_groups))$p.value
  alt_status_vs_age_two_group_chisq_pvalue <- round(alt_status_vs_age_two_group_chisq_pvalue, digits = 3)
  
  # ALT status vs Age (three groups)
  alt_status_vs_age_three_group_chisq_pvalue <- chisq.test(x = factor(x$ALT_status), y = factor(x$age_three_groups))$p.value
  alt_status_vs_age_three_group_chisq_pvalue <- round(alt_status_vs_age_three_group_chisq_pvalue, digits = 3)
  
  # ALT status vs Sex
  alt_status_vs_gender_chisq_pvalue <- chisq.test(x = factor(x$ALT_status), y = factor(x$HARMONY_Gender))$p.value
  alt_status_vs_gender_chisq_pvalue <- round(alt_status_vs_gender_chisq_pvalue, digits = 3)
  
  # ALT status vs Protein cluster
  alt_status_vs_protein_cluster_chisq_pvalue <- chisq.test(x = factor(x$ALT_status), y = factor(x$rdt.cc))$p.value
  alt_status_vs_protein_cluster_chisq_pvalue <- round(alt_status_vs_protein_cluster_chisq_pvalue, digits = 3)
  
  # ALT status vs Protein cluster name
  alt_status_vs_protein_cluster_name_chisq_pvalue <- chisq.test(x = factor(x$ALT_status), y = factor(x$rdt.name))$p.value
  alt_status_vs_protein_cluster_name_chisq_pvalue <- round(alt_status_vs_protein_cluster_name_chisq_pvalue, digits = 3)
  
  # ALT status vs MSI (paired)
  # msi_pvalue <- chisq.test(x = factor(x$ALT_status), y = x$msi_paired)$p.value
  alt_status_vs_msi_paired_chisq_pvalue <- wilcox.test(formula = x$msi_paired ~ factor(x$ALT_status))$p.value
  alt_status_vs_msi_paired_chisq_pvalue <- round(alt_status_vs_msi_paired_chisq_pvalue, digits = 3)
  
  # ALT status vs TMB (paired)
  # tmb_paired_pvalue <- chisq.test(x = factor(x$ALT_status), y = x$tmb_paired)$p.value
  alt_status_vs_tmb_paired_chisq_pvalue <- wilcox.test(formula = x$tmb_paired ~ factor(x$ALT_status))$p.value
  alt_status_vs_tmb_paired_chisq_pvalue <- round(alt_status_vs_tmb_paired_chisq_pvalue, digits = 3)

  # summarize
  output_df <- data.frame(tel_content_vs_msi_paired_correlation,
                          tel_content_vs_msi_paired_pvalue,
                          tel_content_vs_tmb_paired_correlation,
                          tel_content_vs_tmb_paired_pvalue,
                          alt_status_vs_age_two_group_chisq_pvalue, 
                          alt_status_vs_age_three_group_chisq_pvalue, 
                          alt_status_vs_gender_chisq_pvalue, 
                          alt_status_vs_protein_cluster_chisq_pvalue, 
                          alt_status_vs_protein_cluster_name_chisq_pvalue,
                          alt_status_vs_msi_paired_chisq_pvalue, 
                          alt_status_vs_tmb_paired_chisq_pvalue)
  return(output_df)
}

# 1) correlations
df <- compute_corr(x = output_df)
df <- reshape2::melt(df)
write_tsv(x = df, file = file.path(output_dir, "alt_correlations_paired.tsv"))

# 2) ALT telomere content vs TMB (paired) scatter plot
p <- ggplot(output_df, aes(x = t_n_telomere_content, y = tmb_paired)) +
  geom_point(pch = 21) +
  theme_pubr() + 
  xlab("Telomere Content") + ylab("tmb_paired") + ggtitle("Telomere content vs tmb_paired") +
  stat_cor(method = "pearson", color = "red")
p
ggsave(plot = p, filename = file.path(output_dir, "alt_content_vs_tmb_paired.png"))

# 3) ALT status vs TMB (paired) boxplot
plot_data <- output_df %>% 
  group_by(ALT_status) %>% 
  mutate(n = n(), ALT_status = paste0(ALT_status, "\n(n = ", n, ")"))
p <- ggplot(plot_data, aes(x = ALT_status, y = tmb_paired, color = ALT_status)) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.4, outlier.size = 1) +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", 
                     size = 4) +
  xlab("") + 
  ylab("tmb_paired") +
  ggtitle("ALT Status vs tmb_paired") +
  geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
  theme(legend.position = "none") 
ggsave(plot = p, filename = file.path(output_dir, "alt_status_vs_tmb_paired.png"), height = 6, width = 6)

# 4) ALT status vs telomere content
p <- ggplot(plot_data, aes(x = ALT_status, y = t_n_telomere_content, color = ALT_status)) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.4, outlier.size = 1) +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", 
                     size = 4) +
  xlab("") + 
  ylab("tmb_paired") +
  ggtitle("ALT Status vs Telomere Content") +
  geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
  theme(legend.position = "none") 
ggsave(plot = p, filename = file.path(output_dir, "alt_status_vs_telomere_content.png"), height = 6, width = 4)

# 5) ALT telomere content vs Age three groups
plot_df <- output_df %>%
  group_by(age_three_groups) %>%
  mutate(n = n(), age_three_groups = paste0(age_three_groups, "\n(n = ", n, ")"))
plot_df$age_three_groups <- factor(plot_df$age_three_groups,  levels = c("[0,15]\n(n = 42)", "(15,26]\n(n = 20)", "(26,40]\n(n = 8)"))
p <- ggplot(plot_df, aes(x = age_three_groups, y = t_n_telomere_content, color = age_three_groups)) +
  geom_boxplot(coef=0) +
  geom_sina() +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", size = 4, method = "kruskal.test") +
  xlab("\nHARMONY_age_class_derived") + 
  ylab("t_n_telomere_content\n") +
  theme(legend.position = "none") 
p

# 6) ALT telomere content vs Age two groups 
# 2 class
plot_df <- output_df %>%
  group_by(age_two_groups) %>%
  mutate(n = n(), age_two_groups = paste0(age_two_groups, "\n(n = ", n, ")"))
plot_df$age_two_groups <- factor(plot_df$age_two_groups, levels = c("[0,15]\n(n = 42)", "(15,40]\n(n = 28)"))
q <- ggplot(plot_df, aes(x = age_two_groups, y = t_n_telomere_content, color = age_two_groups)) +
  geom_boxplot(coef=0) +
  geom_sina() +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", size = 4, method = "wilcox.test") +
  xlab("\nHARMONY_age_class_derived") + 
  ylab("t_n_telomere_content\n") +
  theme(legend.position = "none") 
q

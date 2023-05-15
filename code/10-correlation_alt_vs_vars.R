# script to correlate ALT status to clinical variables like
# protein clusters, Age (2 and 3 groups), Gender, MSI and TMB
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
})

# output directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
output_dir <- file.path(root_dir, "results", "alt_correlations")
dir.create(output_dir, showWarnings = F, recursive = T)

# function to compute correlations
compute_corr <- function(x){

  # t_n_telomere_content vs MSI
  tel_content_vs_msi_corr = cor.test(x = x$t_n_telomere_content, y = x$MSI_Percent)$estimate
  tel_content_vs_msi_pvalue = cor.test(x = x$t_n_telomere_content, y = x$MSI_Percent)$p.value
  
  # t_n_telomere_content vs TMB
  tel_content_vs_tmb_corr = cor.test(x = x$t_n_telomere_content, y = x$TMB)$estimate
  tel_content_vs_tmb_pvalue = cor.test(x = x$t_n_telomere_content, y = x$TMB)$p.value
  
  # for ALT vs Age (two groups)
  age_two_group_pvalue <- chisq.test(x = factor(x$ALT_status), y = factor(x$age_two_groups))$p.value
  age_two_group_pvalue <- round(age_two_group_pvalue, digits = 3)
  
  # for ALT vs Age (three groups)
  age_three_group_pvalue <- chisq.test(x = factor(x$ALT_status), y = factor(x$age_three_groups))$p.value
  age_three_group_pvalue <- round(age_three_group_pvalue, digits = 3)
  
  # Sex
  gender_pvalue <- chisq.test(x = factor(x$ALT_status), y = factor(x$Gender))$p.value
  gender_pvalue <- round(gender_pvalue, digits = 3)
  
  # Protein cluster
  rdt_cc_pvalue <- chisq.test(x = factor(x$ALT_status), y = factor(x$rdt.cc))$p.value
  rdt_cc_pvalue <- round(rdt_cc_pvalue, digits = 3)
  
  # Protein cluster name
  rdt_name_pvalue <- chisq.test(x = factor(x$ALT_status), y = factor(x$rdt.name))$p.value
  rdt_name_pvalue <- round(rdt_name_pvalue, digits = 3)
  
  # MSI
  # msi_pvalue <- chisq.test(x = factor(x$ALT_status), y = x$MSI_Percent)$p.value
  msi_pvalue <- wilcox.test(formula = x$MSI_Percent ~ factor(x$ALT_status))$p.value
  msi_pvalue <- round(msi_pvalue, digits = 3)
  
  # TMB
  # tmb_pvalue <- chisq.test(x = factor(x$ALT_status), y = x$TMB)$p.value
  tmb_pvalue <- wilcox.test(formula = x$TMB ~ factor(x$ALT_status))$p.value
  tmb_pvalue <- round(tmb_pvalue, digits = 3)

  output_df <- data.frame(tel_content_vs_msi_corr,
                          tel_content_vs_msi_pvalue,
                          tel_content_vs_tmb_corr,
                          tel_content_vs_tmb_pvalue,
                          age_two_group_pvalue, age_three_group_pvalue, 
                          gender_pvalue, 
                          rdt_cc_pvalue, rdt_name_pvalue,
                          msi_pvalue, tmb_pvalue)
  return(output_df)
}

# master histology file
hist_df <- read_tsv(file = file.path(data_dir, "master_histology_hope_cohort.tsv"))

# 1) tumor-normal-paired comparisons
alt_status <- hist_df %>%
  filter(!is.na(msi_paired)) %>%
  dplyr::rename("TMB" = "TMB_paired",
                "MSI_Percent" = "msi_paired")
df <- compute_corr(x = alt_status)
rownames(df) <- "Tumor-Normal"
write_tsv(x = df, file = file.path(output_dir, "alt_correlations_paired.tsv"))

# ALT telomere vs MSI scatter plot
p <- ggplot(alt_status, aes(x = t_n_telomere_content, y = MSI_Percent)) +
  ggpubr::theme_pubr() +
  geom_point(position = "jitter", pch = 21) +
  ggpubr::stat_cor(method = "pearson", color = "red")
ggsave(plot = p, filename = file.path(output_dir, "alt_content_vs_msi_paired.png"))

# ALT status vs MSI boxplot
alt_status <- alt_status %>%
  group_by(ALT_status) %>%
  mutate(n = n(),
         ALT_status_label = paste0(ALT_status, "\n(n = ", n, ")"))
p <- ggplot(alt_status, aes(x = ALT_status_label, y = MSI_Percent, color = ALT_status_label)) +
  ggpubr::theme_pubr() +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.5, outlier.size = 1) +  
  geom_point() +
  ggpubr::stat_compare_means(color = "black", ref.group = ".all.", label = "p.format") + 
  ggpubr::stat_compare_means(color = "red", label.y = 10) +
  theme(legend.position = "none") + xlab("")
p
ggsave(plot = p, filename = file.path(output_dir, "alt_status_vs_msi_paired.png"))

# ALT telomere vs TMB scatter plot
p <- ggplot(alt_status, aes(x = t_n_telomere_content, y = TMB)) +
  ggpubr::theme_pubr() +
  geom_point(position = "jitter", pch = 21) +
  ggpubr::stat_cor(method = "pearson", color = "red")
ggsave(plot = p, filename = file.path(output_dir, "alt_content_vs_tmb_paired.png"))

# ALT status vs TMB boxplot
p <- ggplot(alt_status, aes(x = ALT_status_label, y = TMB, color = ALT_status_label)) +
  ggpubr::theme_pubr() +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.5, outlier.size = 1) +  geom_point() +
  ggpubr::stat_compare_means(color = "black", ref.group = ".all.", label = "p.format") + 
  ggpubr::stat_compare_means(color = "red", label.y = 15) +
  theme(legend.position = "none") + 
  xlab("")
p
ggsave(plot = p, filename = file.path(output_dir, "alt_status_vs_tmb_paired.png"))

# 2) tumor-only 
# alt_status_tumor_only <- hist_df %>%
#   filter(!is.na(msi_tumor_only)) %>%
#   dplyr::rename("TMB" = "TMB_tumor_only",
#                 "MSI_Percent" = "msi_tumor_only")
# df_tumor_only <- compute_corr(x = alt_status_tumor_only)
# rownames(df_tumor_only) <- "Tumor-only"
# 
# # ALT telomere vs MSI scatter plot
# p <- ggplot(alt_status_tumor_only, aes(x = t_n_telomere_content, y = MSI_Percent)) +
#   ggpubr::theme_pubr() +
#   geom_point(position = "jitter", pch = 21) +
#   ggpubr::stat_cor(method = "pearson", color = "red")
# ggsave(plot = p, filename = file.path(output_dir, "alt_content_vs_msi_tumor_only.png"))
# 
# # ALT status vs MSI boxplot
# alt_status_tumor_only <- alt_status_tumor_only %>%
#   group_by(ALT_status) %>%
#   mutate(n = n(),
#          ALT_status_label = paste0(ALT_status, "\n(n = ", n, ")"))
# p <- ggplot(alt_status_tumor_only, aes(x = ALT_status_label, y = MSI_Percent, color = ALT_status_label)) +
#   ggpubr::theme_pubr() +
#   stat_boxplot(geom ='errorbar', width = 0.2) +
#   geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.5, outlier.size = 1) +  geom_point() +
#   ggpubr::stat_compare_means(color = "black", ref.group = ".all.", label = "p.format") + 
#   ggpubr::stat_compare_means(color = "red", label.y = 15) +
#   xlab("") +
#   theme(legend.position = "none")
# ggsave(plot = p, filename = file.path(output_dir, "alt_status_vs_msi_tumor_only.png"))
# 
# # ALT telomere vs TMB scatter plot
# p <- ggplot(alt_status_tumor_only, aes(x = t_n_telomere_content, y = TMB)) +
#   ggpubr::theme_pubr() +
#   geom_point(position = "jitter", pch = 21) +
#   ggpubr::stat_cor(method = "pearson", color = "red")
# ggsave(plot = p, filename = file.path(output_dir, "alt_content_vs_tmb_tumor_only.png"))
# 
# # ALT status vs TMB boxplot
# p <- ggplot(alt_status_tumor_only, aes(x = ALT_status_label, y = TMB, color = ALT_status_label)) +
#   ggpubr::theme_pubr() +
#   stat_boxplot(geom ='errorbar', width = 0.2) +
#   geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.5, outlier.size = 1) +  geom_point() +
#   ggpubr::stat_compare_means(color = "black", ref.group = ".all.", label = "p.format") + 
#   ggpubr::stat_compare_means(color = "red", label.y = 15) +
#   xlab("") +
#   theme(legend.position = "none")
# p
# ggsave(plot = p, filename = file.path(output_dir, "alt_status_vs_tmb_tumor_only.png"))
# 
# # write output
# rbind(df, df_tumor_only) %>%
#   rownames_to_column("Analysis_Type") %>%
#   write_tsv(file = file.path(output_dir, "alt_correlations.tsv"))


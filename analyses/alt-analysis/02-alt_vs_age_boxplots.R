# script to correlate ALT status to clinical variables like
# Protein clusters, Age (2 and 3 groups), Gender and TMB
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(ggforce)
})

# output directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
# input_dir <- file.path(root_dir, "analyses", "master-annotation", "results")
analyses_dir <- file.path(root_dir, "analyses", "alt-analysis")
output_dir <- file.path(analyses_dir, "results")
dir.create(output_dir, showWarnings = F, recursive = T)

# master histology file
hist_df <- read_tsv(file = file.path(root_dir, "data", "Hope-GBM-histologies.tsv"))
alt_status <- read_tsv(file.path(output_dir, "alt_status_aya_hgg.tsv"))
output_df <- hist_df %>%
  filter(!is.na(molecular_subtype)) %>%
  inner_join(alt_status) %>%
  dplyr::select(sample_id, experimental_strategy, ALT_status, t_n_telomere_content, HARMONY_age_class_derived) %>%
  group_by(sample_id) %>%
  mutate(experimental_strategy = toString(experimental_strategy)) %>%
  unique()
output_df <- output_df[grep("WGS", output_df$experimental_strategy),]

# 3 class
plot_df <- output_df %>%
  group_by(HARMONY_age_class_derived) %>%
  mutate(n = n(), HARMONY_age_class_derived = paste0(HARMONY_age_class_derived, "\n(n = ", n, ")"))
plot_df$HARMONY_age_class_derived <- factor(plot_df$HARMONY_age_class_derived, 
                                            levels = c("[0,15]\n(n = 43)", "(15,26]\n(n = 20)", "(26,40]\n(n = 8)"))
p <- ggplot(plot_df, aes(x = HARMONY_age_class_derived, y = t_n_telomere_content, color = HARMONY_age_class_derived)) +
  # geom_violin() +
  geom_boxplot(coef=0) +
  geom_sina() +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", size = 4, method = "kruskal.test") +
  xlab("\nHARMONY_age_class_derived") + 
  ylab("t_n_telomere_content\n") +
  theme(legend.position = "none") 
p

# 2 class
plot_df <- output_df %>%
  mutate(HARMONY_age_class_derived = as.character(HARMONY_age_class_derived)) %>%
  mutate(HARMONY_age_class_derived = ifelse(HARMONY_age_class_derived %in% c("(15,26]", "(26,40]"), "(15,40]", HARMONY_age_class_derived))
plot_df <- plot_df %>%
  group_by(HARMONY_age_class_derived) %>%
  mutate(n = n(), HARMONY_age_class_derived = paste0(HARMONY_age_class_derived, "\n(n = ", n, ")"))
plot_df$HARMONY_age_class_derived <- factor(plot_df$HARMONY_age_class_derived, levels = c("[0,15]\n(n = 43)", "(15,40]\n(n = 28)"))
q <- ggplot(plot_df, aes(x = HARMONY_age_class_derived, y = t_n_telomere_content, color = HARMONY_age_class_derived)) +
  # geom_violin() +
  geom_boxplot(coef=0) +
  geom_sina() +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", size = 4, method = "wilcox.test") +
  xlab("\nHARMONY_age_class_derived") + 
  ylab("t_n_telomere_content\n") +
  theme(legend.position = "none") 
q

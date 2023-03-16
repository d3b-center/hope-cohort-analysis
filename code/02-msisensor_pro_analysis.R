# script to summarise msi-sensor pro output
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
})

# output directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
input_dir <- file.path(data_dir, "msisensor_pro")
results_dir <- file.path(root_dir, "results", "msisensor-pro")
dir.create(results_dir, showWarnings = F, recursive = T)

# read msisensor pro files
lf <- list.files(path = input_dir, pattern = 'msisensor_pro', recursive = T, full.names = T)
df <- lapply(lf, read_tsv)
df <- do.call(rbind,df)
colnames(df)[3] <- "Percent"
df$name <- gsub('.*/', '', lf)
manifest1 <- read_tsv(file.path(input_dir, "manifest_20220927_180139.tsv"))
manifest2 <- read_tsv(file.path(input_dir, "manifest_20230125_101423.tsv"))
manifest <- rbind(manifest1, manifest2)
manifest <- manifest %>% inner_join(df, by = "name")

# updated clinical data from Mateusz
annot <- readr::read_tsv(file.path(data_dir, "hopeonly_clinical_table_011823.tsv"))
manifest <- manifest %>%
  filter(sample_id %in% annot$Sample_ID)

tumor_normal_pair <- read_tsv('data/Tumour_normal_participate.tsv')
manifest <- manifest %>%
  inner_join(tumor_normal_pair, by = c("Kids First Biospecimen ID" = "Tumour")) %>%
  dplyr::select(`Kids First Biospecimen ID`, `Kids First Participant ID`, sample_id, gender, Total_Number_of_Sites, Number_of_Somatic_Sites, Percent) %>%
  unique()
output_df <- manifest %>%
  mutate(Type = ifelse(Percent >= 3.5, "High", "Low")) %>%
  arrange(Percent)
fname <- file.path(results_dir, "hope_cohort_msi_sensor_output.tsv")
write_tsv(output_df, file = fname)

# output_df <- output_df %>%
#   filter(Percent < 1)

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

# modify type
output_df <- output_df %>%
  dplyr::mutate(Type = ifelse(Type == "High", sample_id, ""))

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

# age
plot_data <- output_df %>%
  filter(!is.na(age)) %>%
  mutate(age = as.character(age)) %>%
  group_by(age) %>%
  mutate(n = n()) %>%
  mutate(age = paste0(age, "\n(n = ",n,")")) 
plot_data$age <- factor(plot_data$age, levels = c("[0,15]\n(n = 43)", "(15,26]\n(n = 20)", "(26,40]\n(n = 7)"))
r <- ggplot(plot_data, aes(x = age, y = Percent, color = as.character(age))) +
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
ggsave(filename = file.path(results_dir, "msi_sensor_vs_age.png"), width = 6, height = 6)

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
  geom_text_repel(aes(label = Type), na.rm=TRUE, hjust=0, vjust=0, size = 2, color = "black") +
  ggpubr::theme_pubr(base_size = 10) + ylab("") + 
  stat_compare_means(color = "red", size = 4) +
  xlab("") + 
  ylab("% Microsatellite Instability") +
  ggtitle("% Microsatellite Instability vs. Gender") +
  geom_hline(yintercept = 3.5, linetype = 'dotted', col = 'red') +
  theme(legend.position = "none") 
ggsave(filename = file.path(results_dir, "msi_sensor_vs_gender.png"), width = 6, height = 6)

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


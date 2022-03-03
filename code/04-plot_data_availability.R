# data availability heatmap

suppressPackageStartupMessages({
  library(reshape2)
  library(ggplot2)
  library(tidyverse) 
})

# hope cohort subset of 95 samples to plot
dat <- read.delim('data/hope_cohort_subset.tsv', header = F, col.names = "Sample")
dat$proteomics <- TRUE
dat$phosphoproteomics <- TRUE

# add methylation
methylation <- read.delim('data/methylation_subset.tsv', header = F)
dat$methylation <- dat$Sample %in% methylation$V1

# add RNA-seq and WGS info
annot_info <- read.delim(file.path("results", "annotation.txt"), header = TRUE, check.names = TRUE)
annot_info <- annot_info %>%
  dplyr::select(Sample, Sequencing_Experiment) %>%
  mutate(Sequencing_Experiment = strsplit(as.character(Sequencing_Experiment), ", ")) %>% 
  unnest(Sequencing_Experiment) %>%
  unique()
wgs <- annot_info %>%
  filter(Sequencing_Experiment  == "WGS")
rnaseq <- annot_info %>%
  filter(Sequencing_Experiment  == "RNA-Seq")
dat$WGS <- dat$Sample %in% wgs$Sample
dat$RNASeq <- dat$Sample %in% rnaseq$Sample

# wgs only
wgs_only <- read.delim('data/wgs_tumor_only.tsv', header = F)
dat$WGS_tumor_only <- dat$Sample %in% wgs_only$V1

# order samples
sample_order <- dat %>%
  arrange(proteomics, phosphoproteomics, desc(WGS), desc(WGS_tumor_only), desc(RNASeq), desc(methylation)) %>%
  pull(Sample)

# plot
dat <- melt(dat, id.vars = "Sample", variable.name = "data_type", value.name = "data_availability")
dat$data_type <- factor(dat$data_type, levels=c("methylation", "RNASeq", "WGS_tumor_only", "WGS", "phosphoproteomics", "proteomics"))
dat$Sample <- factor(dat$Sample, levels = sample_order)
dat <- dat %>%
  mutate(label = ifelse(data_availability == TRUE, as.character(data_type), FALSE))
p <- ggplot(dat, aes(Sample, data_type, fill = label)) + 
  geom_tile(colour = "white", aes(height = 1)) + ggpubr::theme_pubr() +
  scale_fill_manual(values = c("FALSE" = "white", 
                               "proteomics" = "#08306B", 
                               "phosphoproteomics" = "#08519C",
                               "WGS" = "#2171B5", 
                               "WGS_tumor_only" = "#4292C6",
                               "RNASeq" = "#6BAED6", 
                               "methylation" = "#9ECAE1")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10)) + 
  xlab("") + ylab("") +
  theme(legend.position = "none")
ggsave(filename = "results/hope_cohort_data_availability.png", width = 15, height = 3)

# updated version
p <- ggplot(dat %>% filter(data_type != "WGS_tumor_only"), aes(Sample, data_type, fill = label)) + 
  geom_tile(colour = "white", aes(height = 1)) + ggpubr::theme_pubr() +
  scale_fill_manual(values = c("FALSE" = "white", 
                               "proteomics" = "#08306B", 
                               "phosphoproteomics" = "#08519C",
                               "WGS" = "#2171B5", 
                               "RNASeq" = "#4292C6", 
                               "methylation" = "#6BAED6")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8)) + 
  xlab("") + ylab("") +
  theme(legend.position = "none")
dat_tmp <- dat %>%
  filter(data_type == "WGS_tumor_only") %>%
  mutate(data_type = "WGS", 
         label = ifelse(label != FALSE | Sample %in% wgs$Sample, "WGS", FALSE))
p + geom_tile(data = dat_tmp, aes(height = 0.5, width = 0.9)) +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.margin = margin(2, 1, 1, 1, "cm"))
ggsave(filename = "results/hope_cohort_data_availability.png", width = 14, height = 3)


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

# order samples
sample_order <- dat %>%
  arrange(proteomics, phosphoproteomics, desc(WGS), desc(RNASeq), desc(methylation)) %>%
  pull(Sample)

# plot
dat <- melt(dat, id.vars = "Sample", variable.name = "data_type", value.name = "data_availability")
dat$data_type <- factor(dat$data_type, levels=c("methylation","RNASeq","WGS", "phosphoproteomics", "proteomics"))
dat$Sample <- factor(dat$Sample, levels = sample_order)
dat <- dat %>%
  mutate(label = ifelse(data_availability == TRUE, as.character(data_type), FALSE))
p <- ggplot(dat, aes(Sample, data_type, fill = label)) + 
  geom_tile(colour = "grey50") + ggpubr::theme_pubr() +
  scale_fill_manual(values = c("FALSE" = "white", 
                               "proteomics" = "#08306B", 
                               "phosphoproteomics" = "#08519C",
                               "WGS" = "#2171B5", 
                               "RNASeq" = "#4292C6", 
                               "methylation" = "#6BAED6")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10)) + 
  xlab("") + ylab("") +
  theme(legend.position = "none")
ggsave(filename = "results/hope_cohort_data_availability.png", width = 15, height = 3)

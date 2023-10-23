# Function: sample availability heatmap

# load libraries
suppressPackageStartupMessages({
  library(reshape2)
  library(ggplot2)
  library(tidyverse) 
})

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data", "v1")
analyses_dir <- file.path(root_dir, "analyses", "data-availability")
output_dir <- file.path(analyses_dir, "plots")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# read histologies
annot <- read_tsv(file.path(data_dir, "Hope-GBM-histologies.tsv"))
annot <- annot %>%
  filter(!is.na(HOPE_diagnosis))
dat <- annot %>% dplyr::select(sample_id) %>% unique()

# proteomics 
dat$Proteomics <- dat$sample_id %in% (annot %>% 
                                        filter(experimental_strategy == "Whole Cell Proteomics") %>% 
                                        pull(sample_id))

# phospho-proteomics
dat$Phosphoproteomics <- dat$sample_id %in% (annot %>% 
                                               filter(experimental_strategy == "Phospho-Proteomics") %>% 
                                               pull(sample_id))

# methylation
dat$Methylation <- dat$sample_id %in% (annot %>% 
                                         filter(experimental_strategy == "Methylation") %>% 
                                         pull(sample_id))

# RNA-seq
dat$RNAseq <- dat$sample_id %in% (annot %>% 
                                    filter(experimental_strategy == "RNA-Seq") %>% 
                                    pull(sample_id))

# add WGS 
cnv = readRDS(file.path(data_dir, "Hope-cnv-controlfreec.rds"))
dat$WGS <- dat$sample_id %in% (annot %>% 
                                 filter(Kids_First_Biospecimen_ID %in% cnv$Kids_First_Biospecimen_ID) %>% 
                                 pull(sample_id))

# add WGS (tumor-only)
cnv_tumor_only = readRDS(file.path(data_dir, "Hope-cnv-controlfreec.rds"))
dat$WGS_tumor_only <- dat$sample_id %in% (annot %>% 
                                            filter(Kids_First_Biospecimen_ID %in% cnv_tumor_only$Kids_First_Biospecimen_ID) %>% 
                                            pull(sample_id))

# snRNASeq
dat$snRNASeq <- dat$sample_id %in% (annot %>% 
                                      filter(experimental_strategy == "snRNA-Seq") %>% 
                                      pull(sample_id))


# order samples
sample_order <- dat %>%
  arrange(desc(Proteomics), 
          desc(Phosphoproteomics), 
          desc(WGS), 
          desc(WGS_tumor_only), 
          desc(RNAseq), 
          desc(Methylation),
          desc(snRNASeq)) %>%
  pull(sample_id)

# plot
dat <- melt(dat, id.vars = "sample_id", variable.name = "data_type", value.name = "data_availability")
dat$data_type <- factor(dat$data_type, levels=c("snRNASeq", "Methylation", "RNAseq", 
                                                "WGS_tumor_only", "WGS", 
                                                "Phosphoproteomics", "Proteomics"))
dat$sample_id <- factor(dat$sample_id, levels = sample_order)
dat <- dat %>%
  mutate(label = ifelse(data_availability == TRUE, as.character(data_type), FALSE))

# generate plot
q <- ggplot(dat %>% filter(!data_type %in% c("WGS_tumor_only")), aes(sample_id, data_type, fill = label)) + 
  geom_tile(colour = "white", aes(height = 1)) + ggpubr::theme_pubr() +
  scale_fill_manual(values = c("FALSE" = "white", 
                               "Proteomics" = "#08306B", 
                               "Phosphoproteomics" = "#08519C",
                               "WGS" = "#2171B5", 
                               "RNAseq" = "#4292C6", 
                               "Methylation" = "#6BAED6",
                               "snRNASeq" = "#9ECAE1")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8)) + 
  xlab("") + ylab("") +
  theme(legend.position = "none")

# add numbers for each datatype
labs = dat %>% 
  filter(label != FALSE,
         label != "WGS_tumor_only") %>%
  group_by(label) %>% 
  summarise(n = n()) %>% 
  mutate(label2 = paste0(label, ' (n = ', n,')')) %>%
  dplyr::select(-c(n))

q <- q + scale_y_discrete(limit = labs$label[match(setdiff(levels(dat$data_type), "WGS_tumor_only"), labs$label)],
                          labels = labs$label2[match(setdiff(levels(dat$data_type), "WGS_tumor_only"), labs$label)])

# now add WGS tumor-only
dat_tmp <- dat %>%
  filter(data_type == "WGS_tumor_only") %>%
  mutate(data_type = "WGS", 
         label = ifelse(label != FALSE | sample_id %in% (annot %>% filter(Kids_First_Biospecimen_ID %in% cnv$Kids_First_Biospecimen_ID) %>% pull(sample_id)), "WGS", FALSE))
q <- q + geom_tile(data = dat_tmp, aes(height = 0.5, width = 0.9), color = "black", fill = "#2171B5") +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.margin = margin(2, 1, 1, 1, "cm"))
q
ggsave(filename = file.path(output_dir, "hope_sample_availability.pdf"), plot = q, width = 12, height = 4)

# data availability heatmap

# load libraries
suppressPackageStartupMessages({
  library(reshape2)
  library(ggplot2)
  library(tidyverse) 
})

# output directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
output_dir <- file.path(root_dir, "results", "data_plots")
dir.create(output_dir, recursive = T, showWarnings = F)

# updated clinical data from Mateusz
dat = readr::read_tsv(file.path(data_dir, "hopeonly_clinical_table_011823.tsv"))
dat <- dat %>%
  dplyr::select(Sample_ID) %>%
  filter(Sample_ID != "7316-3625")

# these two are in 10x but not here
# but Mateusz has asked to add them in
dat <- rbind(dat, data.frame(Sample_ID = c("7316-3158", "7316-1464"))) 

# add proteomics from Nicole's file
proteomics <- read_tsv('data/cluster_data101922.tsv')
dat$Proteomics <- dat$Sample_ID %in% proteomics$id

# add phosphoproteomics from Nicole's file
dat$Phosphoproteomics <- dat$Sample_ID %in% proteomics$id

# add methylation
methylation <- read.delim(file.path(data_dir, "methylation_subset.tsv"), header = F)
dat$Methylation <- dat$Sample_ID %in% methylation$V1

# add RNA-seq from Cavatica manifest
rna_metadata = read_tsv(file.path(data_dir, "manifest", "manifest_20230120_101120_rna.tsv"))
dat$RNAseq <- dat$Sample_ID %in% rna_metadata$sample_id

# add WGS from Cavatica manifest
snv_metadata = read_tsv(file.path(data_dir, "manifest", "manifest_20230120_100627_snv.tsv"))
cnv_metadata = read_tsv(file.path(data_dir, "manifest", "manifest_20230120_095501_cnv.tsv"))
dat$WGS <- dat$Sample_ID %in% cnv_metadata$sample_id

# WGS tumor-only from Cavatica manifest
wgs_tumor_only <- read_tsv(file.path(data_dir, "manifest_tumor_only", "manifest_20230216_162009_snv.tsv"))
dat$WGS_tumor_only <- dat$Sample_ID %in% wgs_tumor_only$sample_id

# add 10x single cell RNAseq (two samples: "7316-1464" and "7316-3158" missing)
single_cell_rnaseq_10x <- readxl::read_xlsx("data/CPTAC_Project Hope cohorts.xlsx", sheet = "Aaron's Manifest")
colnames(single_cell_rnaseq_10x)[10] <- "Sample"
single_cell_rnaseq_10x <- single_cell_rnaseq_10x %>%
  filter(`10X` == "Yes",
         !is.na(`Subject ID`)) %>%
  dplyr::mutate(Type = ifelse(grepl("-A-|-B-|-C-", Sample), "Full", "Half"))
dat$Single_Cell_RNAseq_10x <- dat$Sample_ID %in% (single_cell_rnaseq_10x %>% filter(Type == "Full") %>% pull(`BioSTOR ID`))
dat$Single_Cell_RNAseq_10x_todo <- dat$Sample_ID %in% (single_cell_rnaseq_10x %>% filter(Type == "Half") %>% pull(`BioSTOR ID`))

# add Smart-Seq2 single cell RNAseq (this is from cavatica project)
single_cell_rnaseq_smartseq2 <- readxl::read_xlsx("data/single_cell_smartseq_manifest.xlsx")
single_cell_rnaseq_smartseq2$Name <- gsub("Sample_|-P[1-9]$", "", single_cell_rnaseq_smartseq2$Name)
# this list was given by Mateusz to keep as half-blocks
half_blocks <- c("7316UP-1962", "7316UP-2058", "7316UP-2333", "7316UP-2403", "7316-4842", "7316-4844")
dat$Single_Cell_RNAseq_SmartSeq2 <- dat$Sample_ID %in% setdiff(single_cell_rnaseq_smartseq2$Name, half_blocks)
dat$Single_Cell_RNAseq_SmartSeq2_todo <- dat$Sample_ID %in% half_blocks

# order samples
sample_order <- dat %>%
  arrange(desc(Proteomics), 
          desc(Phosphoproteomics), 
          desc(WGS), 
          desc(WGS_tumor_only), 
          desc(RNAseq), 
          desc(Methylation),
          desc(Single_Cell_RNAseq_SmartSeq2),
          desc(Single_Cell_RNAseq_SmartSeq2_todo),
          desc(Single_Cell_RNAseq_10x),
          desc(Single_Cell_RNAseq_10x_todo)) %>%
  pull(Sample_ID)

# plot
dat <- melt(dat, id.vars = "Sample_ID", variable.name = "data_type", value.name = "data_availability")
dat$data_type <- factor(dat$data_type, levels=c("Single_Cell_RNAseq_10x_todo", "Single_Cell_RNAseq_10x",
                                                "Single_Cell_RNAseq_SmartSeq2_todo", "Single_Cell_RNAseq_SmartSeq2",
                                                "Methylation", "RNAseq", "WGS_tumor_only", "WGS", 
                                                "Phosphoproteomics", "Proteomics"))
dat$Sample_ID <- factor(dat$Sample_ID, levels = sample_order)
dat <- dat %>%
  mutate(label = ifelse(data_availability == TRUE, as.character(data_type), FALSE))
# p <- ggplot(dat, aes(Sample, data_type, fill = label)) +
#   geom_tile(colour = "white", aes(height = 1)) + ggpubr::theme_pubr() +
#   scale_fill_manual(values = c("FALSE" = "white",
#                                "proteomics" = "#08306B",
#                                "phosphoproteomics" = "#08519C",
#                                "WGS" = "#2171B5",
#                                "WGS_tumor_only" = "#4292C6",
#                                "RNASeq" = "#6BAED6",
#                                "methylation" = "#9ECAE1")) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10)) +
#   xlab("") + ylab("") +
#   theme(legend.position = "none")
# ggsave(filename = "results/hope_cohort_data_availability.png", width = 15, height = 3)

# updated version
q <- ggplot(dat %>% filter(!data_type %in% c("WGS_tumor_only", 
                                             "Single_Cell_RNAseq_10x_todo", 
                                             "Single_Cell_RNAseq_SmartSeq2_todo")), 
            aes(Sample_ID, data_type, fill = label)) + 
  geom_tile(colour = "white", aes(height = 1)) + ggpubr::theme_pubr() +
  scale_fill_manual(values = c("FALSE" = "white", 
                               "Proteomics" = "#08306B", 
                               "Phosphoproteomics" = "#08519C",
                               "WGS" = "#2171B5", 
                               "RNAseq" = "#4292C6", 
                               "Methylation" = "#6BAED6",
                               "Single_Cell_RNAseq_SmartSeq2" = "#9ECAE1",
                               "Single_Cell_RNAseq_10x" = "#C6DBEF")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8)) + 
  xlab("") + ylab("") +
  theme(legend.position = "none")

# add WGS tumor only
dat_tmp <- dat %>%
  filter(data_type == "WGS_tumor_only") %>%
  mutate(data_type = "WGS", 
         label = ifelse(label != FALSE | Sample_ID %in% snv_metadata$sample_id, "WGS", FALSE))
q <- q + geom_tile(data = dat_tmp, aes(height = 0.5, width = 0.9)) +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.margin = margin(2, 1, 1, 1, "cm"))

# add single-cell 10x (half-blocks)
dat_tmp <- dat %>%
  filter(data_type == "Single_Cell_RNAseq_10x_todo") %>%
  mutate(data_type = "Single_Cell_RNAseq_10x", 
         label = ifelse(label != FALSE | Sample_ID %in% single_cell_rnaseq_10x$`BioSTOR ID`, "Single_Cell_RNAseq_10x", FALSE))
q <- q + geom_tile(data = dat_tmp, aes(height = 0.5, width = 0.9)) +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.margin = margin(2, 1, 1, 1, "cm"))

# add single-cell smartseq (half-blocks)
dat_tmp <- dat %>%
  filter(data_type == "Single_Cell_RNAseq_SmartSeq2_todo") %>%
  mutate(data_type = "Single_Cell_RNAseq_SmartSeq2", 
         label = ifelse(label != FALSE | Sample_ID %in% single_cell_rnaseq_smartseq2$Name, "Single_Cell_RNAseq_SmartSeq2", FALSE))
q <- q + geom_tile(data = dat_tmp, aes(height = 0.5, width = 0.9)) +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.margin = margin(2, 1, 1, 1, "cm"))
q
ggsave(filename = file.path(output_dir, "hope_cohort_data_availability.pdf"), plot = q, width = 12, height = 4)

# add annotations
# annot <- readxl::read_xlsx('data/clini_m_030722-for_Komal.xlsx')
# annot <- annot %>%
#   dplyr::select(Sample_ID, Diagnosis_demoted,Gender, age.class, Diagnosis.Type_demoted, Sample.annotation, Tumor.Location.condensed3)
# annot <- melt(annot, id.vars = "Sample_ID", variable.name = "data_type", value.name = "data_availability")
# annot <- annot %>%
#   dplyr::rename("Sample" = "Sample_ID") %>%
#   mutate(label = data_availability)
# res <- rbind(annot, dat)
# res$label[res$label == "FALSE"] <- "No_data"
# all_together <- ggplot(res, aes(Sample, data_type, fill = label)) + 
#   geom_tile(colour = "white", aes(height = 1)) + ggpubr::theme_pubr() +
#   scale_fill_manual(values = c("No_data" = "white", 
#                                "proteomics" = "#08306B", 
#                                "phosphoproteomics" = "#08519C",
#                                "WGS" = "#2171B5", 
#                                "WGS_tumor_only" = "#4292C6",
#                                "RNASeq" = "#6BAED6", 
#                                "methylation" = "#9ECAE1",
#                                "Male" = "navy",
#                                "Female" = "deeppink4",
#                                "High-grade glioma/astrocytoma (WHO grade III/IV)"="lightseagreen",
#                                "Astrocytoma;Oligoastrocytoma" = "mediumorchid2",
#                                "Astrocytoma" = "brown2", 
#                                "Glioblastoma" = "orange",
#                                "Ependymoma" = "darkgreen",
#                                "Low-grade glioma/astrocytoma (WHO grade I/II)" = "blue2",
#                                "[0,15]" = "gold",
#                                "(15,40]" = "purple",
#                                "Treatment naive" = "lightgray",
#                                "Post-treatment" = "gray50",
#                                "Post-mortem" = "black",
#                                "Progressive" = "#827397", 
#                                "Primary" = "#d8b9c3",
#                                "Initial CNS Tumor" = "#cee397",
#                                "Recurrence" = "#363062",
#                                "Second Malignancy" = "#005082",
#                                "Cortical" = "magenta",
#                                "Other/Multiple locations/NOS" = "pink",
#                                "Midline" = "purple",
#                                "Cerebellar" = "navy")) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
#         axis.text.y = element_text(size = 8)) + 
#   xlab("") + ylab("")

# linear representation of clinical data
# annot <- readxl::read_xlsx('data/clini_m_030722-for_Komal.xlsx')
# annot <- annot %>%
#   filter(Sample_ID != "7316-4065") %>%
#   dplyr::select(Sample_ID, Diagnosis_demoted, age.class, Gender, Diagnosis.Type_demoted, Sample.annotation, Tumor.Location.condensed3)
# sample_order <- annot %>%
#   arrange(Diagnosis_demoted, age.class, Gender, Diagnosis.Type_demoted, Sample.annotation, Tumor.Location.condensed3) %>%
#   pull(Sample_ID)
# annot <- melt(annot, id.vars = "Sample_ID", variable.name = "data_type", value.name = "data_availability")
# annot <- annot %>%
#   dplyr::rename("Sample" = "Sample_ID") %>%
#   mutate(label = data_availability)
# annot$Sample <- factor(annot$Sample, levels = sample_order)
# annot$data_type <- factor(annot$data_type, levels = rev(levels(annot$data_type)))
# only_clinical <- ggplot(annot, aes(Sample, data_type, fill = label)) + 
#   geom_tile(colour = "white", aes(height = 1)) + ggpubr::theme_pubr() +
#   scale_fill_manual(values = c("High-grade glioma/astrocytoma (WHO grade III/IV)"="lightseagreen",
#                                "Astrocytoma;Oligoastrocytoma" = "mediumorchid2",
#                                "Astrocytoma" = "brown2", 
#                                "Glioblastoma" = "orange",
#                                "Low-grade glioma/astrocytoma (WHO grade I/II)" = "blue2",
#                                "Male" = "navy",
#                                "Female" = "deeppink4",
#                                "[0,15]" = "gold",
#                                "(15,40]" = "purple",
#                                "Progressive" = "#827397", 
#                                "Initial CNS Tumor" = "#cee397",
#                                "Recurrence" = "#363062",
#                                "Second Malignancy" = "#005082",
#                                "Treatment naive" = "lightgray",
#                                "Post-treatment" = "gray50",
#                                "Post-mortem" = "black",
#                                "Cortical" = "magenta",
#                                "Other/Multiple locations/NOS" = "pink",
#                                "Midline" = "purple",
#                                "Cerebellar" = "navy")) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
#         axis.text.y = element_text(size = 8)) + 
#   xlab("") + ylab("") 
# only_clinical
# ggsave(filename = 'results/hope_cohort_data_availability_clinical.png', plot = only_clinical, width = 15, height = 5)

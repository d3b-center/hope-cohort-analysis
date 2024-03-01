suppressPackageStartupMessages({
  library(tidyverse)
})

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analyses_dir <- file.path(root_dir, "analyses", "master_annotation")
output_dir <- file.path(analyses_dir, "results")
dir.create(output_dir, recursive = T, showWarnings = F)

# rna/wgs ids
wgs_ids <- readRDS(file.path(data_dir, "Hope-cnv-controlfreec.rds")) %>%
  pull(Kids_First_Biospecimen_ID) %>%
  unique()
wgs_tumor_only_ids <- readRDS(file.path(data_dir, "Hope-cnv-controlfreec-tumor-only.rds")) %>%
  pull(Kids_First_Biospecimen_ID) %>%
  unique()
wgs_tumor_only_ids <- setdiff(wgs_tumor_only_ids, wgs_ids)
rna_ids <- file.path(data_dir, "Hope-gene-expression-rsem-tpm-collapsed.rds") %>%
  readRDS() %>%
  colnames()

# read histologies
annot <- read_tsv(file.path(data_dir, "Hope-GBM-histologies.tsv"))
annot <- annot %>%
  filter(!is.na(HOPE_diagnosis))

# fix HOPE_diagnosis_type
annot <- annot %>%
  mutate(HOPE_diagnosis_type = case_when(
    HOPE_diagnosis_type %in% c("Recurrent", "recurrent", "Recurrent, residual")  ~ "Recurrence",
    HOPE_diagnosis_type == "Primary" ~ "Initial CNS Tumor",
    TRUE ~ as.character(HOPE_diagnosis_type)))

# fix WHO Grade
annot <- annot %>%
  mutate(HARMONY_WHO.Grade = ifelse(HARMONY_WHO.Grade == "1-2?", NA, HARMONY_WHO.Grade))

# add WHO Grade to Diagnosis 
diagnosis_who_grade_map <- annot %>% 
  filter(!is.na(HOPE_diagnosis)) %>%
  mutate(Diagnosis = gsub(" [(].*", "", HOPE_diagnosis)) %>%
  group_by(Diagnosis, HOPE_diagnosis) %>% 
  summarise(HARMONY_WHO.Grade = toString(as.roman(sort(unique(na.omit(HARMONY_WHO.Grade)))))) %>%
  mutate(HARMONY_WHO.Grade = gsub(", ", "/", HARMONY_WHO.Grade)) %>%
  mutate(HARMONY_WHO.Grade = paste0("(WHO grade ", HARMONY_WHO.Grade, ")")) %>%
  mutate(Diagnosis = paste(Diagnosis, HARMONY_WHO.Grade)) %>%
  dplyr::select(HOPE_diagnosis, Diagnosis)

# add new Diagnosis to annotation 
annot <- annot %>%
  left_join(diagnosis_who_grade_map)

# for tumor-only samples use WGS_tumor_only under experimental_strategy
annot <- annot %>%
  mutate(experimental_strategy = ifelse(Kids_First_Biospecimen_ID %in% wgs_tumor_only_ids, "WGS_tumor_only", experimental_strategy))

# pull TMB (for T/N paired samples)
tmb_paired_output <- read_tsv("analyses/tmb-calculation/results/wgs_paired/snv-mutation-tmb-coding.tsv") %>%
  dplyr::rename("TMB" = "tmb") %>%
  dplyr::select(Tumor_Sample_Barcode, TMB) %>%
  filter(Tumor_Sample_Barcode %in% wgs_ids) %>%
  inner_join(annot, by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::select(sample_id, TMB) %>%
  unique()

# pull TMB (for tumor-only samples)
tmb_tumor_only_output <- read_tsv("analyses/tmb-calculation/results/wgs_tumor_only/snv-mutation-tmb-coding.tsv") %>%
  dplyr::rename("TMB" = "tmb") %>%
  dplyr::select(Tumor_Sample_Barcode, TMB) %>%
  filter(Tumor_Sample_Barcode %in% wgs_tumor_only_ids) %>%
  inner_join(annot, by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::select(sample_id, TMB) %>%
  unique()

# combine both
tmb_output <- rbind(tmb_paired_output, tmb_tumor_only_output)

# add to histology
annot <- annot %>%
  left_join(tmb_output, by = c("sample_id"))

# select columns of interest
annot_info <- annot %>%
  filter(!is.na(CNS_region),
         !is.na(molecular_subtype)) %>%
  group_by(sample_id) %>%
  mutate(Sequencing_Experiment = toString(sort(experimental_strategy))) %>%
  dplyr::rename("Sample" = "sample_id",
                "Annotation" = "HOPE_sample_annotation",
                "Diagnosis_Type" = "HOPE_diagnosis_type",
                "Tumor_Location" = "HOPE_Tumor.Location.condensed",
                "Gender" = "HARMONY_Gender",
                "Age_at_Initial_Diagnosis" = "HOPE_Age_at_Initial_Diagnosis",
                "Age" = "HARMONY_age_class_derived",
                "Molecular_Subtype" = "molecular_subtype",
                "Cancer_Group" = "cancer_group_short") %>%
  dplyr::select(Sample, Sequencing_Experiment, Annotation, Diagnosis, Diagnosis_Type, Tumor_Location, CNS_region, Gender, Age_at_Initial_Diagnosis, Age, Molecular_Subtype, Cancer_Group, TMB) %>%
  unique()

# some samples were removed
removed <- setdiff(annot$sample_id, annot_info$Sample)
removed <- annot %>%
  filter(!is.na(CNS_region)) %>%
  filter(sample_id %in% removed) %>%
  group_by(sample_id) %>%
  mutate(Sequencing_Experiment = toString(sort(experimental_strategy))) %>%
  dplyr::rename("Sample" = "sample_id",
                "Annotation" = "HOPE_sample_annotation",
                "Diagnosis_Type" = "HOPE_diagnosis_type",
                "Tumor_Location" = "HOPE_Tumor.Location.condensed",
                "Gender" = "HARMONY_Gender",
                "Age_at_Initial_Diagnosis" = "HOPE_Age_at_Initial_Diagnosis",
                "Age" = "HARMONY_age_class_derived",
                "Molecular_Subtype" = "molecular_subtype",
                "Cancer_Group" = "cancer_group_short") %>%
  dplyr::select(Sample, Sequencing_Experiment, Annotation, Diagnosis, Diagnosis_Type, Tumor_Location, CNS_region, Gender, Age_at_Initial_Diagnosis, Age, Molecular_Subtype, Cancer_Group, TMB) %>%
  unique()

# combine both
annot_info <- rbind(annot_info, removed)

# Cancer_Group should be PXA where Molecular_Subtype is PXA
annot_info <- annot_info %>%
  mutate(Cancer_Group = ifelse(Molecular_Subtype == "PXA", "PXA", Cancer_Group))

# long names for Cancer_Group
annot_info <- annot_info %>%
  mutate(Cancer_Group = case_when(Cancer_Group == "DHG" ~ "(DHG) Diffuse Hemispheric Glioma",
                                  Cancer_Group == "DMG" ~ "(DMG) Diffuse Midline Glioma",
                                  Cancer_Group == "HGG" ~ "(HGG) High Grade Glioma (not otherwise specified)",
                                  Cancer_Group == "IHG" ~ "(IHG) Infantile Hemispheric Glioma",
                                  Cancer_Group == "PXA" ~ "(PXA) Pleomorphic Xanthoastrocytoma"))

# I think you want a modified master histology file. Can you just tell me what exact columns you/others want so I can create a master 
# data_availability <- read_tsv("analyses/data-availability/results/hope_clinical_data_availability_age_continuous.tsv")
# oncoplot <- read_tsv("analyses/oncoplots/results/annotation_add_tumor_only.txt")
# oncoplot <- oncoplot %>%
#   dplyr::select(-c(Tumor_Location, Diagnosis_Type, Sex))
# final_df <- data_availability %>%
#   left_join(oncoplot)
write_tsv(annot_info, file.path(output_dir, "master_annotation.tsv"))

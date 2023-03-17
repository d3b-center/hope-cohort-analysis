# Author: Komal S. Rathi
# Function: Script to generat Oncogrid matrix/additional files using CNV, SNV, Fusion and Expression data
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
})

# output directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
output_dir <- file.path(root_dir, "results", "oncoplots_two_groups")
dir.create(output_dir, recursive = T, showWarnings = F)

# read driver list from OpenPBTA
brain_goi_list <- read.delim(file = file.path(data_dir, "brain-goi-list-long-gencode-v39.txt"), header = F)

# MMR genes
mmr_genes <- read.delim(file = file.path(data_dir, "mmr_genes.tsv"), header = F)

# combined list
driver_genes <- data.frame(V1 = c(brain_goi_list$V1, mmr_genes$V1)) %>% unique()

# read reference gene lists based on PNOC003
snv <- read.delim(file.path(data_dir, "snv_genes.tsv"), header = F) %>%
  rbind(driver_genes) %>% unique()
fusion <- read.delim(file.path(data_dir, "fusion_genes.tsv"), header = F) %>%
  rbind(driver_genes) %>% unique()
cnv <- read.delim(file.path(data_dir, "cnv_genes.tsv"), header = F) %>%
  rbind(driver_genes) %>% unique()
deg <- read.delim(file.path(data_dir, "deg_genes.tsv"), header = F) %>%
  rbind(driver_genes) %>% unique()

# rna/wgs ids
wgs_ids <- readRDS(file.path(data_dir, "merged_files", "cnv_merged.rds")) %>%
  pull(Kids_First_Biospecimen_ID) %>%
  unique()
rna_ids <- file.path(data_dir, "merged_files", "gene-expression-rsem-tpm-collapsed.rds") %>%
  readRDS() %>%
  colnames()

# histology
manifest <- list.files(path = file.path(data_dir, "manifest"), pattern = "manifest.*.tsv", full.names = T)
manifest <- lapply(manifest, FUN = function(x) readr::read_tsv(x))
manifest <- data.table::rbindlist(manifest)
colnames(manifest) <- gsub(" ", "_", colnames(manifest))
manifest <- manifest %>%
  filter(Kids_First_Biospecimen_ID %in% c(rna_ids, wgs_ids)) %>%
  mutate(experimental_strategy = ifelse(Kids_First_Biospecimen_ID %in% rna_ids, "RNA-Seq", "WGS")) %>%
  dplyr::mutate(Sample = sample_id,
                Sequencing_Experiment = experimental_strategy) %>%
  dplyr::select(Kids_First_Biospecimen_ID, Sample, Sequencing_Experiment) %>%
  unique()

# annotation with sample and experimental strategy
annot_info <- manifest %>%
  dplyr::select(-c(Kids_First_Biospecimen_ID)) %>%
  group_by(Sample) %>%
  mutate(Sequencing_Experiment = toString(sort(Sequencing_Experiment))) %>%
  unique()

# get clinical variables from Mateusz's file
hist <- readr::read_tsv(file.path(data_dir, "hopeonly_clinical_table_011823.tsv"))
hist$diagnosis_type[hist$diagnosis_type == "Recurrent"] = "Recurrence"
hist$diagnosis_type[hist$diagnosis_type == "recurrent"] = "Recurrence"
hist$diagnosis_type[hist$diagnosis_type == "Recurrent, residual"] = "Recurrence"
hist$diagnosis_type[hist$diagnosis_type == "Primary"] = "Initial CNS Tumor"

hist <- hist %>%
  dplyr::select(Sample_ID, diagnosis, diagnosis_type, Gender) %>%
  dplyr::rename("Tumor_Descriptor" = "diagnosis_type",
                "Integrated_Diagnosis" = "diagnosis",
                "Sex" = "Gender",
                "Sample" = "Sample_ID")

# get Age from Nicole's file
cluster_data <- readxl::read_xlsx(file.path(data_dir, "clini_m_030722-for_Komal.xlsx"))
hist <- hist %>%
  inner_join(cluster_data %>% dplyr::select(age.class, id), by = c("Sample" = "id")) %>%
  dplyr::rename("Age" = "age.class")

# hist1 <- readxl::read_xlsx('data/CPTAC Cohort 2 (100) Clinical Data Manifest 2021-06-14-corr.xlsx', sheet = 2) # upenn
# hist2 <- readxl::read_xlsx('data/CPTAC Cohort 2 (100) Clinical Data Manifest 2021-06-14-corr.xlsx', sheet = 3) # cbtn
# hist <- rbind(hist1, hist2)
# hist <- hist %>%
#   filter(`Sequenced CBTN Specimen Group ID` %in% annot_info$Sample) %>%
#   mutate(Sample = `Sequenced CBTN Specimen Group ID`) %>%
#   group_by(Sample) %>%
#   mutate(Tumor_Descriptor = toString(unique(`Diagnosis Type`)),
#          Integrated_Diagnosis = toString(unique(Diagnosis)),
#          Sex = toString(unique(Gender))) %>%
#   dplyr::select(Sample, Tumor_Descriptor, Integrated_Diagnosis, Sex) %>%
#   unique()
# hist$Tumor_Descriptor[hist$Tumor_Descriptor == "recurrent"] <- "Recurrence"

# add Age from Nicole's file
# cluster_data <- read_tsv('data/cluster_data101922.tsv')
# cluster_data <- cluster_data %>%
#   dplyr::select(id, age) 
# hist <- hist %>%
#   left_join(cluster_data, by = c("Sample" = "id")) %>%
#   dplyr::rename("Age" = "age")

# add TMB and MSI-sensor information (TMB from consensus maf file)
tmb_msi <- read_tsv('results/msisensor-pro/msi_output_merged.tsv') %>%
  dplyr::select(sample_id, MSI_Percent, TMB) 

# add to histology
hist <- hist %>%
  left_join(tmb_msi, by = c("Sample" = "sample_id")) %>% 
  distinct(Sample, .keep_all = T)

# add to annot
annot_info <- hist %>%
  inner_join(annot_info, by = "Sample")
annot_info <- annot_info %>%
  dplyr::select(MSI_Percent, TMB, Sample, Tumor_Descriptor, Integrated_Diagnosis, Sex, Age, Sequencing_Experiment)
write.table(annot_info, file = file.path(output_dir, "annotation.txt"), quote = F, sep = "\t", row.names = F)

# 1. get degene info PNOC008 patientss vs GTEx Brain
genes_df <- readRDS(file.path("results/hope_cohort_vs_gtex_brain_degs_edgeR.rds"))
deg_genes <- genes_df %>%
  dplyr::mutate(label = ifelse(diff_expr == "up", "OVE", "UNE"),
         Gene_name = gene_symbol,
         sample_name = sample) %>%
  filter(Gene_name %in% deg$V1) %>%
  inner_join(manifest, by = c("sample_name" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::select(Sample, Gene_name, label) %>%
  unique()

# 2. get cnv info
cnv_genes <- readRDS(file.path("data/merged_files/cnv_merged.rds"))
# tumor_only_cnv_genes <- readRDS(file.path("data/merged_files/cnv_tumor_only_merged.rds"))
# cnv_genes <- plyr::rbind.fill(cnv_genes, tumor_only_cnv_genes)
cnv_genes <- cnv_genes %>%
  dplyr::mutate(label = ifelse(status == "Gain", "GAI", "LOS")) %>%
  filter(hgnc_symbol %in% cnv$V1) %>%
  dplyr::mutate(Gene_name = hgnc_symbol) %>%
  inner_join(manifest, by = c("Kids_First_Biospecimen_ID")) %>%
  dplyr::select(Sample, Gene_name, label) %>%
  unique()

# 3. get snv info
mut_genes <- readRDS(file.path("data/merged_files/snv_merged.rds"))
# tumor_only_mut_genes <- readRDS(file.path("data/merged_files/snv_tumor_only_merged.rds"))
# mut_genes <- plyr::rbind.fill(mut_genes, tumor_only_mut_genes)

# # check if the wgs tumor-only match the snv calls
# wgs_only <- read.delim('data/wgs_tumor_only.tsv', header = F)
# wgs_only_samples <- manifest %>% 
#   filter(Sample %in% wgs_only$V1)
# mut_genes_tumor_only <- mut_genes %>%
#   filter(Kids_First_Biospecimen_ID %in% wgs_only_samples$Kids_First_Biospecimen_ID) %>%
#   nrow()
# if(mut_genes_tumor_only == 0){
#   print("WGS tumor only has no SNV data")
# }

# filter to variant classification of interest
mut_genes <- mut_genes %>%
  filter(Variant_Classification %in% c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "Splice_Site")) %>%
  dplyr::mutate(label = case_when(Variant_Classification %in% "Missense_Mutation" ~ "MIS",
                                  Variant_Classification %in% "Nonsense_Mutation" ~ "NOS",
                                  Variant_Classification %in% "Frame_Shift_Del" ~ "FSD",
                                  Variant_Classification %in% "Frame_Shift_Ins" ~ "FSI",
                                  Variant_Classification %in% "In_Frame_Del" ~ "IFD",
                                  Variant_Classification %in% "Splice_Site" ~ "SPS")) %>%
  filter(Hugo_Symbol %in% snv$V1) %>%
  mutate(Gene_name = Hugo_Symbol,
         sample_name = Kids_First_Biospecimen_ID) %>%
  inner_join(manifest, by = c("sample_name" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::select(Sample, Gene_name, label) %>%
  unique()

# 4. get fusion info
fus_genes <- readRDS(file.path("data/merged_files/fusions_merged.rds"))
fus_genes <- fus_genes %>%
  unite(., col = "Gene_name", Gene1A, Gene1B, Gene2A, Gene2B, na.rm = TRUE, sep = ",") %>%
  mutate(Gene_name = strsplit(as.character(Gene_name), ",")) %>% 
  unnest(Gene_name) %>%
  unique()
fus_genes <- fus_genes %>%
  dplyr::mutate(label = "FUS") %>%
  filter(Gene_name %in% fusion$V1) %>%
  inner_join(manifest, by = c("Kids_First_Biospecimen_ID")) %>%
  dplyr::select(Sample, Gene_name, label) %>%
  unique()

# combine fus + snv
snv_fus <- rbind(mut_genes, fus_genes)

# combine deg + cnv
cnv_deg <- rbind(cnv_genes, deg_genes)

# uniquify rows
snv_fus <- snv_fus %>%
  group_by(Sample, Gene_name) %>%
  dplyr::summarise(label = paste0(label, collapse = ';'))
cnv_deg <- cnv_deg %>%
  group_by(Sample, Gene_name) %>%
  dplyr::summarise(label = paste0(label, collapse = ';'))

# convert to matrix
snv_fus <- snv_fus %>%
  spread(key = Gene_name, value = 'label') %>%
  column_to_rownames('Sample')
cnv_deg <- cnv_deg %>%
  spread(key = Gene_name, value = 'label') %>%
  column_to_rownames('Sample')

# add an * to common genes 
colnames(cnv_deg) <- ifelse(colnames(cnv_deg) %in% colnames(snv_fus), paste0(colnames(cnv_deg),'*'), colnames(cnv_deg))

# merge both matrices
oncogrid_mat <- snv_fus %>%
  rownames_to_column('Sample') %>%
  full_join(cnv_deg %>%
              rownames_to_column('Sample'), by = "Sample")
write.table(oncogrid_mat, file = file.path(output_dir, "oncoprint.txt"), quote = F, sep = "\t", row.names = F)

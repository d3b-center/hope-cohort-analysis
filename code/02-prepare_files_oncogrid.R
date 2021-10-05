# Author: Komal S. Rathi
# Function: Script to generat Oncogrid matrix/additional files using CNV, SNV, Fusion and Expression data
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
})

# read reference gene lists from OMPARE/PNOC003
oncogrid_path_input <- '~/Projects/OMPARE/data/oncogrid/input'
snv <- read.delim(file.path(oncogrid_path_input, "snv-genes"), header = F)
fusion <- read.delim(file.path(oncogrid_path_input, "fusion_genes"), header = F)
cnv <- read.delim(file.path(oncogrid_path_input, "copy_number_gene"), header = F)
deg <- read.delim(file.path(oncogrid_path_input, "all_cnv_tgen_genes"), header = F)

# histology
hist <- read.csv('data/1633446748485-manifest.csv') 
annot_info <- hist %>%
  dplyr::select(sample_id, Kids.First.Biospecimen.ID) %>%
  unique()

# 1. get degene info PNOC008 patientss vs GTEx Brain
genes_df <- readRDS(file.path("results/hope_cohort_vs_gtex_brain_degs.rds"))
deg_genes <- genes_df %>%
  dplyr::mutate(label = ifelse(diff_expr == "up", "OVE", "UNE"),
         Gene_name = gene_symbol,
         sample_name = sample) %>%
  filter(Gene_name %in% deg$V1) %>%
  inner_join(annot_info, by = c("sample_name" = "Kids.First.Biospecimen.ID")) %>%
  dplyr::select(sample_id, Gene_name, label) %>%
  unique()

# 2. get cnv info
cnv_genes <- readRDS(file.path("data/merged_files/cnv_filtered.rds"))
cnv_genes <- cnv_genes %>%
  dplyr::mutate(label = ifelse(status == "Gain", "GAI", "LOS")) %>%
  filter(hgnc_symbol %in% cnv$V1) %>%
  dplyr::mutate(Gene_name = hgnc_symbol) %>%
  inner_join(annot_info, by = c("sample_name" = "Kids.First.Biospecimen.ID")) %>%
  dplyr::select(sample_id, Gene_name, label) %>%
  unique()

# 3. get snv info
mut_genes <- readRDS(file.path("data/merged_files/consensus_mutation_filtered.rds"))
mut_genes <- mut_genes %>%
  filter(!Variant_Classification %in% c("3'Flank", "5'Flank", "3'UTR", "5'UTR", "IGR", "Intron", "RNA")) %>%
  dplyr::mutate(label = case_when(Variant_Classification %in% "Missense_Mutation" ~ "MIS",
                                  Variant_Classification %in% "Nonsense_Mutation" ~ "NOS",
                                  Variant_Classification %in% "Frame_Shift_Del" ~ "FSD",
                                  Variant_Classification %in% "Frame_Shift_Ins" ~ "FSI",
                                  Variant_Classification %in% "In_Frame_Del" ~ "IFD",
                                  Variant_Classification %in% "Splice_Site" ~ "SPS")) %>%
  filter(Hugo_Symbol %in% snv$V1) %>%
  mutate(Gene_name = Hugo_Symbol,
         sample_name = Tumor_Sample_Barcode) %>%
  inner_join(annot_info, by = c("sample_name" = "Kids.First.Biospecimen.ID")) %>%
  dplyr::select(sample_id, Gene_name, label) %>%
  unique()

# 4. get fusion info
fus_genes <- readRDS(file.path("data/merged_files/fusions_filtered.rds"))
fus_genes <- fus_genes %>%
  dplyr::mutate(label = "FUS") %>%
  filter(Gene %in% fusion$V1) %>%
  dplyr::mutate(Gene_name = Gene,
         sample_name = Kids.First.Biospecimen.ID) %>%
  inner_join(annot_info, by = c("sample_name" = "Kids.First.Biospecimen.ID")) %>%
  dplyr::select(sample_id, Gene_name, label) %>%
  unique()

# combine fus + snv
snv_fus <- rbind(mut_genes, fus_genes)

# combine deg + cnv
cnv_deg <- rbind(cnv_genes, deg_genes)

# uniquify rows
snv_fus <- snv_fus %>%
  group_by(sample_id, Gene_name) %>%
  dplyr::summarise(label = paste0(label, collapse = ';'))
cnv_deg <- cnv_deg %>%
  group_by(sample_id, Gene_name) %>%
  dplyr::summarise(label = paste0(label, collapse = ';'))

# convert to matrix
snv_fus <- snv_fus %>%
  spread(key = Gene_name, value = 'label') %>%
  column_to_rownames('sample_id')
cnv_deg <- cnv_deg %>%
  spread(key = Gene_name, value = 'label') %>%
  column_to_rownames('sample_id')

# add an * to common genes 
colnames(cnv_deg) <- ifelse(colnames(cnv_deg) %in% colnames(snv_fus), paste0(colnames(cnv_deg),'*'), colnames(cnv_deg))

# merge both matrices
oncogrid_mat <- snv_fus %>%
  rownames_to_column('Sample') %>%
  full_join(cnv_deg %>%
              rownames_to_column('Sample'), by = "Sample")
write.table(oncogrid_mat, file = file.path("results", "oncoprint.txt"), quote = F, sep = "\t", row.names = F)

# annotation
annot_info <- hist %>%
  group_by(sample_id) %>%
  dplyr::mutate(Sample = sample_id,
                Sequencing_Experiment = toString(unique(sort(experimental_strategy))),
                Tumor_Descriptor = toString(unique(Tumor.Descriptor)),
                Integrated_Diagnosis = toString(unique(disease_type)),
                Gender = toString(unique(gender))) %>%
  dplyr::select(Sample, Sequencing_Experiment, Tumor_Descriptor, Integrated_Diagnosis, Gender) %>%
  unique()
annot_info$sample_id <- NULL
write.table(annot_info, file = file.path("results", "annotation.txt"), quote = F, sep = "\t", row.names = F)

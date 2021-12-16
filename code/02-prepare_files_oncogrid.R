# Author: Komal S. Rathi
# Function: Script to generat Oncogrid matrix/additional files using CNV, SNV, Fusion and Expression data
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
})

# read driver list from OpenPBTA
drivers <- read.delim(file.path("data", "brain-goi-list-long.txt"), header = F)

# read reference gene lists based on PNOC003
snv <- read.delim(file.path("data", "snv_genes.tsv"), header = F) %>%
  rbind(drivers) %>% unique()
fusion <- read.delim(file.path("data", "fusion_genes.tsv"), header = F) %>%
  rbind(drivers) %>% unique()
cnv <- read.delim(file.path("data", "cnv_genes.tsv"), header = F) %>%
  rbind(drivers) %>% unique()
deg <- read.delim(file.path("data", "deg_genes.tsv"), header = F) %>%
  rbind(drivers) %>% unique()

# rna/wgs ids
wgs_ids <- readRDS('data/merged_files/cnv_filtered.rds') %>%
  pull(sample_name) %>%
  unique()
rna_ids <- unique(colnames(readRDS('data/merged_files/gene-expression-rsem-tpm-collapsed.rds')))

# histology
manifest <- list.files(path = 'data/', pattern = "csv", full.names = T)
manifest <- lapply(manifest, FUN = function(x) read.csv(x))
manifest <- data.table::rbindlist(manifest)
manifest$experimental_strategy[manifest$experimental_strategy == ""] <- "WGS"
manifest <- manifest %>%
  filter(Kids.First.Biospecimen.ID %in% c(rna_ids, wgs_ids)) %>%
  mutate(experimental_strategy = ifelse(Kids.First.Biospecimen.ID %in% rna_ids, "RNA-Seq", "WGS")) %>%
  dplyr::mutate(Sample = sample_id,
                Sequencing_Experiment = experimental_strategy) %>%
  dplyr::select(Kids.First.Biospecimen.ID, Sample, Sequencing_Experiment) %>%
  unique()

# annotation with sample and experimental strategy
annot_info <- manifest %>%
  dplyr::select(-c(Kids.First.Biospecimen.ID)) %>%
  group_by(Sample) %>%
  mutate(Sequencing_Experiment = toString(sort(Sequencing_Experiment))) %>%
  unique()

# from Mateusz
hist1 <- readxl::read_xlsx('data/CPTAC Cohort 2 (100) Clinical Data Manifest 2021-06-14-corr.xlsx', sheet = 2) # upenn
hist2 <- readxl::read_xlsx('data/CPTAC Cohort 2 (100) Clinical Data Manifest 2021-06-14-corr.xlsx', sheet = 3) # cbtn
hist <- rbind(hist1, hist2)
hist <- hist %>%
  filter(`Sequenced CBTN Specimen Group ID` %in% annot_info$Sample) %>%
  mutate(Sample = `Sequenced CBTN Specimen Group ID`) %>%
  group_by(Sample) %>%
  mutate(Age = toString(unique(sort(`Age at Specimen Diagnosis`))),
         Tumor_Descriptor = toString(unique(`Diagnosis Type`)),
         Integrated_Diagnosis = toString(unique(Diagnosis)),
         Gender = toString(unique(Gender))) %>%
  dplyr::select(Sample, Tumor_Descriptor, Integrated_Diagnosis, Age, Gender) %>%
  unique()
hist$Tumor_Descriptor[hist$Tumor_Descriptor == "recurrent"] <- "Recurrence"
hist$Age <- as.numeric(hist$Age)/365
hist$Age <- ifelse(hist$Age > 0 & hist$Age <= 15, "0-15", "15-35")

annot_info <- hist %>%
  inner_join(annot_info, by = "Sample")
annot_info <- annot_info %>%
  dplyr::select(Sample, Tumor_Descriptor, Integrated_Diagnosis, Gender, Age, Sequencing_Experiment)
write.table(annot_info, file = file.path("results", "annotation.txt"), quote = F, sep = "\t", row.names = F)

# 1. get degene info PNOC008 patientss vs GTEx Brain
genes_df <- readRDS(file.path("results/hope_cohort_vs_gtex_brain_degs_edgeR.rds"))
deg_genes <- genes_df %>%
  dplyr::mutate(label = ifelse(diff_expr == "up", "OVE", "UNE"),
         Gene_name = gene_symbol,
         sample_name = sample) %>%
  filter(Gene_name %in% deg$V1) %>%
  inner_join(manifest, by = c("sample_name" = "Kids.First.Biospecimen.ID")) %>%
  dplyr::select(Sample, Gene_name, label) %>%
  unique()

# 2. get cnv info
cnv_genes <- readRDS(file.path("data/merged_files/cnv_filtered.rds"))
cnv_genes <- cnv_genes %>%
  dplyr::mutate(label = ifelse(status == "Gain", "GAI", "LOS")) %>%
  filter(hgnc_symbol %in% cnv$V1) %>%
  dplyr::mutate(Gene_name = hgnc_symbol) %>%
  inner_join(manifest, by = c("sample_name" = "Kids.First.Biospecimen.ID")) %>%
  dplyr::select(Sample, Gene_name, label) %>%
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
  inner_join(manifest, by = c("sample_name" = "Kids.First.Biospecimen.ID")) %>%
  dplyr::select(Sample, Gene_name, label) %>%
  unique()

# 4. get fusion info
fus_genes <- readRDS(file.path("data/merged_files/fusions_filtered.rds"))
fus_genes <- fus_genes %>%
  dplyr::mutate(label = "FUS") %>%
  filter(Gene %in% fusion$V1) %>%
  dplyr::mutate(Gene_name = Gene,
         sample_name = Kids.First.Biospecimen.ID) %>%
  inner_join(manifest, by = c("sample_name" = "Kids.First.Biospecimen.ID")) %>%
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
write.table(oncogrid_mat, file = file.path("results", "oncoprint.txt"), quote = F, sep = "\t", row.names = F)



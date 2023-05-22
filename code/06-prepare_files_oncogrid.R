# Author: Komal S. Rathi
# Function: Script to generat Oncogrid matrix/additional files using CNV, SNV, Fusion and Expression data
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(msigdbr)
})

# output directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
output_dir <- file.path(root_dir, "results")
dir.create(output_dir, recursive = T, showWarnings = F)

# read driver list from OT
brain_goi_list <- read.delim(file = file.path(data_dir, "hgat_goi_list.tsv"), header = T)

# MMR genes
geneset_db <- msigdbr::msigdbr(category = "C2", subcategory = "KEGG")
geneset_db <- geneset_db[grep("MISMATCH_REPAIR|KEGG_BASE_EXCISION_REPAIR|KEGG_HOMOLOGOUS_RECOMBINATION", geneset_db$gs_name),]
mmr_genes <- unique(geneset_db$human_gene_symbol)

# combined list
driver_genes <- data.frame(V1 = c(brain_goi_list$HGAT, mmr_genes)) %>% unique()

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

# read master histology
hist_df <- read_tsv("data/master_histology_hope_cohort.tsv")

# format diagnosis type 
hist_df$diagnosis_type[hist_df$diagnosis_type == "Recurrent"] = "Recurrence"
hist_df$diagnosis_type[hist_df$diagnosis_type == "recurrent"] = "Recurrence"
hist_df$diagnosis_type[hist_df$diagnosis_type == "Recurrent, residual"] = "Recurrence"
hist_df$diagnosis_type[hist_df$diagnosis_type == "Primary"] = "Initial CNS Tumor"

# pull sequencing experiment information from cavatica manifests
manifest <- list.files(path = file.path(data_dir, "manifest"), pattern = "manifest.*.tsv", full.names = T)
manifest <- grep("msi|cram|fusion", manifest, invert = T, value = T)
manifest <- lapply(manifest, FUN = function(x) readr::read_tsv(x))
manifest <- data.table::rbindlist(manifest)
colnames(manifest) <- gsub(" ", "_", colnames(manifest))
manifest <- manifest %>%
  filter(Kids_First_Biospecimen_ID %in% c(rna_ids, wgs_ids),
         sample_id %in% hist_df$Sample_ID) %>%
  mutate(experimental_strategy = ifelse(Kids_First_Biospecimen_ID %in% rna_ids, "RNA-Seq", "WGS")) %>%
  dplyr::mutate(Sample_ID = sample_id,
                Sequencing_Experiment = experimental_strategy) %>%
  dplyr::select(Kids_First_Biospecimen_ID, Sample_ID, Sequencing_Experiment) %>%
  unique()

# annotation with sample and experimental strategy
seq_info <- manifest %>%
  dplyr::select(-c(Kids_First_Biospecimen_ID)) %>%
  group_by(Sample_ID) %>%
  mutate(Sequencing_Experiment = toString(sort(Sequencing_Experiment))) %>%
  filter(Sample_ID %in% hist_df$Sample_ID) %>%
  unique()

# select columns of interest
annot_info <- hist_df %>%
  left_join(seq_info, by = "Sample_ID") %>%
  dplyr::rename("Tumor_Descriptor" = "diagnosis_type",
                "Integrated_Diagnosis" = "diagnosis",
                "Sex" = "Gender",
                "Sample" = "Sample_ID",
                "MSI" = "msi_paired",
                "TMB" = "TMB_paired",
                "Age_Two_Groups" = "age_two_groups",
                "Age_Three_Groups" = "age_three_groups") %>%
  dplyr::select(Sample, Tumor_Descriptor, Integrated_Diagnosis, Sex, Sequencing_Experiment, MSI, TMB, Age_Two_Groups, Age_Three_Groups)

# save annotation file
write.table(annot_info, file = file.path(output_dir, "annotation.txt"), quote = F, sep = "\t", row.names = F)

# 1. get degene info vs GTEx Brain
genes_df <- readRDS(file.path("results/hope_cohort_vs_gtex_brain_degs.rds"))
deg_genes <- genes_df %>%
  dplyr::mutate(label = ifelse(diff_expr == "up", "OVE", "UNE")) %>%
  dplyr::rename("Gene_name" = "genes") %>%
  filter(Gene_name %in% deg$V1) %>%
  inner_join(manifest, by = c("sample" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::select(Sample_ID, Gene_name, label) %>%
  unique()

# 2. get cnv info
cnv_genes <- readRDS(file.path("data/merged_files/cnv_merged.rds"))
cnv_genes <- cnv_genes %>%
  dplyr::mutate(label = ifelse(status == "Gain", "GAI", "LOS")) %>%
  filter(hgnc_symbol %in% cnv$V1) %>%
  dplyr::rename("Gene_name" = "hgnc_symbol") %>%
  inner_join(manifest, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::select(Sample_ID, Gene_name, label) %>%
  unique()

# 3. get snv info
mut_genes <- readRDS(file.path("data/merged_files/snv_merged.rds"))

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
  dplyr::rename("Gene_name" = "Hugo_Symbol") %>%
  inner_join(manifest, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::select(Sample_ID, Gene_name, label) %>%
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
  dplyr::select(Sample_ID, Gene_name, label) %>%
  unique()

# combine fus + snv
snv_fus <- rbind(mut_genes, fus_genes)

# combine deg + cnv
cnv_deg <- rbind(cnv_genes, deg_genes)

# uniquify rows
snv_fus <- snv_fus %>%
  group_by(Sample_ID, Gene_name) %>%
  dplyr::summarise(label = paste0(label, collapse = ';'))
cnv_deg <- cnv_deg %>%
  group_by(Sample_ID, Gene_name) %>%
  dplyr::summarise(label = paste0(label, collapse = ';'))

# convert to matrix
snv_fus <- snv_fus %>%
  spread(key = Gene_name, value = 'label') %>%
  column_to_rownames('Sample_ID')
cnv_deg <- cnv_deg %>%
  spread(key = Gene_name, value = 'label') %>%
  column_to_rownames('Sample_ID')

# add an * to common genes 
colnames(cnv_deg) <- ifelse(colnames(cnv_deg) %in% colnames(snv_fus), paste0(colnames(cnv_deg),'*'), colnames(cnv_deg))

# merge both matrices
oncogrid_mat <- snv_fus %>%
  rownames_to_column('Sample') %>%
  full_join(cnv_deg %>%
              rownames_to_column('Sample'), by = "Sample")

print("Samples missing from oncogrid:")
print(setdiff(hist_df$Sample_ID, oncogrid_mat$Sample))
# "7316-1455"  "7316-942"   "7316-24"    "7316-255"   "7316-1889"  "7316-212"   "7316UP-904"

# save output matrix
oncogrid_mat <- oncogrid_mat %>%
  filter(Sample %in% annot_info$Sample)
write.table(oncogrid_mat, file = file.path(output_dir, "oncoprint.txt"), quote = F, sep = "\t", row.names = F)

# Author: Komal S. Rathi
# Function: Script to generate Oncogrid matrix/additional files using CNV, SNV, Fusion and Expression data
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(msigdbr)
})

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analyses_dir <- file.path(root_dir, "analyses", "oncoplots")
input_dir <- file.path(analyses_dir, "input")
output_dir <- file.path(analyses_dir, "results")
dir.create(output_dir, recursive = T, showWarnings = F)

# read driver list from OT
brain_goi_list <- read.delim(file = file.path(input_dir, "hgat_goi_list.tsv"), header = T)

# MMR genes
geneset_db <- msigdbr::msigdbr(category = "C2", subcategory = "KEGG")
geneset_db <- geneset_db[grep("MISMATCH_REPAIR|KEGG_BASE_EXCISION_REPAIR|KEGG_HOMOLOGOUS_RECOMBINATION", geneset_db$gs_name),]
mmr_genes <- unique(geneset_db$human_gene_symbol)

# combined list
driver_genes <- data.frame(V1 = c(brain_goi_list$HGAT, mmr_genes)) %>% unique()

# read reference gene lists based on PNOC003
snv <- read.delim(file.path(input_dir, "snv_genes.tsv"), header = F) %>%
  rbind(driver_genes) %>% unique()
fusion <- read.delim(file.path(input_dir, "fusion_genes.tsv"), header = F) %>%
  rbind(driver_genes) %>% unique()
cnv <- read.delim(file.path(input_dir, "cnv_genes.tsv"), header = F) %>%
  rbind(driver_genes) %>% unique()
deg <- read.delim(file.path(input_dir, "deg_genes.tsv"), header = F) %>%
  rbind(driver_genes) %>% unique()

# rna/wgs ids
wgs_ids <- readRDS(file.path(data_dir, "Hope-cnv-controlfreec.rds")) %>%
  pull(biospecimen_id) %>%
  unique()
wgs_tumor_only_ids <- readRDS(file.path(data_dir, "Hope-cnv-controlfreec-tumor-only.rds")) %>%
  pull(biospecimen_id) %>%
  unique()
wgs_tumor_only_ids <- setdiff(wgs_tumor_only_ids, wgs_ids)
rna_ids <- file.path(data_dir, "Hope-gene-expression-rsem-tpm-collapsed.rds") %>%
  readRDS() %>%
  colnames()

# read master histology
hist_df <- read_tsv(file.path(data_dir, "master_histology_hope_cohort.tsv"))

# format diagnosis type 
hist_df$diagnosis_type[hist_df$diagnosis_type == "Recurrent"] = "Recurrence"
hist_df$diagnosis_type[hist_df$diagnosis_type == "recurrent"] = "Recurrence"
hist_df$diagnosis_type[hist_df$diagnosis_type == "Recurrent, residual"] = "Recurrence"
hist_df$diagnosis_type[hist_df$diagnosis_type == "Primary"] = "Initial CNS Tumor"

# pull sequencing experiment information from cavatica manifests
manifest <- list.files(path = file.path("analyses", "merge-files", "input", "manifest"), pattern = "manifest.*.tsv", full.names = T)
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

# add tumor-only for remaining ids
manifest_tumor_only <- list.files(path = file.path("analyses", "merge-files", "input", "manifest_tumor_only"), pattern = "manifest.*.tsv", full.names = T)
manifest_tumor_only <- grep("msi", manifest_tumor_only, invert = T, value = T)
manifest_tumor_only <- lapply(manifest_tumor_only, FUN = function(x) readr::read_tsv(x))
manifest_tumor_only <- data.table::rbindlist(manifest_tumor_only)
colnames(manifest_tumor_only) <- gsub(" ", "_", colnames(manifest_tumor_only))
manifest_tumor_only <- manifest_tumor_only %>%
  filter(Kids_First_Biospecimen_ID %in% c(wgs_tumor_only_ids),
         sample_id %in% hist_df$Sample_ID) %>%
  dplyr::mutate(Sample_ID = sample_id,
                Sequencing_Experiment = "WGS_Tumor_Only",
                experimental_strategy = "WGS_Tumor_Only") %>%
  dplyr::select(Kids_First_Biospecimen_ID, Sample_ID, Sequencing_Experiment) %>%
  unique()

# combine both into manifest
manifest <- rbind(manifest, manifest_tumor_only)

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
  dplyr::rename("Sample" = "Sample_ID",
                "Integrated_Diagnosis" = "diagnosis",
                "Diagnosis_Type" = "diagnosis_type",
                "Tumor_Location" = "Tumor.Location.condensed",
                "Sex" = "Gender",
                "Age" = "age_two_groups",
                # "Age_Three_Groups" = "age_three_groups",
                # "MSI" = "msi_paired",
                "TMB" = "TMB_paired") %>%
  dplyr::select(Sample, Sequencing_Experiment, Integrated_Diagnosis, Diagnosis_Type, Tumor_Location, Sex, Age, TMB)

# save annotation file
write.table(annot_info, file = file.path(output_dir, "annotation_add_tumor_only.txt"), quote = F, sep = "\t", row.names = F)

# 1. get degene info vs GTEx Brain
genes_df <- readRDS(file.path(output_dir, "hope_cohort_vs_gtex_brain_degs.rds"))
deg_genes <- genes_df %>%
  dplyr::mutate(label = ifelse(diff_expr == "up", "OVE", "UNE")) %>%
  dplyr::rename("Gene_name" = "genes") %>%
  filter(Gene_name %in% deg$V1) %>%
  inner_join(manifest, by = c("sample" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::select(Sample_ID, Gene_name, label) %>%
  unique()

# 2. get cnv info
cnv_genes <- readRDS(file.path(data_dir, "Hope-cnv-controlfreec.rds"))
cnv_genes <- cnv_genes %>%
  dplyr::mutate(label = ifelse(status == "Gain", "GAI", "LOS")) %>%
  filter(gene_symbol %in% cnv$V1) %>%
  dplyr::rename("Gene_name" = "gene_symbol",
                "Kids_First_Biospecimen_ID" = "biospecimen_id") %>%
  inner_join(manifest, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::select(Sample_ID, Gene_name, label) %>%
  unique()

# get cnv info (tumor-only)
cnv_genes_tumor_only <- readRDS(file.path(data_dir, "Hope-cnv-controlfreec-tumor-only.rds"))
cnv_genes_tumor_only <- cnv_genes_tumor_only %>%
  dplyr::mutate(label = ifelse(status == "Gain", "GAI", "LOS")) %>%
  filter(gene_symbol %in% cnv$V1,
         biospecimen_id %in% wgs_tumor_only_ids) %>%
  dplyr::rename("Gene_name" = "gene_symbol",
                "Kids_First_Biospecimen_ID" = "biospecimen_id") %>%
  inner_join(manifest, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::select(Sample_ID, Gene_name, label) %>%
  unique()

# combine both
cnv_genes <- rbind(cnv_genes, cnv_genes_tumor_only)

# 3. get snv info
mut_genes <- readRDS(file.path(data_dir, "Hope-consensus-mutation.maf.rds"))

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
  dplyr::rename("Gene_name" = "Hugo_Symbol",
                "Kids_First_Biospecimen_ID" = "Tumor_Sample_Barcode") %>%
  inner_join(manifest, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::select(Sample_ID, Gene_name, label) %>%
  unique()

# get snv info (tumor-only)
mut_genes_tumor_only <- readRDS(file.path(data_dir, "Hope-mutect2-mutation-tumor-only.maf.rds"))
mut_genes_tumor_only <- mut_genes_tumor_only %>%
  filter(Variant_Classification %in% c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "Splice_Site"),
         Tumor_Sample_Barcode %in% wgs_tumor_only_ids) %>%
  dplyr::mutate(label = case_when(Variant_Classification %in% "Missense_Mutation" ~ "MIS",
                                  Variant_Classification %in% "Nonsense_Mutation" ~ "NOS",
                                  Variant_Classification %in% "Frame_Shift_Del" ~ "FSD",
                                  Variant_Classification %in% "Frame_Shift_Ins" ~ "FSI",
                                  Variant_Classification %in% "In_Frame_Del" ~ "IFD",
                                  Variant_Classification %in% "Splice_Site" ~ "SPS")) %>%
  filter(Hugo_Symbol %in% snv$V1) %>%
  dplyr::rename("Gene_name" = "Hugo_Symbol",
                "Kids_First_Biospecimen_ID" = "Tumor_Sample_Barcode") %>%
  inner_join(manifest, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::select(Sample_ID, Gene_name, label) %>%
  unique()

# combine both
mut_genes <- rbind(mut_genes, mut_genes_tumor_only)

# 4. get fusion info
fus_genes <- readRDS(file.path(data_dir, "Hope-fusion-putative-oncogenic.rds")) %>%
  as.data.frame()
fus_genes <- fus_genes %>%
  unite(., col = "Gene_name", Gene1A, Gene1B, Gene2A, Gene2B, na.rm = TRUE, sep = ",") %>%
  mutate(Gene_name = strsplit(as.character(Gene_name), ",")) %>% 
  unnest(Gene_name) %>%
  unique()
fus_genes <- fus_genes %>%
  dplyr::mutate(label = "FUS") %>%
  dplyr::rename("Kids_First_Biospecimen_ID" = "Sample") %>%
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
# "7316-942", "7316-24", "7316-212", "7316UP-904"

# save output matrix
oncogrid_mat <- oncogrid_mat %>%
  filter(Sample %in% annot_info$Sample)
write.table(oncogrid_mat, file = file.path(output_dir, "oncoprint_add_tumor_only.txt"), quote = F, sep = "\t", row.names = F)

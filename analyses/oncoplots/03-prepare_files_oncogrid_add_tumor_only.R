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

# read HGAT driver list from OT
brain_goi_list <- read.delim(file = file.path(input_dir, "hgat_goi_list.tsv"), header = T)

# MMR genes
geneset_db <- msigdbr::msigdbr(category = "C2", subcategory = "KEGG")
geneset_db <- geneset_db[grep("MISMATCH_REPAIR|KEGG_BASE_EXCISION_REPAIR|KEGG_HOMOLOGOUS_RECOMBINATION", geneset_db$gs_name),]
mmr_genes <- unique(geneset_db$human_gene_symbol)

# combined list
driver_genes <- data.frame(V1 = c(brain_goi_list$HGAT, mmr_genes)) %>% unique()

# read reference gene lists based on PNOC003
snv <- read.delim(file.path(input_dir, "snv_genes.tsv")) %>%
  dplyr::select(hg38) %>%
  dplyr::rename("V1" = "hg38") %>%
  rbind(driver_genes) %>% unique()
fusion <- read.delim(file.path(input_dir, "fusion_genes.tsv")) %>%
  dplyr::select(hg38) %>%
  dplyr::rename("V1" = "hg38") %>%
  rbind(driver_genes) %>% unique()
cnv <- read.delim(file.path(input_dir, "cnv_genes.tsv")) %>%
  dplyr::select(hg38) %>%
  dplyr::rename("V1" = "hg38") %>%
  rbind(driver_genes) %>% unique()
deg <- read.delim(file.path(input_dir, "deg_genes.tsv")) %>%
  dplyr::select(hg38) %>%
  dplyr::rename("V1" = "hg38") %>%
  rbind(driver_genes) %>% unique()

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

# read master histology
hist_df <- read_tsv(file.path(data_dir, "Hope-GBM-histologies.tsv"))
hist_df <- hist_df %>%
  filter(!is.na(HOPE_diagnosis))

# fix HOPE_diagnosis_type
hist_df <- hist_df %>%
  mutate(HOPE_diagnosis_type = case_when(
    HOPE_diagnosis_type %in% c("Recurrent", "recurrent", "Recurrent, residual")  ~ "Recurrence",
    HOPE_diagnosis_type == "Primary" ~ "Initial CNS Tumor",
    TRUE ~ as.character(HOPE_diagnosis_type)))

# fix WHO Grade
hist_df <- hist_df %>%
  mutate(HARMONY_WHO.Grade = ifelse(HARMONY_WHO.Grade == "1-2?", NA, HARMONY_WHO.Grade))

# add WHO Grade to Diagnosis 
diagnosis_who_grade_map <- hist_df %>% 
  filter(!is.na(HOPE_diagnosis)) %>%
  mutate(Diagnosis = gsub(" [(].*", "", HOPE_diagnosis)) %>%
  group_by(Diagnosis, HOPE_diagnosis) %>% 
  summarise(HARMONY_WHO.Grade = toString(as.roman(sort(unique(na.omit(HARMONY_WHO.Grade)))))) %>%
  mutate(HARMONY_WHO.Grade = gsub(", ", "/", HARMONY_WHO.Grade)) %>%
  mutate(HARMONY_WHO.Grade = paste0("(WHO grade ", HARMONY_WHO.Grade, ")")) %>%
  mutate(Diagnosis = paste(Diagnosis, HARMONY_WHO.Grade)) %>%
  dplyr::select(HOPE_diagnosis, Diagnosis)

# add new Diagnosis to annotation 
hist_df <- hist_df %>%
  left_join(diagnosis_who_grade_map)

# for tumor-only samples use WGS_tumor_only
hist_df <- hist_df %>%
  mutate(experimental_strategy = ifelse(Kids_First_Biospecimen_ID %in% wgs_tumor_only_ids, "WGS_tumor_only", experimental_strategy))

# pull TMB (for T/N paired samples)
tmb_paired_output <- read_tsv("../tmb-calculation/results/wgs_paired/snv-mutation-tmb-coding.tsv") %>%
  dplyr::rename("TMB" = "tmb") %>%
  dplyr::select(Tumor_Sample_Barcode, TMB) %>%
  filter(Tumor_Sample_Barcode %in% wgs_ids) %>%
  inner_join(hist_df, by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::select(sample_id, TMB) %>%
  unique()

# pull TMB (for tumor-only samples)
tmb_tumor_only_output <- read_tsv("../tmb-calculation/results/wgs_tumor_only/snv-mutation-tmb-coding.tsv") %>%
  dplyr::rename("TMB" = "tmb") %>%
  dplyr::select(Tumor_Sample_Barcode, TMB) %>%
  filter(Tumor_Sample_Barcode %in% wgs_tumor_only_ids) %>%
  inner_join(hist_df, by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::select(sample_id, TMB) %>%
  unique()

# combine both
tmb_output <- rbind(tmb_paired_output, tmb_tumor_only_output)

# add to histology
hist_df <- hist_df %>%
  left_join(tmb_output, by = c("sample_id"))

# select columns of interest
annot_info <- hist_df %>%
  filter(Kids_First_Biospecimen_ID %in% c(rna_ids, wgs_ids, wgs_tumor_only_ids)) %>%
  group_by(sample_id) %>%
  mutate(Sequencing_Experiment = toString(sort(experimental_strategy))) %>%
  dplyr::rename("Sample" = "sample_id",
                "Diagnosis_Type" = "HOPE_diagnosis_type",
                "Tumor_Location" = "HOPE_Tumor.Location.condensed",
                "Sex" = "HARMONY_Gender",
                "Age" = "HARMONY_age_class_derived",
                "Molecular_Subtype" = "molecular_subtype",
                "Cancer_Group" = "cancer_group_short") %>%
  dplyr::select(Sample, Sequencing_Experiment, Diagnosis, Molecular_Subtype, Diagnosis_Type, Tumor_Location, CNS_region, Cancer_Group, Sex, Age, TMB) %>%
  unique()

# save annotation file
write.table(annot_info, file = file.path(output_dir, "annotation_add_tumor_only.txt"), quote = F, sep = "\t", row.names = F)

# 1. get degene info vs GTEx Brain
genes_df <- readRDS(file.path(output_dir, "hope_cohort_vs_gtex_brain_degs.rds"))
deg_genes <- genes_df %>%
  dplyr::mutate(label = ifelse(diff_expr == "up", "OVE", "UNE")) %>%
  dplyr::rename("Gene_name" = "genes",
                "Sample" = "sample") %>%
  filter(Gene_name %in% deg$V1) %>%
  inner_join(hist_df, by = c("Sample" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::select(sample_id, Gene_name, label) %>%
  unique()

# 2. get cnv info
cnv_genes <- readRDS(file.path(data_dir, "Hope-cnv-controlfreec.rds"))
cnv_genes <- cnv_genes %>%
  dplyr::mutate(label = ifelse(status == "Gain", "GAI", "LOS")) %>%
  filter(gene_symbol %in% cnv$V1) %>%
  dplyr::rename("Gene_name" = "gene_symbol") %>%
  inner_join(hist_df, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::select(sample_id, Gene_name, label) %>%
  unique()

# get cnv info (tumor-only)
cnv_genes_tumor_only <- readRDS(file.path(data_dir, "Hope-cnv-controlfreec-tumor-only.rds"))
cnv_genes_tumor_only <- cnv_genes_tumor_only %>%
  dplyr::mutate(label = ifelse(status == "Gain", "GAI", "LOS")) %>%
  filter(gene_symbol %in% cnv$V1,
         Kids_First_Biospecimen_ID %in% wgs_tumor_only_ids) %>%
  dplyr::rename("Gene_name" = "gene_symbol") %>%
  inner_join(hist_df, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::select(sample_id, Gene_name, label) %>%
  unique()

# combine both
cnv_genes <- rbind(cnv_genes, cnv_genes_tumor_only)

# 3. get snv info
mut_genes <- data.table::fread(file.path(data_dir, "Hope-snv-consensus-plus-hotspots.maf.tsv.gz"))

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
  inner_join(hist_df, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::select(sample_id, Gene_name, label) %>%
  unique()

# get snv info (tumor-only)
mut_genes_tumor_only <- data.table::fread(file.path(data_dir, "Hope-tumor-only-snv-mutect2.maf.tsv.gz"))
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
  inner_join(hist_df, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::select(sample_id, Gene_name, label) %>%
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
  inner_join(hist_df, by = c("Kids_First_Biospecimen_ID")) %>%
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
  rownames_to_column('sample_id') %>%
  full_join(cnv_deg %>%
              rownames_to_column('sample_id'), by = "sample_id")

print("Samples missing from oncogrid:")
print(setdiff(hist_df$sample_id, oncogrid_mat$sample_id))
# [1] "7316-212"   "7316-24"    "7316-942"   "7316UP-904"

# save output matrix
oncogrid_mat <- oncogrid_mat %>%
  filter(sample_id %in% annot_info$Sample)
write.table(oncogrid_mat, file = file.path(output_dir, "oncoprint_add_tumor_only.txt"), quote = F, sep = "\t", row.names = F)

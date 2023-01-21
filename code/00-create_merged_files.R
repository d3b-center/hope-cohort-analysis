# script to generate master genomics files
suppressPackageStartupMessages({
  library(reshape2)
  library(GenomicRanges)
  library(tidyverse)
})

# source functions
# source('~/Projects/d3b-patient-report-analysis/code/utils/filter_mutations.R')
# source('~/Projects/d3b-patient-report-analysis/code/utils/filter_cnv.R')
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")

# load reference
chr_map <- read.delim(file.path(root_dir, "../OMPARE", "data", "mart_export_genechr_mapping.txt"), stringsAsFactors = F, check.names = F)
colnames(chr_map) <- c("hgnc_symbol", "gene_start", "gene_end", "chromosome")
chr_map <- chr_map %>%
  filter(hgnc_symbol != "")
# cancer_genes <- readRDS('~/Projects/d3b-patient-report-analysis/data/cancer_gene_list.rds')

# histologies
manifest <- list.files(path = file.path(data_dir, "manifest"), pattern = "manifest.*.tsv", full.names = T)
manifest <- lapply(manifest, FUN = function(x) readr::read_tsv(x))
manifest <- plyr::rbind.fill(manifest)
colnames(manifest) <- gsub(" ", "_", colnames(manifest))
hist <- manifest %>%
  dplyr::select(name, Kids_First_Biospecimen_ID, case_id, sample_id) %>%
  mutate(file_name = gsub('.*/', '', name)) %>%
  unique()

# TPM
fname <- file.path(data_dir, "merged_files", "gene-expression-rsem-tpm-collapsed.rds")
if(!file.exists(fname)){
  # merge
  cmd <- paste('Rscript', '~/Projects/Utils/merge_rsem.R', 
               '--sourcedir', 'data/gene_expression', 
               '--prefix', 'gene-expression-rsem-tpm',
               '--type', 'TPM',
               '--outdir', 'data/merged_files')
  system(cmd)
  
  # collapse
  cmd <- paste('Rscript', '~/Projects/Utils/collapse_rnaseq.R', 
               '--mat', 'data/merged_files/gene-expression-rsem-tpm.rds', 
               '--gene_sym', 'FALSE',
               '--outfile', 'data/merged_files/gene-expression-rsem-tpm-collapsed.rds')
  system(cmd)
  
  # add BS id to TPM (n = 82)
  fname <- file.path("data", "merged_files", "gene-expression-rsem-tpm-collapsed.rds")
  tpm <- readRDS(fname)
  rna_hist <- hist %>%
    mutate(file_name = gsub(".rsem.genes.results.gz", "", file_name)) %>%
    dplyr::filter(file_name %in% colnames(tpm))
  tpm <- tpm %>%
    dplyr::select(rna_hist$file_name)
  identical(rna_hist$file_name, colnames(tpm))
  colnames(tpm) <- rna_hist$Kids_First_Biospecimen_ID
  saveRDS(tpm, file = fname)
}

# expected counts
fname <- file.path("data", "merged_files", "gene-counts-rsem-expected_count-collapsed.rds")
if(!file.exists(fname)){
  # merge
  cmd <- paste('Rscript', '~/Projects/Utils/merge_rsem.R', 
               '--sourcedir', 'data/gene_expression', 
               '--prefix', 'gene-counts-rsem-expected_count',
               '--type', 'expected_count',
               '--outdir', 'data/merged_files')
  system(cmd)
  
  # collapse
  cmd <- paste('Rscript', '~/Projects/Utils/collapse_rnaseq.R', 
               '--mat', 'data/merged_files/gene-counts-rsem-expected_count.rds', 
               '--gene_sym', 'FALSE',
               '--outfile', 'data/merged_files/gene-counts-rsem-expected_count-collapsed.rds')
  system(cmd)
  
  # add BS id to expected counts (n = 82)
  fname <- file.path("data", "merged_files", "gene-counts-rsem-expected_count-collapsed.rds")
  exp_count <- readRDS(fname)
  rna_hist <- hist %>%
    mutate(file_name = gsub(".rsem.genes.results.gz", "", file_name)) %>%
    filter(file_name %in% colnames(exp_count))
  exp_count <- exp_count %>%
    dplyr::select(rna_hist$file_name)
  identical(rna_hist$file_name, colnames(exp_count))
  colnames(exp_count) <- rna_hist$Kids_First_Biospecimen_ID
  saveRDS(exp_count, file = fname)
}

# merge fusions (n = 82)
hope_cohort_fusions <- list.files(path = file.path("data", "gene_fusions"), pattern = "*annoFuse_filter.tsv", recursive = TRUE, full.names = T)
hope_cohort_fusions <- lapply(hope_cohort_fusions, FUN = function(x) readr::read_tsv(x))
hope_cohort_fusions <- plyr::rbind.fill(hope_cohort_fusions)
hope_cohort_fusions <- hope_cohort_fusions %>%
  dplyr::rename("Kids_First_Biospecimen_ID" = "Sample")
length(unique(hope_cohort_fusions$Kids_First_Biospecimen_ID))
saveRDS(hope_cohort_fusions, file = file.path("data", "merged_files", "fusions_merged.rds"))

# merge mutations (n = 71)
hope_cohort_mutations <- list.files(path = file.path("data", "consenus_maf"), recursive = T, full.names = T)
hope_cohort_mutations <- lapply(hope_cohort_mutations, FUN = function(x) readr::read_tsv(x, skip = 1))
hope_cohort_mutations <- plyr::rbind.fill(hope_cohort_mutations)
hope_cohort_mutations <- hope_cohort_mutations %>%
  dplyr::rename("Kids_First_Biospecimen_ID" = "Tumor_Sample_Barcode")
length(unique(hope_cohort_mutations$Kids_First_Biospecimen_ID))
saveRDS(hope_cohort_mutations, file = file.path("data", "merged_files", "snv_merged.rds"))

# function to merge cnv
merge_cnv <- function(nm){
  print(nm)
  sample_name <- gsub(".*/", "", nm)
  # sample_name <- gsub(".controlfreec.CNVs.p.value.txt", "", sample_name)
  sample_name <- hist %>%
    filter(file_name == sample_name) %>%
    pull(Kids_First_Biospecimen_ID) %>%
    unique()
  x <- data.table::fread(nm)

  # map to gene symbols
  query <- with(x, GRanges(chr, IRanges(start = start, end = end)))
  subject <- with(chr_map, GRanges(chromosome, IRanges(start = gene_start, end = gene_end, names = hgnc_symbol)))
  output <- findOverlaps(query = query, subject = subject, type = "within")
  output <- data.frame(x[queryHits(output),], chr_map[subjectHits(output),])

  # modify
  output$status <- stringr::str_to_title(output$status)
  # output <- filter_cnv(myCNVData = output, myCancerGenes = cancer_genes)
  if(nrow(output) > 1){
    output$Kids_First_Biospecimen_ID <- sample_name
    return(output)
  }
}

# merge cnv (n = 71)
hope_cohort_cnv <- list.files(path = 'data/copy_number', pattern = "*.txt", recursive = TRUE, full.names = T)
hope_cohort_cnv <- lapply(hope_cohort_cnv, FUN = function(x) merge_cnv(nm = x))
hope_cohort_cnv <- data.table::rbindlist(hope_cohort_cnv)
print(length(unique(hope_cohort_cnv$Kids_First_Biospecimen_ID)))
saveRDS(hope_cohort_cnv, file = file.path("data/merged_files", "cnv_merged.rds"))

# merge cnv (n = 71)
# hope_cohort_cnv <- list.files(path = file.path("data", "copy_number"), pattern = "*.cns", recursive = TRUE, full.names = T)
# hope_cohort_cnv <- sapply(hope_cohort_cnv, FUN = function(x) read_tsv(x), simplify = FALSE)
# hope_cohort_cnv <- data.table::rbindlist(hope_cohort_cnv, idcol = TRUE)
# hope_cohort_cnv <- hope_cohort_cnv %>%
#   dplyr::rename("file_name" = ".id") %>%
#   mutate(file_name = gsub(".*/", "", file_name)) %>%
#   inner_join(hist %>% dplyr::select(file_name, Kids_First_Biospecimen_ID), by = "file_name")  %>%
#   dplyr::select(-c(file_name)) 
# print(length(unique(hope_cohort_cnv$Kids_First_Biospecimen_ID)))
# saveRDS(hope_cohort_cnv, file = file.path("data", "merged_files", "cnv_merged.rds"))

# merge tumor-only mutations (n = 85)
# hope_cohort_mutations <- list.files(path = file.path("data", "tumor_only_maf"), pattern = "public", recursive = T, full.names = T)
# hope_cohort_mutations <- lapply(hope_cohort_mutations, FUN = function(x) merge_files(x))
# hope_cohort_mutations <- data.table::rbindlist(hope_cohort_mutations)
# hope_cohort_mutations$sample_name <- NULL
# length(unique(hope_cohort_mutations$Tumor_Sample_Barcode))
# saveRDS(hope_cohort_mutations, file = file.path("data", "merged_files", "snv_tumor_only_merged.rds"))

# # merge tumor only cnv (n = 85)
# hope_cohort_cnv <- list.files(path = file.path("data", "tumor_only_copy_number"), pattern = "*.txt", recursive = TRUE, full.names = T)
# hope_cohort_cnv <- lapply(hope_cohort_cnv, FUN = function(x) merge_cnv(nm = x))
# hope_cohort_cnv <- data.table::rbindlist(hope_cohort_cnv)
# print(length(unique(hope_cohort_cnv$sample_name)))
# saveRDS(hope_cohort_cnv, file = file.path("data", "merged_files", "cnv_tumor_only_merged.rds"))

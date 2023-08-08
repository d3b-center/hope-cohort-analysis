# script to generate master genomics files
suppressPackageStartupMessages({
  library(reshape2)
  library(tidyverse)
  library(GenomicRanges)
})

# source functions
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analyses_dir <- file.path(root_dir, "analyses", "merge-files")
input_dir <- file.path(analyses_dir, "input")
output_dir <- file.path(analyses_dir, "results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# scripts
merge_rsem <- file.path(analyses_dir, "utils", "merge_rsem.R")
collapse_rsem <- file.path(analyses_dir, "utils", "collapse_rnaseq.R")

# manifest for gene expression files
rna_manifest <- read_tsv(file.path(input_dir, "manifest", "manifest_20230427_144425_rna.tsv"))
colnames(rna_manifest) <- gsub(" ", "_", colnames(rna_manifest))
rna_manifest <- rna_manifest %>%
  dplyr::select(name, Kids_First_Biospecimen_ID, case_id, sample_id) %>%
  mutate(file_name = gsub('.*/', '', name)) %>%
  unique()

# TPM
fname <- file.path(output_dir, "Hope-gene-expression-rsem-tpm-collapsed.rds")
if(!file.exists(fname)){
  # merge files from cavatica
  cmd <- paste('Rscript', merge_rsem, 
               '--sourcedir', file.path(input_dir, 'gene_expression'), 
               '--output_file', file.path(output_dir, 'Hope-gene-expression-rsem-tpm.rds'),
               '--type', 'TPM')
  system(cmd)
  
  # add BS id to TPM (n = 90)
  tpm <- file.path(output_dir, "Hope-gene-expression-rsem-tpm.rds") %>%
    readRDS() %>%
    column_to_rownames("gene_id")
  rna_hist <- rna_manifest %>%
    mutate(file_name = gsub(".rsem.genes.results.gz", "", file_name)) %>%
    dplyr::filter(file_name %in% colnames(tpm))
  tpm <- tpm %>%
    dplyr::select(rna_hist$file_name)
  identical(rna_hist$file_name, colnames(tpm))
  colnames(tpm) <- rna_hist$Kids_First_Biospecimen_ID
  tpm %>%
    rownames_to_column("gene_id") %>%
    saveRDS(file = file.path(output_dir, "Hope-gene-expression-rsem-tpm.rds"))
  
  # collapse to unique gene symbols
  cmd <- paste('Rscript', collapse_rsem, 
               '--mat', file.path(output_dir, 'Hope-gene-expression-rsem-tpm.rds'), 
               '--gene_sym', 'FALSE',
               '--outfile', file.path(output_dir, 'Hope-gene-expression-rsem-tpm-collapsed.rds'))
  system(cmd)
}

# expected counts
fname <- file.path(output_dir, "Hope-gene-counts-rsem-expected_count-collapsed.rds")
if(!file.exists(fname)){
  # merge files from cavatica
  cmd <- paste('Rscript', merge_rsem, 
               '--sourcedir', file.path(input_dir, 'gene_expression'), 
               '--output_file', file.path(output_dir, 'Hope-gene-counts-rsem-expected_count.rds'),
               '--type', 'TPM')
  system(cmd)
  
  # add BS id to counts (n = 90)
  counts <- file.path(output_dir, "Hope-gene-counts-rsem-expected_count.rds") %>%
    readRDS() %>%
    column_to_rownames("gene_id")
  rna_hist <- rna_manifest %>%
    mutate(file_name = gsub(".rsem.genes.results.gz", "", file_name)) %>%
    dplyr::filter(file_name %in% colnames(counts))
  counts <- counts %>%
    dplyr::select(rna_hist$file_name)
  identical(rna_hist$file_name, colnames(counts))
  colnames(counts) <- rna_hist$Kids_First_Biospecimen_ID
  counts %>%
    rownames_to_column("gene_id") %>%
    saveRDS(file = file.path(output_dir, "Hope-gene-counts-rsem-expected_count.rds"))

  # collapse to unique gene symbols
  cmd <- paste('Rscript', collapse_rsem, 
               '--mat', file.path(output_dir, 'Hope-gene-counts-rsem-expected_count.rds'), 
               '--gene_sym', 'FALSE',
               '--outfile', file.path(output_dir, 'Hope-gene-counts-rsem-expected_count-collapsed.rds'))
  system(cmd)
}

# merge mutations (n = 76)
hope_cohort_mutations <- list.files(path = file.path(input_dir, "consenus_maf"), recursive = T, full.names = T)
hope_cohort_mutations <- lapply(hope_cohort_mutations, FUN = function(x) readr::read_tsv(x, skip = 1))
hope_cohort_mutations <- plyr::rbind.fill(hope_cohort_mutations)
length(unique(hope_cohort_mutations$Tumor_Sample_Barcode))
saveRDS(hope_cohort_mutations, file = file.path(output_dir, "Hope-consensus-mutation.maf.rds"))

# manifest for cnv files
cnv_manifest <- read_tsv(file.path(input_dir, "manifest", "manifest_20230427_150037_cnv.tsv"))
colnames(cnv_manifest) <- gsub(" ", "_", colnames(cnv_manifest))
cnv_manifest <- cnv_manifest %>%
  dplyr::select(name, Kids_First_Biospecimen_ID, case_id, sample_id) %>%
  mutate(file_name = gsub('.*/', '', name)) %>%
  unique()

# function to merge cnv
merge_cnv <- function(nm){
  print(nm)
  sample_name <- gsub(".*/", "", nm)
  sample_name <- cnv_manifest %>%
    filter(file_name == sample_name) %>%
    pull(Kids_First_Biospecimen_ID) %>%
    unique()
  x <- data.table::fread(nm)

  # map to gene symbols
  query <- with(x, GRanges(chr, IRanges(start = start, end = end)))
  subject <- with(gencode_gtf, GRanges(seqnames, IRanges(start = start, end = end, names = gene_symbol)))
  output <- findOverlaps(query = query, subject = subject, type = "within")
  output <- data.frame(x[queryHits(output),], gencode_gtf[subjectHits(output),])

  # modify
  output$status <- stringr::str_to_title(output$status)
  if(nrow(output) > 1){
    output$Kids_First_Biospecimen_ID <- sample_name
    return(output)
  }
}

# get coordinates of genes from gencode v39
gencode_gtf <- rtracklayer::import(con = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.primary_assembly.annotation.gtf.gz")
gencode_gtf <- as.data.frame(gencode_gtf)
gencode_gtf <- gencode_gtf %>%
  dplyr::select(seqnames, start, end, gene_name) %>%
  mutate(seqnames = gsub("^chr", "", seqnames)) %>%
  dplyr::rename("gene_symbol" = "gene_name") %>%
  unique()

# merge cnv (n = 76)
hope_cohort_cnv <- list.files(path = file.path(input_dir, "copy_number"), pattern = "*.txt", recursive = TRUE, full.names = T)
hope_cohort_cnv <- lapply(hope_cohort_cnv, FUN = function(x) merge_cnv(nm = x))
hope_cohort_cnv <- data.table::rbindlist(hope_cohort_cnv)
hope_cohort_cnv <- hope_cohort_cnv %>%
  dplyr::rename("biospecimen_id" = "Kids_First_Biospecimen_ID",
                "copy_number" = "copy.number") %>%
  dplyr::select(biospecimen_id, status,  copy_number, gene_symbol) %>%
  unique()
print(length(unique(hope_cohort_cnv$biospecimen_id)))
saveRDS(hope_cohort_cnv, file = file.path(output_dir, "Hope-cnv-controlfreec.rds"))

# manifest for msi files
msi_manifest <- read_tsv(file.path(input_dir, "manifest", "manifest_20230504_084539_msi.tsv"))
colnames(msi_manifest) <- gsub(" ", "_", colnames(msi_manifest))
msi_manifest <- msi_manifest %>%
  dplyr::select(name, Kids_First_Biospecimen_ID, case_id, sample_id) %>%
  mutate(file_name = gsub('.*/', '', name)) %>%
  unique()

# merge msi (n = 76)
lf <- list.files(path = file.path(input_dir, "msisensor_pro"), pattern = 'msisensor_pro', recursive = T, full.names = T)
hope_cohort_msi <- lapply(lf, read_tsv)
hope_cohort_msi <- do.call(rbind, hope_cohort_msi)
colnames(hope_cohort_msi)[3] <- "Percent"
hope_cohort_msi$name <- gsub('.*/', '', lf)

# add identifiers 
hope_cohort_msi <- hope_cohort_msi %>%
  inner_join(msi_manifest %>% dplyr::select(-c(name)), by = c("name" = "file_name")) %>%
  dplyr::select(Kids_First_Biospecimen_ID, case_id, sample_id, Total_Number_of_Sites, Number_of_Somatic_Sites, Percent)
hope_cohort_msi <- unique(hope_cohort_msi)
print(length(unique(hope_cohort_msi$Kids_First_Biospecimen_ID)))
saveRDS(hope_cohort_msi, file = file.path(output_dir, "Hope-msi-paired.rds"))

# merge cnv (n = 72)
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

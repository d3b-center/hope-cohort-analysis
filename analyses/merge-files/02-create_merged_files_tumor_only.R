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

# histologies
manifest <- list.files(path = file.path(input_dir, "manifest_tumor_only"), pattern = "manifest.*.tsv", full.names = T)
manifest <- lapply(manifest, FUN = function(x) readr::read_tsv(x))
manifest <- plyr::rbind.fill(manifest)
colnames(manifest) <- gsub(" ", "_", colnames(manifest))
hist <- manifest %>%
  dplyr::select(name, Kids_First_Biospecimen_ID, case_id, sample_id) %>%
  mutate(file_name = gsub('.*/', '', name)) %>%
  unique()

# merge tumor-only mutations (n = 93)
hope_cohort_mutations <- list.files(path = file.path(input_dir, "mutect_maf_tumor_only"), pattern = "public", recursive = T, full.names = T)
hope_cohort_mutations <- lapply(hope_cohort_mutations, FUN = function(x) readr::read_tsv(x, skip = 1))
hope_cohort_mutations <- plyr::rbind.fill(hope_cohort_mutations)
length(unique(hope_cohort_mutations$Tumor_Sample_Barcode))
saveRDS(hope_cohort_mutations, file = file.path(output_dir, "Hope-mutect2-mutation-tumor-only.maf.rds"))

# function to merge cnv
merge_cnv <- function(nm){
  print(nm)
  sample_name <- gsub(".*/", "", nm)
  sample_name <- hist %>%
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

# merge tumor only cnv (n = 93)
hope_cohort_cnv <- list.files(path = file.path(input_dir, "copy_number_tumor_only"), pattern = "*.txt", recursive = TRUE, full.names = T)
hope_cohort_cnv <- lapply(hope_cohort_cnv, FUN = function(x) merge_cnv(nm = x))
hope_cohort_cnv <- data.table::rbindlist(hope_cohort_cnv)
hope_cohort_cnv <- hope_cohort_cnv %>%
  dplyr::rename("biospecimen_id" = "Kids_First_Biospecimen_ID",
                "copy_number" = "copy.number") %>%
  dplyr::select(biospecimen_id, status,  copy_number, gene_symbol) %>%
  unique()
print(length(unique(hope_cohort_cnv$biospecimen_id)))
saveRDS(hope_cohort_cnv, file = file.path(output_dir, "Hope-cnv-controlfreec-tumor-only.rds"))

# merge msi (n = 93)
lf <- list.files(path = file.path(input_dir, "msisensor_pro_tumor_only"), pattern = 'msisensor_pro', recursive = T, full.names = T)
hope_cohort_msi <- lapply(lf, read_tsv)
hope_cohort_msi <- do.call(rbind, hope_cohort_msi)
colnames(hope_cohort_msi)[3] <- "Percent"
hope_cohort_msi$name <- gsub('.*/', '', lf)

# add identifiers 
hope_cohort_msi <- hope_cohort_msi %>%
  inner_join(hist %>% dplyr::select(-c(name)), by = c("name" = "file_name")) %>%
  dplyr::select(Kids_First_Biospecimen_ID, case_id, sample_id, Total_Number_of_Sites, Number_of_Somatic_Sites, Percent)
hope_cohort_msi <- unique(hope_cohort_msi)
print(length(unique(hope_cohort_msi$Kids_First_Biospecimen_ID)))
saveRDS(hope_cohort_msi, file = file.path(output_dir, "Hope-msi-tumor_only.rds"))

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

# histology file
hist_df <- readr::read_tsv(file.path(root_dir, "data", "Hope-GBM-histologies-base.tsv"))

# merge tumor-only mutations (n = 90)
hope_cohort_mutations <- list.files(path = file.path(input_dir, "mutect_maf_tumor_only"), pattern = "public", recursive = T, full.names = T)
hope_cohort_mutations <- lapply(hope_cohort_mutations, FUN = function(x) readr::read_tsv(x, skip = 1))
hope_cohort_mutations <- plyr::rbind.fill(hope_cohort_mutations)
hope_cohort_mutations <- hope_cohort_mutations %>%
  filter(t_alt_count != 0, 
         t_depth > 5)
hope_cohort_mutations <- hope_cohort_mutations %>%
  filter(Tumor_Sample_Barcode %in% hist_df$Kids_First_Biospecimen_ID)
length(unique(hope_cohort_mutations$Tumor_Sample_Barcode))
data.table::fwrite(x = hope_cohort_mutations, file = file.path(output_dir, "Hope-tumor-only-snv-mutect2.maf.tsv.gz"), sep = "\t")

# manifest for cnv files
cnv_manifest <- read_tsv(file.path(input_dir, "manifest_tumor_only", "manifest_20230830_152155_cnv.tsv"))
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
  subject <- with(x, GRanges(chr, IRanges(start = start, end = end)))
  query <- with(gencode_gtf, GRanges(seqnames, IRanges(start = start, end = end, names = gene_symbol)))
  output <- findOverlaps(query = query, subject = subject, type = "within")
  output <- data.frame(x[subjectHits(output),], gencode_gtf[queryHits(output),])
  output <- output %>%
    dplyr::select(chr, start, end, gene_symbol, copy.number, status, 
                  genotype, uncertainty,WilcoxonRankSumTestPvalue, KolmogorovSmirnovPvalue) %>%
    unique()
  
  # modify
  # output$status <- stringr::str_to_title(output$status)
  if(nrow(output) > 1){
    output$Kids_First_Biospecimen_ID <- sample_name
    return(output)
  }
}

# get coordinates of genes from gencode v39
gencode_gtf <- rtracklayer::import(con = file.path(root_dir, "data", "gencode.v39.primary_assembly.annotation.gtf.gz"))
gencode_gtf <- as.data.frame(gencode_gtf)
gencode_gtf <- gencode_gtf %>%
  dplyr::select(seqnames, start, end, gene_name) %>%
  mutate(seqnames = gsub("^chr", "", seqnames)) %>%
  dplyr::rename("gene_symbol" = "gene_name") %>%
  unique()

# merge tumor only cnv (n = 90)
hope_cohort_cnv <- list.files(path = file.path(input_dir, "copy_number_tumor_only"), pattern = "*.txt", recursive = TRUE, full.names = T)
hope_cohort_cnv <- lapply(hope_cohort_cnv, FUN = function(x) merge_cnv(nm = x))
hope_cohort_cnv <- data.table::rbindlist(hope_cohort_cnv)
hope_cohort_cnv <- hope_cohort_cnv %>%
  dplyr::rename("copy number" = "copy.number") %>%
  dplyr::select(Kids_First_Biospecimen_ID, chr, start, end, gene_symbol, `copy number`, status, genotype, uncertainty,
                WilcoxonRankSumTestPvalue, KolmogorovSmirnovPvalue) %>%
  unique()
hope_cohort_cnv <- hope_cohort_cnv %>%
  filter(Kids_First_Biospecimen_ID %in% hist_df$Kids_First_Biospecimen_ID)
print(length(unique(hope_cohort_cnv$Kids_First_Biospecimen_ID)))
saveRDS(hope_cohort_cnv, file = file.path(output_dir, "Hope-cnv-controlfreec-tumor-only.rds"))

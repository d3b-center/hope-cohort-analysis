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

# histologies
manifest <- list.files(path = file.path(data_dir, "manifest_tumor_only"), pattern = "manifest.*.tsv", full.names = T)
manifest <- lapply(manifest, FUN = function(x) readr::read_tsv(x))
manifest <- plyr::rbind.fill(manifest)
colnames(manifest) <- gsub(" ", "_", colnames(manifest))
hist <- manifest %>%
  dplyr::select(name, Kids_First_Biospecimen_ID, case_id, sample_id) %>%
  mutate(file_name = gsub('.*/', '', name)) %>%
  unique()

# merge tumor-only mutations (n = 92)
hope_cohort_mutations <- list.files(path = file.path("data", "mutect_maf_tumor_only"), pattern = "public", recursive = T, full.names = T)
hope_cohort_mutations <- lapply(hope_cohort_mutations, FUN = function(x) readr::read_tsv(x, skip = 1))
hope_cohort_mutations <- plyr::rbind.fill(hope_cohort_mutations)
hope_cohort_mutations <- hope_cohort_mutations %>%
  dplyr::rename("Kids_First_Biospecimen_ID" = "Tumor_Sample_Barcode")
length(unique(hope_cohort_mutations$Kids_First_Biospecimen_ID))
saveRDS(hope_cohort_mutations, file = file.path("data", "merged_files", "snv_merged_tumor_only.rds"))

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

# merge tumor only cnv (n = 85)
hope_cohort_cnv <- list.files(path = file.path("data", "copy_number_tumor_only"), pattern = "*.txt", recursive = TRUE, full.names = T)
hope_cohort_cnv <- lapply(hope_cohort_cnv, FUN = function(x) merge_cnv(nm = x))
hope_cohort_cnv <- data.table::rbindlist(hope_cohort_cnv)
print(length(unique(hope_cohort_cnv$Kids_First_Biospecimen_ID)))
saveRDS(hope_cohort_cnv, file = file.path("data", "merged_files", "cnv_merged_tumor_only.rds"))

# merge msi (n = 93)
lf <- list.files(path = "data/msisensor_pro_tumor_only", pattern = 'msisensor_pro', recursive = T, full.names = T)
hope_cohort_msi <- lapply(lf, read_tsv)
hope_cohort_msi <- do.call(rbind, hope_cohort_msi)
colnames(hope_cohort_msi)[3] <- "Percent"
hope_cohort_msi$name <- gsub('.*/', '', lf)

# add sample_id 
hist$case_id[hist$name == "25bb3e6c-40cd-420a-a080-8bf49f1157a0_tumor_only_msisensor_pro"] <- "C1063212"
hist$sample_id[hist$name == "25bb3e6c-40cd-420a-a080-8bf49f1157a0_tumor_only_msisensor_pro"] <- "7316-3636"
hope_cohort_msi <- hope_cohort_msi %>%
  inner_join(hist %>% dplyr::select(-c(name)), by = c("name" = "file_name")) %>%
  dplyr::select(case_id, sample_id, Total_Number_of_Sites, Number_of_Somatic_Sites, Percent)
hope_cohort_msi <- unique(hope_cohort_msi)
saveRDS(hope_cohort_msi, file = file.path(output_dir, "msi_merged.rds"))

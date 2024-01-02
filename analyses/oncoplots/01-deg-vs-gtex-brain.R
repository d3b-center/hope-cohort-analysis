# script to find degs in HOPE cohort vs GTEx brain
suppressPackageStartupMessages({
  library(tidyverse)
  library(reshape2)
})

# data directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analyses_dir <- file.path(root_dir, "analyses", "oncoplots")
output_dir <- file.path(analyses_dir, "results")
dir.create(output_dir, showWarnings = F, recursive = T)
data_dir <- file.path(root_dir, "data")

# source function
source(file.path(analyses_dir, "utils", "noiseq_ss_diffexpr.R"))

output_file <- file.path(output_dir, "hope_cohort_vs_gtex_brain_degs.rds")
# only do this if file does not exists because this takes a lot of time
if(!file.exists(output_file)){
  # gencode reference (v39)
  gencode_gtf <- rtracklayer::import(con = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.primary_assembly.annotation.gtf.gz")
  gencode_gtf <- as.data.frame(gencode_gtf)
  gencode_gtf <- gencode_gtf %>%
    dplyr::select(gene_id, gene_name, gene_type) %>%
    filter(gene_type == "protein_coding") %>%
    unique()
  
  # pull GTEx Brain samples from OT (v12)
  openpedcan_dir <- "~/Projects/OpenPedCan-analysis"
  gtex_sample_info <- read_tsv(file.path(openpedcan_dir, "data", "histologies.tsv"))
  gtex_sample_info <- gtex_sample_info %>%
    filter(cohort == "GTEx", gtex_group == "Brain")
  
  # Dataset 1: GTEx Brain counts
  gtex_brain_counts <- readRDS(file.path(openpedcan_dir, "data", "gene-counts-rsem-expected_count-collapsed.rds"))
  gtex_brain_counts <- gtex_brain_counts %>%
    dplyr::select(gtex_sample_info$Kids_First_Biospecimen_ID)
  gtex_brain_counts <- gtex_brain_counts[grep("^HIST", rownames(gtex_brain_counts), invert = T),]
  gtex_brain_counts <- gtex_brain_counts[rownames(gtex_brain_counts) %in% gencode_gtf$gene_name,]
  
  # Dataset2: HOPE cohort counts
  hope_cohort_counts <- readRDS(file.path(data_dir, 'Hope-gene-counts-rsem-expected_count-collapsed.rds'))
  hope_cohort_counts <- hope_cohort_counts[grep("^HIST", rownames(hope_cohort_counts), invert = T),]
  hope_cohort_counts <- hope_cohort_counts[rownames(hope_cohort_counts) %in% gencode_gtf$gene_name,]
  
  # histology file to get samples of interest
  hist_df <- read_tsv(file = file.path(data_dir, "Hope-GBM-histologies.tsv"))
  hist_df <- hist_df %>% 
    filter(!is.na(HOPE_diagnosis),
           experimental_strategy == "RNA-Seq") %>%
    dplyr::select(Kids_First_Biospecimen_ID, experimental_strategy, RNA_library)
  hope_cohort_counts <- hope_cohort_counts %>%
    dplyr::select(hist_df$Kids_First_Biospecimen_ID)
  
  output_df <- data.frame()
  for(i in 1:ncol(hope_cohort_counts)){
    # extract counts for patient of interest
    patient_of_interest <- colnames(hope_cohort_counts)[i]
    poi_counts <- hope_cohort_counts %>%
      dplyr::select(patient_of_interest)
    poi_sample_info <- hist_df %>%
      filter(Kids_First_Biospecimen_ID %in% patient_of_interest)
    tmp <- run_rnaseq_analysis_noiseq(poi_counts = poi_counts, 
                                      poi_sample_info = poi_sample_info,
                                      ref_counts = gtex_brain_counts, 
                                      ref_sample_info = gtex_sample_info,
                                      sample_name = patient_of_interest)
    # bind to final output
    output_df <- rbind(output_df, tmp)
  }
  
  # save output
  saveRDS(output_df, file = output_file)
}

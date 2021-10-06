# script to find degs in HOPE cohort vs GTEx brain
suppressPackageStartupMessages({
  library(tidyverse)
  library(reshape2)
})

# data directory
data_dir <- 'data/'

# source function
source('utils/ss_diffexpr.R')

# gencode reference
gencode_v27 <- read.delim(file.path(data_dir, 'gencode.v27.primary_assembly.annotation.txt'))
gencode_v27_pc <- gencode_v27 %>%
  filter(biotype == "protein_coding")

# Dataset1: GTex Brain
gtex_brain_tpm <- readRDS(file.path(data_dir, "gtex", "gtex_brain_tpm.rds"))
gtex_brain_tpm <- gtex_brain_tpm[grep("^HIST", rownames(gtex_brain_tpm), invert = T),]
gtex_brain_tpm <- gtex_brain_tpm[rownames(gtex_brain_tpm) %in% gencode_v27_pc$gene_symbol,]

gtex_brain_counts <- readRDS(file.path(data_dir, "gtex", "gtex_brain_counts.rds"))
gtex_brain_counts <- gtex_brain_counts[grep("^HIST", rownames(gtex_brain_counts), invert = T),]
gtex_brain_counts <- gtex_brain_counts[rownames(gtex_brain_counts) %in% gencode_v27_pc$gene_symbol,]

# Dataset2: HOPE cohort
hope_cohort_tpm <- readRDS(file.path(data_dir, 'merged_files', 'gene-expression-rsem-tpm-collapsed.rds'))
hope_cohort_tpm <- hope_cohort_tpm[grep("^HIST", rownames(hope_cohort_tpm), invert = T),]
hope_cohort_tpm <- hope_cohort_tpm[rownames(hope_cohort_tpm) %in% gencode_v27_pc$gene_symbol,]

hope_cohort_counts <- readRDS(file.path(data_dir, 'merged_files', 'gene-counts-rsem-expected_count-collapsed.rds'))
hope_cohort_counts <- hope_cohort_counts[grep("^HIST", rownames(hope_cohort_counts), invert = T),]
hope_cohort_counts <- hope_cohort_counts[rownames(hope_cohort_counts) %in% gencode_v27_pc$gene_symbol,]

# cancer Genes
# cancer_genes <- readRDS(file.path(data_dir, 'cancer_gene_list.rds'))

# input data
hope_cohort_tpm_melt <- hope_cohort_tpm %>% rownames_to_column("gene_symbol") %>%  gather('sample', "tpm", -gene_symbol) 
hope_cohort_counts_melt <- hope_cohort_counts %>% rownames_to_column("gene_symbol") %>% gather('sample', "tpm", -gene_symbol) 

# z-score and return only patient's value (i.e. last column)
get_zscore <- function(x) {
  x <- log2(x+1)
  out <- (x-mean(x))/sd(x)
  return(out[length(out)])
}


calc_degs <- function(expData, gtexData) {
  
  sample_of_interest <- unique(expData$sample)
  print(sample_of_interest)
  
  # Merge GTEx and Patient data on common genes
  intGenesTmp <- intersect(rownames(gtexData), expData$gene_symbol)
  mergeDF <- gtexData %>%
    rownames_to_column("gene_symbol") %>%
    inner_join(expData %>%
                 dplyr::select(-c(sample)), by = "gene_symbol") %>%
    column_to_rownames("gene_symbol")
  colnames(mergeDF)[ncol(mergeDF)] <- sample_of_interest 
  
  # Filter in Patient: TPM > 10
  mergeDF <- mergeDF[mergeDF[,sample_of_interest] > 10,] 
  
  # z-score
  output <- apply(mergeDF, FUN = get_zscore, MARGIN = 1) 
  
  # full data
  # combine tpm and z-score for sample of interest
  genes_df <- data.frame(z_score = output, tpm = mergeDF[names(output),sample_of_interest], sample = sample_of_interest)
  thresh <- 1.5
  genes_df <- genes_df %>%
    mutate(diff_expr = ifelse(z_score < (-1*thresh), "down", 
                              ifelse(z_score > thresh, "up", NA))) %>%
    rownames_to_column("gene_symbol")
  return(genes_df)
}

res <- plyr::ddply(hope_cohort_tpm_melt, 
            .variables = "sample", 
            .fun = function(x) calc_degs(expData = x, gtexData = gtex_brain_tpm))
res <- res %>%
  filter(!is.na(diff_expr))
saveRDS(res, file = 'results/hope_cohort_vs_gtex_brain_degs.rds')

# ssexpr
calc_degs_ssexpr <- function(expData_counts, gtexData_counts) {
  
  sample_of_interest <- unique(expData_counts$sample)
  print(sample_of_interest)
  
  # Merge GTEx and Patient data on common genes
  intGenesTmp <- intersect(rownames(gtexData_counts), expData_counts$gene_symbol)
  mergeDF_counts <- gtexData_counts %>%
    rownames_to_column("gene_symbol") %>%
    inner_join(expData_counts %>%
                 dplyr::select(-c(sample)), by = "gene_symbol") %>%
    column_to_rownames("gene_symbol")
  colnames(mergeDF_counts)[ncol(mergeDF_counts)] <- "sample_of_interest" 
  
  # apply single sample differential expression
  genes_df <- ss_diffexpr(expr = mergeDF_counts, norm_method = "tmm", housekeeping_genes = NULL)
  genes_df$sample <- sample_of_interest
  genes_df <- genes_df %>%
    rownames_to_column("gene_symbol")
  return(genes_df)
}

res <- plyr::ddply(hope_cohort_counts_melt, 
                   .variables = "sample", 
                   .fun = function(x) calc_degs_ssexpr(expData_counts = x, gtexData_counts = gtex_brain_counts))
saveRDS(res, file = 'results/hope_cohort_vs_gtex_brain_degs_edgeR.rds')

# Author: Komal S. Rathi
# NOISeq based differential expression (modified version from d3b-patient-report-analysis)

suppressPackageStartupMessages({
  library(tidyverse)
  library(DGCA)
  library(NOISeq)
  library(sva)
})

run_rnaseq_analysis_noiseq <- function(poi_counts, poi_sample_info, ref_counts, ref_sample_info, sample_name){
  
  # merge ref_counts and sample of interest data on common genes
  merged_counts <- ref_counts %>% 
    rownames_to_column('gene_symbol') %>% 
    inner_join(poi_counts %>% 
                 rownames_to_column('gene_symbol'), by = 'gene_symbol') %>% 
    column_to_rownames('gene_symbol')
  
  # combined sample info
  ref_sample_info <- ref_sample_info %>%
    filter(experimental_strategy == "RNA-Seq",
           Kids_First_Biospecimen_ID %in% colnames(ref_counts)) %>%
    dplyr::select(Kids_First_Biospecimen_ID, RNA_library)
  poi_sample_info <- poi_sample_info %>%
    filter(experimental_strategy == "RNA-Seq") %>%
    dplyr::select(Kids_First_Biospecimen_ID, RNA_library)
  combined_sample_info <- plyr::rbind.fill(ref_sample_info, poi_sample_info)
  combined_sample_info <- combined_sample_info %>%
    filter(!is.na(RNA_library))
  
  # counts columns and sample info should match in order
  merged_counts <- merged_counts %>%
    dplyr::select(combined_sample_info$Kids_First_Biospecimen_ID)
  combined_sample_info <- combined_sample_info %>%
    mutate(group = ifelse(Kids_First_Biospecimen_ID %in% sample_name, "POI", "Others")) %>%
    column_to_rownames('Kids_First_Biospecimen_ID')
  combined_sample_info$group <- factor(combined_sample_info$group, levels = c("POI", "Others"))
  
  # print samples of interest
  print("samples of interest")
  combined_sample_info %>% 
    filter(group == "POI") %>% 
    rownames_to_column("Kids_First_Biospecimen_ID") %>% 
    pull(Kids_First_Biospecimen_ID) %>% 
    print()
  
  # use DCGA to filter out low count, low variance features
  merged_counts <- DGCA::filterGenes(inputMat = merged_counts, 
                                     filterTypes = c("central", "dispersion"),
                                     filterDispersionType = "cv", 
                                     filterDispersionPercentile = 0.2,
                                     sequential = TRUE)
  
  # convert to ExpressionSet
  merged_counts <- readData(merged_counts, factors = combined_sample_info)
  
  # using ComBat batch correction only if more than 1 type of RNA library
  if(length(unique(merged_counts$RNA_library)) > 1){
    print("batch correct")
    merged_counts <- ComBat(dat = log2(merged_counts@assayData$exprs + 1), batch = merged_counts$RNA_library, mean.only = FALSE)
    merged_counts <- 2^merged_counts
    merged_counts <- readData(merged_counts, factors = combined_sample_info)
  }
  
  # run noiseq
  # using corrected counts so no need to normalize again i.e. use norm = "n"
  set.seed(42)
  myresults <- noiseq(input = merged_counts, 
                      factor = "group", 
                      norm = "tmm", 
                      pnr = 0.2, 
                      nss = 20, 
                      v = 0.02, 
                      replicates = "no")
  dge_output <- myresults@results[[1]]
  dge_output <- dge_output %>%
    rownames_to_column('genes') %>%
    dplyr::rename("logFC" = "M") %>%
    filter(prob > 0.9) %>%
    mutate(diff_expr = ifelse(logFC < 0, "down", "up"),
           sample = sample_name)
  dge_output <- dge_output %>%
    dplyr::select(genes, diff_expr, logFC, sample) %>%
    unique()
  return(dge_output)
}
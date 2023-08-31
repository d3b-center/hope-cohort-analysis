# modified version adapted from annoFuse::expression_filter_fusion to retain all input columns 
expression_filter_fusion <- function (standardFusioncalls, expressionMatrix, expressionFilter){
  # standardFusioncalls <- .check_annoFuse_calls(standardFusioncalls)
  stopifnot(is.data.frame(expressionMatrix))
  stopifnot(is.numeric(expressionFilter))
  fusion_sample_gene_df <- standardFusioncalls %>% 
    dplyr::select(Sample, FusionName, Gene1A, Gene1B, Gene2A, Gene2B) %>% 
    tidyr::gather(Gene1A, Gene1B, Gene2A, Gene2B, key = gene_position, value = GeneSymbol) %>% 
    dplyr::select(-gene_position) %>% 
    dplyr::filter(GeneSymbol != "") %>% 
    dplyr::arrange(Sample, FusionName) %>% 
    dplyr::distinct()
  
  expression_long_df <- expressionMatrix %>% 
    dplyr::select(-one_of("gene_id", "EnsembleID")) %>% 
    reshape2::melt(variable.name = "Sample", value.name = "expression_value")
  
  if (!all(fusion_sample_gene_df$Sample %in% expression_long_df$Sample)) {
    warning("Not all samples in expression file. Only returning fusions for samples in expressionMatrix.")
  }
  
  if (length(unique(setdiff(fusion_sample_gene_df$GeneSymbol, expression_long_df$GeneSymbol))) > 0) {
    warning("Not all genes in expression file. Adding genes unique to fusion in output.")
  }
  
  expression_long_df <- expression_long_df %>% 
    dplyr::filter(expression_value > expressionFilter)
  
  expression_filtered_fusions <- fusion_sample_gene_df %>% 
    dplyr::left_join(expression_long_df, by = c("Sample", "GeneSymbol")) %>% 
    dplyr::group_by(FusionName, Sample) %>% 
    dplyr::mutate(all_low_expression = all(is.na(expression_value))) %>% 
    dplyr::filter(!all_low_expression) %>% 
    dplyr::select(FusionName, Sample) %>% 
    dplyr::distinct() %>% 
    dplyr::inner_join(standardFusioncalls, by = c("FusionName", "Sample")) 
  # %>% 
  #   dplyr::select(c("id","LeftBreakpoint", 
  #                   "RightBreakpoint", "FusionName", "Sample", "Caller", 
  #                   "Fusion_Type", "JunctionReadCount", "SpanningFragCount", 
  #                   "Confidence", "annots", "Gene1A", "Gene2A", "Gene1B", 
  #                   "Gene2B", "BreakpointLocation", "SpanningDelta"))
  return(expression_filtered_fusions)
}

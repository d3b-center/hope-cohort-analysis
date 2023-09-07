# collapse rnaseq function
library(optparse)
library(tidyverse)
library(dplyr)

# function
collapse_rnaseq <- function(expr.mat){
  # separate gene_id and gene_symbol
  expr.mat <- expr.mat %>% 
    mutate(gene_id = str_replace(gene_id, "_PAR_Y_", "_"))  %>%
    separate(gene_id, c("gene_id", "gene_symbol"), sep = "\\_", extra = "merge") %>%
    unique()
  
  # uniquify gene_symbol
  expr.collapsed <- expr.mat %>% 
    mutate(means = rowMeans(dplyr::select(.,-gene_id, -gene_symbol))) %>% # take rowMeans
    arrange(desc(means)) %>% # arrange decreasing by means
    distinct(gene_symbol, .keep_all = TRUE) %>% # keep the ones with greatest mean value. If ties occur, keep the first occurencce
    dplyr::select(-c(means, gene_id)) %>%
    unique() %>%
    remove_rownames() %>%
    column_to_rownames('gene_symbol')
  
  return(expr.collapsed)
}

# parameters
option_list <- list(
  make_option(c("--mat"), type = "character",
              help = "Expression Matrix (RSEM TPM, FPKM or Expected counts) (.rds)"),
  make_option(c("--gene_sym"),  type = "logical",
              help = "Is gene symbol present?"),
  make_option(c("--outfile"), type = "character",
              help = "Output filename (.rds)"))

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
mat <- opt$mat
gene_sym <- opt$gene_sym
outfile <- opt$outfile

# read data
mat <- readRDS(mat)

# collapse to unique gene symbols
if(gene_sym == TRUE){
  mat_collapsed <- mat %>%
    column_to_rownames('gene_symbol')
} else {
  mat_collapsed <- collapse_rnaseq(mat)
}

# save output
print(dim(mat_collapsed))
saveRDS(mat_collapsed, file = outfile)





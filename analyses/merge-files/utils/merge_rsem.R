# Script to merge RSEM files

library(reshape2)
library(optparse)

option_list <- list(
  make_option(c("--sourcedir"), type = "character",
              help = "Source directory with all RSEM files"),
  make_option(c("--output_file"), type = "character",
              help = "output file name"),
  make_option(c("--type"), type = "character",
              help = "Type of data"))

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
sourcedir <- opt$sourcedir
output_file <- opt$output_file
type <- opt$type

# universal function to merge expr and fusion results
merge.res <- function(nm){
  print(nm)
  sample_name <- gsub('.*[/]|[.].*','',nm)
  x <- data.table::fread(nm)
  if(nrow(x) > 1){
    x <- as.data.frame(x)
    x$sample_name <- sample_name
    return(x)
  } 
}

# expression (FPKM)
lf <- list.files(path = sourcedir, pattern = "*.genes.results", recursive = T, full.names = T)
expr.mat <- lapply(lf, FUN = function(x) merge.res(x))
expr.mat <- data.table::rbindlist(expr.mat)
expr.mat <- dcast(expr.mat, gene_id~sample_name, value.var = type)
saveRDS(expr.mat, file = output_file)

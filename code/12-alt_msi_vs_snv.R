# script to correlate ALT and MSI to major SNV > 6% 
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(reshape2)
})

# output directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
output_dir <- file.path(root_dir, "results", "alt_correlations")
dir.create(output_dir, showWarnings = F, recursive = T)

# list of major SNV >= 6% mutation
major_snv <- c("TP53", "H3-3A", "NF1", "ATRX", "BRAF", "IDH1", "PDGFRA",
               "TSC2", "ATM", "BCOR", "PIK3CA", "ARID1A", "MSH6", "CIC", "PTPN11", 
               "SETD2", "SMARCA4")

# read oncoprint matrix
mat = read.table(file.path("results", "oncoprint.txt"),  header = TRUE, stringsAsFactors=FALSE, sep = "\t",check.names = FALSE)
mat[is.na(mat)] = ""
rownames(mat) = mat[, 1]
mat = mat[, -1]
mat = mat[, colnames(mat) %in% major_snv]
mat <- melt(as.matrix(mat), value.name = "mutation", varnames = c("sample_id", "gene"))

# only keep SNV annotations
mat <- mat %>%
  dplyr::mutate(type = ifelse(grepl("MIS|NOS|FSD|FSI|IFD|SPS", mutation), 1, 0))
mat <- dcast(mat, sample_id~gene, value.var = "type")
# mat[rowSums(mat) > 0,] %>% dim()

# ALT status, telomere content and MSI percent correlation with the binary matrix
# read combined file for paired MSI output (n = 69)
merged_output <- read_tsv(file.path("results", "msisensor-pro", "msi_output_merged.tsv"))

# combine the file with oncoprint matrix (n = 68)
output_df <- mat %>%
  filter(sample_id != "7316-2810") %>% # this sample does not have any SNV data so remove it
  inner_join(merged_output, by = "sample_id")
output_df$ALT_status <- ifelse(output_df$ALT_status == "ALT+", 1, ifelse(output_df$ALT_status == "ALT-", 0, NA))
vars <- c("t_n_telomere_content", "ALT_status", "MSI_Percent")
final_df <- data.frame()
for(j in 1:length(vars)){
  for(i in 1:length(major_snv)){
    var_val <- vars[j]
    col_val <- major_snv[i]
    pval <- cor.test(output_df[,var_val], output_df[,col_val])$p.value
    estimate <- cor.test(output_df[,var_val], output_df[,col_val])$estimate
    df1 <- data.frame(variable1 = var_val, variable2 = col_val, corr_pvalue = pval, corr_value = estimate)
    final_df <- rbind(final_df, df1)
  }
}
final_df <- final_df %>%
  mutate(Significant = ifelse(corr_pvalue < 0.05, TRUE, FALSE))
write_tsv(final_df, file = file.path(output_dir, "alt_msi_correlation_with_snv.tsv"))

# BRAF and TP53 mutations are mutually exclusive?
tp53_samples = output_df[which(output_df$TP53 == 1),"sample_id"]
braf_samples = output_df[which(output_df$BRAF == 1),"sample_id"]
length(tp53_samples) # 28
length(braf_samples) # 8
intersect(tp53_samples, braf_samples)
# [1] "7316-4215"

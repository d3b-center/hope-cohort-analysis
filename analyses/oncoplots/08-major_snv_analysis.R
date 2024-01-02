# script to correlate ALT and MSI to major SNV > 6% (later when Mateusz responds*)
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(reshape2)
})

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analyses_dir <- file.path(root_dir, "analyses", "oncoplots")
input_dir <- file.path(analyses_dir, "results")
output_dir <- file.path(analyses_dir, "results", "major_snv")
dir.create(output_dir, recursive = T, showWarnings = F)

# list of major SNV >= 6% mutation
mat = read.table(file.path(input_dir, "oncoprint.txt"),  header = TRUE, stringsAsFactors=FALSE, sep = "\t",check.names = FALSE)
mat[is.na(mat)] = ""
rownames(mat) = mat[, 1]
mat = mat[, -1]
mat = t(as.matrix(mat))

# read annotation 
annot_info <- read.delim(file.path(input_dir, "annotation.txt"), header = TRUE, check.names = TRUE)
annot_info <- annot_info %>%
  filter(Sample %in% colnames(mat),
         Sequencing_Experiment != "RNA-Seq") %>%
  remove_rownames() %>%
  column_to_rownames('Sample') %>%
  as.data.frame()
samples_to_use <- intersect(rownames(annot_info), colnames(mat))
mat <- mat[,samples_to_use]
annot_info <- annot_info[samples_to_use,]

# top 20 genes 
major_snv <- apply(mat, 1, FUN = function(x) (length(grep("MIS|NOS|FSD|FSI|NOT|SPS|IFD", x))/ncol(mat))*100)
major_snv <- names(sort(major_snv, decreasing = TRUE)[1:20])

# read oncoprint matrix
mat <- mat[rownames(mat) %in% major_snv,]
mat <- melt(as.matrix(t(mat)), value.name = "mutation", varnames = c("sample_id", "gene"))

# only keep SNV annotations
mat <- mat %>%
  dplyr::mutate(type = ifelse(grepl("MIS|NOS|FSD|FSI|IFD|SPS", mutation), 1, 0))
mat <- dcast(mat, sample_id~gene, value.var = "type")

# ALT status, telomere content and MSI percent correlation with the binary matrix

# read alt status and telomere content
alt_status_output <- read_tsv(file.path("../alt-analysis/results/alt_status_aya_hgg.tsv")) %>%
  dplyr::select(sample_id, t_n_telomere_content, ALT_status)

# read msi output
msi_paired_output <- read_tsv(file.path("../msi-sensor-analysis/results/msisensor-pro-paired/Hope-msi-paired.tsv")) %>%
  dplyr::mutate(Type = ifelse(Percent > 3.5, sample_id, "")) %>%
  dplyr::rename("msi_paired" = "Percent") %>%
  dplyr::select(sample_id, msi_paired)

# combine with oncoprint matrix 
output_df <- mat %>%
  inner_join(alt_status_output) %>%
  inner_join(msi_paired_output)
output_df$ALT_status <- ifelse(output_df$ALT_status == "ALT+", 1, ifelse(output_df$ALT_status == "ALT-", 0, NA))
vars <- c("t_n_telomere_content", "ALT_status", "msi_paired")
final_df <- data.frame()
for(j in 1:length(vars)){
  for(i in 1:length(major_snv)){
    var_val <- vars[j]
    col_val <- major_snv[i]
    if(length(which(colnames(output_df) == col_val)) == 0){
      print("Not present")
    } else {
      pval <- cor.test(output_df[,var_val], output_df[,col_val])$p.value
      estimate <- cor.test(output_df[,var_val], output_df[,col_val])$estimate
      df1 <- data.frame(variable1 = var_val, variable2 = col_val, corr_pvalue = pval, corr_value = estimate)
      final_df <- rbind(final_df, df1)
    }
  }
}
final_df <- final_df %>%
  mutate(Significant = ifelse(corr_pvalue < 0.05, TRUE, FALSE)) %>%
  arrange(desc(Significant))
write_tsv(final_df, file = file.path(output_dir, "alt_msi_correlation_with_snv.tsv"))

# BRAF and TP53 mutations are mutually exclusive?
tp53_samples = output_df[which(output_df$TP53 == 1),"sample_id"]
braf_samples = output_df[which(output_df$BRAF == 1),"sample_id"]
length(tp53_samples) # 30
length(braf_samples) # 8
intersect(tp53_samples, braf_samples)
# [1] "7316-4215"

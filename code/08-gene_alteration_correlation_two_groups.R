# correlation of gene alterations with clinical variables
suppressPackageStartupMessages({
  library(tidyverse)
  library(reshape2)
})

# output directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
input_dir <- file.path(root_dir, "results", "oncoplots_two_groups")
output_dir <- file.path(root_dir, "results", "correlation_analysis_two_groups")
dir.create(output_dir, showWarnings = F, recursive = T)

# matrix
mat <- readr::read_tsv(file = file.path(input_dir, "oncoprint.txt"))

# subset to 95 samples 
hope_cohort_subset <- read.delim(file.path(data_dir, "hope_cohort_subset.tsv"), header = F)
mat <- mat %>%
  filter(Sample %in% hope_cohort_subset$V1)
mat <- melt(mat, id.vars = "Sample", variable.name = "gene", value.name = "alteration_type")
mat$gene <- gsub('[*]', '', mat$gene)
mat <- unique(mat)
mat$alteration_type <- gsub(";OVE", NA, mat$alteration_type)
mat$alteration_type <- gsub("OVE", NA, mat$alteration_type)
mat$alteration_type <- gsub(";UNE", NA, mat$alteration_type)
mat$alteration_type <- gsub("UNE", NA, mat$alteration_type)
mat <- mat %>%
  group_by(Sample, gene) %>%
  summarise(alteration_type = toString(na.omit(alteration_type)))
mat <- mat %>%
  mutate(alteration = ifelse(alteration_type == "", "Absent", "Present"))

# remove samples that do not have DNA info  
annot_info <- read.delim(file.path(input_dir, "annotation.txt"), header = TRUE, check.names = TRUE)
annot_info <- annot_info %>%
  filter(Sample %in% mat$Sample, 
         Sequencing_Experiment != "RNA-Seq") 
mat <- mat %>%
  inner_join(annot_info, by = "Sample") %>%
  as.data.frame()

# do correlation with all genes
cols_to_use <- c("Age", "Sex")
compute_corr <- function(dat, cols){
  print(unique(dat$gene))
  print(nrow(dat))
  # stopifnot(nrow(dat) == 69)
  for(i in 1:length(cols)){
    if(length(unique(dat$alteration)) > 1){
      kw_test <- broom::tidy(kruskal.test(formula = factor(alteration) ~ factor(get(cols[i])), data = dat))
      kruskal_wallis_pval <- kw_test$p.value 
      chisq_test <- chisq.test(x = factor(dat$alteration), y = factor(dat[,cols[i]]))
      chisq_test_pval <- chisq_test$p.value
    } else {
      kruskal_wallis_pval <- NA
      chisq_test_pval <- NA
    }
    
    df <- data.frame(variable = cols[i], 
                     # kruskal_wallis_pval = kruskal_wallis_pval,
                     chisq_test_pval = chisq_test_pval)
    if(i == 1){
      total <- df
    } else {
      total <- rbind(total, df)
    }
  }
  return(total)
}

genes_of_interest <- c("NF1", "TP53", "IDH1", "ATM", "PDGFRA", "CDKN2A", "TSC2", "ASXL1")
output <- plyr::ddply(.data = mat, 
                          .variables = "gene", 
                          .fun = function(x) compute_corr(x, cols = cols_to_use))
output <- output %>%
  dplyr::arrange(chisq_test_pval) %>%
  dplyr::mutate(signif = ifelse(chisq_test_pval < 0.05, TRUE, FALSE))
write.table(output, file = file.path(output_dir, "gene_correlation_with_age_or_sex.tsv"), quote = F, row.names = F, sep = "\t")
output <- output %>%
  filter(gene %in% genes_of_interest)

# 1) coexistence correlation between 2 genes
coexistence_analysis <- function(mat, gene1, gene2){
  dat <- mat %>%
    filter(gene %in% c(gene1, gene2))
  dat <- dcast(dat, Sample + Age + Sex ~ gene, value.var = "alteration")
  dat <- dat %>%
    mutate(coexistence = ifelse(get(gene1) == "Present" & get(gene2) == "Present", "Present", "Absent"))
  chisq_test_coexistence <- chisq.test(x = dat[,gene1], y = dat[,gene2])
  chisq_test_coexistence_age <- chisq.test(x = factor(dat$coexistence), y = factor(dat$Age))
  chisq_test_coexistence_sex <- chisq.test(x = factor(dat$coexistence), y = factor(dat$Sex))
  df <- data.frame(Comparison = paste0(gene1, "-", gene2, " Co-existence"),
                   chisq_pval = round(chisq_test_coexistence$p.value, 5),
                   chisq_pval_Age = round(chisq_test_coexistence_age$p.value, 2),
                   chisq_pval_Sex = round(chisq_test_coexistence_sex$p.value, 2))
  write_tsv(df, file = file.path(output_dir, paste0(gene1, "-", gene2, "-coexistence.tsv")))
  
  # stratified by age
  all_age <- unique(dat$Age)
  age_total <- data.frame()
  for(i in 1:length(all_age)){
    all_age[i]
    gene1_vals <- dat %>%
      filter(Age == all_age[i]) %>% 
      pull(gene1)
    gene2_vals <- dat %>%
      filter(Age == all_age[i]) %>% 
      pull(gene2)
    if(length(unique(gene1_vals)) > 1 & length(unique(gene2_vals)) > 1){
      pval <- chisq.test(x = gene1_vals, y = gene2_vals) 
      pval <- pval$p.value
    } else {
      pval <- NA
    }
    df <- data.frame(comparison = paste(gene1, "vs", gene2), 
                     variable = "Age",
                     value = all_age[i], 
                     chisq_pval = round(pval, 5))
    age_total <- rbind(age_total, df)
  }
  # write_tsv(age_total, file = file.path(output_dir, paste0(gene1, "-vs-", gene2, "-by-age.tsv")))
  
  # stratified by sex
  all_sex <- unique(dat$Sex)
  sex_total <- data.frame()
  for(i in 1:length(all_sex)){
    all_sex[i]
    gene1_vals <- dat %>%
      filter(Sex == all_sex[i]) %>% 
      pull(gene1)
    gene2_vals <- dat %>%
      filter(Sex == all_sex[i]) %>% 
      pull(gene2)
    if(length(unique(gene1_vals)) > 1 & length(unique(gene2_vals)) > 1){
      pval <- chisq.test(x = gene1_vals, y = gene2_vals) 
      pval <- pval$p.value
    } else {
      pval <- NA
    }
    df <- data.frame(comparison = paste(gene1, "vs", gene2), 
                     variable = "Sex",
                     value = all_sex[i], 
                     chisq_pval = round(pval, 5))
    sex_total <- rbind(sex_total, df)
  }
  # write_tsv(sex_total, file = file.path(output_dir, paste0(gene1, "-vs-", gene2, "-by-sex.tsv")))
  
  # stratified by age and sex
  dat$Age_Sex <- paste0(dat$Age, "_", dat$Sex)
  all_age_sex <- unique(dat$Age_Sex)
  age_sex_total <- data.frame()
  for(i in 1:length(all_age_sex)){
    all_age_sex[i]
    gene1_vals <- dat %>%
      filter(Age_Sex == all_age_sex[i]) %>% 
      pull(gene1)
    gene2_vals <- dat %>%
      filter(Age_Sex == all_age_sex[i]) %>% 
      pull(gene2)
    if(length(unique(gene1_vals)) > 1 & length(unique(gene2_vals)) > 1){
      pval <- chisq.test(x = gene1_vals, y = gene2_vals) 
      pval <- pval$p.value
    } else {
      pval <- NA
    }
    df <- data.frame(comparison = paste(gene1, "vs", gene2), 
                     variable = "Age + Sex",
                     value = all_age_sex[i],
                     # Age = gsub("_.*", "", all_age_sex[i]), 
                     # Sex = gsub(".*_", "", all_age_sex[i]), 
                     chisq_pval = round(pval, 5))
    age_sex_total <- rbind(age_sex_total, df)
  }
  # write_tsv(age_sex_total, file = file.path(output_dir, paste0(gene1, "-vs-", gene2, "-by-age-and-sex.tsv")))
  
  # combine all results
  df <- rbind(age_total, sex_total, age_sex_total)
  write_tsv(df, file = file.path(output_dir, paste0(gene1, "-vs-", gene2, "-by-age-and-sex.tsv")))
}

# run co-existence analysis
coexistence_analysis(mat = mat, gene1 = "ATRX", gene2 = "TP53")
coexistence_analysis(mat = mat, gene1 = "H3-3A", gene2 = "TP53")
coexistence_analysis(mat = mat, gene1 = "NF1", gene2 = "ATRX")
coexistence_analysis(mat = mat, gene1 = "LRP1B", gene2 = "NTRK3")


# correlation of gene alterations with clinical variables
suppressPackageStartupMessages({
  library(tidyverse)
  library(reshape2)
})

# matrix
mat <- read.table(file.path("results", "oncoprint.txt"),  header = TRUE, stringsAsFactors=FALSE, sep = "\t",check.names = FALSE)

# subset to 95 samples 
hope_cohort_subset <- read.delim('data/hope_cohort_subset.tsv', header = F)
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
annot_info <- read.delim(file.path("results", "annotation.txt"), header = TRUE, check.names = TRUE)
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
  stopifnot(nrow(dat) == 63)
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
write.table(output, file = 'results/correlation.tsv', quote = F, row.names = F, sep = "\t")
output <- output %>%
  filter(gene %in% genes_of_interest)

# coexistence correlation ATRX and TP53
mat_atrx_tp53 <- mat %>%
  filter(gene %in% c("ATRX", "TP53"))
mat_atrx_tp53 <- dcast(mat_atrx_tp53, Sample + Age + Sex ~ gene, value.var = "alteration")
mat_atrx_tp53 <- mat_atrx_tp53 %>%
  mutate(coexistence = ifelse(ATRX == "Present" & TP53 == "Present", "Present", "Absent"))
chisq_test_atrx_tp53 <- chisq.test(x = mat_atrx_tp53$ATRX, y = mat_atrx_tp53$TP53)
round(chisq_test_atrx_tp53$p.value, 5) # 0.00015

# stratified by age and sex
chisq.test(x = mat_atrx_tp53 %>%
             filter(Age == "14-33.5") %>% pull(ATRX), 
           y = mat_atrx_tp53 %>%
             filter(Age == "14-33.5") %>% pull(TP53)) # 0.007049
chisq.test(x = mat_atrx_tp53 %>%
             filter(Age == "0-14") %>% pull(ATRX), 
           y = mat_atrx_tp53 %>%
             filter(Age == "0-14") %>% pull(TP53)) # 0.0298
chisq.test(x = mat_atrx_tp53 %>%
             filter(Sex == "Male") %>% pull(ATRX), 
           y = mat_atrx_tp53 %>%
             filter(Sex == "Male") %>% pull(TP53)) # 0.004222
chisq.test(x = mat_atrx_tp53 %>%
             filter(Sex == "Female") %>% pull(ATRX), 
           y = mat_atrx_tp53 %>%
             filter(Sex == "Female") %>% pull(TP53)) # 0.06625
chisq.test(x = mat_atrx_tp53 %>%
             filter(Age == "0-14", 
                    Sex == "Female") %>% pull(ATRX), 
           y = mat_atrx_tp53 %>%
             filter(Age == "0-14", 
                    Sex == "Female") %>% pull(TP53)) # 0.4545
chisq.test(x = mat_atrx_tp53 %>%
             filter(Age == "0-14", 
                    Sex == "Male") %>% pull(ATRX), 
           y = mat_atrx_tp53 %>%
             filter(Age == "0-14", 
                    Sex == "Male") %>% pull(TP53)) # 0.1138

# coexistence is correlated with age group or sex 
chisq_test_atrx_tp53_age <- chisq.test(x = factor(mat_atrx_tp53$coexistence), y = factor(mat_atrx_tp53$Age))
round(chisq_test_atrx_tp53_age$p.value, 3) # 0.321
chisq_test_atrx_tp53_sex <- chisq.test(x = factor(mat_atrx_tp53$coexistence), y = factor(mat_atrx_tp53$Sex))
round(chisq_test_atrx_tp53_sex$p.value, 3) # 0.354

# coexistence correlation H3F3A and TP53
mat_h3f3a_tp53 <- mat %>%
  filter(gene %in% c("H3F3A", "TP53"))
mat_h3f3a_tp53 <- dcast(mat_h3f3a_tp53, Sample + Age + Sex ~ gene, value.var = "alteration")
mat_h3f3a_tp53 <- mat_h3f3a_tp53 %>%
  mutate(coexistence = ifelse(H3F3A == "Present" & TP53 == "Present", "Present", "Absent"))
chisq_test_h3f3a_tp53 <- chisq.test(x = mat_h3f3a_tp53$H3F3A, y = mat_h3f3a_tp53$TP53)
round(chisq_test_h3f3a_tp53$p.value, 5) # 0.00939

# stratified by age and sex
chisq.test(x = mat_h3f3a_tp53 %>%
             filter(Age == "14-33.5") %>% pull(H3F3A), 
           y = mat_h3f3a_tp53 %>%
             filter(Age == "14-33.5") %>% pull(TP53)) # 0.3691
chisq.test(x = mat_h3f3a_tp53 %>%
             filter(Age == "0-14") %>% pull(H3F3A), 
           y = mat_h3f3a_tp53 %>%
             filter(Age == "0-14") %>% pull(TP53)) # 0.02095
chisq.test(x = mat_h3f3a_tp53 %>%
             filter(Sex == "Male") %>% pull(H3F3A), 
           y = mat_h3f3a_tp53 %>%
             filter(Sex == "Male") %>% pull(TP53)) # 0.5493
chisq.test(x = mat_h3f3a_tp53 %>%
             filter(Sex == "Female") %>% pull(H3F3A), 
           y = mat_h3f3a_tp53 %>%
             filter(Sex == "Female") %>% pull(TP53)) # 0.004127
chisq.test(x = mat_h3f3a_tp53 %>%
             filter(Age == "0-14", 
                    Sex == "Female") %>% pull(H3F3A), 
           y = mat_h3f3a_tp53 %>%
             filter(Age == "0-14", 
                    Sex == "Female") %>% pull(TP53)) # 0.03131
chisq.test(x = mat_h3f3a_tp53 %>%
             filter(Age == "0-14", 
                    Sex == "Male") %>% pull(H3F3A), 
           y = mat_h3f3a_tp53 %>%
             filter(Age == "0-14", 
                    Sex == "Male") %>% pull(TP53)) # 0.5029


# coexistence is correlated with age group or sex 
chisq_test_h3f3a_tp53_age <- chisq.test(x = factor(mat_h3f3a_tp53$coexistence), y = factor(mat_h3f3a_tp53$Age))
round(chisq_test_h3f3a_tp53_age$p.value, 3) # 1
chisq_test_h3f3a_tp53_sex <- chisq.test(x = factor(mat_h3f3a_tp53$coexistence), y = factor(mat_h3f3a_tp53$Sex))
round(chisq_test_h3f3a_tp53_sex$p.value, 3) # 0.747

# coexistence correlation NF1 and ATRX
mat_nf1_atrx <- mat %>%
  filter(gene %in% c("NF1", "ATRX"))
mat_nf1_atrx <- dcast(mat_nf1_atrx, Sample + Age + Sex ~ gene, value.var = "alteration")
mat_nf1_atrx <- mat_nf1_atrx %>%
  mutate(coexistence = ifelse(NF1 == "Present" & ATRX == "Present", "Present", "Absent"))
chisq_test_nf1_atrx <- chisq.test(x = mat_nf1_atrx$NF1, y = mat_nf1_atrx$ATRX)
round(chisq_test_nf1_atrx$p.value, 5) # 0.04196

# stratified by age and sex
chisq.test(x = mat_nf1_atrx %>%
             filter(Age == "14-33.5") %>% pull(NF1), 
           y = mat_nf1_atrx %>%
             filter(Age == "14-33.5") %>% pull(ATRX)) # 1
chisq.test(x = mat_nf1_atrx %>%
             filter(Age == "0-14") %>% pull(NF1), 
           y = mat_nf1_atrx %>%
             filter(Age == "0-14") %>% pull(ATRX)) # 0.001736
chisq.test(x = mat_nf1_atrx %>%
             filter(Sex == "Male") %>% pull(NF1), 
           y = mat_nf1_atrx %>%
             filter(Sex == "Male") %>% pull(ATRX)) # 0.4446
chisq.test(x = mat_nf1_atrx %>%
             filter(Sex == "Female") %>% pull(NF1), 
           y = mat_nf1_atrx %>%
             filter(Sex == "Female") %>% pull(ATRX)) # 0.1127
chisq.test(x = mat_nf1_atrx %>%
             filter(Age == "0-14", 
                    Sex == "Female") %>% pull(NF1), 
           y = mat_nf1_atrx %>%
             filter(Age == "0-14", 
                    Sex == "Female") %>% pull(ATRX)) # 0.01519
chisq.test(x = mat_nf1_atrx %>%
             filter(Age == "0-14", 
                    Sex == "Male") %>% pull(NF1), 
           y = mat_nf1_atrx %>%
             filter(Age == "0-14", 
                    Sex == "Male") %>% pull(ATRX)) # 0.1138

# coexistence is correlated with age group or sex
chisq_test_nf1_atrx_age <- chisq.test(x = factor(mat_nf1_atrx$coexistence), y = factor(mat_nf1_atrx$Age))
round(chisq_test_nf1_atrx_age$p.value, 3) # 0.335
chisq_test_nf1_atrx_sex <- chisq.test(x = factor(mat_nf1_atrx$coexistence), y = factor(mat_nf1_atrx$Sex))
round(chisq_test_nf1_atrx_sex$p.value, 3) # 0.561



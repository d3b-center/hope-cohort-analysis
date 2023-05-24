# script to do survival analysis with ALT and MSI  
suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(survminer)
  library(data.table)
})

# output directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
output_dir <- file.path(root_dir, "results", "survival_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

# source function
source(file.path(root_dir, "utils", "plotForest.R"))

# master histology
merged_output <- read_tsv(file.path(data_dir, "master_histology_hope_cohort.tsv"))
merged_output <- merged_output %>%
  mutate(OS_days = as.numeric(Age_at_Initial_Diagnosis),
         OS_status = as.numeric(ifelse(last_known_clinical_status == "Alive", 0, 1)))

# add biospecimen identifiers
wgs_ids <- readRDS(file.path(data_dir, "merged_files", "cnv_merged.rds")) %>%
  pull(Kids_First_Biospecimen_ID) %>%
  unique()
rna_ids <- file.path(data_dir, "merged_files", "gene-expression-rsem-tpm-collapsed.rds") %>%
  readRDS() %>%
  colnames()

# add molecular subtype from OT v12
hist_df <- read_tsv(file.path(data_dir, "histologies.tsv"))
hist_df <- hist_df %>%
  filter(sample_id %in% merged_output$Sample_ID,
         Kids_First_Biospecimen_ID %in% c(wgs_ids, rna_ids)) %>%
  dplyr::select(sample_id, molecular_subtype) %>%
  unique()
merged_output <- merged_output %>%
  left_join(hist_df, by = c("Sample_ID" = "sample_id"))

# survival curves stratified by molecular subtype + ALT status
merged_output <- merged_output %>%
  mutate(molecular_subtype = as.factor(molecular_subtype),
         ALT_status = as.factor(ALT_status),
         t_n_telomere_content = as.numeric(t_n_telomere_content),
         MSI_status = as.factor(ifelse(msi_paired > (msi_paired %>% mean(na.rm = T)), "MSI-high", "MSI-low")))
merged_output$molecular_subtype <- relevel(merged_output$molecular_subtype, ref = "HGG, H3 wildtype")
merged_output$ALT_status <- relevel(merged_output$ALT_status, ref = "ALT-")
merged_output$MSI_status <- relevel(merged_output$MSI_status, ref = "MSI-low")
res.cox <- coxph(Surv(OS_days, OS_status) ~ molecular_subtype + ALT_status, data = merged_output)
# survminer::ggforest(res.cox, data = as.data.table(merged_output))
pdf(file = file.path(output_dir, "alt_status_vs_survival_multivariate.pdf"), width = 8, height = 6)
plotForest(model = res.cox)
dev.off()

# kaplan meier
pdf(file = file.path(output_dir, "alt_vs_survival_km.pdf"), height = 8, width = 10, onefile = FALSE)
fit <- survival::survfit(Surv(as.numeric(OS_days), OS_status) ~ ALT_status, data = merged_output)
sdf <- survival::survdiff(Surv(as.numeric(OS_days), OS_status) ~ ALT_status, data = merged_output)
pvalue <- round(1 - pchisq(sdf$chisq, length(sdf$n) - 1), 3)
p <- ggsurvplot(fit, 
                title = "ALT Status vs Survival",
                data = merged_output, 
                font.x = c(20, face = "bold"),
                font.tickslab = c(20),
                font.y = c(20, face = "bold"),
                font.legend = c(20, "bold"),
                pval = FALSE, 
                pval.method = FALSE,
                ggtheme = theme_minimal(),
                linetype = "strata",
                legend = "bottom")  +
  guides(colour = guide_legend(ncol = 2)) 
p <- p$plot +
  annotate("text", x = 12000, y = 0.75, label = paste0("Log-rank\nP-value: ", pvalue), cex=6, col="black", vjust=0, hjust = 1.1, fontface=1)
print(p)
dev.off()

# survival curves stratified by molecular subtype + MSI (High > 3.5 and Low) 
res.cox <- coxph(Surv(OS_days, OS_status) ~ molecular_subtype + MSI_status, data = merged_output)
# survminer::ggforest(res.cox, data = as.data.table(merged_output))
pdf(file = file.path(output_dir, "msi_status_vs_survival_multivariate.pdf"), width = 8, height = 6)
plotForest(model = res.cox)
dev.off()

# kaplan meier
pdf(file = file.path(output_dir, "msi_status_vs_survival_km.pdf"), height = 8, width = 10, onefile = FALSE)
fit <- survival::survfit(Surv(as.numeric(OS_days), OS_status) ~ MSI_status, data = merged_output)
sdf <- survival::survdiff(Surv(as.numeric(OS_days), OS_status) ~ MSI_status, data = merged_output)
pvalue <- round(1 - pchisq(sdf$chisq, length(sdf$n) - 1), 3)
p <- ggsurvplot(fit, 
                title = "MSI Status vs Survival",
                data = merged_output, 
                font.x = c(20, face = "bold"),
                font.tickslab = c(20),
                font.y = c(20, face = "bold"),
                font.legend = c(20, "bold"),
                pval = FALSE, 
                pval.method = FALSE,
                ggtheme = theme_minimal(),
                linetype = "strata",
                legend = "bottom")  +
  guides(colour = guide_legend(ncol = 2)) 
p <- p$plot +
  annotate("text", x = 12000, y = 0.75, label = paste0("Log-rank\nP-value: ", pvalue), cex=6, col="black", vjust=0, hjust = 1.1, fontface=1)
print(p)
dev.off()


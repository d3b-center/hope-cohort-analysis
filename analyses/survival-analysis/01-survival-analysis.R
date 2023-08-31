# script to do survival analysis with ALT and MSI  
suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(survminer)
  library(data.table)
})

# output directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
input_dir <- file.path(root_dir, "analyses", "master-annotation", "results")
analyses_dir <- file.path(root_dir, "analyses", "survival-analysis")
output_dir <- file.path(analyses_dir, "results")
dir.create(output_dir, showWarnings = F, recursive = T)

# source function
source(file.path(analyses_dir, "utils", "plotForest.R"))

# master histology
merged_output <- read_tsv(file.path(input_dir, "master_histology_hope_cohort.tsv"))
merged_output <- merged_output %>%
  mutate(OS_days = as.numeric(Age_at_Initial_Diagnosis),
         OS_status = as.numeric(ifelse(last_known_clinical_status == "Alive", 0, 1)))

# update columns
merged_output <- merged_output %>%
  mutate(molecular_subtype = as.factor(molecular_subtype),
         ALT_status = as.factor(ALT_status),
         t_n_telomere_content = as.numeric(t_n_telomere_content),
         MSI_status = as.factor(ifelse(msi_paired > (msi_paired %>% mean(na.rm = T)), "MSI-high", "MSI-low")))
merged_output$molecular_subtype <- relevel(merged_output$molecular_subtype, ref = "HGG, H3 wildtype")

# keep molecular subtypes with >= 3 samples
merged_output_filtered <- merged_output %>%
  group_by(molecular_subtype) %>%
  mutate(n = n()) %>%
  filter(n >= 3)

# survival curves stratified by molecular subtype + ALT status
merged_output$ALT_status <- relevel(merged_output$ALT_status, ref = "ALT-")
pdf(file = file.path(output_dir, "alt_status_vs_survival_multivariate.pdf"), width = 8, height = 6)
res.cox <- coxph(Surv(OS_days, OS_status) ~ molecular_subtype + ALT_status, data = merged_output)
plotForest(model = res.cox)

# filter subtypes with < 3 samples
res.cox <- coxph(Surv(OS_days, OS_status) ~ molecular_subtype + ALT_status, data = merged_output_filtered)
plotForest(model = res.cox, title_extra = "Molecular subtypes >= 3 samples")
dev.off()

# molecular subtype and telomere content 
pdf(file = file.path(output_dir, "t_n_telomere_content_vs_survival_multivariate.pdf"), width = 8, height = 6)
res.cox <- coxph(Surv(OS_days, OS_status) ~ molecular_subtype + t_n_telomere_content, data = merged_output)
plotForest(model = res.cox)

# filter subtypes with < 3 samples
res.cox <- coxph(Surv(OS_days, OS_status) ~ molecular_subtype + t_n_telomere_content, data = merged_output_filtered)
plotForest(model = res.cox, title_extra = "Molecular subtypes >= 3 samples")
dev.off()

# kaplan meier
pdf(file = file.path(output_dir, "alt_vs_survival_km.pdf"), height = 8, width = 10)
merged_output$ALT_status <- relevel(merged_output$ALT_status, ref = "ALT+")
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

# filter subtypes with < 3 samples
merged_output_filtered$ALT_status <- relevel(merged_output_filtered$ALT_status, ref = "ALT+")
fit <- survival::survfit(Surv(as.numeric(OS_days), OS_status) ~ ALT_status, data = merged_output_filtered)
sdf <- survival::survdiff(Surv(as.numeric(OS_days), OS_status) ~ ALT_status, data = merged_output_filtered)
pvalue <- round(1 - pchisq(sdf$chisq, length(sdf$n) - 1), 3)
p <- ggsurvplot(fit, 
                title = "ALT Status vs Survival\nMolecular subtypes >= 3 samples",
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
pdf(file = file.path(output_dir, "msi_status_vs_survival_multivariate.pdf"), width = 8, height = 6)
merged_output$MSI_status <- relevel(merged_output$MSI_status, ref = "MSI-low")
res.cox <- coxph(Surv(OS_days, OS_status) ~ molecular_subtype + MSI_status, data = merged_output)
plotForest(model = res.cox)

# filter subtypes with < 3 samples
merged_output_filtered$MSI_status <- relevel(merged_output_filtered$MSI_status, ref = "MSI-low")
res.cox <- coxph(Surv(OS_days, OS_status) ~ molecular_subtype + MSI_status, data = merged_output_filtered)
plotForest(model = res.cox, title_extra = "Molecular subtypes >= 3 samples")
dev.off()

# molecular subtype and MSI value 
pdf(file = file.path(output_dir, "msi_paired_vs_survival_multivariate.pdf"), width = 8, height = 6)
res.cox <- coxph(Surv(OS_days, OS_status) ~ molecular_subtype + msi_paired, data = merged_output)
plotForest(model = res.cox)

# filter subtypes with < 3 samples
res.cox <- coxph(Surv(OS_days, OS_status) ~ molecular_subtype + msi_paired, data = merged_output_filtered)
plotForest(model = res.cox, title_extra = "Molecular subtypes >= 3 samples")
dev.off()

# kaplan meier
pdf(file = file.path(output_dir, "msi_status_vs_survival_km.pdf"), height = 8, width = 10)
merged_output$MSI_status <- relevel(merged_output$MSI_status, ref = "MSI-high")
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

# filter subtypes with < 3 samples
merged_output_filtered$MSI_status <- relevel(merged_output_filtered$MSI_status, ref = "MSI-high")
fit <- survival::survfit(Surv(as.numeric(OS_days), OS_status) ~ MSI_status, data = merged_output_filtered)
sdf <- survival::survdiff(Surv(as.numeric(OS_days), OS_status) ~ MSI_status, data = merged_output_filtered)
pvalue <- round(1 - pchisq(sdf$chisq, length(sdf$n) - 1), 3)
p <- ggsurvplot(fit, 
                title = "MSI Status vs Survival\nMolecular subtypes >= 3 samples",
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

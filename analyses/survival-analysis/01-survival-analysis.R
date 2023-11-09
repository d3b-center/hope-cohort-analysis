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
analyses_dir <- file.path(root_dir, "analyses", "survival-analysis")
input_dir <- file.path(analyses_dir, "input")
output_dir <- file.path(analyses_dir, "plots")
dir.create(output_dir, showWarnings = F, recursive = T)

# source function
source(file.path(analyses_dir, "utils", "plotForest.R"))

# read histologies 
hist_df <- read_tsv(file = file.path(data_dir, "Hope-GBM-histologies.tsv"))
hist_df <- hist_df %>% 
  filter(!is.na(molecular_subtype)) %>%
  mutate(OS_days = as.numeric(HOPE_Age_at_Initial_Diagnosis),
         OS_status = as.numeric(ifelse(HOPE_last_known_clinical_status == "Alive", 0, 1)))

# filter to HOPE annotation binary matrix
binary_matrix <- read_tsv(file.path(input_dir, "compare_HOPE_v2plot.annotation.txt"))
binary_matrix <- binary_matrix %>%
  filter(Remove == 0)
hist_df <- hist_df %>%
  filter(sample_id %in% binary_matrix$sample_id)

# 1) add ALT status
alt_status_output <- read_tsv("../alt-analysis/results/alt_status_aya_hgg.tsv") %>%
  dplyr::select(sample_id, t_n_telomere_content, ALT_status)
annot <- hist_df %>%
  inner_join(alt_status_output) %>%
  mutate(molecular_subtype = as.factor(molecular_subtype),
         ALT_status = as.factor(ALT_status),
         t_n_telomere_content = as.numeric(t_n_telomere_content))
annot$molecular_subtype <- relevel(annot$molecular_subtype, ref = "HGG, H3 wildtype")

# keep molecular subtypes with >= 3 samples
annot_filtered <- annot %>%
  group_by(molecular_subtype) %>%
  mutate(n = n()) %>%
  filter(n >= 3)

# survival curves stratified by molecular subtype + ALT status
annot$ALT_status <- relevel(annot$ALT_status, ref = "ALT-")
pdf(file = file.path(output_dir, "alt_status_vs_survival_multivariate.pdf"), width = 8, height = 6)
res.cox <- coxph(Surv(OS_days, OS_status) ~ molecular_subtype + ALT_status, data = annot)
plotForest(model = res.cox)

# # filter subtypes with < 3 samples
# res.cox <- coxph(Surv(OS_days, OS_status) ~ molecular_subtype + ALT_status, data = annot_filtered)
# plotForest(model = res.cox, title_extra = "Molecular subtypes >= 3 samples")
dev.off()

# molecular subtype and telomere content 
pdf(file = file.path(output_dir, "t_n_telomere_content_vs_survival_multivariate.pdf"), width = 8, height = 6)
res.cox <- coxph(Surv(OS_days, OS_status) ~ molecular_subtype + t_n_telomere_content, data = annot)
plotForest(model = res.cox)

# # filter subtypes with < 3 samples
# res.cox <- coxph(Surv(OS_days, OS_status) ~ molecular_subtype + t_n_telomere_content, data = annot_filtered)
# plotForest(model = res.cox, title_extra = "Molecular subtypes >= 3 samples")
dev.off()

# kaplan meier
pdf(file = file.path(output_dir, "alt_vs_survival_km.pdf"), height = 8, width = 10)
annot$ALT_status <- relevel(annot$ALT_status, ref = "ALT+")
fit <- survival::survfit(Surv(as.numeric(OS_days), OS_status) ~ ALT_status, data = annot)
sdf <- survival::survdiff(Surv(as.numeric(OS_days), OS_status) ~ ALT_status, data = annot)
pvalue <- round(1 - pchisq(sdf$chisq, length(sdf$n) - 1), 3)
p <- ggsurvplot(fit, 
                title = "ALT Status vs Survival",
                data = annot, 
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

# # filter subtypes with < 3 samples
# annot_filtered$ALT_status <- relevel(annot_filtered$ALT_status, ref = "ALT+")
# fit <- survival::survfit(Surv(as.numeric(OS_days), OS_status) ~ ALT_status, data = annot_filtered)
# sdf <- survival::survdiff(Surv(as.numeric(OS_days), OS_status) ~ ALT_status, data = annot_filtered)
# pvalue <- round(1 - pchisq(sdf$chisq, length(sdf$n) - 1), 3)
# p <- ggsurvplot(fit, 
#                 title = "ALT Status vs Survival\nMolecular subtypes >= 3 samples",
#                 data = annot, 
#                 font.x = c(20, face = "bold"),
#                 font.tickslab = c(20),
#                 font.y = c(20, face = "bold"),
#                 font.legend = c(20, "bold"),
#                 pval = FALSE, 
#                 pval.method = FALSE,
#                 ggtheme = theme_minimal(),
#                 linetype = "strata",
#                 legend = "bottom")  +
#   guides(colour = guide_legend(ncol = 2)) 
# p <- p$plot +
#   annotate("text", x = 12000, y = 0.75, label = paste0("Log-rank\nP-value: ", pvalue), cex=6, col="black", vjust=0, hjust = 1.1, fontface=1)
# print(p)
dev.off()

# 2) add MSI paired output 
msi_paired_output <- read_tsv(file.path("../msi-sensor-analysis/results/msisensor-pro-paired/Hope-msi-paired.tsv")) %>%
  dplyr::mutate(Type = ifelse(Percent > 3.5, sample_id, "")) %>%
  dplyr::rename("msi_paired" = "Percent")
annot <- hist_df %>%
  inner_join(msi_paired_output) %>%
  mutate(molecular_subtype = as.factor(molecular_subtype),
         MSI_status = as.factor(ifelse(msi_paired > (msi_paired %>% mean(na.rm = T)), "MSI-high", "MSI-low")))
annot$molecular_subtype <- relevel(annot$molecular_subtype, ref = "HGG, H3 wildtype")

# keep molecular subtypes with >= 3 samples
annot_filtered <- annot %>%
  group_by(molecular_subtype) %>%
  mutate(n = n()) %>%
  filter(n >= 3)

# survival curves stratified by molecular subtype + MSI (High > 3.5 and Low) 
pdf(file = file.path(output_dir, "msi_status_vs_survival_multivariate.pdf"), width = 8, height = 6)
annot$MSI_status <- relevel(annot$MSI_status, ref = "MSI-low")
res.cox <- coxph(Surv(OS_days, OS_status) ~ molecular_subtype + MSI_status, data = annot)
plotForest(model = res.cox)

# filter subtypes with < 3 samples
annot_filtered$MSI_status <- relevel(annot_filtered$MSI_status, ref = "MSI-low")
res.cox <- coxph(Surv(OS_days, OS_status) ~ molecular_subtype + MSI_status, data = annot_filtered)
plotForest(model = res.cox, title_extra = "Molecular subtypes >= 3 samples")
dev.off()

# molecular subtype and MSI value 
pdf(file = file.path(output_dir, "msi_paired_vs_survival_multivariate.pdf"), width = 8, height = 6)
res.cox <- coxph(Surv(OS_days, OS_status) ~ molecular_subtype + msi_paired, data = annot)
plotForest(model = res.cox)

# filter subtypes with < 3 samples
res.cox <- coxph(Surv(OS_days, OS_status) ~ molecular_subtype + msi_paired, data = annot_filtered)
plotForest(model = res.cox, title_extra = "Molecular subtypes >= 3 samples")
dev.off()

# kaplan meier
pdf(file = file.path(output_dir, "msi_status_vs_survival_km.pdf"), height = 8, width = 10)
annot$MSI_status <- relevel(annot$MSI_status, ref = "MSI-high")
fit <- survival::survfit(Surv(as.numeric(OS_days), OS_status) ~ MSI_status, data = annot)
sdf <- survival::survdiff(Surv(as.numeric(OS_days), OS_status) ~ MSI_status, data = annot)
pvalue <- round(1 - pchisq(sdf$chisq, length(sdf$n) - 1), 3)
p <- ggsurvplot(fit, 
                title = "MSI Status vs Survival",
                data = annot, 
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
annot_filtered$MSI_status <- relevel(annot_filtered$MSI_status, ref = "MSI-high")
fit <- survival::survfit(Surv(as.numeric(OS_days), OS_status) ~ MSI_status, data = annot_filtered)
sdf <- survival::survdiff(Surv(as.numeric(OS_days), OS_status) ~ MSI_status, data = annot_filtered)
pvalue <- round(1 - pchisq(sdf$chisq, length(sdf$n) - 1), 3)
p <- ggsurvplot(fit, 
                title = "MSI Status vs Survival\nMolecular subtypes >= 3 samples",
                data = annot, 
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

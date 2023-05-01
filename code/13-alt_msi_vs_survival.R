# script to correlate ALT and MSI to major SNV > 6% 
suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(survminer)
})

# output directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
output_dir <- file.path(root_dir, "results", "alt_correlations")
dir.create(output_dir, showWarnings = F, recursive = T)

# read combined file for paired MSI output (n = 69)
merged_output <- read_tsv(file.path("results", "msisensor-pro", "msi_output_merged.tsv"))

# hope clinical file
clin_df <- read_tsv("data/hopeonly_clinical_table_011823.tsv")
clin_df <- clin_df %>%
  dplyr::select(Sample_ID, Age_at_Initial_Diagnosis, last_known_clinical_status) %>%
  mutate(OS_days = Age_at_Initial_Diagnosis,
         OS_status = ifelse(last_known_clinical_status == "Alive", 0, 1))

# add to ALT and MSI output
surv_data <- merged_output %>%
  inner_join(clin_df, by = c("sample_id" = "Sample_ID"))

# survival curves stratified by ALT status 
surv_data$ALT_status <- factor(surv_data$ALT_status, levels = c("ALT+", "ALT-"))
fit <- survival::survfit(Surv(as.numeric(OS_days), OS_status) ~ ALT_status, data = surv_data)
diff <- survdiff(Surv(as.numeric(OS_days), OS_status) ~ ALT_status, data = surv_data)
pvalue <- round(pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE), 3)

pdf(file = file.path(output_dir, "alt_vs_survival.pdf"), height = 8, width = 10, onefile = FALSE)
p <- ggsurvplot(fit, 
                title = "ALT Status vs Survival",
                data = surv_data, 
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

# survival curves stratified by MSI (High > 3.5 and Low) 
surv_data$Type <- ifelse(surv_data$MSI_Percent > (surv_data$MSI_Percent %>% mean), TRUE, FALSE)
fit <- survival::survfit(Surv(as.numeric(OS_days), OS_status) ~ Type, data = surv_data)
pdf(file = file.path(output_dir, "msi_vs_survival.pdf"), height = 8, width = 10, onefile = FALSE)
p <- ggsurvplot(fit,
                title = "MSI vs Survival",
                data = surv_data,
                font.x = c(14, face = "bold"),
                font.tickslab = c(10, face = "bold"),
                font.y = c(14, face = "bold"),
                font.legend = c(12, "bold"),
                pval = TRUE, 
                pval.method = TRUE,
                ggtheme = theme_minimal(),
                linetype = "strata",
                legend = "bottom")  +
  guides(colour = guide_legend(ncol = 2)) 
print(p)
dev.off()

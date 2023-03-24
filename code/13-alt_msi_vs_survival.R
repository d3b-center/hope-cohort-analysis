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
fit <- survival::survfit(Surv(as.numeric(OS_days), OS_status) ~ ALT_status, data = surv_data)
pdf(file = file.path(output_dir, "alt_vs_survival.pdf"), height = 8, width = 10, onefile = FALSE)
p <- ggsurvplot(fit,
                title = "ALT Status",
                data = surv_data,
                font.x = c(12, face = "bold"),
                font.tickslab = c(10, face = "bold"),
                font.y = c(12, face = "bold"),
                font.legend = c(10, "bold"),
                pval = TRUE, 
                pval.method = TRUE,
                ggtheme = theme_minimal(),
                linetype = "strata",
                legend = "bottom")  +
  guides(colour = guide_legend(ncol = 2)) 
print(p)
dev.off()

# survival curves stratified by MSI (High > 3.5 and Low) 
fit <- survival::survfit(Surv(as.numeric(OS_days), OS_status) ~ Type, data = surv_data)
pdf(file = file.path(output_dir, "msi_vs_survival.pdf"), height = 8, width = 10, onefile = FALSE)
p <- ggsurvplot(fit,
                title = "MSI",
                data = surv_data,
                font.x = c(12, face = "bold"),
                font.tickslab = c(10, face = "bold"),
                font.y = c(12, face = "bold"),
                font.legend = c(10, "bold"),
                pval = TRUE, 
                pval.method = TRUE,
                ggtheme = theme_minimal(),
                linetype = "strata",
                legend = "bottom")  +
  guides(colour = guide_legend(ncol = 2)) 
print(p)
dev.off()

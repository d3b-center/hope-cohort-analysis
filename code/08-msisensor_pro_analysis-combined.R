# script to summarise msi-sensor pro output
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
})

# read msisensor pro files
fname <- 'results/msisensor-pro-tumor-only/hope_cohort_msi_sensor_output.tsv'
if(!file.exists(fname)){
  lf <- list.files(path = 'data/msisensor_pro_tumor_only/', pattern = 'msisensor_pro', recursive = T, full.names = T)
  df <- lapply(lf, read_tsv)
  df <- do.call(rbind,df)
  colnames(df)[3] <- "Percent"
  df$name <- gsub('.*/', '', lf)
  manifest <- read_tsv('data/msisensor_pro_tumor_only/manifest_20230131_010544.tsv')
  manifest <- manifest %>% inner_join(df, by = "name")
  manifest <- manifest %>%
    dplyr::select(`Kids First Biospecimen ID`, `Kids First Participant ID`, sample_id, gender, Total_Number_of_Sites, Number_of_Somatic_Sites, Percent) %>%
    unique()
  output_df <- manifest %>%
    mutate(Type = ifelse(Percent >= 3.5, "High", "Low")) %>%
    arrange(Percent)
  write_tsv(output_df, file = fname)
} else {
  output_df <- read_tsv(fname)
}

paired_results <- read_tsv('results/msisensor-pro/hope_cohort_msi_sensor_output.tsv')
paired_results$analysis_type <- "paired"
tumor_only_results <- read_tsv('results/msisensor-pro-tumor-only/hope_cohort_msi_sensor_output.tsv')
tumor_only_results$analysis_type <- "tumor_only"
combined_df <- rbind(paired_results, tumor_only_results)

# barplot
unique_ids <- combined_df[which(duplicated(combined_df$`Kids First Biospecimen ID`)),"Kids First Biospecimen ID"] %>% 
  unique() %>%
  pull(`Kids First Biospecimen ID`)
pdf(file = "results/msisensor-pro-combined/tumor_only_vs_paired_analysis.pdf", width = 20, height = 6)
ggplot(combined_df %>% filter(`Kids First Biospecimen ID` %in% unique_ids), aes(x = `Kids First Biospecimen ID`, y = Percent, fill = analysis_type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  ggpubr::theme_pubr(base_size = 10, legend = "bottom") + 
  ggtitle("Comparison of samples processed using both Paired analysis and Tumor-only analysis") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()

# script to generate telomere ratio
suppressPackageStartupMessages({
  library(tidyverse)
})

# output directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analyses_dir <- file.path(root_dir, "analyses", "alt-analysis")
input_dir <- file.path(analyses_dir, "input")
output_dir <- file.path(analyses_dir, "results")
dir.create(output_dir, showWarnings = F, recursive = T)

telhunt_file <- file.path(input_dir, "76_pairs_tumor_normal_summary.tsv")

telhunt_df <- read_tsv(telhunt_file) %>%
  separate(PID, into = c("cohort_partipant_id", "sample_id"), sep = "_") %>%
  dplyr::rename(sample_type = sample) %>%
  mutate(sample_type = str_to_sentence(sample_type),
         sample_type = case_when(sample_type == "Control" ~ "Normal",
                                 TRUE ~ sample_type))

spec1 <- telhunt_df %>%
  build_wider_spec(names_from = sample_type, 
                   values_from = c(total_reads, intratel_reads, gc_bins_for_correction, total_reads_with_tel_gc, tel_content))

alt_ratio <- telhunt_df %>%
  pivot_wider_spec(spec = spec1, names_repair = "unique") %>%
  mutate(t_n_telomere_content = tel_content_Tumor / tel_content_Normal,
         ALT_status = case_when(t_n_telomere_content >= 1.068 ~ "ALT+",
                                t_n_telomere_content < 1.068 ~ "ALT-")) %>%
  write_tsv(file.path(output_dir, "alt_status_aya_hgg.tsv"))

# how many samples with each?
table(alt_ratio$ALT_status)





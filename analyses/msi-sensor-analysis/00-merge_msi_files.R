# Function: merge sample level MSI output files from cavatica into cohort level TSV files

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analyses_dir <- file.path(root_dir, "analyses", "msi-sensor-analysis")
input_dir <- file.path(analyses_dir, "input")
output_dir <- file.path(analyses_dir, "results")
dir.create(output_dir, showWarnings = F, recursive = T)

# read histologies
hist_df <- read_tsv(file.path(data_dir, "Hope-GBM-histologies.tsv"))

# manifest for msi files
msi_manifest <- read_tsv(file.path(input_dir, "manifest", "manifest_20230504_084539_msi.tsv"))
colnames(msi_manifest) <- gsub(" ", "_", colnames(msi_manifest))
msi_manifest <- msi_manifest %>%
  dplyr::select(name, Kids_First_Biospecimen_ID, case_id, sample_id) %>%
  mutate(file_name = gsub('.*/', '', name)) %>%
  unique()

# filter to BS-ids in histology file only
msi_manifest <- msi_manifest %>%
  filter(Kids_First_Biospecimen_ID %in% hist_df$Kids_First_Biospecimen_ID)

# merge msi (n = 73)
lf <- list.files(path = file.path(input_dir, "msisensor_pro"), pattern = 'msisensor_pro', recursive = T, full.names = T)
hope_cohort_msi <- lapply(lf, read_tsv)
hope_cohort_msi <- do.call(rbind, hope_cohort_msi)
colnames(hope_cohort_msi)[3] <- "Percent"
hope_cohort_msi$name <- gsub('.*/', '', lf)

# add identifiers 
hope_cohort_msi <- hope_cohort_msi %>%
  inner_join(msi_manifest %>% dplyr::select(-c(name)), by = c("name" = "file_name")) %>%
  dplyr::select(Kids_First_Biospecimen_ID, case_id, sample_id, Total_Number_of_Sites, Number_of_Somatic_Sites, Percent)
hope_cohort_msi <- unique(hope_cohort_msi)
print(length(unique(hope_cohort_msi$Kids_First_Biospecimen_ID)))
write_tsv(hope_cohort_msi, file = file.path(output_dir, "msisensor-pro-paired", "Hope-msi-paired.tsv"))

# manifest for msi files
msi_manifest <- read_tsv(file.path(input_dir, "manifest", "manifest_20230504_083829_msi_tumor_only.tsv"))
colnames(msi_manifest) <- gsub(" ", "_", colnames(msi_manifest))
msi_manifest <- msi_manifest %>%
  dplyr::select(name, Kids_First_Biospecimen_ID, case_id, sample_id) %>%
  mutate(file_name = gsub('.*/', '', name)) %>%
  unique()

# filter to BS-ids in histology file only
msi_manifest <- msi_manifest %>%
  filter(Kids_First_Biospecimen_ID %in% hist_df$Kids_First_Biospecimen_ID)

# merge msi (n = 90)
lf <- list.files(path = file.path(input_dir, "msisensor_pro_tumor_only"), pattern = 'msisensor_pro', recursive = T, full.names = T)
hope_cohort_msi <- lapply(lf, read_tsv)
hope_cohort_msi <- do.call(rbind, hope_cohort_msi)
colnames(hope_cohort_msi)[3] <- "Percent"
hope_cohort_msi$name <- gsub('.*/', '', lf)

# add identifiers 
hope_cohort_msi <- hope_cohort_msi %>%
  inner_join(msi_manifest %>% dplyr::select(-c(name)), by = c("name" = "file_name")) %>%
  dplyr::select(Kids_First_Biospecimen_ID, case_id, sample_id, Total_Number_of_Sites, Number_of_Somatic_Sites, Percent)
hope_cohort_msi <- unique(hope_cohort_msi)
print(length(unique(hope_cohort_msi$Kids_First_Biospecimen_ID)))
write_tsv(hope_cohort_msi, file = file.path(output_dir, "msisensor-pro-tumor-only", "Hope-msi-tumor_only.tsv"))

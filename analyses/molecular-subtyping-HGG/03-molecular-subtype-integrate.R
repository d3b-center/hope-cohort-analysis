library(tidyverse)

## set directories and read histology file
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "molecular-subtyping-HGG")
results_dir <- file.path(analysis_dir, "results")
data_dir <- file.path(root_dir, "data")

hist <- readr::read_tsv(file.path(data_dir, "Hope-GBM-histologies-base.tsv")) %>% 
  select(-c(molecular_subtype, gtex_group, gtex_subgroup))

HGG_mol_subtype <- readr::read_tsv(file.path(results_dir, "Hope_subtype.tsv")) %>% 
  select(Kids_First_Biospecimen_ID, molecular_subtype)

## add molecular subtype to histology-base file
hist_with_subtype <- hist %>%
  left_join(HGG_mol_subtype) %>%
  ## Change the molecular subtype of these seven samples with wrong pathology diagnosis and add Note
  ## 7316-1723, 7316-1746, 7316-194, 7316-212, 7316-2857
  mutate(molecular_subtype = case_when(sample_id %in% c("7316-1723", "7316-1746", "7316-194", "7316-212") ~ NA_character_, 
                                       TRUE ~ molecular_subtype), 
         Notes = case_when(sample_id %in% c("7316-1723", "7316-1746", "7316-194", "7316-212", "7316-2151", "7316-2857", "7316-4844") ~ "Final diagnoses updated based on genomic and pathology review", 
                          TRUE ~ NA_character_)) %>% 
  mutate(integrated_diagnosis = case_when(grepl("DMG, H3 K28", molecular_subtype) ~ "Diffuse midline glioma, H3 K28-mutant",
                                          grepl("DHG, H3 G35", molecular_subtype) ~ "Diffuse hemispheric glioma, H3 G35-mutant",
                                          grepl("HGG, H3 wildtype", molecular_subtype) ~ "High-grade glioma, IDH-wildtype and H3-wildtype",
                                          grepl("HGG, IDH", molecular_subtype) ~ "High-grade glioma, IDH-mutant",
                                          molecular_subtype == "HGG, To be classified" ~ "High-grade glioma",
                                          # account for IGH subtypes
                                          grepl("IHG, NTRK-altered", molecular_subtype) ~ "Infant-type hemispheric glioma, NTRK-altered",
                                          grepl("IHG, ALK-altered", molecular_subtype) ~ "Infant-type hemispheric glioma, ALK-altered",
                                          grepl("IHG, ROS1-altered", molecular_subtype) ~ "Infant-type hemispheric glioma, ROS1-altered",
                                          grepl("IHG, MET-altered", molecular_subtype) ~ "Infant-type hemispheric glioma, MET-altered",
                                          molecular_subtype == "IHG, To be classified" ~ "Infant-type hemispheric glioma",
                                          molecular_subtype == "PXA" ~ "Pleomorphic xanthoastrocytoma",
                                          sample_id == "7316-1723" ~ "CIC-rearranged sarcoma",
                                          sample_id == "7316-1746" ~ "Primary intracranial sarcoma, DICER1-mutant",
                                          sample_id == "7316-194" ~ "Rosette-forming glioneuronal tumor",
                                          sample_id == "7316-2151" ~ "Neuroepithelial tumor with PATZ1 fusion",
                                          sample_id == "7316-2857" ~ "CNS tumor with BCOR internal tandem duplication",
                                          TRUE~ NA_character_),
         broad_histology = case_when(molecular_subtype == "PXA" ~ "Pleomorphic xanthoastrocytoma", 
                                     grepl(paste(c("IHG", "HGG", "DHG", "DMG"), collapse = "|"), molecular_subtype) ~ "Diffuse astrocytic and oligodendroglial tumor", 
                                     sample_id %in% c("7316-1723", "7316-1746") ~ "Mesenchymal non‐meningothelial tumor",
                                     sample_id == "7316-194" ~ "Low-grade glial/glioneuronal tumors",
                                     sample_id == "7316-2151" ~ "High-grade neuroepithelial tumor",
                                     sample_id == "7316-2857" ~ "High-grade neuroepithelial tumor",
                                     TRUE ~ NA_character_)) %>% 
  mutate(cancer_group = str_extract(integrated_diagnosis, "[^,]*")) %>%
  select(colnames(.)[!grepl(paste(c("^HARMONY_", "^HOPE_"), collapse = "|"), colnames(.))], 
         starts_with("HARMONY_"), starts_with("HOPE_")) %>%
  select(-short_histology) %>%
  write_tsv(file.path(results_dir, "Hope-GBM-histologies.tsv"))



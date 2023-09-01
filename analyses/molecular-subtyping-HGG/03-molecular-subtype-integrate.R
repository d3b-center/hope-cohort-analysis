library(tidyverse)

## set directories and read histology file
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "molecular-subtyping-HGG")
results_dir <- file.path(analysis_dir, "results")
data_dir <- file.path(root_dir, "data")

hist <- readr::read_tsv(file.path(data_dir, "Hope-GBM-histologies-base.tsv")) %>% 
  select(-molecular_subtype)

HGG_mol_subtype <- readr::read_tsv(file.path(results_dir, "Hope_subtype.tsv")) %>% 
  select(Kids_First_Biospecimen_ID, molecular_subtype)

## add molecular subtype to histology-base file
hist_with_subtype <- hist %>%
  left_join(HGG_mol_subtype) %>%
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
                                          TRUE~ NA_character_),
         broad_histology = case_when(molecular_subtype == "PXA" ~ "Pleomorphic xanthoastrocytoma", 
                                     grepl(paste(c("IHG", "HGG", "DHG", "DMG"), collapse = "|"), molecular_subtype) ~ "Diffuse astrocytic and oligodendroglial tumor", 
                                     TRUE ~ NA_character_), 
         short_histology = case_when(molecular_subtype == "PXA" ~ "LGAT", 
                                     grepl(paste(c("IHG", "HGG", "DHG", "DMG"), collapse = "|"), molecular_subtype) ~ "HGAT", 
                                     TRUE ~ NA_character_)) %>% 
  mutate(cancer_group = str_extract(integrated_diagnosis, "[^,]*")) %>%
  select(colnames(.)[!grepl(paste(c("^HARMONY_", "^HOPE_"), collapse = "|"), colnames(.))], 
         starts_with("HARMONY_"), starts_with("HOPE_")) %>%
  write_tsv(file.path(results_dir, "Hope-GBM-histologies.tsv"))



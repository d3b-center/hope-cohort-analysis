library(tidyverse)

## set directoies
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
input_dir <- file.path(root_dir, "analyses", "proteomic-data", "input")
output_dir <- file.path(root_dir, "analyses", "proteomic-data", "output")
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

## read in files
hist <- readr::read_tsv(file.path(data_dir, "Hope-GBM-histologies-base.tsv"))
whole_cell_protein <- data.table::fread(file.path(input_dir, "harmonized_proteome_imputed_data_09152023.tsv"))
pho <- data.table::fread(file.path(input_dir, "harmonized_phosphosite_imputed_data_09152023.tsv"))
pho_motif <- data.table::fread(file.path(input_dir, "harmonized_phosphosite_imputed_data_ischemia_removed_motif_09152023.tsv"))
prot_ID <- read_tsv(file.path(input_dir, "proteinID_annotation.tsv"))

# Add bs_id to whole cell proteomic
whole_cell_id <- data.frame(sample_id = colnames(whole_cell_protein)[-c(1:3)]) %>% 
  left_join(hist) %>% 
  filter(experimental_strategy == "Whole Cell Proteomics" | grepl("^C3", sample_id)) %>%
  select(Kids_First_Biospecimen_ID, sample_id) %>% 
  arrange(sample_id)

## generate a named vector with sample_id and bs_id
col_mapping <- setNames(whole_cell_id$Kids_First_Biospecimen_ID, whole_cell_id$sample_id)

## assign bs_id to whole cell proteomic matrix
colnames(whole_cell_protein) <- c(colnames(whole_cell_protein)[1:3], 
                                  col_mapping[colnames(whole_cell_protein)[-c(1:3)]])

## remove OldSymbol column
whole_cell_protein <- whole_cell_protein %>% 
  select(-c("OldSymbol", "Symbol.V5")) %>% 
  ## add protein symbols (NP id) %>% 
  left_join(prot_ID[, c("ApprovedGeneSymbol", "protein_HOPE")]) %>% 
  select(protein_HOPE, everything())


## save the file
write_tsv(whole_cell_protein, file.path(output_dir, "Hope_proteome_imputed_data.tsv"))


## add bs_id for phospho data
pho_id <- data.frame(sample_id = colnames(pho)[-c(1:8)]) %>% 
  left_join(hist) %>% 
  filter(experimental_strategy == "Phospho-Proteomics" | grepl("^C3", sample_id)) %>%
  select(Kids_First_Biospecimen_ID, sample_id) %>% 
  arrange(sample_id)

## generate a named vector with sample_id and bs_id
col_mapping <- setNames(pho_id$Kids_First_Biospecimen_ID, pho_id$sample_id)

## assign bs_id to whole cell proteomic matrix
colnames(pho) <- c(colnames(pho)[1:8], 
                   col_mapping[colnames(pho)[-c(1:8)]])

## remove OldSymbol column
pho <- pho %>% 
  select(-OldSymbol)

## save the file
write_tsv(pho, file.path(output_dir, "Hope_phosphosite_imputed_data.tsv"))

## add bs_id for phospho motif data
### in this file, sample_id is in different format, change the format to make it consistant with our sample_id

colnames(pho_motif) <- colnames(pho_motif) %>% 
  gsub("^X", "", .) %>% 
  gsub("\\.", "-", .)

## add bs_id for phospho_motif data
pho_motif_id <- data.frame(sample_id = colnames(pho_motif)[-c(1:10)]) %>% 
  left_join(hist) %>% 
  filter(experimental_strategy == "Phospho-Proteomics" | grepl("^C3", sample_id)) %>%
  select(Kids_First_Biospecimen_ID, sample_id) %>% 
  arrange(sample_id)

## generate a named vector with sample_id and bs_id
col_mapping <- setNames(pho_motif_id$Kids_First_Biospecimen_ID, pho_motif_id$sample_id)

## assign bs_id to whole cell proteomic matrix
colnames(pho_motif) <- c(colnames(pho_motif)[1:10], 
                   col_mapping[colnames(pho_motif)[-c(1:10)]])

## remove OldSymbol column
pho_motif <- pho_motif %>% 
  select(-OldSymbol)

## save the file
write_tsv(pho_motif, file.path(output_dir, "Hope_phosphosite_imputed_data_ischemia_removed_motif.tsv"))


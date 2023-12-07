library(tidyverse)

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
output_dir <- file.path(root_dir, "analyses", "proteomic-data", "output")

## read the gene liftover file 
HGNC_symbol <- readr::read_tsv(file.path(root_dir, "analyses", "proteomic-data", "/input", "hgnc_complete_set.txt"))
prot_df <- data.table::fread(file.path(output_dir, "Hope_proteome_imputed_data_liftover.tsv.gz")) 
Pho_df <- data.table::fread(file.path(output_dir, "Hope_phosphosite_imputed_data_liftover.tsv.gz"))
motif_df <- data.table::fread(file.path(output_dir, "Hope_phosphosite_imputed_data_ischemia_removed_motif_liftover.tsv.gz"))

## identify the old gene symbols that have multiple genes after liftover
dup_gene_symbol <- HGNC_symbol %>% 
  group_by(prev_symbol) %>% 
  summarise(count_gene = n()) %>% 
  filter(count_gene > 1, 
         !is.na(prev_symbol)) %>% 
  select(-count_gene) %>%
  left_join(HGNC_symbol[, c("prev_symbol", "symbol")], by = "prev_symbol")

## set a function to label the gene symbols that has multiple gene symbols after liftover
label_dup <- function(input, file_dir){
  dup_df <- dup_gene_symbol %>% 
    filter(symbol %in% input$ApprovedGeneSymbol) %>% 
    group_by(prev_symbol) %>% 
    summarise(count = n()) %>% 
    filter(count > 1) %>%
    left_join(HGNC_symbol[, c("prev_symbol", "symbol")], by = "prev_symbol") %>% 
    select(-count)
  
  input <- input %>%
    left_join(dup_df, by = c("ApprovedGeneSymbol" = "symbol")) %>% 
    dplyr::rename("Note" = "prev_symbol") %>% 
    mutate(Note = replace_na(Note, "non-duplicate"))
  
  write_tsv(input, file_dir)
}

## label the genes that share the same previous gene symbol and save the file
label_dup(prot_df, file_dir = file.path(output_dir, "Hope_proteome_imputed_data_liftover.tsv.gz"))
label_dup(Pho_df, file_dir = file.path(output_dir, "Hope_phosphosite_imputed_data_liftover.tsv.gz"))
label_dup(motif_df, file_dir = file.path(output_dir, "Hope_phosphosite_imputed_data_ischemia_removed_motif_liftover.tsv.gz"))



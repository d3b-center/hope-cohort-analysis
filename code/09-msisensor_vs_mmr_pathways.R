# 1. KEGG_MISMATCH_REPAIR GSVA
# script to get pathways of interest
suppressPackageStartupMessages({
  library(msigdbr)
  library(tidyverse)
})

# set working directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
results_dir <- file.path(root_dir, "results")
input_dir <- file.path(root_dir, "data")

# tpm_dat
tpm_dat <- file.path(input_dir, "merged_files/gene-expression-rsem-tpm-collapsed.rds")
tpm_dat <- readRDS(tpm_dat)
rna_ids <- unique(colnames(tpm_dat))

# histology file
manifest <- list.files(path = 'data/', pattern = "^manifest.*.csv", full.names = T)
manifest <- lapply(manifest, FUN = function(x) read.csv(x))
manifest <- data.table::rbindlist(manifest)
manifest <- manifest %>%
  filter(Kids.First.Biospecimen.ID %in% c(rna_ids)) %>%
  mutate(experimental_strategy = ifelse(Kids.First.Biospecimen.ID %in% rna_ids, "RNA-Seq", "WGS")) %>%
  dplyr::mutate(Sample = sample_id,
                Sequencing_Experiment = experimental_strategy) %>%
  dplyr::select(Kids.First.Biospecimen.ID, Sample, Sequencing_Experiment) %>%
  unique()

# kegg pathways 
geneset_db <- msigdbr::msigdbr(category = "C2", subcategory = "KEGG")
geneset_db <- geneset_db[grep("MISMATCH_REPAIR|KEGG_BASE_EXCISION_REPAIR|KEGG_HOMOLOGOUS_RECOMBINATION", geneset_db$gs_name),]
write.table(unique(geneset_db$human_gene_symbol), file = 'data/mmr_genes.tsv', col.names = F, row.names = F, quote = F)
geneset_db <- base::split(geneset_db$human_gene_symbol, list(geneset_db$gs_name))

# log2
tpm_dat <- as.matrix(log2(tpm_dat + 1))

# then calculate the Gaussian-distributed scores
ssgsea_scores_each <- GSVA::gsva(expr = tpm_dat,
                                 gset.idx.list = geneset_db,
                                 method = "ssgsea",
                                 min.sz = 1, 
                                 max.sz = 500,
                                 mx.diff = TRUE) %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("pathway_name")

ssgsea_scores_each <- ssgsea_scores_each %>%
  tidyr::gather(Kids_First_Biospecimen_ID, gsea_score, -c('pathway_name')) %>%
  dplyr::select(Kids_First_Biospecimen_ID, pathway_name, gsea_score)

# merge with sample ids
ssgsea_scores_each <- ssgsea_scores_each %>%
  inner_join(manifest, by = c("Kids_First_Biospecimen_ID" = "Kids.First.Biospecimen.ID"))

# merge with msi
fname <- 'results/msisensor-pro/hope_cohort_msi_sensor_output.tsv'
output_df <- read_tsv(fname)
output_df <- ssgsea_scores_each %>%
  inner_join(output_df, by = c("Sample" = "sample_id")) %>%
  dplyr::select(Sample, pathway_name, gsea_score, Percent, Type) %>%
  unique()
ggplot(output_df, aes(x = Percent, y = gsea_score, color = Type)) +
  geom_point(position = "jitter", pch = 21) +
  facet_wrap(~pathway_name)

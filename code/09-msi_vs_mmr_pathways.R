# script to find correlation between MSI and MMR pathways from KEGG
suppressPackageStartupMessages({
  library(msigdbr)
  library(tidyverse)
  library(ggpubr)
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
manifest <- readr::read_tsv(file = 'data/manifest/manifest_20230120_101120_rna.tsv')
manifest <- manifest %>%
  filter(`Kids First Biospecimen ID` %in% c(rna_ids)) %>%
  dplyr::mutate(Sample = sample_id,
                Sequencing_Experiment = experimental_strategy) %>%
  dplyr::select(`Kids First Biospecimen ID`, Sample, Sequencing_Experiment) %>%
  unique()

# kegg pathways 
geneset_db <- msigdbr::msigdbr(category = "C2", subcategory = "KEGG")
geneset_db <- geneset_db[grep("MISMATCH_REPAIR|KEGG_BASE_EXCISION_REPAIR|KEGG_HOMOLOGOUS_RECOMBINATION", geneset_db$gs_name),]
# write.table(unique(geneset_db$human_gene_symbol), file = 'data/mmr_genes.tsv', col.names = F, row.names = F, quote = F)
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
  inner_join(manifest, by = c("Kids_First_Biospecimen_ID" = "Kids First Biospecimen ID"))

# 1) merge with MSI (tumor only)
output_df <- read_tsv(file.path("results", "msisensor-pro-tumor-only", "hope_cohort_msi_sensor_output.tsv"))
output_df <- ssgsea_scores_each %>%
  inner_join(output_df, by = c("Sample" = "sample_id")) %>%
  dplyr::select(Sample, pathway_name, gsea_score, Percent, Type) %>%
  unique()
pdf(file = file.path("results", "msisensor-pro-tumor-only", "msi_vs_mmr_pathways.pdf"), height = 6, width = 10)
ggplot(output_df, aes(x = Percent, y = gsea_score)) +
  xlab("MSI Percent") + 
  ylab("GSVA score") +
  ggpubr::theme_pubr() +
  geom_point(position = "jitter", pch = 21) +
  facet_wrap(~pathway_name, scales = "free") +
  stat_cor(method = "pearson", color = "red")
dev.off()

# 2) merge with MSI (paired)
output_df <- read_tsv(file.path("results", "msisensor-pro", "hope_cohort_msi_sensor_output.tsv"))
output_df <- ssgsea_scores_each %>%
  inner_join(output_df, by = c("Sample" = "sample_id")) %>%
  dplyr::select(Sample, pathway_name, gsea_score, Percent, Type) %>%
  unique()

pdf(file = file.path("results", "msisensor-pro", "msi_vs_mmr_pathways.pdf"), height = 6, width = 10)
ggplot(output_df, aes(x = Percent, y = gsea_score)) +
  xlab("MSI Percent") + 
  ylab("GSVA score") +
  ggpubr::theme_pubr() +
  geom_point(position = "jitter", pch = 21) +
  facet_wrap(~pathway_name, scales = "free") +
  stat_cor(method = "pearson", color = "red")
dev.off()

# script to find correlation between MSI and MMR pathways from KEGG
suppressPackageStartupMessages({
  library(msigdbr)
  library(tidyverse)
  library(ggpubr)
})

# set working directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
input_dir <- file.path(root_dir, "data")
analyses_dir <- file.path(root_dir, "analyses", "msi-sensor-analysis")
output_dir <- file.path(analyses_dir, "results")

# get coordinates of genes from gencode v39
gencode_gtf <- rtracklayer::import(con = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.primary_assembly.annotation.gtf.gz")
gencode_gtf <- as.data.frame(gencode_gtf)
gencode_gtf <- gencode_gtf %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  filter(gene_type == "protein_coding") %>%
  unique()

# tpm_dat
tpm_dat <- file.path(input_dir, "Hope-gene-expression-rsem-tpm-collapsed.rds") %>%
  readRDS()
tpm_dat <- tpm_dat[rownames(tpm_dat) %in% gencode_gtf$gene_name,]

# master histology file
hist_df <- file.path(input_dir, "master_histology_hope_cohort.tsv") %>% read_tsv()

# filter to biospecimens of interest
rna_manifest <- list.files(path = file.path(root_dir, "analyses", "merge-files", "input", "manifest"), pattern = "rna", full.names = T) %>%
  read_tsv()
rna_manifest <- rna_manifest %>%
  filter(sample_id %in% hist_df$Sample_ID) 
tpm_dat <- tpm_dat %>%
  dplyr::select(rna_manifest$`Kids First Biospecimen ID`)

# kegg pathways 
geneset_db <- msigdbr::msigdbr(category = "C2", subcategory = "KEGG")
geneset_db <- geneset_db[grep("MISMATCH_REPAIR|KEGG_BASE_EXCISION_REPAIR|KEGG_HOMOLOGOUS_RECOMBINATION", geneset_db$gs_name),]
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
  inner_join(rna_manifest, by = c("Kids_First_Biospecimen_ID" = "Kids First Biospecimen ID"))
ssgsea_scores_each <- ssgsea_scores_each %>%
  dplyr::select(Kids_First_Biospecimen_ID, pathway_name, gsea_score, sample_id)

# add MSI (paired and tumor-only)
ssgsea_scores_each <- ssgsea_scores_each %>%
  inner_join(hist_df %>%
  dplyr::select(Sample_ID, msi_paired, msi_tumor_only), by = c("sample_id" = "Sample_ID"))

# 1) merge with MSI (tumor only)
pdf(file = file.path(output_dir, "msisensor-pro-tumor-only", "msi_vs_mmr_pathways.pdf"), height = 6, width = 12)
ggplot(ssgsea_scores_each, aes(x = msi_tumor_only, y = gsea_score)) +
  xlab("MSI Percent") + 
  ylab("GSVA score") +
  ggpubr::theme_pubr() +
  geom_point(position = "jitter", pch = 21) +
  facet_wrap(~pathway_name, scales = "free") +
  stat_cor(method = "pearson", color = "red")
dev.off()

# 2) merge with MSI (paired)
pdf(file = file.path(output_dir, "msisensor-pro-paired", "msi_vs_mmr_pathways.pdf"), height = 6, width = 12)
ggplot(ssgsea_scores_each, aes(x = msi_paired, y = gsea_score)) +
  xlab("MSI Percent") + 
  ylab("GSVA score") +
  ggpubr::theme_pubr() +
  geom_point(position = "jitter", pch = 21) +
  facet_wrap(~pathway_name, scales = "free") +
  stat_cor(method = "pearson", color = "red")
dev.off()

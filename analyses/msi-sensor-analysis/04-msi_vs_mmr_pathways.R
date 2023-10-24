# Function: correlation between MSI and MMR pathways from KEGG

# load libraries
suppressPackageStartupMessages({
  library(msigdbr)
  library(tidyverse)
  library(ggpubr)
})

# set working directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data", "v1")
analyses_dir <- file.path(root_dir, "analyses", "msi-sensor-analysis")
output_dir <- file.path(analyses_dir, "results")

# get coordinates of genes from gencode v39
gencode_gtf <- rtracklayer::import(con = file.path(root_dir, "data", "gencode.v39.primary_assembly.annotation.gtf.gz"))
gencode_gtf <- as.data.frame(gencode_gtf)
gencode_gtf <- gencode_gtf %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  filter(gene_type == "protein_coding") %>%
  unique()

# read TPM data
tpm_dat <- readRDS(file = file.path(data_dir, "Hope-gene-expression-rsem-tpm-collapsed.rds"))
tpm_dat <- tpm_dat %>%
  filter(rownames(tpm_dat) %in% gencode_gtf$gene_name)

# read histology file
annot <- read_tsv(file = file.path(data_dir, "Hope-GBM-histologies.tsv"))
annot <- annot %>% 
  filter(!is.na(HOPE_diagnosis),
         experimental_strategy == "RNA-Seq")

# filter to biospecimens of interest
tpm_dat <- tpm_dat %>%
  dplyr::select(annot$Kids_First_Biospecimen_ID)

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
  inner_join(annot)
ssgsea_scores_each <- ssgsea_scores_each %>%
  dplyr::select(Kids_First_Biospecimen_ID, pathway_name, gsea_score, sample_id)

# read MSI paired output 
msi_paired_output <- read_tsv(file.path(output_dir, "msisensor-pro-paired", "Hope-msi-paired.tsv")) %>%
  dplyr::mutate(Type = ifelse(Percent > 3.5, sample_id, "")) %>%
  dplyr::rename("msi_paired" = "Percent") %>%
  dplyr::select(sample_id, msi_paired)

# read MSI tumor-only output 
msi_tumor_only_output <- read_tsv(file.path(output_dir, "msisensor-pro-tumor-only", "Hope-msi-tumor_only.tsv")) %>%
  dplyr::mutate(Type = ifelse(Percent > 3.5, sample_id, "")) %>%
  dplyr::rename("msi_tumor_only" = "Percent") %>%
  dplyr::select(sample_id, msi_tumor_only)

# add MSI (paired and tumor-only)
ssgsea_scores_each <- ssgsea_scores_each %>%
  inner_join(msi_tumor_only_output) %>%
  left_join(msi_paired_output)

# 1) MSI (tumor only) vs MMR pathways
pdf(file = file.path(output_dir, "msisensor-pro-tumor-only", "msi_vs_mmr_pathways.pdf"), height = 6, width = 12)
ggplot(ssgsea_scores_each, aes(x = msi_tumor_only, y = gsea_score)) +
  xlab("% Microsatellite Instability") + 
  ylab("GSVA score") +
  ggpubr::theme_pubr() +
  geom_point(position = "jitter", pch = 21) +
  facet_wrap(~pathway_name, scales = "free") +
  stat_cor(method = "pearson", color = "red")
dev.off()

# 2) MSI (paired) vs MMR pathways
pdf(file = file.path(output_dir, "msisensor-pro-paired", "msi_vs_mmr_pathways.pdf"), height = 6, width = 12)
ggplot(ssgsea_scores_each, aes(x = msi_paired, y = gsea_score)) +
  xlab("% Microsatellite Instability") + 
  ylab("GSVA score") +
  ggpubr::theme_pubr() +
  geom_point(position = "jitter", pch = 21) +
  facet_wrap(~pathway_name, scales = "free") +
  stat_cor(method = "pearson", color = "red")
dev.off()

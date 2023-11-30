# script to liftover PNOC003 genes from hg19 to hg38
suppressPackageStartupMessages({
  library(tidyverse)
})

# read hgnc genes for old and new names (https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt)
hgnc_genes <- read_tsv("../proteomic-data/input/hgnc_complete_set.txt")
hgnc_genes <- hgnc_genes %>%
  mutate(prev_symbol = strsplit(as.character(prev_symbol), "[|]")) %>% 
  unnest(prev_symbol) %>%
  dplyr::select(prev_symbol, symbol) %>%
  unique()

# read gencode v39
gencode_v39 <- rtracklayer::import("../../data/gencode.v39.primary_assembly.annotation.gtf.gz")
gencode_v39 <- gencode_v39 %>%
  as.data.frame()

# 1) snv genes
snv_genes = read_tsv("input/snv_genes.tsv", col_names = "gene_name")
print(dim(snv_genes))
snv_genes <- snv_genes %>%
  mutate(hg19 = ifelse(!gene_name %in% gencode_v39$gene_name, gene_name, ""))
snv_genes <- snv_genes %>%
  left_join(hgnc_genes, by = c("hg19" = "prev_symbol"))
snv_genes <- snv_genes %>%
  mutate(hg38 = ifelse(!is.na(symbol), symbol, gene_name))
snv_genes$symbol <- NULL
print(dim(snv_genes))
write_tsv(snv_genes, file = "input/snv_genes.tsv")

# 2) cnv genes
cnv_genes = read_tsv("input/cnv_genes.tsv", col_names = "gene_name")
cnv_genes <- cnv_genes %>%
  mutate(hg19 = ifelse(!gene_name %in% gencode_v39$gene_name, gene_name, ""))
cnv_genes <- cnv_genes %>%
  left_join(hgnc_genes, by = c("hg19" = "prev_symbol"))
cnv_genes <- cnv_genes %>%
  mutate(hg38 = ifelse(!is.na(symbol), symbol, gene_name))
cnv_genes$symbol <- NULL
write_tsv(cnv_genes, file = "input/cnv_genes.tsv")

# 3) fusion genes
fusion_genes = read_tsv("input/fusion_genes.tsv", col_names = "gene_name")
fusion_genes <- fusion_genes %>%
  mutate(hg19 = ifelse(!gene_name %in% gencode_v39$gene_name, gene_name, ""))
fusion_genes <- fusion_genes %>%
  left_join(hgnc_genes, by = c("hg19" = "prev_symbol"))
fusion_genes <- fusion_genes %>%
  mutate(hg38 = ifelse(!is.na(symbol), symbol, gene_name))
fusion_genes$symbol <- NULL
write_tsv(fusion_genes, file = "input/fusion_genes.tsv")

# 4) deg genes
deg_genes = read_tsv("input/deg_genes.tsv", col_names = "gene_name")
deg_genes <- deg_genes %>%
  mutate(hg19 = ifelse(!gene_name %in% gencode_v39$gene_name, gene_name, ""))
deg_genes <- deg_genes %>%
  left_join(hgnc_genes, by = c("hg19" = "prev_symbol"))
deg_genes <- deg_genes %>%
  mutate(hg38 = ifelse(!is.na(symbol), symbol, gene_name))
deg_genes$symbol <- NULL
write_tsv(deg_genes, file = "input/deg_genes.tsv")

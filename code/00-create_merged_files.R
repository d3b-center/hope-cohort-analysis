# script to generate master genomics files
suppressPackageStartupMessages({
  library(reshape2)
  library(GenomicRanges)
})

# source functions
source('~/Projects/OMPARE/code/utils/filter_mutations.R')
source('~/Projects/OMPARE/code/utils/filter_cnv.R')

# load reference
chr_map <- read.delim(file.path('~/Projects/OMPARE/data/', 'mart_export_genechr_mapping.txt'), stringsAsFactors = F, check.names = F)
colnames(chr_map) <- c("hgnc_symbol", "gene_start", "gene_end", "chromosome")
chr_map <- chr_map %>%
  filter(hgnc_symbol != "")
cancer_genes <- readRDS('~/Projects/OMPARE/data/cancer_gene_list.rds')

# histologies
hist <- read.csv('data/1633446748485-manifest.csv')
hist$name <- gsub('.call.seg|.kallisto.abundance.tsv.gz', '', hist$name)
hist <- hist %>%
  dplyr::select(name, Kids.First.Biospecimen.ID) %>%
  unique()

# TPM
# merge
cmd <- paste('Rscript', '~/Projects/Utils/merge_rsem.R', 
      '--sourcedir', 'data/gene_expression', 
      '--prefix', 'gene-expression-rsem-tpm',
      '--type', 'TPM',
      '--outdir', 'data/merged_files')
system(cmd)

# collapse
cmd <- paste('Rscript', '~/Projects/Utils/collapse_rnaseq.R', 
             '--mat', 'data/merged_files/gene-expression-rsem-tpm.rds', 
             '--gene_sym', 'FALSE',
             '--outfile', 'data/merged_files/gene-expression-rsem-tpm-collapsed.rds')
system(cmd)

# expected counts
# merge
cmd <- paste('Rscript', '~/Projects/Utils/merge_rsem.R', 
             '--sourcedir', 'data/gene_expression', 
             '--prefix', 'gene-counts-rsem-expected_count',
             '--type', 'expected_count',
             '--outdir', 'data/merged_files')
system(cmd)

# collapse
cmd <- paste('Rscript', '~/Projects/Utils/collapse_rnaseq.R', 
             '--mat', 'data/merged_files/gene-counts-rsem-expected_count.rds', 
             '--gene_sym', 'FALSE',
             '--outfile', 'data/merged_files/gene-counts-rsem-expected_count-collapsed.rds')
system(cmd)

# add BS id to TPM
tpm <- readRDS('data/merged_files/gene-expression-rsem-tpm-collapsed.rds')
rna_hist <- hist %>%
  filter(name %in% colnames(tpm))
tpm <- tpm %>%
  dplyr::select(rna_hist$name)
identical(rna_hist$name, colnames(tpm))
colnames(tpm) <- rna_hist$Kids.First.Biospecimen.ID
saveRDS(tpm, file = 'data/merged_files/gene-expression-rsem-tpm-collapsed.rds')

# add BS id to expected counts
exp_count <- readRDS('data/merged_files/gene-counts-rsem-expected_count-collapsed.rds')
rna_hist <- hist %>%
  filter(name %in% colnames(exp_count))
exp_count <- exp_count %>%
  dplyr::select(rna_hist$name)
identical(rna_hist$name, colnames(exp_count))
colnames(exp_count) <- rna_hist$Kids.First.Biospecimen.ID
saveRDS(exp_count, file = 'data/merged_files/gene-counts-rsem-expected_count-collapsed.rds')

# merge cnv
merge_cnv <- function(nm){
  print(nm)
  sample_name <- gsub(".*/", "", nm)
  sample_name <- gsub(".controlfreec.CNVs.p.value.txt", "", sample_name)
  sample_name <- hist %>%
    filter(name == sample_name) %>%
    .$Kids.First.Biospecimen.ID
  x <- data.table::fread(nm)
  
  # map to gene symbols
  query <- with(x, GRanges(chr, IRanges(start = start, end = end)))
  subject <- with(chr_map, GRanges(chromosome, IRanges(start = gene_start, end = gene_end, names = hgnc_symbol)))
  output <- findOverlaps(query = query, subject = subject, type = "within")
  output <- data.frame(x[queryHits(output),], chr_map[subjectHits(output),])
  
  # modify
  output$status <- stringr::str_to_title(output$status)
  output <- filter_cnv(myCNVData = output, myCancerGenes = cancer_genes)
  if(nrow(output) > 1){
    output$sample_name <- sample_name
    return(output)
  }
}

# merge mutations/fusions
merge_files <- function(nm){
  x <- data.table::fread(nm)
  sample_name <- gsub(".*/", "", nm)
  sample_name <- gsub('.arriba.fusions.tsv', '', sample_name)
  if(nrow(x) > 1){
    x$sample_name <- sample_name
    x <- as.data.frame(x)
    return(x)
  } 
}

# merge mutations
hope_cohort_mutations <- list.files(path = 'data/consensus_maf/', pattern = "*.maf", recursive = TRUE, full.names = T)
hope_cohort_mutations <- lapply(hope_cohort_mutations, FUN = function(x) merge_files(x))
hope_cohort_mutations <- data.table::rbindlist(hope_cohort_mutations)

# filter mutations
hope_cohort_mutations <- filter_mutations(myMutData = hope_cohort_mutations, myCancerGenes = cancer_genes)
hope_cohort_mutations$sample_name <- NULL
saveRDS(hope_cohort_mutations, file = file.path("data/merged_files", "consensus_mutation_filtered.rds"))

# merge and filter cnv
hope_cohort_cnv <- list.files(path = 'data/copy_number', pattern = "*.txt", recursive = TRUE, full.names = T)
hope_cohort_cnv <- lapply(hope_cohort_cnv, FUN = function(x) merge_cnv(nm = x))
hope_cohort_cnv <- data.table::rbindlist(hope_cohort_cnv)
saveRDS(hope_cohort_cnv, file = file.path("data/merged_files", "cnv_filtered.rds"))

# fusions
hope_cohort_fusions <- list.files(path = 'data/arriba_fusions/', pattern = "*.arriba.fusions.tsv", recursive = TRUE, full.names = T)
hope_cohort_fusions <- lapply(hope_cohort_fusions, FUN = function(x) merge_files(x))
hope_cohort_fusions <- data.table::rbindlist(hope_cohort_fusions)
colnames(hope_cohort_fusions)[1] <- "gene1"
hope_cohort_fusions <- hope_cohort_fusions %>%
  inner_join(hist, by = c("sample_name" = "name")) %>%
  separate_rows(gene1, gene2, sep = ",", convert = TRUE) %>%
  mutate(gene1 = gsub('[(].*', '', gene1),
         gene2 = gsub('[(].*',' ', gene2),
         reading_frame = ifelse(reading_frame == ".", "other", reading_frame)) %>%
  unite(col = "Gene", gene1, gene2, sep = ", ", na.rm = T) %>%
  separate_rows(Gene, convert = TRUE) %>%
  filter(Gene %in% cancer_genes$Gene_Symbol) %>%
  unique()
saveRDS(hope_cohort_fusions, file = file.path('data/merged_files', "fusions_filtered.rds"))


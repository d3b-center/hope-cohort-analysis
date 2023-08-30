# script to filter by expression and annotate fusion data with PFAM domain information using annoFuse
suppressPackageStartupMessages({
  library(tidyverse)
  library(annoFuse)
})

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analyses_dir <- file.path(root_dir, "analyses", "merge-files")
input_dir <- file.path(analyses_dir, "input")
output_dir <- file.path(analyses_dir, "results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# merge fusions (n = 90)
hope_cohort_fusions <- list.files(path = file.path(input_dir, "gene_fusions"), pattern = "*annoFuse_filter.tsv", recursive = TRUE, full.names = T)
hope_cohort_fusions <- lapply(hope_cohort_fusions, FUN = function(x) readr::read_tsv(x))
hope_cohort_fusions <- plyr::rbind.fill(hope_cohort_fusions)
length(unique(hope_cohort_fusions$Sample))

# get merged expression matrix
expr_mat <- file.path(output_dir, "Hope-gene-expression-rsem-tpm.rds") %>%
  readRDS()

# add expression filter TPM < 1
# modified version adapted from annoFuse::expression_filter_fusion to retain all input columns
# the original function selects and returns only a subset of columns
source(file.path(analyses_dir, "utils", "expression_filter_fusion.R"))
hope_cohort_fusions <- expression_filter_fusion(standardFusioncalls = hope_cohort_fusions, 
                                                                   expressionMatrix = expr_mat %>% 
                                                                     mutate(gene_id = str_replace(gene_id, "_PAR_Y_", "_"))  %>%
                                                                     separate(gene_id, c("gene_id", "GeneSymbol"), sep = "\\_", extra = "merge") %>%
                                                                     unique(), 
                                                                   expressionFilter = 1)
hope_cohort_fusions <- unique(hope_cohort_fusions)

# get pfam data from annoFuse
bioMartDataPfam <- readRDS(system.file("extdata", "pfamDataBioMart.RDS", package = "annoFuseData"))

# add unique id to join domain info obtained from gene1a and gene1b
hope_cohort_fusions <- hope_cohort_fusions %>%
  rownames_to_column("id") 

# get domain info on standard fusion calls
df <- annoFuse::get_Pfam_domain(standardFusioncalls = hope_cohort_fusions, 
                                bioMartDataPfam = bioMartDataPfam)

# add gene1A domain info
gene1a <- df$Gene1A %>% 
  dplyr::select(id, LeftBreakpoint, RightBreakpoint, FusionName, Sample, Caller, Gene1A, Gene2A, Gene1B, Gene2B, Gene1A_DOMAIN_RETAINED_IN_FUSION) %>%
  unique()
hope_cohort_fusions <- hope_cohort_fusions %>%
  left_join(gene1a)

# add gene1B domain info
gene1b <- df$Gene1B %>% 
  dplyr::select(id, LeftBreakpoint, RightBreakpoint, FusionName, Sample, Caller, Gene1A, Gene2A, Gene1B, Gene2B, Gene1B_DOMAIN_RETAINED_IN_FUSION) %>%
  unique()
hope_cohort_fusions <- hope_cohort_fusions %>%
  left_join(gene1b)

# add caller count
caller_count <- annoFuse::called_by_n_callers(standardFusioncalls = hope_cohort_fusions, numCaller = 1)
caller_count <- caller_count %>%
  dplyr::select(-c(note)) %>%
  unique()
hope_cohort_fusions <- hope_cohort_fusions %>%
  left_join(caller_count)

# update column names to match those in OpenPedCan repo
hope_cohort_fusions <- hope_cohort_fusions %>%
  dplyr::rename("DomainRetainedGene1A" = "Gene1A_DOMAIN_RETAINED_IN_FUSION",
                "DomainRetainedGene1B" = "Gene1B_DOMAIN_RETAINED_IN_FUSION",
                "caller.count" = "caller_count") %>%
  dplyr::select(-c(id, SpanningFragCount, JunctionReadCount, Confidence, SpanningDelta, Caller))

# add Kids_First_Participant_ID using manifest for fusion files
fusion_manifest <- readr::read_tsv(file.path(input_dir, "manifest", "manifest_20230830_150316_fusion.tsv"))
colnames(fusion_manifest) <- gsub(" ", "_", colnames(fusion_manifest))
fusion_manifest <- fusion_manifest %>%
  dplyr::select(Kids_First_Biospecimen_ID, Kids_First_Participant_ID) %>%
  unique()
hope_cohort_fusions <- hope_cohort_fusions %>%
  inner_join(fusion_manifest, by = c("Sample" = "Kids_First_Biospecimen_ID")) 

# write output
saveRDS(hope_cohort_fusions, file = file.path(output_dir, "Hope-fusion-putative-oncogenic.rds"))

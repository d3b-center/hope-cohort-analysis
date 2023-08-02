# script to annotate fusion data with PFAM domain information using annoFuse
suppressPackageStartupMessages({
  library(tidyverse)
  library(annoFuse)
})

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
output_dir <- file.path(data_dir, "merged_files")

# get pfam data from annoFuse
bioMartDataPfam <- readRDS(system.file("extdata", "pfamDataBioMart.RDS", package = "annoFuseData"))

# pull hope cohort fusion files
hope_cohort_fusions <- readRDS(file.path(output_dir, "fusions_merged.rds"))
hope_cohort_fusions <- hope_cohort_fusions %>%
  dplyr::rename("Sample" = "Kids_First_Biospecimen_ID") # to keep column names consistent with OpenPedCan repo

# add unique id to join domain info obtained from gene1a and gene1b
hope_cohort_fusions <- hope_cohort_fusions %>%
  rownames_to_column("id") 

# get domain info on standard fusion calls
df <- annoFuse::get_Pfam_domain(standardFusioncalls = hope_cohort_fusions, bioMartDataPfam = bioMartDataPfam)

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

# add Kids_First_Participant_ID
manifest <- list.files(path = file.path(data_dir, "manifest"), pattern = "manifest.*fusion.tsv", full.names = T)
manifest <- readr::read_tsv(manifest)
manifest <- manifest %>%
  dplyr::select(`Kids First Biospecimen ID`, `Kids First Participant ID`) %>%
  unique()
colnames(manifest) <- gsub(" ", "_", colnames(manifest))
hope_cohort_fusions <- hope_cohort_fusions %>%
  inner_join(manifest, by = c("Sample" = "Kids_First_Biospecimen_ID")) 

# write output
saveRDS(hope_cohort_fusions, file = file.path(output_dir, "Hope-fusion-putative-oncogenic.rds"))

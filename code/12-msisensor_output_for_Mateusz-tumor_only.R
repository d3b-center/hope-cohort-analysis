# msisensor pro
fname <- 'results/msisensor-pro-tumor-only/hope_cohort_msi_sensor_output.tsv'
msi_output <- read_tsv(fname)
msi_output <- msi_output %>%
  dplyr::select(-c(Type)) %>%
  dplyr::rename("MSI_Percent" = "Percent")

# get TMB from OT
tmb_output <- read_tsv('~/Projects/OpenPedCan-analysis/analyses/tmb-calculation/results/snv-mutation-tmb-all.tsv')
tmb_output <- tmb_output %>%
  dplyr::rename("TMB" = "tmb") %>%
  dplyr::select(Tumor_Sample_Barcode, TMB)
msi_output <- msi_output %>%
  left_join(tmb_output, by = c("Kids First Biospecimen ID" = "Tumor_Sample_Barcode"))

# add positive controls from JL paper
# positive controls
pos_controls <- readxl::read_xlsx('data/media-4.xlsx')
pos_controls <- pos_controls %>%
  dplyr::select(`SAMPLE ID`, MMR_SOMATIC, MMR_GERMLINE) %>%
  dplyr::rename("sample_id" = "SAMPLE ID",
                "MMR_Somatic" = "MMR_SOMATIC",
                "MMR_Germline" = "MMR_GERMLINE")
msi_output <- msi_output %>%
  left_join(pos_controls, by = "sample_id")

# add proteomics cluster 
proteomics <- read_tsv('../d3b-miRNA-analysis/analyses/actmir-analysis/input/cluster_data_090722.tsv')
proteomics <- proteomics %>%
  dplyr::select(id, rdt.cc) 
dev_clusters <- read_tsv('data/cluster_data101922.tsv')
dev_clusters <- dev_clusters %>%
  dplyr::select(id, dtt.cc, age, rdt.name) 
cluster_meta <- proteomics %>%
  inner_join(dev_clusters)
msi_output <- msi_output %>%
  left_join(cluster_meta, by = c("sample_id" = "id")) 
write_tsv(msi_output, file = "results/msisensor-pro-tumor-only/msi_output_merged_tumor_only.tsv")

## Purpose

The goal of this module is to generate downstream plots for  `% microsatellite instability` generated using [MSIsensor-pro](https://github.com/xjtu-omics/msisensor-pro).

## Input files 

1. Input files from `data/` directory:

```
../../data/v1
├── Hope-GBM-histologies.tsv
└── Hope-gene-expression-rsem-tpm-collapsed.rds # gene expression data

../../data
└── gencode.v39.primary_assembly.annotation.gtf.gz
```

2. Module specific input files:

```
input
├── cluster_data_090722.tsv # proteomics cluster data
├── cluster_data_101922.tsv # developmental cluster data
├── manifest
│	├── manifest_20230504_083829_msi_tumor_only.tsv
│	└── manifest_20230504_084539_msi.tsv
├── msisensor_pro # input files from cavatica for T/N pipeline
│	├── ...
│	└── {uuid}_msisensor_pro
└── msisensor_pro_tumor_only # input files from cavatica for tumor-only pipeline
    ├── ...
    └── {uuid}_tumor_only_msisensor_pro
```

## Run script

This module can be run using the following bash script:

```
bash run_msisensor_analysis.sh
```

## Description of scripts

### 1) Merge MSI sample-level files to cohort-level files

`00-merge_msi_files.R`: This script merges individual sample-level files generated using MSIsensor-pro `T/N paired pipeline` and `Tumor-only pipeline` into a single cohort-level file.

```
# T/N paired pipeline
results/msisensor-pro-paired
└── Hope-msi-paired.tsv

# tumor-only pipeline
results/msisensor-pro-tumor-only
└── Hope-msi-tumor_only.tsv
```

### 2) MSI (T/N paired) vs MSI (tumor-only)

`01-msi_paired_vs_tumor_only.R`: This script generates a scatter plot of `% microsatellite instability` from T/N paired vs Tumor-only MSIsensor-pro pipeline run on common samples. 

Following are the output files for this script:

```
msi-sensor-analysis/results/msisensor-pro-combined
└── tumor_only_vs_paired_analysis.pdf
```

### 3) MSI (T/N paired) vs clinical variables

`02-msi_paired_comparisons.R`:  This script generates a correlations of  `% microsatellite instability` from T/N MSisensor-pro pipeline vs. various clinical variables like `TMB`, `Proteomics clusters (rdt_cc)`, `Developmental clusters (dtt_cc)`, `Age (three groups)`, `Age (two groups)`, `Developmental cluster name (rdt_name)`, `Sex`, `ALT status` and `Telomere content`. 

Following are the output files for this script:

```
results/msisensor-pro-paired
├── msi_vs_age_three_groups.png
├── msi_vs_age_two_groups.png
├── msi_vs_alt_status.png
├── msi_vs_dev_cluster_name.png
├── msi_vs_dev_clusters.png
├── msi_vs_gender.png
├── msi_vs_proteomics_clusters.png
├── msi_vs_telomere_content.png
└── msi_vs_tmb.png
```

### 4) MSI (tumor-only) vs clinical variables

`03-msi_tumor_only_comparisons.R`:  This script generates a correlations of  `% microsatellite instability` from Tumor-only MSisensor-pro pipeline vs. various clinical variables like `TMB`, `Proteomics clusters (rdt_cc)`, `Developmental clusters (dtt_cc)`, `Age (three groups)`, `Age (two groups)`, `Developmental cluster name (rdt_name)`, and `Sex`. 

Following are the output files for this script:

```
results/msisensor-pro-tumor-only
├── msi_vs_age_three_groups.png
├── msi_vs_age_two_groups.png
├── msi_vs_dev_cluster_name.png
├── msi_vs_dev_clusters.png
├── msi_vs_gender.png
├── msi_vs_proteomics_clusters.png
└── msi_vs_tmb.png
```

### 5) MSI (T/N paired) and MSI (tumor-only) vs MMR pathways

`04-msi_vs_mmr_pathways.R`: This script generates a scatter plot of `% microsatellite instability` from T/N paired and Tumor-only MSIsensor-pro pipeline vs GSVA scores obtained on TPM data for MMR pathways like `KEGG_MISMATCH_REPAIR`, `KEGG_BASE_EXCISION_REPAIR`, `KEGG_HOMOLOGOUS_RECOMBINATION`.

Following are the output files for this script:

```
# tumor-only
results/msisensor-pro-tumor-only
└── msi_vs_mmr_pathways.pdf

# paired
results/msisensor-paired
└── msi_vs_mmr_pathways.pdf
```

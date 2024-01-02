#### Authors: 

Komal S. Rathi

### Purpose

The goal of this module is to generate downstream plots for  `% microsatellite instability` generated using [MSIsensor-pro](https://github.com/xjtu-omics/msisensor-pro).

### Description of scripts

***
`00-merge_msi_files.R`: This script merges individual sample-level files generated using MSIsensor-pro `T/N paired pipeline` and `Tumor-only pipeline` into a single cohort-level file.

#### Inputs

```
# data files
../../data
└── Hope-GBM-histologies.tsv

# input files
input
├── manifest
│   ├── manifest_20230504_083829_msi_tumor_only.tsv
│   └── manifest_20230504_084539_msi.tsv
├── msisensor_pro # input files from cavatica for T/N pipeline
│   ├── ...
│   └── {uuid}_msisensor_pro
└── msisensor_pro_tumor_only # input files from cavatica for tumor-only pipeline
    ├── ...
    └── {uuid}_tumor_only_msisensor_pro
```

#### Outputs

```
# T/N paired pipeline
results/msisensor-pro-paired
└── Hope-msi-paired.tsv

# Tumor-only pipeline
results/msisensor-pro-tumor-only
└── Hope-msi-tumor_only.tsv
```
***

`01-msi_paired_vs_tumor_only.R`: This script generates a scatter plot of `% microsatellite instability` from T/N paired vs Tumor-only MSIsensor-pro pipeline run on common samples. 

#### Inputs

```
# T/N paired pipeline
results/msisensor-pro-paired
└── Hope-msi-paired.tsv

# Tumor-only pipeline
results/msisensor-pro-tumor-only
└── Hope-msi-tumor_only.tsv
```

#### Outputs

```
msi-sensor-analysis/results/msisensor-pro-combined
└── tumor_only_vs_paired_analysis.pdf
```
***
`02-msi_paired_comparisons.R`:  This script generates a correlations of  `% microsatellite instability` from T/N MSisensor-pro pipeline vs. various clinical variables like `TMB`, `Proteomics clusters (rdt_cc)`, `Developmental clusters (dtt_cc)`, `Age (three groups)`, `Age (two groups)`, `Developmental cluster name (rdt_name)`, `Sex`, `ALT status` and `Telomere content`. 

#### Inputs

```
../../data
└── Hope-GBM-histologies.tsv

# % MSI from T/N pipeline
results/msisensor-pro-paired
└── Hope-msi-paired.tsv

# TMB values for T/N samples
../tmb-calculation/results/wgs_paired
└── snv-mutation-tmb-coding.tsv

# ALT-status output
../alt-analysis/results
└── alt_status_aya_hgg.tsv

# developmental cluster data
input
└── cluster_data_101922.tsv
```

#### Outputs

```
results/msisensor-pro-paired
├── msi_vs_age_three_groups.pdf
├── msi_vs_age_two_groups.pdf
├── msi_vs_alt_status.pdf
├── msi_vs_dev_cluster_name.pdf
├── msi_vs_gender.pdf
├── msi_vs_telomere_content.pdf
└── msi_vs_tmb.pdf
```
***
`03-msi_tumor_only_comparisons.R`:  This script generates a correlations of  `% microsatellite instability` from Tumor-only MSisensor-pro pipeline vs. various clinical variables like `TMB`, `Proteomics clusters (rdt_cc)`, `Developmental clusters (dtt_cc)`, `Age (three groups)`, `Age (two groups)`, `Developmental cluster name (rdt_name)`, and `Sex`. 

#### Inputs

```
../../data
└── Hope-GBM-histologies.tsv

# % MSI for tumor-only pipeline
results/msisensor-pro-tumor-only
└── Hope-msi-tumor_only.tsv

# TMB values for tumor-only samples
../tmb-calculation/results/wgs_tumor_only
└── snv-mutation-tmb-coding.tsv

# developmental cluster data
input
└── cluster_data_101922.tsv
```

#### Outputs

```
results/msisensor-pro-tumor-only
├── msi_vs_age_three_groups.pdf
├── msi_vs_age_two_groups.pdf
├── msi_vs_dev_cluster_name.pdf
├── msi_vs_gender.pdf
└── msi_vs_tmb.pdf
```
***
`04-msi_vs_mmr_pathways.R`: This script generates a scatter plot of `% microsatellite instability` from T/N paired and Tumor-only MSIsensor-pro pipeline vs GSVA scores obtained on TPM data for MMR pathways like `KEGG_MISMATCH_REPAIR`, `KEGG_BASE_EXCISION_REPAIR`, `KEGG_HOMOLOGOUS_RECOMBINATION`.

#### Inputs

```
../../data
├── gencode.v39.primary_assembly.annotation.gtf.gz
├── Hope-GBM-histologies.tsv
└── Hope-gene-expression-rsem-tpm-collapsed.rds 

# results from MSI paired pipeline
results/msisensor-pro-paired
└── Hope-msi-paired.tsv

# results from MSI tumor-only pipeline
results/msisensor-pro-tumor-only
└── Hope-msi-tumor_only.tsv
```

#### Outputs

```
# tumor-only
results/msisensor-pro-tumor-only
└── msi_vs_mmr_pathways.pdf

# paired
results/msisensor-paired
└── msi_vs_mmr_pathways.pdf
```

### Run script

This module can be run using the following bash script:

```
bash run_msisensor_analysis.sh
```



### Introduction

The goal of this module is to generate plots for  `% microsatellite instability` generated using [MSIsensor-pro](https://github.com/xjtu-omics/msisensor-pro).

### Input files 

The input files for this analysis are MSIsensor-pro output files generated on [cavatica](https://cavatica.sbgenomics.com/u/cavatica/project-hope/files).  The files are then merged into single `rds` files:

```
../merge-files/results
├── Hope-msi-paired.rds
└── Hope-msi-tumor_only.rds
```

### Scripts

`01-msi_paired_vs_tumor_only.R`: This script generates a scatter plot of `% microsatellite instability` from T/N paired vs Tumor-only MSIsensor-pro pipeline run on common samples. 

Following are the output files for this script:

```
msi-sensor-analysis/results/msisensor-pro-combined
└── tumor_only_vs_paired_analysis.pdf
```

`02-msi_paired_comparisons.R`:  This script generates a correlations of  `% microsatellite instability` from T/N MSisensor-pro pipeline vs. various clinical variables like `TMB`, `Proteomics clusters (rdt_cc)`, `Developmental clusters (dtt_cc)`, `Age (three groups)`, `Age (two groups)`, `Developmental cluster name (rdt_name)`, `Sex`, `ALT-status` and `Telomere content`. 

Following are the output files for this script:

```
msi-sensor-analysis/results/msisensor-pro-paired
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

`03-msi_tumor_only_comparisons.R`:  This script generates a correlations of  `% microsatellite instability` from Tumor-only MSisensor-pro pipeline vs. various clinical variables like `TMB`, `Proteomics clusters (rdt_cc)`, `Developmental clusters (dtt_cc)`, `Age (three groups)`, `Age (two groups)`, `Developmental cluster name (rdt_name)`, `Sex`, `ALT-status` and `Telomere content`.. 

Following are the output files for this script:

```
msi-sensor-analysis/results/msisensor-pro-tumor-only
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

`04-msi_vs_mmr_pathways.R`: This script generates a scatter plot of `% microsatellite instability` from T/N paired and Tumor-only MSIsensor-pro pipeline vs GSVA scores obtained on TPM data for MMR pathways like `KEGG_MISMATCH_REPAIR`, `KEGG_BASE_EXCISION_REPAIR`, `KEGG_HOMOLOGOUS_RECOMBINATION`.

Following are the output files for this script:

```
# tumor-only
msi-sensor-analysis/results/msisensor-pro-tumor-only
└── msi_vs_mmr_pathways.pdf

# paired
msi-sensor-analysis/results/msisensor-paired
└── msi_vs_mmr_pathways.pdf
```

This module can be run using the following bash script:

```
bash run_msisensor_analysis.sh
```

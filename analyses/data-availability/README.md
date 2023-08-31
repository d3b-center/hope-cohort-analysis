
### Introduction

The goal of this module is to generate circular and linear heatmaps representing clinical variables and sample availability.

### Input files 

The input files for this module are files obtained from various sources: 

```
data-availability/input
├── CPTAC_Project Hope cohorts.xlsx
├── methylation_subset.tsv
└── single_cell_smartseq_manifest.xlsx
```

### Scripts

`01-clinical_data_availability_diagnosis.R`: This script generates a circular heatmap with clinical variables like `Age`, `Sex`, `Diagnosis`, `Diagnosis Type`, `Annotation`, `Tumor Location`. There are two plots generated, one with Age divided into two categorical groups and another with Age divided into three categorical groups. 

Following are the output files for this script:

```
data-availability/results
├── hope_cohort_data_availability_clinical_three_groups.pdf
└── hope_cohort_data_availability_clinical_two_groups.pdf
```

`02-clinical_data_availability_diagnosis_v2.R`:  This script generates a circular heatmap with clinical variables like `Age`, `Sex`, `Diagnosis`, `Diagnosis Type`, `Annotation`, `Tumor Location`. In this heatmap, `Age` is represented as a continuous variable.

Following are the output files for this script:

```
data-availability/results
└── hope_cohort_data_availability_clinical_v2.pdf
```

`03-clinical_data_availability_WHO.R`:  This script generates a circular heatmap with clinical variables like `Age`, `Sex`, `WHO Grade`, `Diagnosis Type`, `Annotation`, `Tumor Location`. There are two plots generated, one with Age divided into two categorical groups and another with Age divided into three categorical groups. 

Following are the output files for this script:

```
data-availability/results
├── hope_cohort_data_availability_clinical_with_WHOGrade_three_groups.pdf
└── hope_cohort_data_availability_clinical_with_WHOGrade_two_groups.pdf
```

`04-clinical_data_availability_WHO_v2.R`: This script generates a circular heatmap with clinical variables like `Age`, `Sex`, `WHO Grade`, `Diagnosis Type`, `Annotation`, `Tumor Location`. In this heatmap, `Age` is represented as a continuous variable.

Following are the output files for this script:

```
data-availability/results
└── hope_cohort_data_availability_clinical_with_WHOGrade_v2.pdf
```

`05-plot_data_availability.R`:  This script generates a linear heatmap with sample availability of various data types like `Proteomics`, `Phosphoproteomics`, `WGS`, `RNAseq`, `Methylation`, `Single Cell RNAseq (SmartSeq2)` and `Single Cell RNAseq (10x)`.

Following are the output files for this script:

```
data-availability/results
└── hope_cohort_data_availability.pdf
```

This module can be run using the following bash script:

```
bash run_data_plots.sh
```

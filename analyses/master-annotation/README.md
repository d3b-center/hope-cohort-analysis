
### Introduction

The goal of this module is to generate a master annotation file with clinical and sample variables from various resources. 

### Input files 

The input files for this module are files obtained from various sources: 

```
master-annotation/input
├── alt_status_aya_hgg.tsv
├── clini_m_030722-for_Komal.xlsx
├── cluster_data101922.tsv
├── cluster_data_090722.tsv
└── hopeonly_clinical_table_011823.tsv
```

### Scripts

`01-create-master-annotation.R`: This scripts pulls variables of interest from the input files and generates a master annotation file.  

Following is the output of this script:
```
master-annotation/results
└── master_histology_hope_cohort.tsv
```

This module can be run using the following R script:

```
Rscript 01-create-master-annotation.R
```



### Authors: 

Komal S. Rathi
 
### Purpose

The goal of this module is to generate circular and linear heatmaps representing clinical and sample data availability.

### Inputs  

```
../../data
└── Hope-GBM-histologies.tsv
```

### Description of scripts

***

`01-clinical_data_availability_age_three_groups.R`: This script generates a circular heatmap with clinical variables like `Age` (three groups), `Sex`, `Diagnosis`, `Diagnosis Type`, `Annotation`, `Tumor Location`. 

#### Output

```
results
└── hope_clinical_data_availability_age_three_groups.pdf
```

***
`02-clinical_data_availability_age_two_groups.R`: This script generates a circular heatmap with clinical variables like `Age` (two groups), `Sex`, `Diagnosis`, `Diagnosis Type`, `Annotation`, `Tumor Location`. 

#### Output

```
results
└── hope_clinical_data_availability_age_two_groups.pdf
```

***

`03-clinical_data_availability_age_continuous.R`:  This script generates a circular heatmap with clinical variables like `Age`, `Sex`, `Diagnosis`, `Diagnosis Type`, `Annotation`, `Tumor Location`. In this heatmap, `Age` is represented as a continuous variable.

#### Output

```
results
└── hope_clinical_data_availability_age_continuous.pdf
```
***
`03-clinical_data_availability_age_continuous_weipingma.R`:  This script generates a circular heatmap with clinical variables like `Age`, `Sex`, `Diagnosis`, `Diagnosis Type`, `Annotation`, `Tumor Location`. In this heatmap, `Age` is represented as a continuous variable. This is an updated script by Weiping Ma (cc: weiping.ma@mssm.edu) which adds label to each annotation layer.

#### Output

```
results
└── hope_clinical_data_availability_age_continuous_weipingma.pdf
```
***
`04-sample_data_availability.R`:  This script generates a linear heatmap with sample availability of various data types like `Proteomics`, `Phosphoproteomics`, `WGS` (`WGS tumor-only` as smaller blocks within), `RNAseq`, `Methylation`, and `snRNASeq`.

#### Output

```
results
└── hope_sample_availability.pdf
```

### Run module

This module can be run using the following bash script:

```
bash run_data_plots.sh
```

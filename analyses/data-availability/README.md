## Author: Komal S. Rathi
 
### Purpose

The goal of this module is to generate circular and linear heatmaps representing clinical and sample data availability.

### Data version  

```
data/v1
└── Hope-GBM-histologies.tsv
```

### 1. Clinical data with diagnosis and Age as category

`01-clinical_data_availability_diagnosis_age_category.R`: This script generates a circular heatmap with clinical variables like `Age`, `Sex`, `Diagnosis`, `Diagnosis Type`, `Annotation`, `Tumor Location`. There are two plots generated, one with Age divided into two categorical groups and another with Age divided into three categorical groups. 

#### Output

```
results
├── hope_clinical_data_availability_diagnosis_age_three_groups.pdf
└── hope_clinical_data_availability_diagnosis_age_two_groups.pdf
```

### 2. Clinical data with diagnosis and Age as continuous variable

`02-clinical_data_availability_diagnosis_age_continuous.R`:  This script generates a circular heatmap with clinical variables like `Age`, `Sex`, `Diagnosis`, `Diagnosis Type`, `Annotation`, `Tumor Location`. In this heatmap, `Age` is represented as a continuous variable.

#### Output

```
results
└── hope_clinical_data_availability_diagnosis_age_continuous.pdf
```

### 3. Clinical data with WHO grade and Age as category

`03-clinical_data_availability_WHO_age_category.R`:  This script generates a circular heatmap with clinical variables like `Age`, `Sex`, `WHO Grade`, `Diagnosis Type`, `Annotation`, `Tumor Location`. There are two plots generated, one with Age divided into two categorical groups and another with Age divided into three categorical groups. 

#### Output

```
results
├── hope_clinical_data_availability_who_grade_age_three_groups.pdf
└── hope_clinical_data_availability_who_grade_age_two_groups.pdf
```

### 4. Clinical data with WHO grade and Age as continuous variable


`04-clinical_data_availability_WHO_age_continuous.R`: This script generates a circular heatmap with clinical variables like `Age`, `Sex`, `WHO Grade`, `Diagnosis Type`, `Annotation`, `Tumor Location`. In this heatmap, `Age` is represented as a continuous variable.

#### Output

```
results
└── hope_clinical_data_availability_who_grade_age_continuous.pdf
```

### 5. Sample data availability


`05-sample_data_availability.R`:  This script generates a linear heatmap with sample availability of various data types like `Proteomics`, `Phosphoproteomics`, `WGS` (`WGS tumor-only` as smaller blocks within), `RNAseq`, `Methylation`, and `snRNASeq`.

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

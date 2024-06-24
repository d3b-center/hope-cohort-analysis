
### Authors: 

Jo Lynne Rokita
Komal S. Rathi

### Purpose

The goal of this module is to compute correlations between ALT-status/Telomere content and various sample/clinical variables. 

### Description of scripts

***
`00-generate-telomere-ratio.R`: This script computes the telomere ratio using telomere content for tumor and normal samples output from Telomere Hunter.

This script can be run using the following R script:

```
Rscript 00-generate-telomere-ratio.R
```

#### Inputs

```
input
└── 76_pairs_tumor_normal_summary.tsv
```

#### Outputs

```
results
└── alt_status_aya_hgg.tsv
```

***

`01-correlation_alt_vs_vars.R`: This scripts computes the correlation between ALT-status vs various clinical variables and T/N telomere content vs various clinical variables like Age (two and three groups), Gender, Protein clusters, and TMB.

This script can be run using the following R script:

```
Rscript 01-correlation_alt_vs_vars.R
```

#### Inputs

```
../../data
└── Hope-GBM-histologies.tsv

input
└── cluster_data_101922.tsv
```

#### Outputs 

```
results
├── alt_status_chisq_output.tsv # chisq output of ALT-status vs Age (2 and 3 groups), gender and protein clusters
├── alt_status_vs_tmb.pdf # ALT status vs TMB boxplot
├── telomere_content_vs_age_three_groups.pdf # Telomere content vs Age (3. groups) boxplot
├── telomere_content_vs_age_two_groups.pdf # Telomere content vs Age (2 groups) boxplot
├── telomere_content_vs_alt_status.pdf # Telomere content vs ALT-status boxplot
├── telomere_content_vs_dev_cluster_name.pdf # Telomere content vs Protein clusters boxplot
├── telomere_content_vs_gender.pdf # Telomere content vs Gender boxplot
└── telomere_content_vs_tmb.pdf # Telomere content vs TMB scatter plot
```

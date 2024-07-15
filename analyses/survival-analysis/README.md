
### Authors

Komal S. Rathi

### Purpose

The goal of this module is to generate univariate and multivariate survival analysis.

### Description of scripts
***
`01-survival-analysis.R`: This script generates survival curves for all molecular subtypes with `>= 3` samples as well as all samples. Univariate analysis was performed using `Kaplan-Meier` and multivariate analysis were performed using `Cox proportional hazards regression` model.

Command to run this script:

```
Rscript 01-survival-analysis.R
```

#### Inputs 

```
../../data
└── Hope-GBM-histologies.tsv

# ALT-status and telomere content
../alt-analysis/results
└── alt_status_aya_hgg.tsv
```

#### Outputs

```
results
├── alt_status_vs_survival_multivariate.pdf # multivariate analysis with ALT status and molecular subtype
├── alt_vs_survival_km.pdf # univariate analysis with ALT status only
├── msi_paired_vs_survival_multivariate.pdf # multivariate analysis with % microsatellite instability and molecular subtype 
├── msi_status_vs_survival_km.pdf # univariate analysis with microsatellite instability status only
├── msi_status_vs_survival_multivariate.pdf # multivariate analysis with microsatellite instability status and molecular subtype 
└── t_n_telomere_content_vs_survival_multivariate.pdf # univariate analysis with telomere content only
```
***
`02-summary-statics-TP53.Rmd`: This script generates the summary statistics for `TP53 score` in HOPE cohort. Analyse is designed to identify the relationship between `TP53 score` and other clinical features, such as `gender`, `race`, and `molecular subtype`. The plots are save in `plot` folder. Data was subset into two files, `HGG_survial_tp53.tsv` and `GBM_survial_tp53.tsv`, based on the cancer group, which can be found in `input` folder. 

This script can be run using the following bash script:

```
Rscript -e "rmarkdown::render('02-summary-statics-TP53.Rmd')"
```
***
`03-survival-tp53.Rmd`: This script generates the survival analysis in HOPE cohort, HGG subtype. `Kaplan-Meier` and `Hazard ratio` are plotted to show the impact of factors, such as `molecular subtyping`, `TP53 score`, `TP53 status` and other features, on overall survival. The plots are save in `plot` folder. 

This script can be run using the following bash script:

```
Rscript -e "rmarkdown::render('03-survival-tp53.Rmd')"
```

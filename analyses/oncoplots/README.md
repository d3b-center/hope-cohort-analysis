### Authors

Komal S. Rathi

### Purpose

The goal of this module is to generate oncoplots and cascade plots.

### Description of scripts

***
`01-deg-vs-gtex-brain.R`: The function of this script is to perform single-sample differential analysis for each sample in the HOPE cohort vs the GTEx normal brain tissue cohort using the R package [NOISeq](https://www.bioconductor.org/packages/release/bioc/html/NOISeq.html). 

Command to run this script:

```
Rscript --vanilla 01-deg-vs-gtex-brain.R
```

#### Inputs 

```
# GTEx brain tissue data from https://github.com/d3b-center/OpenPedCan-analysis
data/gene-counts-rsem-expected_count-collapsed.rds 

# data 
../../data
├── Hope-gene-counts-rsem-expected_count-collapsed.rds
└── Hope-GBM-histologies.tsv
```

#### Outputs

```
results
└── hope_cohort_vs_gtex_brain_degs.rds
```
***
`02-prepare_files_oncogrid.R`: This script generates the input matrix and annotation file required to generate various oncoplots and cascade plots.

Command to run this script:

```
Rscript --vanilla 02-prepare_files_oncogrid.R
```

#### Inputs 

```
# specific gene lists
input
├── cnv_genes.tsv
├── deg_genes.tsv
├── fusion_genes.tsv
├── hgat_goi_list.tsv -> OpenPedCan-analysis/analyses/oncoprint-landscape/data/hgat_goi_list.tsv
└── snv_genes.tsv

# data files
../../data
├── Hope-GBM-histologies.tsv
├── Hope-cnv-controlfreec.rds
├── Hope-snv-consensus-plus-hotspots.maf.tsv.gz
└── Hope-fusion-putative-oncogenic.rds

# DEG output
results
└── hope_cohort_vs_gtex_brain_degs.rds

# TMB information
../tmb-calculation/results/wgs_paired
└── snv-mutation-tmb-coding.tsv
```

#### Outputs

```
results
├── annotation.txt
└── oncoprint.txt
```
***
`03-prepare_files_oncogrid_add_tumor_only.R`: This script generates the input matrix and annotation file required to generate various oncoplots and cascade plots. Here we add SNVs and CNVs for patients that have only tumor-only samples available and T/N paired samples are absent. 

Command to run this script:

```
Rscript --vanilla 03-prepare_files_oncogrid_add_tumor_only.R
```

#### Inputs 

```
# specific gene lists
input
├── cnv_genes.tsv
├── deg_genes.tsv
├── fusion_genes.tsv
├── hgat_goi_list.tsv -> OpenPedCan-analysis/analyses/oncoprint-landscape/data/hgat_goi_list.tsv
└── snv_genes.tsv

# data files
../../data
├── Hope-GBM-histologies.tsv
├── Hope-cnv-controlfreec.rds
├── Hope-cnv-controlfreec-tumor-only.rds
├── Hope-snv-consensus-plus-hotspots.maf.tsv.gz
├── Hope-tumor-only-snv-mutect2.maf.tsv.gz
└── Hope-fusion-putative-oncogenic.rds

# DEG output
results
└── hope_cohort_vs_gtex_brain_degs.rds

# TMB information
../tmb-calculation/results/wgs_paired
└── snv-mutation-tmb-coding.tsv
../tmb-calculation/results/wgs_tumor_only
└── snv-mutation-tmb-coding.tsv
```

#### Outputs

```
results
├── annotation_add_tumor_only.txt
└── oncoprint_add_tumor_only.txt
```
***
`04-plot_oncogrid.R`: Script to generate oncoplots.

Command to run this script:

```
Rscript --vanilla 04-plot_oncogrid.R
```

#### Inputs 

```
results
├── annotation.txt
└── oncoprint.txt
```

#### Outputs

```
results/oncoplots
├── oncoplot.pdf
├── oncoplot_orderby_H3F3A.pdf
├── oncoplot_orderby_sex.pdf
├── oncoplot_orderby_sex_H3F3A_status.pdf
├── oncoplot_orderby_sex_age.pdf
└── oncoplot_orderby_sex_age_H3F3A_status.pdf
```
***
`05-plot_oncogrid_no_RNA.R`: Script to generate oncoplots without gene expression information.

Command to run this script:

```
Rscript --vanilla 05-plot_oncogrid_no_RNA.R
```

#### Inputs 

```
results
├── annotation.txt
└── oncoprint.txt
```

#### Outputs

```
results/oncoplots
├── oncoplot_norna.pdf
├── oncoplot_orderby_H3F3A_norna.pdf
├── oncoplot_orderby_age_H3F3A_status_norna.pdf
├── oncoplot_orderby_age_norna.pdf
├── oncoplot_orderby_age_sex_H3F3A_status_norna.pdf
├── oncoplot_orderby_sex_H3F3A_status_norna.pdf
├── oncoplot_orderby_sex_age_H3F3A_status_norna.pdf
├── oncoplot_orderby_sex_age_norna.pdf
└── oncoplot_orderby_sex_norna.pdf
```
***
`06-cascade_plots.R`: Script to generate cascade plots.

Command to run this script:

```
Rscript --vanilla 06-cascade_plots.R
```

#### Inputs 

```
results
├── annotation.txt
└── oncoprint.txt
```

#### Outputs

```
results/cascade_plots
├── cascade_orderby_age.pdf
├── cascade_orderby_sex.pdf
└── cascade_plot.pdf
```
***
`07-cascade_plots_add_tumor_only.R`: Script to generate cascade-style oncoplots with tumor-only samples added for patients where T/N samples are unavailable.

Command to run this script:

```
Rscript --vanilla 07-cascade_plots_add_tumor_only.R
```

#### Inputs 

```
results
├── annotation_add_tumor_only.txt
└── oncoprint_add_tumor_only.txt
```

#### Outputs

```
results/cascade_plots_add_tumor_only
├── cascade_orderby_age.pdf
├── cascade_orderby_sex.pdf
└── cascade_plot.pdf
```
***
`08-major_snv_analysis.R`: This script computes correlations between ALT-status and % MSI to major SNV i.e. `> 6%`. 

Command to run this script:

```
Rscript --vanilla 08-major_snv_analysis.R
```

#### Inputs 

```
results
├── annotation.txt
└── oncoprint.txt

# ALT-status
../alt-analysis/results
└── alt_status_aya_hgg.tsv

# % MSI 
../msi-sensor-analysis/results/msisensor-pro-paired
└── Hope-msi-paired.tsv
```

#### Outputs

```
results/major_snv
└── alt_msi_correlation_with_snv.tsv
```
***
`09-gene_alteration_correlation.R`: This script computes coexistence analysis between pairs of gene alterations and with various clinical variables.

Command to run this script:

```
Rscript --vanilla 09-gene_alteration_correlation.R
```

#### Inputs 

```
results
├── annotation.txt
└── oncoprint.txt
```

#### Outputs

```
results/correlation_analysis
├── ATRX-TP53-coexistence.tsv
├── ATRX-vs-TP53-by-age-and-sex.tsv
├── H3-3A-TP53-coexistence.tsv
├── H3-3A-vs-TP53-by-age-and-sex.tsv
├── NF1-ATRX-coexistence.tsv
├── NF1-vs-ATRX-by-age-and-sex.tsv
└── gene_correlation_with_age_or_sex.tsv
```


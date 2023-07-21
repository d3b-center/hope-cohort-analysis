# HOPE cohort analysis

Module authors: [Komal S. Rathi](https://github.com/komalsrathi/)

## Run the script on docker

To pull the docker image, run the command line:
```
docker pull pgc-images.sbgenomics.com/zhuangzhuanggeng/d3b_hope_analysis:last
```

To run the docker, run the command line below. For mac M1 user, add `--platform=linux/arm64`.
```
docker run -e PASSWORD=pass -p 8787:8787 --name <CONTAINER_NAME> -v $PWD:/home/rstudio/HOPE pgc-images.sbgenomics.com/zhuangzhuanggeng/d3b_hope_analysis:last

```

## Structure

code
├── 00-create_merged_files.R 
├── 01-deg-vs-gtex-brain.R
├── 02-prepare_files_oncogrid.R
└── 03-plot_oncogrid.R
```

1. 00-create_merged_files.R: Script to merge files obtained from cavatica. Merge fusions, copy number (controlfreec), consensus mutations and RSEM gene expressions.
2. 01-deg-vs-gtex-brain.R: Script to identify differentially expressed genes in HOPE cohort vs GTEx Brain tissue samples. We used TPM from both datasets, normalized to z-score, and used a cut-off of 1.5 to get up/down genes in individual samples from HOPE cohort.
3. 02-prepare_files_oncogrid.R: Script to prepare input files for oncoplot generation. Reference files and genelists were obtained from https://github.com/d3b-center/d3b-pnoc003-HGG-DMG-omics/tree/master/analyses/Oncoplot.
4. 03-plot_oncogrid.R: Script to generate oncoplot.

## Results

Please refer to the below folder descriptions to find the corresponding results:

```
results
├── alt_correlations.tsv # ALT correlations with Age (two and three groups), Sex, MSI, TMB and Protein clusters
├── correlation_analysis_three_groups # Gene alteration correlations with clinical variables
├── correlation_analysis_two_groups # Gene alteration correlations with clinical variables
├── data_plots # Data availability and Clinical data availability plots
├── msisensor-pro # MSI sensor pro plots/results for tumor-normal paired data
├── msisensor-pro-tumor-only # MSI sensor pro plots/results for tumor-only data
├── msisensor-pro-combined # Comparison of MSI sensor pro tumor-normal paired vs tumor-only data
├── oncoplots_three_groups # oncoplots using three age groups
└── oncoplots_two_groups # oncoplots using two age groups
```


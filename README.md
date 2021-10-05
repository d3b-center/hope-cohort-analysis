# HOPE cohort analysis

Module authors: [Komal S. Rathi](https://github.com/komalsrathi/)

## Introduction

RNA-sequencing samples: n = 74
WGS samples: n = 57

## Structure

```
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
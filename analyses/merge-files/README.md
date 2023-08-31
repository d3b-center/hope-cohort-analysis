
### Introduction

The goal of this module is to merge sample-level files into project-level files. 

### Input files 

The input files are downloaded from Project HOPE on Cavatica: https://cavatica.sbgenomics.com/u/cavatica/project-hope. 

```
merge-files/input
├── consenus_maf
│	└── {id}.consensus_somatic.norm.annot.public.maf
├──mutect_maf_tumor_only
│   └── {id}.mutect2.norm.annot.public.maf
├── copy_number
│	└── {id}.controlfreec.CNVs.p.value.txt
├── copy_number_tumor_only
│	└── {id}.controlfreec.CNVs.p.value.txt
├── gene_expression
│	└── {id}.rsem.genes.results.gz
├── gene_fusions
│	└── {id}.annoFuse_filter.tsv
├── msisensor_pro
│	└── {id}_msisensor_pro
├── msisensor_pro_tumor_only
│	└── {id}_tumor_only_msisensor_pro
├── manifest
│	├── manifest_20230504_084539_msi.tsv
│	├── manifest_20230830_145506_snv.tsv
│	├── manifest_20230830_150316_fusion.tsv
│	├── manifest_20230830_150931_rna.tsv
│	└── manifest_20230830_151211_cnv.tsv
└── manifest_tumor_only
	├── manifest_20230504_083829_msi.tsv
	├── manifest_20230830_151430_snv.tsv
	└── manifest_20230830_152155_cnv.tsv
```

### Scripts

`01-create_merged_files.R`: This scripts merges all files generated using the RSEM gene-expression quantification and T/N paired pipeline. Following are the outputs of this script:

```
merge-files/results
├── Hope-cnv-controlfreec.rds
├── Hope-gene-counts-rsem-expected_count-collapsed.rds
├── Hope-gene-counts-rsem-expected_count.rds
├── Hope-gene-expression-rsem-tpm-collapsed.rds
├── Hope-gene-expression-rsem-tpm.rds
├── Hope-msi-paired.rds
├── Hope-snv-consensus-plus-hotspots.maf.tsv.gz
```

`02-create_merged_files_tumor_only.R`: This scripts merges all files generated using Tumor-only pipeline. Following are the outputs of this script:
```
merge-files/results
├── Hope-cnv-controlfreec-tumor-only.rds
├── Hope-msi-tumor_only.rds
├── Hope-tumor-only-snv-mutect2.maf.tsv.gz
```

`03-filter-annotate-fusions.R`: This scripts merges all fusion files and applies [annoFuse](https://github.com/d3b-center/annoFuse) expression, artifact filters and kinase domain annotations.

. Following are the outputs of this script:
```
merge-files/results
├── Hope-fusion-putative-oncogenic.rds
```

The entire module can be run using the following bash script:

```
bash run_merge_files.sh
```

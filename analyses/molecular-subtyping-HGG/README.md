## Molecular subtyping for high-grade glioma (HGG) 

This module is adapted from: [`d3b-center/OpenPedCan-analysis`](https://github.com/d3b-center/OpenPedCan-analysis).


#### Running the analysis

The analysis can be run with the following (assuming you are in the analysis directory):

```
bash run-molecular-subtyping.sh
```

### Order of analysis

`00-fusion-summary.Rmd`: This script is adapted from `analyses/fusion-summary` and designed to generate fusion files specifically for HGG molecular subtyping. Input files needed for this script include `Hope-GBM-histologies-base.tsv`, `Hope-fusion-putative-oncogenic.rds`.  Output of this script is `input/fusion_summary_hgg_foi.tsv`, which contains the fusion that is needed for HGG molecular subtyping. 

`01-HOPE-HGG-subtyping-subset-file.Rmd`: This script is used to subset the files that are needed for HGG subtyping, including CNV, SV, SNV, and fusion data. 

####  output
  -  `results/Hope_HGG_defining_lesions.tsv`: This file contains the target lesions (H3 K28M/K28I/G35R/G35V) for HGG molecular subtyping.
  -  `results/HGG_cleaned_cnv.tsv`: This file contains the copy number variants (EGFR, PTEN, PDGFRA and MYCN) for subtyping. 
  -  `results/HGG_cleaned_mutation.tsv`: This file contains mutations related to `H3 K28 mutant`, `H3 G35 mutant`, `IDH mutant` and `H3.3 and IDH wildtype`. 
  -  `results/HGG_cleaned_fusion.tsv`: This file contains the fusion of interests (FGFR1, NTRK, MET, ROS1, and ALK). 
  -  `results/HGG_cleaned_expression.tsv`: This file contains the scaled expression of target genes (FOXG1, OLIG2, TP73_AS1, EGFR and EZHIP).

`02-HOPE-molecular-subtyping.Rmd`: This script follows WHO guidebook to assign the HGG molecular subtype. It takes the output from `01-HOPE-HGG-subtyping-subset-file.Rmd` and `tp53_nf1_score` module and assigns the molecular subtype to each HOPE cohort sample.

#### output
  -  `results/Hope_subtype.tsv`: This contains the subtyping results of each sample, with molecular subtype, methylation molecular subtype and TP53 status.

`03-molecular-subtype-integrate.R`: This script integrates the molecular subtyping results of samples from Hope cohort with `Hope-GBM-histologies-base.tsv` and generate `Hope-GBM-histologies.tsv`. 


#### Inputs files

All these files can be downloaded from s3, using `hope-cohort-analysis/download-data.sh` script. 

`Hope-GBM-histologies-base.tsv`
`Hope-fusion-putative-oncogenic.rds`
`Hope-tumor-only-snv-mutect2.maf.tsv.gz`
`Hope-snv-consensus-plus-hotspots.maf.tsv.gz`
`Hope-cnv-controlfreec-tumor-only.rds`
`Hope-cnv-controlfreec.rds`

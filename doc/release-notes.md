# Release note

## Current release (V2)
- release data: 2023-11-21
- Overview:
  - Sequence sample:
    - 184 RNA-seq
    - 87 WGS
    - 78 DNA methylation
    - 29 snRNA-seq
    - 91 Proteomic
- Updates: 
  - Added `cancer_group_short` column in final histology file (https://github.com/d3b-center/hope-cohort-analysis/pull/65)
  - Genomics data from two samples (`7316-1106` and `7316-3000`) were removed because they do not have proteomic data.
  - DNA methylation data of sample `7316-212` was removed from the dataset due to a sample swap issue. 
  - Proteomic data was not included in this release.  


## Previous release (V1)
- release date: 2023-09-01
- Overview: 
  - Sequence sample:
    - 186 RNA-seq
    - 91 WGS
    - 81 DNA methylation
    - 29 snRNA-seq
    - 91 Proteomic
    
  - Files added
    - Hope-GBM-histologies-base.tsv
    - Hope-GBM-histologies.tsv
    - Hope-gene-expression-rsem-tpm-collapsed.rds
    - Hope-and-CPTAC-GBM-gene-expression-rsem-tpm-collapsed.rds  
    - Hope-gene-expression-rsem-tpm.rds
    - Hope-and-CPTAC-GBM.splice-events-rmats.tsv.gz		   
    - Hope-snv-consensus-plus-hotspots.maf.tsv.gz
    - Hope-cnv-controlfreec-tumor-only.rds			   
    - Hope-sv-manta.tsv.gz
    - Hope-cnv-controlfreec.rds				   
    - Hope-tumor-only-snv-mutect2.maf.tsv.gz (filtered to remove `t_alt_count ==0` and `t_depth < 4` )
    - Hope-fusion-putative-oncogenic.rds
    - Hope-gene-counts-rsem-expected_count-collapsed.rds	   
    - Hope-gene-counts-rsem-expected_count.rds
    - Hope-methyl-beta-values.rds
    - Hope-methyl-m-values.rds

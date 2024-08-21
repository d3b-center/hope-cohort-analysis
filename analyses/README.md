
Order of analyses:

1. [merge-files](https://github.com/d3b-center/hope-cohort-analysis/tree/master/analyses/merge-files): This module has scripts to merge files obtained from cavatica i.e. RSEM gene expression, Consensus MAF, ControlFREEC, Fusions which are then filtered and annotated. The merged files are then used for data releases.
2. [tp53_nf1_score](https://github.com/d3b-center/hope-cohort-analysis/tree/master/analyses/tp53_nf1_score): The purpose of this module is to apply the TP53 inactivation, NF1 inactivation, and Ras activation classifiers trained on TCGA PanCan data to the HOPE cohort dataset.
3. [molecular-subtyping-HGG](https://github.com/d3b-center/hope-cohort-analysis/tree/master/analyses/molecular-subtyping-HGG): Adds `Cancer_Group`, `Cancer_Group_Short` and `molecular_subtype` to `../data/Hope-GBM-histologies-base.tsv` and generates `Hope-GBM-histologies.tsv`.
4. [tmb-calculation](https://github.com/d3b-center/hope-cohort-analysis/tree/master/analyses/tmb-calculation): TMB calculation for `T/N paired` and `Tumor-only` WGS samples. This module is adapted from [OpenPedCan-anaysis](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/tmb-calculation)
5. [master-annotation](https://github.com/d3b-center/hope-cohort-analysis/tree/master/analyses/master-annotation): This module generates a master annotation file by combining various sources of information from the HOPE group into one single tsv file for downstream analyses.
6. [data-availability](https://github.com/d3b-center/hope-cohort-analysis/tree/master/analyses/data-availability): This module has scripts to create clinical and sample data availability plots.
7. [alt-analysis](https://github.com/d3b-center/hope-cohort-analysis/tree/master/analyses/alt-analysis): Downstream analysis of `ALT status` and `telomere content` for `T/N paired` WGS samples. 
8. [msi-sensor-analysis](https://github.com/d3b-center/hope-cohort-analysis/tree/master/analyses/msi-sensor-analysis): Downstream analysis of MSI outputs for `T/N paired` and `Tumor-only` WGS samples. 
9. [survival-analysis](https://github.com/d3b-center/hope-cohort-analysis/tree/master/analyses/survival-analysis): This module has scripts to do survival analysis with `ALT status` and `Molecular subtypes`.  
10. [oncoplots](https://github.com/d3b-center/hope-cohort-analysis/tree/master/analyses/oncoplots): This module has scripts to create oncoplots and cascade style plots. Reference files and genelists were obtained from [PNOC003](https://github.com/d3b-center/d3b-pnoc003-HGG-DMG-omics/tree/master/analyses/Oncoplot)

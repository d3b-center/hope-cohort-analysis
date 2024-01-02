
# Tumor Mutation Burden calculation

This analysis utilizes the TMB calculation module from [Pediatric Open Targets, OPenPedCan-analysis](https://github.com/PediatricOpenTargets/OpenPedCan-analysis) to calculate Tumor Mutation Burden (TMB) for each tumor sample with SNV calls. 

The TMB calculation for T/N paired samples uses SNV calls from Mutect2, Strelka2, Lancet, and Vardict callers. On the other hand, TMB calculation for Tumor-only samples using SNV calls from Mutect2 caller only.

## Input files 

1. Input files from `data/` directory:

```
../../data/v1
├── Hope-GBM-histologies.tsv
├── Hope-snv-consensus-plus-hotspots.maf.tsv.gz # SNV consensus MAF file for T/N paired samples
└── Hope-tumor-only-snv-mutect2.maf.tsv.gz # SNV mutect2 MAF file for Tumor-only samples

../../data
└── gencode.v39.primary_assembly.annotation.gtf.gz
```

2. Module specific input files:

```
input
├── biospecimen_id_to_bed_map.tsv # generated per OpenPedCan-analysis repo
├── gencode.v39.primary_assembly.annotation.bed # generated in bash script
├── hg38_strelka.bed # copied from OpenPedCan-analysis repo
├── intersect_strelka_mutect2_vardict_WGS.bed # generated in bash script
└── wgs_canonical_calling_regions.hg38.bed # copied from OpenPedCan-analysis repo
```

## TMB Calculation

For each experimental strategy and TMB calculation, the intersection of the genomic regions effectively being surveyed are used. For TMB calculation using consensus MAF, all calls whether consensus or not are used for TMB calculations. These genomic regions are used for first filtering mutations to these regions and then for using the size in bp of the genomic regions surveyed as the TMB denominator.

### All mutations TMB

For all mutation TMBs, all callers are used. For WGS samples, the size of the genome covered by the intersection of Strelka2 and Mutect2's surveyed areas which are considered representative of all callers is used for the denominator.

```
WGS_all_mutations_TMB = (total # consensus and non-consensus snvs called by all callers) / intersection_strelka_mutect_genome_size
```

### Coding only TMB

Coding only TMB uses all callers as well and the intersection demoninators are calculated by using coding sequence ranges in the gtf from Gencode 39.
This file is included in the project's data download.
SNVs outside of these coding sequences are filtered out before being summed and used for TMB calculations as follows:

```
WGS_coding_only_TMB = (total # consensus and non-consensus snvs called by all callers) / intersection_wgs_strelka_mutect_CDS_genome_size
```

## Description of scripts

### 1) `run_tmb_calculation.sh`

This is a bash script wrapper for setting input file paths for the main analysis script, `01-calculate_tmb.R` and creating additional intermediate input files `input/intersect_strelka_mutect2_vardict_WGS.bed` and `input/gencode.v27.primary_assembly.annotation.bed`.  All file paths set in this script relative to the module directory. Therefore, this script should always run as if it were being called from the directory it lives in, the module directory (`d3b-hope-cohort-analysis/analyses/tmb-calculation`).

### 2) 01-calculate_tmb.R

Uses the snv consensus file for T/N paired samples and mutect2 file for tumor-only samples to calculate TMB for all WGS samples. Two TMB files are created per analysis, one including *all snv* and a *coding snvs only*. For T/N paired samples, outputs are generated using both `consensus` and `non-consensus mutations` called by all callers (`Mutect2`, `Strelka2`, `Lancet`, and `Vardict`) while for Tumor-only samples, outputs are generated using `Mutect2` calls.

#### Usage:

```
Rscript 01-calculate_tmb.R --help

Options:
    --maf_file=MAF_FILE
        path to maf file (tsv.gz)

    --bed_files=BED_FILES
        Input samples to target BED mapping file (.bed)

    --histologies_file=HISTOLOGIES_FILE
        histologies file

    --coding_regions=CODING_REGIONS
        BED file for coding regions to use for coding only TMB

    --nonsynfilter_maf
        If TRUE, filter out synonymous mutations, keep 
              non-synonymous mutations, according to maftools definition.
              Default is FALSE

    --nonsynfilter_focr
        If TRUE, filter out synonymous mutations, keep 
              non-synonymous mutations, according to Friends of Cancer 
              Research definition. Default is FALSE

    --results_dir=RESULTS_DIR
        path to output directory

    -h, --help
        Show this help message and exit
```

#### Outputs

```
results
├── wgs_paired
│   ├── snv-mutation-tmb-all.tsv
│   └── snv-mutation-tmb-coding.tsv
└── wgs_tumor_only
    ├── snv-mutation-tmb-all.tsv
    └── snv-mutation-tmb-coding.tsv
```

### 3) 02-compare_tmb.R

This is a QC script that compares the TMB of samples overlapping between OpenPedCan-analysis repo and d3b-hope-cohort-analysis repo. Following is the correlation value obtained by comparing 44 common samples between the two projects:

```
0.9999734
``` 

### 4) `Util/split_mnv.R`

Contains a function to split multinucleotide variants (MNVs) into single nucleotide variants (SNVs).

### 5) `Util/tmb_functions.R`

Contains functions for calculating Tumor Mutation Burden (TMB).

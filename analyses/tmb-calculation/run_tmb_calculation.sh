#!/bin/bash
# OpenPedCan 2021
# Eric Wafula
set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Set up paths to data files consumed by analysis
data_dir="../../data"
input_dir="input"

# BED and GTF file paths
cds_file="${input_dir}/gencode.v39.primary_assembly.annotation.bed"
wgs_bed="${input_dir}/intersect_strelka_mutect2_vardict_WGS.bed"

# sample to BED mapping file
mapping_file="${input_dir}/biospecimen_id_to_bed_map.tsv"

# Histology file
histology_file="${data_dir}/Hope-GBM-histologies.tsv"

############# Create intersection BED files for TMB calculations ###############
# Make All mutations BED files
bedtools intersect \
  -a ${input_dir}/hg38_strelka.bed \
  -b ${input_dir}/wgs_canonical_calling_regions.hg38.bed \
  > $wgs_bed

#################### Make coding regions file
# Convert GTF to BED file for use in bedtools
# Here we are only extracting lines with as a CDS i.e. are coded in protein
gunzip -c ${data_dir}/gencode.v39.primary_assembly.annotation.gtf.gz \
  | awk '$3 ~ /CDS/' \
  | convert2bed --do-not-sort --input=gtf - \
  | sort -k 1,1 -k 2,2n \
  | bedtools merge  \
  > $cds_file

# ######################### Calculate consensus TMB ##############################
Rscript 01-calculate_tmb.R \
--maf_file "${data_dir}/Hope-snv-consensus-plus-hotspots.maf.tsv.gz" \
--bed_files $mapping_file \
--histologies_file $histology_file \
--coding_regions $cds_file \
--nonsynfilter_maf \
--results_dir "results/wgs_paired"

# ######################### Calculate mutect2 TMB ##############################
Rscript 01-calculate_tmb.R \
--maf_file "${data_dir}/Hope-tumor-only-snv-mutect2.maf.tsv.gz" \
--bed_files $mapping_file \
--histologies_file $histology_file \
--coding_regions $cds_file \
--nonsynfilter_maf \
--results_dir "results/wgs_tumor_only"

#!/bin/bash
# OpenPedCan 2021
# Eric Wafula, updated by Jo Lynne Rokita, 2023
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
scratch_dir="../../scratch"

# BED and GTF file paths
cds_file="${scratch_dir}/gencode.v39.primary_assembly.annotation.bed"

#################### Make coding regions file
# Convert GTF to BED file for use in bedtools
# Here we are only extracting lines with as a CDS i.e. are coded in protein
gunzip -c ${data_dir}/gencode.v39.primary_assembly.annotation.gtf.gz \
  | awk '$3 ~ /CDS/' \
  | convert2bed --do-not-sort --input=gtf - \
  | sort -k 1,1 -k 2,2n \
  | bedtools merge  \
  > $cds_file

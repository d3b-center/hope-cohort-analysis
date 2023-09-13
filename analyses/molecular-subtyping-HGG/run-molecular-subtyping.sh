#!/bin/bash

set -e
set -o pipefail

if [ -f ../../scratch/gencode.v39.primary_assembly.annotation.bed ]; then
  echo "gencode.v39.primary_assembly.annotation.bed is already in scratch folder"
else
  gunzip -c ../../data/gencode.v39.primary_assembly.annotation.gtf.gz \
  | awk '$3 ~ /CDS/' \
  | convert2bed --do-not-sort --input=gtf - \
  > ../../scratch/gencode.v39.primary_assembly.annotation.bed
fi
  

## run fusion summary script
Rscript -e "rmarkdown::render('00-fusion-summary.Rmd', clean = TRUE)"

## Run 01 to subset files
Rscript -e "rmarkdown::render('01-HOPE-HGG-subtyping-subset-file.Rmd', clean = TRUE)"

## Run 02 script for molecular subtyping
Rscript -e "rmarkdown::render('02-HOPE-molecular-subtyping.Rmd', clean = TRUE)"

## Run 03 script to add pathology based on molecular subtype
Rscript 03-molecular-subtype-integrate.R



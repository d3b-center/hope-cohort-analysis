#!/bin/bash

set -e
set -o pipefail

## run fusion summary script
Rscript -e "rmarkdown::render('00-fusion-summary.Rmd', clean = TRUE)"

## Run 01 to subset files
Rscript -e "rmarkdown::render('01-HOPE-HGG-subtyping-subset-file.Rmd', clean = TRUE)"

## Run 02 script for molecular subtyping
Rscript -e "rmarkdown::render('02-HOPE-molecular-subtyping.Rmd', clean = TRUE)"

## Run 03 script to add pathology based on molecular subtype
Rscript 03-molecular-subtype-integrate.R



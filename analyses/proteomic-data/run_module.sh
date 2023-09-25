#!/bin/bash 

set -e
set -o pipefail

## Run R script to substitue the sample_id in proteomic/phospho data to bs_id
Rscript --vanilla 01-add-bsid-to-proteomic-data.R

## Run python script to liftover the gene symbols in proteomic/phospho data 

### whole cell proteomic data
python 01-update_gene_symbols.py \
  -g input/hgnc_complete_set.txt \
  -f output/Hope_proteome_imputed_data.tsv \
  -o output/Hope_proteome_imputed_data_liftover.tsv.gz \
  -u "ApprovedGeneSymbol"
  
### phospho data
python 01-update_gene_symbols.py \
  -g input/hgnc_complete_set.txt \
  -f output/Hope_phosphosite_imputed_data_ischemia_removed.tsv \
  -o output/Hope_phosphosite_imputed_data_ischemia_removed_liftover.tsv.gz \
  -u "ApprovedGeneSymbol" 

### phospho with motif data
python 01-update_gene_symbols.py \
  -g input/hgnc_complete_set.txt \
  -f output/Hope_phosphosite_imputed_data_ischemia_removed_motif.tsv \
  -o output/Hope_phosphosite_imputed_data_ischemia_removed_motif_liftover.tsv.gz \
  -u "ApprovedGeneSymbol" 
  
  
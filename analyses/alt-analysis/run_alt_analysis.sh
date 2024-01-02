# generate telomere ratio
Rscript --vanilla 00-generate-telomere-ratio.R

# compute correlations of alt status/telomere content with other clinical variables
Rscript --vanilla 01-correlation_alt_vs_vars.R

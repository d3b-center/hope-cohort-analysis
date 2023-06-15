# create data availability plots
bash run_data_plots.sh

# create master annotation
Rscript --vanilla code/02-create-master-annotation.R

# MSI analysis
bash run_msisensor_analysis.sh

# oncoplots
bash run_oncoplots.sh

# misc analysis
# gene alterations correlation
Rscript --vanilla code/08-gene_alteration_correlation.R

# MSI vs MMR pathways
Rscript --vanilla code/09-msi_vs_mmr_pathways.R

# ALT vs other variables
Rscript --vanilla code/10-correlation_alt_vs_vars.R

# cascade plots
Rscript --vanilla code/11-cascade_plots.R

# major SNV correlations
Rscript --vanilla code/12-major_snv_analysis.R

# survival analysis
Rscript --vanilla code/13-survival-analysis.R

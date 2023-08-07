# create data availability plots
bash analyses/data-availability/run_data_plots.sh

# create merged files
bash analyses/merge-files/run_merge_files.sh

# create oncoplots
bash analyses/oncoplots/run_oncoplots.sh

# survival analysis
Rscript --vanilla analyses/survival-analysis/01-survival-analysis.R

# ALT correlations
Rscript --vanilla analyses/alt-analysis/01-correlation_alt_vs_vars.R

# MSIsensor-pro analysis
bash analyses/msi-sensor-analysis/run_msisensor_analysis.sh

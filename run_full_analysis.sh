# create merged files
cd analyses/merge-files && bash run_merge_files.sh
cd ../..

# create data availability plots
cd analyses/data-availability && bash run_data_plots.sh
cd ../..

# create oncoplots
cd analyses/oncoplots && bash run_oncoplots.sh
cd ../..

# survival analysis
Rscript --vanilla analyses/survival-analysis/01-survival-analysis.R

# ALT correlations
Rscript --vanilla analyses/alt-analysis/01-correlation_alt_vs_vars.R

# MSIsensor-pro analysis
cd analyses/msi-sensor-analysis && bash run_msisensor_analysis.sh
cd ../..

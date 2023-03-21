# merge MSI sensor files into a single file
Rscript code/02-merge_msi_paired.R
Rscript code/02-merge_msi_tumor_only.R

# compare MSI paired vs tumor-only analysis
Rscript code/03-msi_paired_vs_tumor_only.R

# merged output file for all downstream analyses
Rscript code/04-msi_paired_merged_output.R
Rscript code/04-msi_tumor_only_merged_output.R

# summary files with correlation of MSI with other variables
Rscript code/05-msi_paired_summary.R
Rscript code/05-msi_tumor_only_summary.R

# MSI sensor vs MMR pathways
Rscript code/09-msi_vs_mmr_pathways.R

# compare ALT status with other variables
Rscript code/10-correlation_alt_vs_vars.R

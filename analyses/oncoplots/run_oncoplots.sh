# tumor vs normal differential expression (single-sample)
Rscript --vanilla 01-deg-vs-gtex-brain.R

# prepare files for oncoplots                  
Rscript --vanilla 02-prepare_files_oncogrid.R

# prepare files for oncoplots with tumor-only data
Rscript --vanilla 03-prepare_files_oncogrid_add_tumor_only.R

# plot oncoplot with all modalities
Rscript --vanilla 04-plot_oncogrid.R

# plot oncoplot without RNA-only data                     
Rscript --vanilla 05-plot_oncogrid_no_RNA.R

# plot cascade plot with all modalities
Rscript --vanilla 06-cascade_plots.R

# plot cascade plot with all modalities with tumor-only data
Rscript --vanilla 07-cascade_plots_add_tumor_only.R

# correlate ALT and MSI to major SNV > 6% from oncoplot_orderby_sex_age_H3F3A_status_norna.pdf
Rscript --vanilla 08-major_snv_analysis.R

# correlation of gene alterations with clinical variables
Rscript --vanilla 09-gene_alteration_correlation.R

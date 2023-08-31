# merge files (RNA, CNV, SNV, MSI)
Rscript --vanilla 01-create_merged_files.R

# merge tumor-only files (RNA, CNV, SNV, MSI)
Rscript --vanilla 02-create_merged_files_tumor_only.R

# merge, filter and annotate fusions
Rscript --vanilla 03-filter-annotate-fusions.R

# Function: Script to generate Oncogrid plot (add tumor-only information)

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(tidyverse)
  library(circlize)
})

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analyses_dir <- file.path(root_dir, "analyses", "oncoplots")
input_dir <- file.path(analyses_dir, "results")
output_dir <- file.path(analyses_dir, "results", "cascade_plots_add_tumor_only")
dir.create(output_dir, recursive = T, showWarnings = F)

# matrix
mat = read.table(file.path(input_dir, "oncoprint_add_tumor_only.txt"),  header = TRUE, stringsAsFactors=FALSE, sep = "\t",check.names = FALSE)
mat[is.na(mat)] = ""
rownames(mat) = mat[, 1]
mat = mat[, -1]
mat = t(as.matrix(mat))

# remove RNA related annotations
mat[mat %in% c("OVE", "UNE", "FUS", "GAI", "LOS")] <- ""
mat = gsub(";UNE|;OVE|;FUS", "", mat)
mat = gsub("GAI|LOS|;LOS|GAI;LOS|LOS;GAI", "", mat)
# mat[mat == ""] <- "NO"
mat[mat != ""] <- "Mutation"

keep = apply(mat, 1, FUN = function(x) length(unique(x)) == 2)
mat = mat[keep,]

# get unique values
apply(mat, 1, unique) %>%  unlist %>% unique()

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w, h, gp = gpar(fill = "#ffffff",col= "#595959"))
  },
  Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = "black", col = NA))
  }
)
col = c("Mutation" = "black")

# read annotation and TMB info
annot_info <- read.delim(file.path(input_dir, "annotation_add_tumor_only.txt"), header = TRUE, check.names = TRUE)
annot_info <- annot_info %>%
  filter(Sample %in% colnames(mat),
         Sequencing_Experiment != "RNA-Seq") %>%
  remove_rownames() %>%
  column_to_rownames('Sample') %>%
  as.data.frame()
samples_to_use <- intersect(rownames(annot_info), colnames(mat))
mat <- mat[,samples_to_use]
annot_info <- annot_info[samples_to_use,]

# annotation 1
col_fun_tmb = colorRamp2(c(0, max(annot_info$TMB, na.rm = T)), c("white", "magenta3"))
# col_fun_msi = colorRamp2(c(0, max(annot_info$MSI, na.rm = T)), c("white", "purple3"))
annot_info$Integrated_Diagnosis[is.na(annot_info$Integrated_Diagnosis)] <- "N/A"
annot_info$Diagnosis_Type[is.na(annot_info$Diagnosis_Type)] <- "N/A"
annot_info$Tumor_Location[is.na(annot_info$Tumor_Location)] <- "N/A"
annot_info$Sex[is.na(annot_info$Sex)] <- "N/A"
annot_info$Age[is.na(annot_info$Age)] <- "N/A"
annot_info$Age <- factor(annot_info$Age, levels = c("[0,15]", "(15,40]", "N/A"))
# annot_info$Age_Three_Groups[is.na(annot_info$Age_Three_Groups)] <- "N/A"
# annot_info$Age_Three_Groups <- factor(annot_info$Age_Three_Groups, levels = c("[0,15]", "(15,26]", "(26,40]", "N/A"))
ha = HeatmapAnnotation(df = annot_info %>% dplyr::select(-c(Sequencing_Experiment)), col = list(
  TMB = col_fun_tmb,
  Sequencing_Experiment = c("WGS" = "red",
                            "RNA-Seq, WGS" = "orangered",
                            "RNA-Seq, WGS_Tumor_Only" = "orange"),
  Integrated_Diagnosis = c("High-grade glioma/astrocytoma (WHO grade III/IV)"="lightseagreen",
                           "Astrocytoma;Oligoastrocytoma" = "mediumorchid2",
                           "Astrocytoma" = "brown2", 
                           "Glioblastoma" = "orange",
                           "Pleomorphic xanthoastrocytoma" = "darkgreen",
                           "Diffuse Midline Glioma" = "blue2"),
  Diagnosis_Type = c("Initial CNS Tumor" = "#cee397",
                     "Progressive" = "#827397",
                     "Recurrence" = "#363062",
                     "Second Malignancy" = "#005082"),
  Tumor_Location = c("Cortical" = "magenta",
                     "Other/Multiple locations/NOS" = "pink",
                     "Midline" = "purple",
                     "Cerebellar" = "navy"),
  Sex = c("Female" = "deeppink4",
          "Male" = "navy"),
  Age = c("[0,15]" = "gold",
                     "(15,40]" = "purple",
                     "N/A" = "gray")),
  # Age_Three_Groups = c("[0,15]" = "gold",
  #                      "(15,26]" = "purple",
  #                      "(26,40]" = "darkgreen",
  #                      "N/A" = "gray")),
  annotation_name_gp = gpar(fontsize = 9),
  gp = gpar(col = "#595959"), 
  simple_anno_size = unit(4, "mm"), 
  annotation_name_side = "left",
  annotation_legend_param = list(
    # Sequencing_Experiment = list(nrow = 3)
  ))

# oncoprint
ht = oncoPrint(mat, get_type = function(x)strsplit(x, ";")[[1]],
               alter_fun = alter_fun, 
               col = col, 
               show_column_names = TRUE, 
               column_names_gp = gpar(fontsize = 9),
               row_names_gp = gpar(fontsize = 9),
               column_names_side = "top",
               top_annotation = ha,
               right_annotation = NULL,
               row_names_side = "left",
               pct_side = "right",
               remove_empty_rows = TRUE,
               heatmap_legend_param = list(title = "Alteration", nrow = 9, title_position = "topleft", direction = "horizontal",
                                           at = c("Mutation"),
                                           labels = c("SNV")
               ))

pdf(file = file.path(output_dir, "cascade_plot.pdf"), width = 15, height = 10) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# order by Age (2 groups)
ht = oncoPrint(mat, get_type = function(x)strsplit(x, ";")[[1]],
               alter_fun = alter_fun, 
               column_split = annot_info$Age,
               col = col, 
               show_column_names = TRUE, 
               column_names_gp = gpar(fontsize = 9),
               row_names_gp = gpar(fontsize = 9),
               column_names_side = "top",
               top_annotation = ha,
               right_annotation = NULL,
               row_names_side = "left",
               pct_side = "right",
               remove_empty_rows = TRUE,
               heatmap_legend_param = list(title = "Alteration", nrow = 9, title_position = "topleft", direction = "horizontal",
                                           at = c("Mutation"),
                                           labels = c("SNV")
               ))
pdf(file = file.path(output_dir, "cascade_orderby_age.pdf"), width = 15, height = 10) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# # order by Age (3 groups)
# ht = oncoPrint(mat, get_type = function(x)strsplit(x, ";")[[1]],
#                alter_fun = alter_fun, 
#                column_split = annot_info$Age_Three_Groups,
#                col = col, 
#                show_column_names = TRUE, 
#                column_names_gp = gpar(fontsize = 9),
#                row_names_gp = gpar(fontsize = 9),
#                column_names_side = "top",
#                top_annotation = ha,
#                right_annotation = NULL,
#                row_names_side = "left",
#                pct_side = "right",
#                remove_empty_rows = TRUE,
#                heatmap_legend_param = list(title = "Alteration", nrow = 9, title_position = "topleft", direction = "horizontal",
#                                            at = c("Mutation"),
#                                            labels = c("SNV")
#                ))
# pdf(file = file.path(output_dir, "cascade_orderby_age_three_groups.pdf"), width = 15, height = 10) 
# draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
# dev.off()

# order by Sex
ht = oncoPrint(mat, get_type = function(x)strsplit(x, ";")[[1]],
               alter_fun = alter_fun, 
               column_split = annot_info$Sex,
               col = col, 
               show_column_names = TRUE, 
               column_names_gp = gpar(fontsize = 9),
               row_names_gp = gpar(fontsize = 9),
               column_names_side = "top",
               top_annotation = ha,
               right_annotation = NULL,
               row_names_side = "left",
               pct_side = "right",
               remove_empty_rows = TRUE,
               heatmap_legend_param = list(title = "Alteration", nrow = 9, title_position = "topleft", direction = "horizontal",
                                           at = c("Mutation"),
                                           labels = c("SNV")
               ))
pdf(file = file.path(output_dir, "cascade_orderby_sex.pdf"), width = 15, height = 10) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

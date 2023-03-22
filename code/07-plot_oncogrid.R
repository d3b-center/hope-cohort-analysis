# Function: Script to generate Oncogrid plot

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(tidyverse)
  library(circlize)
})

# output directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
output_dir <- file.path(root_dir, "results", "oncoplots_three_groups")
dir.create(output_dir, recursive = T, showWarnings = F)

# matrix
mat = read.table(file.path("results", "oncoprint.txt"),  header = TRUE, stringsAsFactors=FALSE, sep = "\t",check.names = FALSE)
mat[is.na(mat)] = ""
rownames(mat) = mat[, 1]
mat = mat[, -1]
mat = t(as.matrix(mat))

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w, h, gp = gpar(fill = "#ffffff",col= "#595959"))
  },
  GAI = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = "#ff4d4d", col = NA))
  },
  LOS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = "#0D47A1", col = NA))
  },
  MIS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = "#77b300", col = NA))
  },
  NOS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = "#80bfff", col = NA))
  },
  FSD = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = "#1a53ff", col = NA))
  },
  FSI = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = "#8D6E63", col = NA))
  },
  NOT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = "#9966ff", col = NA))
  },
  SPS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = "#E69F00", col = NA))
  },
  IFD = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = "#827717", col = NA))
  },
  OVE = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h*0.66, gp = gpar(fill = "#dbc6eb", col = NA))
  },
  UNE = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h*0.66, gp = gpar(fill = "#709fb0", col = NA))
  },
  FUS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h*0.66, gp = gpar(fill = "#AB47BC", col = NA))
  }
)
col = c("GAI" = "#ff4d4d", "LOS" = "#0D47A1", 
        "FUS" = "#AB47BC",
        "MIS" = "#77b300", "NOS" ="#80bfff", "FSD" = "#1a53ff", "FSI" ="#8D6E63", "NOT" = "#9966ff","SPS" = "#E69F00","IFD" = "#827717","OVE" = "#dbc6eb","UNE" = "#709fb0")

# read annotation and TMB info
annot_info <- read.delim(file.path("results", "annotation.txt"), header = TRUE, check.names = TRUE)
annot_info <- annot_info %>%
  dplyr::select(-c(age_two_groups)) %>%
  dplyr::rename("Age" = "age_three_groups") %>%
  filter(Sample %in% colnames(mat)) %>%
  remove_rownames() %>%
  column_to_rownames('Sample') %>%
  as.data.frame()
samples_to_use <- intersect(rownames(annot_info), colnames(mat))
mat <- mat[,samples_to_use]
annot_info <- annot_info[samples_to_use,]

# annotation 1
col_fun_tmb = colorRamp2(c(0, max(annot_info$TMB, na.rm = T)), c("white", "magenta3"))
col_fun_msi = colorRamp2(c(0, max(annot_info$MSI_Percent, na.rm = T)), c("white", "purple3"))
annot_info$Tumor_Descriptor[is.na(annot_info$Tumor_Descriptor)] <- "N/A"
annot_info$Integrated_Diagnosis[is.na(annot_info$Integrated_Diagnosis)] <- "N/A"
annot_info$Sex[is.na(annot_info$Sex)] <- "N/A"
annot_info$Age[is.na(annot_info$Age)] <- "N/A"
annot_info$Age <- factor(annot_info$Age, levels = c("[0,15]", "(15,26]", "(26,40]"))
ha = HeatmapAnnotation(df = annot_info , col = list(
  TMB = col_fun_tmb,
  MSI_Percent = col_fun_msi,
  Sequencing_Experiment = c("WGS" = "red",
                            "RNA-Seq" = "yellow",
                            "RNA-Seq, WGS" = "orange",
                            "N/A" = "gray"),
  Tumor_Descriptor = c("Initial CNS Tumor" = "#cee397",
                       "Progressive" = "#827397",
                       "Recurrence" = "#363062",
                       "Second Malignancy" = "#005082",
                       "N/A" = "gray"),
  Integrated_Diagnosis = c("High-grade glioma/astrocytoma (WHO grade III/IV)"="lightseagreen",
                           "Astrocytoma;Oligoastrocytoma" = "mediumorchid2",
                           "Astrocytoma" = "brown2", 
                           "Glioblastoma" = "orange",
                           "Pleomorphic xanthoastrocytoma" = "darkgreen",
                           "Diffuse Midline Glioma" = "blue2",
                           "N/A" = "gray"),
  Sex = c("Female" = "deeppink4",
          "Male" = "navy",
          "N/A" = "gray"),
  Age = c("[0,15]" = "gold",
          "(15,26]" = "purple",
          "(26,40]" = "darkgreen",
          "N/A" = "gray")),
  annotation_name_gp = gpar(fontsize = 9),
  gp = gpar(col = "#595959"), simple_anno_size = unit(4, "mm"), annotation_name_side = "left",
  annotation_legend_param = list(
    Sequencing_Experiment = list(nrow = 3)
  ))

# annotation 2
amp = ifelse(apply(mat, 1, function(x) sum(grepl("GAI", x) + grepl("LOS",x) + grepl("OVE",x) + grepl("UNE",x))/length(x) > 0), "Copy number alteration and RNAseq Expression", "Genetic and Fusion alteration")
amp = factor(amp, levels = c("Genetic and Fusion alteration", "Copy number alteration and RNAseq Expression"))

# oncoprint
ht = oncoPrint(mat, get_type = function(x)strsplit(x, ";")[[1]],
               alter_fun = alter_fun, col = col, show_column_names = TRUE, 
               column_names_gp = gpar(fontsize = 9),
               row_names_gp = gpar(fontsize = 9),
               column_names_side = "top",
               top_annotation = ha,
               right_annotation = NULL,
               row_names_side = "left",
               pct_side = "right",
               row_split = amp,
               remove_empty_rows = TRUE,
               heatmap_legend_param = list(title = "Alterations", nrow = 9, title_position = "topleft", direction = "horizontal",
                                           at = c("GAI", "LOS", "FUS", "MIS", "NOS", "FSD", "FSI", "SPS", "IFD", "OVE", "UNE"),
                                           labels = c("Copy gain", "Copy loss", "Gene Fusion", "Missense", "Nonsense", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice site", "In_Frame_Del", "Over Expression", "Under Expression")
               ))

pdf(file = file.path(output_dir, "oncoplot.pdf"), width = 15, height = 27) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# column order by sex
sex_ordered <- annot_info %>% 
  dplyr::arrange(Sex)
ht = oncoPrint(mat, get_type = function(x)strsplit(x, ";")[[1]],
               alter_fun = alter_fun, col = col, show_column_names = TRUE, 
               column_names_gp = gpar(fontsize = 9),
               row_names_gp = gpar(fontsize = 9),
               column_order = match(rownames(sex_ordered), colnames(mat)),
               column_names_side = "top",
               top_annotation = ha,
               right_annotation = NULL,
               row_names_side = "left",
               pct_side = "right",
               row_split = amp,
               remove_empty_rows = TRUE,
               heatmap_legend_param = list(title = "Alterations", nrow = 9, title_position = "topleft", direction = "horizontal",
                                           at = c("GAI", "LOS", "FUS", "MIS", "NOS", "FSD", "FSI", "SPS", "IFD", "OVE", "UNE"),
                                           labels = c("Copy gain", "Copy loss", "Gene Fusion", "Missense", "Nonsense", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice site", "In_Frame_Del", "Over Expression", "Under Expression")
               ))
pdf(file = file.path(output_dir, "oncoplot_orderby_sex.pdf"), width = 15, height = 27) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# column order by H3F3A status
n <- grep("H3-3A$", rownames(mat))
ht = oncoPrint(mat, get_type = function(x)strsplit(x, ";")[[1]],
               alter_fun = alter_fun, col = col, show_column_names = TRUE, 
               column_names_gp = gpar(fontsize = 9),
               row_names_gp = gpar(fontsize = 9),
               column_order = order(mat[n,], decreasing = T),
               column_names_side = "top",
               top_annotation = ha,
               right_annotation = NULL,
               row_names_side = "left",
               pct_side = "right",
               row_split = amp,
               remove_empty_rows = TRUE,
               heatmap_legend_param = list(title = "Alterations", nrow = 9, title_position = "topleft", direction = "horizontal",
                                           at = c("GAI", "LOS", "FUS", "MIS", "NOS", "FSD", "FSI", "SPS", "IFD", "OVE", "UNE"),
                                           labels = c("Copy gain", "Copy loss", "Gene Fusion", "Missense", "Nonsense", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice site", "In_Frame_Del", "Over Expression", "Under Expression")
               ))
pdf(file = file.path(output_dir, "oncoplot_orderby_H3F3A.pdf"), width = 15, height = 27) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()


# column order by sex + H3F3A
tmp <- reshape2::melt(mat) %>%
  filter(Var1 == "H3-3A")
tmp <- annot_info %>% 
  rownames_to_column("sample_id") %>%
  inner_join(tmp, by = c("sample_id" = "Var2")) %>%
  arrange(Sex, desc(value)) %>%
  column_to_rownames("sample_id")
ht = oncoPrint(mat, get_type = function(x)strsplit(x, ";")[[1]],
               alter_fun = alter_fun, col = col, show_column_names = TRUE, 
               column_names_gp = gpar(fontsize = 9),
               row_names_gp = gpar(fontsize = 9),
               column_order = match(rownames(tmp), colnames(mat)),
               column_names_side = "top",
               top_annotation = ha,
               right_annotation = NULL,
               row_names_side = "left",
               pct_side = "right",
               row_split = amp,
               remove_empty_rows = TRUE,
               heatmap_legend_param = list(title = "Alterations", nrow = 9, title_position = "topleft", direction = "horizontal",
                                           at = c("GAI", "LOS", "FUS", "MIS", "NOS", "FSD", "FSI", "SPS", "IFD", "OVE", "UNE"),
                                           labels = c("Copy gain", "Copy loss", "Gene Fusion", "Missense", "Nonsense", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice site", "In_Frame_Del", "Over Expression", "Under Expression")
               ))
pdf(file = file.path(output_dir, "oncoplot_orderby_sex_H3F3A_status.pdf"), width = 15, height = 27) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# column order by sex + age
sex_age_ordered <- annot_info %>% 
  dplyr::arrange(Sex, Age)
ht = oncoPrint(mat, get_type = function(x)strsplit(x, ";")[[1]],
               alter_fun = alter_fun, col = col, show_column_names = TRUE, 
               column_names_gp = gpar(fontsize = 9),
               row_names_gp = gpar(fontsize = 9),
               column_order = match(rownames(sex_age_ordered), colnames(mat)),
               column_names_side = "top",
               top_annotation = ha,
               right_annotation = NULL,
               row_names_side = "left",
               pct_side = "right",
               row_split = amp,
               remove_empty_rows = TRUE,
               heatmap_legend_param = list(title = "Alterations", nrow = 9, title_position = "topleft", direction = "horizontal",
                                           at = c("GAI", "LOS", "FUS", "MIS", "NOS", "FSD", "FSI", "SPS", "IFD", "OVE", "UNE"),
                                           labels = c("Copy gain", "Copy loss", "Gene Fusion", "Missense", "Nonsense", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice site", "In_Frame_Del", "Over Expression", "Under Expression")
               ))
pdf(file = file.path(output_dir, "oncoplot_orderby_sex_age.pdf"), width = 15, height = 27) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# column order by sex + age + histone
tmp <- reshape2::melt(mat) %>%
  filter(Var1 == "H3-3A")
tmp <- annot_info %>% 
  rownames_to_column("sample_id") %>%
  inner_join(tmp, by = c("sample_id" = "Var2")) %>%
  arrange(Sex, Age, desc(value)) %>%
  column_to_rownames("sample_id")
ht = oncoPrint(mat, get_type = function(x)strsplit(x, ";")[[1]],
               alter_fun = alter_fun, col = col, show_column_names = TRUE, 
               column_names_gp = gpar(fontsize = 9),
               row_names_gp = gpar(fontsize = 9),
               column_order = match(rownames(tmp), colnames(mat)),
               column_names_side = "top",
               top_annotation = ha,
               right_annotation = NULL,
               row_names_side = "left",
               pct_side = "right",
               row_split = amp,
               remove_empty_rows = TRUE,
               heatmap_legend_param = list(title = "Alterations", nrow = 9, title_position = "topleft", direction = "horizontal",
                                           at = c("GAI", "LOS", "FUS", "MIS", "NOS", "FSD", "FSI", "SPS", "IFD", "OVE", "UNE"),
                                           labels = c("Copy gain", "Copy loss", "Gene Fusion", "Missense", "Nonsense", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice site", "In_Frame_Del", "Over Expression", "Under Expression")
               ))
pdf(file = file.path(output_dir, "oncoplot_orderby_sex_age_H3F3A_status.pdf"), width = 15, height = 27) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

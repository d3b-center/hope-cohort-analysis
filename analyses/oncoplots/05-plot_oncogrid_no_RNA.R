# Function: Script to generate Oncogrid plot

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(tidyverse)
  library(circlize)
})

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analyses_dir <- file.path(root_dir, "analyses", "oncoplots")
input_dir <- file.path(analyses_dir, "results")
output_dir <- file.path(analyses_dir, "results", "oncoplots")
dir.create(output_dir, recursive = T, showWarnings = F)

# matrix
mat = read.table(file.path(input_dir, "oncoprint.txt"),  header = TRUE, stringsAsFactors=FALSE, sep = "\t",check.names = FALSE)
mat[is.na(mat)] = ""
rownames(mat) = mat[, 1]
mat = mat[, -1]
mat = t(as.matrix(mat))

# remove RNA related annotations
mat[mat == "OVE"] <- ""
mat[mat == "UNE"] <- ""
mat = gsub(";UNE", "", mat)
mat = gsub(";OVE", "", mat)

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
  FUS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h*0.66, gp = gpar(fill = "#AB47BC", col = NA))
  }
)
col = c("GAI" = "#ff4d4d", "LOS" = "#0D47A1", 
        "FUS" = "#AB47BC",
        "MIS" = "#77b300", "NOS" ="#80bfff", "FSD" = "#1a53ff", "FSI" ="#8D6E63", "NOT" = "#9966ff","SPS" = "#E69F00","IFD" = "#827717")

# read annotation and TMB info
annot_info <- read.delim(file.path(input_dir, "annotation.txt"), header = TRUE, check.names = TRUE)
annot_info <- annot_info %>%
  filter(Sample %in% colnames(mat),
         Sequencing_Experiment != "RNA-Seq") %>%
  remove_rownames() %>%
  column_to_rownames('Sample') %>%
  as.data.frame()
samples_to_use <- intersect(rownames(annot_info), colnames(mat))
mat <- mat[,samples_to_use]
annot_info <- annot_info[samples_to_use,]

# only plot top 20 genes 
snv_fus_genes_to_keep = apply(mat, 1, FUN = function(x) 
  (length(grep("MIS|NOS|FSD|FSI|NOT|SPS|IFD|FUS", x))/ncol(mat))*100)
snv_fus_genes_to_keep <- names(sort(snv_fus_genes_to_keep, decreasing = TRUE)[1:20])
cnv_dge_genes_to_keep = apply(mat, 1, FUN = function(x) 
  (length(grep("GAI|LOS|OVE|UNE", x))/ncol(mat))*100)
cnv_dge_genes_to_keep <- names(sort(cnv_dge_genes_to_keep, decreasing = TRUE)[1:20])
mat <- mat[which(rownames(mat) %in% c(snv_fus_genes_to_keep, cnv_dge_genes_to_keep)),]

# annotation 1
col_fun_tmb = colorRamp2(c(0, max(annot_info$TMB, na.rm = T)), c("white", "magenta3"))
annot_info$Age <- factor(annot_info$Age, levels = c("[0,15]", "(15,26]", "(26,40]"))
annot_info$Molecular_Subtype <- as.character(annot_info$Molecular_Subtype)
annot_info$Cancer_Group <- as.character(annot_info$Cancer_Group)

ha = HeatmapAnnotation(df = annot_info %>% dplyr::select(-c(Sequencing_Experiment)), 
                       col = list(TMB = col_fun_tmb,
                                  Diagnosis = c("High-grade glioma/astrocytoma (WHO grade III/IV)" = "lightseagreen",
                                                "Diffuse Midline Glioma (WHO grade III/IV)" = "darkgreen",
                                                "Astrocytoma;Oligoastrocytoma (WHO grade III)" = "mediumorchid2",
                                                "Astrocytoma (WHO grade III/IV)" = "#5fff57", 
                                                "Glioblastoma (WHO grade IV)" = "#f268d6",
                                                "Pleomorphic xanthoastrocytoma (WHO grade II/III)" = "#005082"),
                                  Molecular_Subtype = c("DMG, H3 K28" = "#053061",
                                                        "DHG, H3 G35, TP53" = "#A6761D",
                                                        "HGG, H3 wildtype" = "#4393c3",
                                                        "HGG, H3 wildtype, TP53" = "darkgreen",
                                                        "DMG, H3 K28, TP53" = "#BC80BD",
                                                        "HGG, IDH, TP53" = "#FFFF99",
                                                        "IHG, NTRK-altered, TP53"  = "#E7298A",
                                                        "IHG, NTRK-altered" = "#f4a582",
                                                        "IHG, ROS1-altered" = "#d6604d",
                                                        "IHG, ALK-altered" = "#E31A1C",
                                                        "PXA" = "#67001f",
                                                        "HGG, IDH" = "#B3DE69",
                                                        "NA" = "#f1f1f1"),
                                  Diagnosis_Type = c("Initial CNS Tumor" = "#cee397",
                                                     "Progressive" = "#827397",
                                                     "Recurrence" = "#363062",
                                                     "Second Malignancy" = "#005082"),
                                  Tumor_Location = c("Cortical" = "#D4806C",
                                                     "Other/Multiple locations/NOS" = "#7C8F97",
                                                     "Midline" = "#344C68",
                                                     "Cerebellar" = "#94004C"),
                                  CNS_region = c("Posterior fossa" = "#D4806C",
                                                 "Other" = "#7C8F97",
                                                 "Midline" = "#344C68",
                                                 "Hemispheric" = "#94004C",
                                                 "Mixed" = "darkgreen"),
                                  Cancer_Group = c("DMG" = "#053061",
                                                   "DHG" = "#A6761D",
                                                   "HGG" = "#4393c3",
                                                   "IHG" = "#E7298A",
                                                   "NA" = "#f1f1f1"),
                                  Sex = c("Male" = "#0707CF",
                                          "Female" = "#CC0303"),
                                  Age = c("[0,15]" = "#C7E9C0",
                                          "(15,26]" = "#74C476",
                                          "(26,40]" = "#238B45")),
                       annotation_name_gp = gpar(fontsize = 9),
                       gp = gpar(col = "#595959"), simple_anno_size = unit(4, "mm"), annotation_name_side = "left")

# annotation 2
amp = ifelse(apply(mat, 1, function(x) sum(grepl("GAI", x) + grepl("LOS",x))/length(x) > 0), "Copy number alteration", "Genetic and Fusion alteration")
amp = factor(amp, levels = c("Genetic and Fusion alteration", "Copy number alteration"))

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
                                           at = c("GAI", "LOS", "FUS", "MIS", "NOS", "FSD", "FSI", "SPS", "IFD"),
                                           labels = c("Copy gain", "Copy loss", "Gene Fusion", "Misense", "Nonsense", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice site", "In_Frame_Del")
               ))


pdf(file = file.path(output_dir, "oncoplot_norna.pdf"), width = 16, height = 8) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# column order by sex
sex_ordered <- annot_info %>% 
  arrange(Sex)
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
                                           at = c("GAI", "LOS", "FUS", "MIS", "NOS", "FSD", "FSI", "SPS", "IFD"),
                                           labels = c("Copy gain", "Copy loss", "Gene Fusion", "Misense", "Nonsense", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice site", "In_Frame_Del")
               ))
pdf(file = file.path(output_dir, "oncoplot_orderby_sex_norna.pdf"), width = 16, height = 8) 
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
                                           at = c("GAI", "LOS", "FUS", "MIS", "NOS", "FSD", "FSI", "SPS", "IFD"),
                                           labels = c("Copy gain", "Copy loss", "Gene Fusion", "Misense", "Nonsense", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice site", "In_Frame_Del")
               ))
pdf(file = file.path(output_dir, "oncoplot_orderby_H3F3A_norna.pdf"), width = 16, height = 8) 
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
                                           at = c("GAI", "LOS", "FUS", "MIS", "NOS", "FSD", "FSI", "SPS", "IFD"),
                                           labels = c("Copy gain", "Copy loss", "Gene Fusion", "Misense", "Nonsense", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice site", "In_Frame_Del")
               ))
pdf(file = file.path(output_dir, "oncoplot_orderby_sex_H3F3A_status_norna.pdf"), width = 16, height = 8) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# column order by sex + age
sex_age_ordered <- annot_info %>% 
  arrange(Sex, Age)
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
                                           at = c("GAI", "LOS", "FUS", "MIS", "NOS", "FSD", "FSI", "SPS", "IFD"),
                                           labels = c("Copy gain", "Copy loss", "Gene Fusion", "Misense", "Nonsense", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice site", "In_Frame_Del")
               ))
pdf(file = file.path(output_dir, "oncoplot_orderby_sex_age_norna.pdf"), width = 16, height = 8) 
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
                                           at = c("GAI", "LOS", "FUS", "MIS", "NOS", "FSD", "FSI", "SPS", "IFD"),
                                           labels = c("Copy gain", "Copy loss", "Gene Fusion", "Misense", "Nonsense", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice site", "In_Frame_Del")
               ))
pdf(file = file.path(output_dir, "oncoplot_orderby_sex_age_H3F3A_status_norna.pdf"), width = 16, height = 8) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# annotation 1
ha = HeatmapAnnotation(df = annot_info %>% dplyr::select(-c(Sequencing_Experiment)), 
                       col = list(TMB = col_fun_tmb,
                                  Diagnosis = c("High-grade glioma/astrocytoma (WHO grade III/IV)" = "lightseagreen",
                                                "Diffuse Midline Glioma (WHO grade III/IV)" = "darkgreen",
                                                "Astrocytoma;Oligoastrocytoma (WHO grade III)" = "mediumorchid2",
                                                "Astrocytoma (WHO grade III/IV)" = "#5fff57", 
                                                "Glioblastoma (WHO grade IV)" = "#f268d6",
                                                "Pleomorphic xanthoastrocytoma (WHO grade II/III)" = "#005082"),
                                  Molecular_Subtype = c("DMG, H3 K28" = "#053061",
                                                        "DHG, H3 G35, TP53" = "#A6761D",
                                                        "HGG, H3 wildtype" = "#4393c3",
                                                        "HGG, H3 wildtype, TP53" = "darkgreen",
                                                        "DMG, H3 K28, TP53" = "#BC80BD",
                                                        "HGG, IDH, TP53" = "#FFFF99",
                                                        "IHG, NTRK-altered, TP53"  = "#E7298A",
                                                        "IHG, NTRK-altered" = "#f4a582",
                                                        "IHG, ROS1-altered" = "#d6604d",
                                                        "IHG, ALK-altered" = "#E31A1C",
                                                        "PXA" = "#67001f",
                                                        "HGG, IDH" = "#B3DE69",
                                                        "NA" = "#f1f1f1"),
                                  Diagnosis_Type = c("Initial CNS Tumor" = "#cee397",
                                                     "Progressive" = "#827397",
                                                     "Recurrence" = "#363062",
                                                     "Second Malignancy" = "#005082"),
                                  Tumor_Location = c("Cortical" = "#D4806C",
                                                     "Other/Multiple locations/NOS" = "#7C8F97",
                                                     "Midline" = "#344C68",
                                                     "Cerebellar" = "#94004C"),
                                  CNS_region = c("Posterior fossa" = "#D4806C",
                                                 "Other" = "#7C8F97",
                                                 "Midline" = "#344C68",
                                                 "Hemispheric" = "#94004C",
                                                 "Mixed" = "darkgreen"),
                                  Cancer_Group = c("DMG" = "#053061",
                                                   "DHG" = "#A6761D",
                                                   "HGG" = "#4393c3",
                                                   "IHG" = "#E7298A",
                                                   "NA" = "#f1f1f1"),
                                  Sex = c("Male" = "#0707CF",
                                          "Female" = "#CC0303"),
                                  Age = c("[0,15]" = "#C7E9C0",
                                          "(15,26]" = "#74C476",
                                          "(26,40]" = "#238B45")),
                       annotation_name_gp = gpar(fontsize = 9),
                       gp = gpar(col = "#595959"), simple_anno_size = unit(4, "mm"), annotation_name_side = "left")

# annotation 2
amp = ifelse(apply(mat, 1, function(x) sum(grepl("GAI", x) + grepl("LOS",x))/length(x) > 0), "Copy number alteration", "Genetic and Fusion alteration")
amp = factor(amp, levels = c("Genetic and Fusion alteration", "Copy number alteration"))

# column order age
age_ordered <- annot_info %>% 
  arrange(Age)
ht = oncoPrint(mat, get_type = function(x)strsplit(x, ";")[[1]],
               alter_fun = alter_fun, col = col, show_column_names = TRUE, 
               column_names_gp = gpar(fontsize = 9),
               row_names_gp = gpar(fontsize = 9),
               column_order = match(rownames(age_ordered), colnames(mat)),
               column_names_side = "top",
               top_annotation = ha,
               right_annotation = NULL,
               row_names_side = "left",
               pct_side = "right",
               row_split = amp,
               remove_empty_rows = TRUE,
               heatmap_legend_param = list(title = "Alterations", nrow = 9, title_position = "topleft", direction = "horizontal",
                                           at = c("GAI", "LOS", "FUS", "MIS", "NOS", "FSD", "FSI", "SPS", "IFD"),
                                           labels = c("Copy gain", "Copy loss", "Gene Fusion", "Misense", "Nonsense", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice site", "In_Frame_Del")
               ))
pdf(file = file.path(output_dir, "oncoplot_orderby_age_norna.pdf"), width = 16, height = 8) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# column order by age + histone
tmp <- reshape2::melt(mat) %>%
  filter(Var1 == "H3-3A")
tmp <- annot_info %>% 
  rownames_to_column("sample_id") %>%
  inner_join(tmp, by = c("sample_id" = "Var2")) %>%
  arrange(Age, desc(value)) %>%
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
                                           at = c("GAI", "LOS", "FUS", "MIS", "NOS", "FSD", "FSI", "SPS", "IFD"),
                                           labels = c("Copy gain", "Copy loss", "Gene Fusion", "Misense", "Nonsense", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice site", "In_Frame_Del")
               ))
pdf(file = file.path(output_dir, "oncoplot_orderby_age_H3F3A_status_norna.pdf"), width = 16, height = 8) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# column order by age + sex + histone
tmp <- reshape2::melt(mat) %>%
  filter(Var1 == "H3-3A")
tmp <- annot_info %>% 
  rownames_to_column("sample_id") %>%
  inner_join(tmp, by = c("sample_id" = "Var2")) %>%
  arrange(Age, Sex, desc(value)) %>%
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
                                           at = c("GAI", "LOS", "FUS", "MIS", "NOS", "FSD", "FSI", "SPS", "IFD"),
                                           labels = c("Copy gain", "Copy loss", "Gene Fusion", "Misense", "Nonsense", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice site", "In_Frame_Del")
               ))
pdf(file = file.path(output_dir, "oncoplot_orderby_age_sex_H3F3A_status_norna.pdf"), width = 16, height = 8) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

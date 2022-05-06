# Function: Script to generate Oncogrid plot

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(tidyverse)
})

# matrix
mat = read.table(file.path("results", "oncoprint.txt"),  header = TRUE, stringsAsFactors=FALSE, sep = "\t",check.names = FALSE)

# subset to 94 samples 
# hope_cohort_subset <- read.delim('data/hope_cohort_subset.tsv', header = F)
hope_cohort_subset <- readxl::read_xlsx('data/hope_clinical_table_033022.xlsx', sheet = 1)
mat <- mat %>%
  filter(Sample %in% hope_cohort_subset$Sample_ID)

mat[is.na(mat)] = ""
rownames(mat) = mat[, 1]
mat = mat[, -1]
mat = t(as.matrix(mat))

# code
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
col = c("GAI" = "#ff4d4d", "LOS" = "#0D47A1" , "MIS" = "#77b300", "FUS" = "#AB47BC","NOS" ="#80bfff", "FSD" = "#1a53ff", "FSI" ="#8D6E63", "NOT" = "#9966ff","SPS" = "#E69F00","IFD" = "#827717")

# read annotation and TMB info
annot_info <- read.delim(file.path("results", "annotation.txt"), header = TRUE, check.names = TRUE)
annot_info <- annot_info %>%
  filter(Sample %in% colnames(mat),
         Sequencing_Experiment != "RNA-Seq") %>%
  remove_rownames() %>%
  column_to_rownames('Sample') %>%
  as.data.frame()
mat <- mat[,rownames(annot_info)]
annot_info <- annot_info[colnames(mat),]

# annotation 1
ha = HeatmapAnnotation(df = annot_info %>% dplyr::select(-c(Sequencing_Experiment)), col = list(
  Sequencing_Experiment = c("WGS" = "red",
                            "RNA-Seq, WGS" = "orange"),
  Tumor_Descriptor = c("Progressive" = "#827397", 
                       "Primary" = "#d8b9c3",
                       "Initial CNS Tumor" = "#cee397",
                       "Recurrence" = "#363062",
                       "Second Malignancy" = "#005082"),
  Integrated_Diagnosis = c("High-grade glioma/astrocytoma (WHO grade III/IV)"="lightseagreen",
                           "Astrocytoma;Oligoastrocytoma" = "mediumorchid2",
                           "Astrocytoma" = "brown2", 
                           "Glioblastoma" = "orange",
                           "Ependymoma" = "darkgreen",
                           "Low-grade glioma/astrocytoma (WHO grade I/II)" = "blue2"),
  Sex = c("Female" = "deeppink4",
             "Male" = "navy"),
  Age = c("0-14" = "gold",
          "14-33.5" = "purple")),
  gp = gpar(col = "#595959"), simple_anno_size = unit(4, "mm"), annotation_name_side = "left",
  annotation_legend_param = list(
    # Sequencing_Experiment = list(nrow =3)
  ))

# annotation 2
amp = ifelse(apply(mat, 1, function(x) sum(grepl("GAI", x) + grepl("LOS",x))/length(x) > 0),"Copy number alteration","Genetic and Fusion alteration")
amp = factor(amp, levels = c("Genetic and Fusion alteration", "Copy number alteration"))

# oncoprint
ht = oncoPrint(mat, get_type = function(x)strsplit(x, ";")[[1]],
               alter_fun = alter_fun, col = col, show_column_names =TRUE, column_names_gp = gpar(fontsize = 9),
               column_names_side = "top",
               top_annotation = ha,
               right_annotation = NULL,
               row_names_side = "left",
               pct_side = "right",
               row_split = amp,
               remove_empty_rows = TRUE,
               heatmap_legend_param = list(title = "Alterations", nrow = 9, title_position = "topleft", direction = "horizontal",
                                           at = c("GAI","LOS","MIS","FUS","NOS","FSD","FSI","SPS","IFD"),
                                           labels = c("Copy gain", "Copy loss", "Misense","Gene Fusion","Nonsense","Frame_Shift_Del","Frame_Shift_Ins","Splice site","In_Frame_Del")
               ))


pdf(file = 'results/oncoplot_norna.pdf', width = 15, height = 12) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# column order by sex
sex_ordered <- annot_info %>% 
  arrange(Sex)
ht = oncoPrint(mat, get_type = function(x)strsplit(x, ";")[[1]],
               alter_fun = alter_fun, col = col, show_column_names =TRUE, column_names_gp = gpar(fontsize = 9),
               column_order = match(rownames(sex_ordered), colnames(mat)),
               column_names_side = "top",
               top_annotation = ha,
               right_annotation = NULL,
               row_names_side = "left",
               pct_side = "right",
               row_split = amp,
               remove_empty_rows = TRUE,
               heatmap_legend_param = list(title = "Alterations", nrow = 9, title_position = "topleft", direction = "horizontal",
                                           at = c("GAI","LOS","MIS","FUS","NOS","FSD","FSI","SPS","IFD"),
                                           labels = c("Copy gain", "Copy loss", "Misense","Gene Fusion","Nonsense","Frame_Shift_Del","Frame_Shift_Ins","Splice site","In_Frame_Del")
               ))
pdf(file = 'results/oncoplot_orderby_sex_norna.pdf', width = 15, height = 12) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# column order by H3F3A status
n <- grep("H3F3A$", rownames(mat))
ht = oncoPrint(mat, get_type = function(x)strsplit(x, ";")[[1]],
               alter_fun = alter_fun, col = col, show_column_names =TRUE, column_names_gp = gpar(fontsize = 9),
               column_order = order(mat[n,], decreasing = T),
               column_names_side = "top",
               top_annotation = ha,
               right_annotation = NULL,
               row_names_side = "left",
               pct_side = "right",
               row_split = amp,
               remove_empty_rows = TRUE,
               heatmap_legend_param = list(title = "Alterations", nrow = 9, title_position = "topleft", direction = "horizontal",
                                           at = c("GAI","LOS","MIS","FUS","NOS","FSD","FSI","SPS","IFD"),
                                           labels = c("Copy gain", "Copy loss", "Misense","Gene Fusion","Nonsense","Frame_Shift_Del","Frame_Shift_Ins","Splice site","In_Frame_Del")
               ))
pdf(file = 'results/oncoplot_orderby_H3F3A_norna.pdf', width = 15, height = 12) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()


# column order by sex + H3F3A
tmp <- reshape2::melt(mat) %>%
  filter(Var1 == "H3F3A")
tmp <- annot_info %>% 
  rownames_to_column("sample_id") %>%
  inner_join(tmp, by = c("sample_id" = "Var2")) %>%
  arrange(Sex, desc(value)) %>%
  column_to_rownames("sample_id")
ht = oncoPrint(mat, get_type = function(x)strsplit(x, ";")[[1]],
               alter_fun = alter_fun, col = col, show_column_names =TRUE, column_names_gp = gpar(fontsize = 9),
               column_order = match(rownames(tmp), colnames(mat)),
               column_names_side = "top",
               top_annotation = ha,
               right_annotation = NULL,
               row_names_side = "left",
               pct_side = "right",
               row_split = amp,
               remove_empty_rows = TRUE,
               heatmap_legend_param = list(title = "Alterations", nrow = 9, title_position = "topleft", direction = "horizontal",
                                           at = c("GAI","LOS","MIS","FUS","NOS","FSD","FSI","SPS","IFD"),
                                           labels = c("Copy gain", "Copy loss", "Misense","Gene Fusion","Nonsense","Frame_Shift_Del","Frame_Shift_Ins","Splice site","In_Frame_Del")
               ))
pdf(file = 'results/oncoplot_orderby_sex_H3F3A_status_norna.pdf', width = 15, height = 12) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# column order by sex + age
sex_age_ordered <- annot_info %>% 
  arrange(Sex, Age)
ht = oncoPrint(mat, get_type = function(x)strsplit(x, ";")[[1]],
               alter_fun = alter_fun, col = col, show_column_names =TRUE, column_names_gp = gpar(fontsize = 9),
               column_order = match(rownames(sex_age_ordered), colnames(mat)),
               column_names_side = "top",
               top_annotation = ha,
               right_annotation = NULL,
               row_names_side = "left",
               pct_side = "right",
               row_split = amp,
               remove_empty_rows = TRUE,
               heatmap_legend_param = list(title = "Alterations", nrow = 9, title_position = "topleft", direction = "horizontal",
                                           at = c("GAI","LOS","MIS","FUS","NOS","FSD","FSI","SPS","IFD"),
                                           labels = c("Copy gain", "Copy loss", "Misense","Gene Fusion","Nonsense","Frame_Shift_Del","Frame_Shift_Ins","Splice site","In_Frame_Del")
               ))
pdf(file = 'results/oncoplot_orderby_sex_age_norna.pdf', width = 15, height = 12) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# column order by sex + age + histone
tmp <- reshape2::melt(mat) %>%
  filter(Var1 == "H3F3A")
tmp <- annot_info %>% 
  rownames_to_column("sample_id") %>%
  inner_join(tmp, by = c("sample_id" = "Var2")) %>%
  arrange(Sex, Age, desc(value)) %>%
  column_to_rownames("sample_id")
ht = oncoPrint(mat, get_type = function(x)strsplit(x, ";")[[1]],
               alter_fun = alter_fun, col = col, show_column_names =TRUE, column_names_gp = gpar(fontsize = 9),
               column_order = match(rownames(tmp), colnames(mat)),
               column_names_side = "top",
               top_annotation = ha,
               right_annotation = NULL,
               row_names_side = "left",
               pct_side = "right",
               row_split = amp,
               remove_empty_rows = TRUE,
               heatmap_legend_param = list(title = "Alterations", nrow = 9, title_position = "topleft", direction = "horizontal",
                                           at = c("GAI","LOS","MIS","FUS","NOS","FSD","FSI","SPS","IFD"),
                                           labels = c("Copy gain", "Copy loss", "Misense","Gene Fusion","Nonsense","Frame_Shift_Del","Frame_Shift_Ins","Splice site","In_Frame_Del")
               ))
pdf(file = 'results/oncoplot_orderby_sex_age_H3F3A_status_norna.pdf', width = 15, height = 12) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# annotation 1
ha = HeatmapAnnotation(df = annot_info %>% dplyr::select(Tumor_Descriptor, Integrated_Diagnosis, Age, Sex), 
                       col = list(
  Sequencing_Experiment = c("WGS" = "red",
                            "RNA-Seq, WGS" = "orange"),
  Tumor_Descriptor = c("Progressive" = "#827397", 
                       "Primary" = "#d8b9c3",
                       "Initial CNS Tumor" = "#cee397",
                       "Recurrence" = "#363062",
                       "Second Malignancy" = "#005082"),
  Integrated_Diagnosis = c("High-grade glioma/astrocytoma (WHO grade III/IV)"="lightseagreen",
                           "Astrocytoma;Oligoastrocytoma" = "mediumorchid2",
                           "Astrocytoma" = "brown2", 
                           "Glioblastoma" = "orange",
                           "Ependymoma" = "darkgreen",
                           "Low-grade glioma/astrocytoma (WHO grade I/II)" = "blue2"),
  Age = c("0-14" = "gold",
          "14-33.5" = "purple"),
  Sex = c("Female" = "deeppink4",
          "Male" = "navy")),
  gp = gpar(col = "#595959"), simple_anno_size = unit(4, "mm"), annotation_name_side = "left",
  annotation_legend_param = list(
    # Sequencing_Experiment = list(nrow =3)
  ))

# annotation 2
amp = ifelse(apply(mat, 1, function(x) sum(grepl("GAI", x) + grepl("LOS",x))/length(x) > 0),"Copy number alteration","Genetic and Fusion alteration")
amp = factor(amp, levels = c("Genetic and Fusion alteration", "Copy number alteration"))

# column order age
age_ordered <- annot_info %>% 
  arrange(Age)
ht = oncoPrint(mat, get_type = function(x)strsplit(x, ";")[[1]],
               alter_fun = alter_fun, col = col, show_column_names =TRUE, column_names_gp = gpar(fontsize = 9),
               column_order = match(rownames(age_ordered), colnames(mat)),
               column_names_side = "top",
               top_annotation = ha,
               right_annotation = NULL,
               row_names_side = "left",
               pct_side = "right",
               row_split = amp,
               remove_empty_rows = TRUE,
               heatmap_legend_param = list(title = "Alterations", nrow = 9, title_position = "topleft", direction = "horizontal",
                                           at = c("GAI","LOS","MIS","FUS","NOS","FSD","FSI","SPS","IFD"),
                                           labels = c("Copy gain", "Copy loss", "Misense","Gene Fusion","Nonsense","Frame_Shift_Del","Frame_Shift_Ins","Splice site","In_Frame_Del")
               ))
pdf(file = 'results/oncoplot_orderby_age_norna.pdf', width = 15, height = 12) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# column order by age + histone
tmp <- reshape2::melt(mat) %>%
  filter(Var1 == "H3F3A")
tmp <- annot_info %>% 
  rownames_to_column("sample_id") %>%
  inner_join(tmp, by = c("sample_id" = "Var2")) %>%
  arrange(Age, desc(value)) %>%
  column_to_rownames("sample_id")
ht = oncoPrint(mat, get_type = function(x)strsplit(x, ";")[[1]],
               alter_fun = alter_fun, col = col, show_column_names =TRUE, column_names_gp = gpar(fontsize = 9),
               column_order = match(rownames(tmp), colnames(mat)),
               column_names_side = "top",
               top_annotation = ha,
               right_annotation = NULL,
               row_names_side = "left",
               pct_side = "right",
               row_split = amp,
               remove_empty_rows = TRUE,
               heatmap_legend_param = list(title = "Alterations", nrow = 9, title_position = "topleft", direction = "horizontal",
                                           at = c("GAI","LOS","MIS","FUS","NOS","FSD","FSI","SPS","IFD"),
                                           labels = c("Copy gain", "Copy loss", "Misense","Gene Fusion","Nonsense","Frame_Shift_Del","Frame_Shift_Ins","Splice site","In_Frame_Del")
               ))
pdf(file = 'results/oncoplot_orderby_age_H3F3A_status_norna.pdf', width = 15, height = 12) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# column order by age + sex + histone
tmp <- reshape2::melt(mat) %>%
  filter(Var1 == "H3F3A")
tmp <- annot_info %>% 
  rownames_to_column("sample_id") %>%
  inner_join(tmp, by = c("sample_id" = "Var2")) %>%
  arrange(Age, Sex, desc(value)) %>%
  column_to_rownames("sample_id")
ht = oncoPrint(mat, get_type = function(x)strsplit(x, ";")[[1]],
               alter_fun = alter_fun, col = col, show_column_names =TRUE, column_names_gp = gpar(fontsize = 9),
               column_order = match(rownames(tmp), colnames(mat)),
               column_names_side = "top",
               top_annotation = ha,
               right_annotation = NULL,
               row_names_side = "left",
               pct_side = "right",
               row_split = amp,
               remove_empty_rows = TRUE,
               heatmap_legend_param = list(title = "Alterations", nrow = 9, title_position = "topleft", direction = "horizontal",
                                           at = c("GAI","LOS","MIS","FUS","NOS","FSD","FSI","SPS","IFD"),
                                           labels = c("Copy gain", "Copy loss", "Misense","Gene Fusion","Nonsense","Frame_Shift_Del","Frame_Shift_Ins","Splice site","In_Frame_Del")
               ))
pdf(file = 'results/oncoplot_orderby_age_sex_H3F3A_status_norna.pdf', width = 15, height = 12) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()





# Function: data availability heatmap with diagnosis as top level annotation

# load libraries
suppressPackageStartupMessages({
  library(reshape2)
  library(tidyverse) 
  library(ComplexHeatmap)
  library(circlize) 
})

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data", "v1")
analyses_dir <- file.path(root_dir, "analyses", "data-availability")
output_dir <- file.path(analyses_dir, "results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# read histologies
annot <- read_tsv(file.path(data_dir, "Hope-GBM-histologies.tsv"))

# remove redundancy
annot$HOPE_diagnosis_type[annot$HOPE_diagnosis_type == "Recurrent"] = "Recurrence"
annot$HOPE_diagnosis_type[annot$HOPE_diagnosis_type == "recurrent"] = "Recurrence"
annot$HOPE_diagnosis_type[annot$HOPE_diagnosis_type == "Recurrent, residual"] = "Recurrence"
annot$HOPE_diagnosis_type[annot$HOPE_diagnosis_type == "Primary"] = "Initial CNS Tumor"

# subset to columns of interest
annot <- annot %>%
  filter(!is.na(molecular_subtype),
         !is.na(HARMONY_age_class_derived)) %>%
  dplyr::rename("Age" = "HARMONY_age_class_derived",
                "Gender" = "HARMONY_Gender",
                "Diagnosis" = "HOPE_diagnosis",
                "Diagnosis Type" = "HOPE_diagnosis_type",
                "Annotation" = "HOPE_sample_annotation",
                "Tumor Location" = "HOPE_Tumor.Location.condensed") %>%
  dplyr::select(sample_id, Age, Gender, Diagnosis, `Diagnosis Type`, Annotation, `Tumor Location`) %>%
  unique() %>%
  column_to_rownames('sample_id') %>%
  dplyr::arrange(Age, Gender, Diagnosis, `Diagnosis Type`, Annotation, `Tumor Location`)

# 1) Two age groups

# merge age into 2 groups
plot_df <- annot %>%
  mutate(Age = as.character(Age)) %>%
  mutate(Age = ifelse(Age %in% c("(15,26]", "(26,40]"), "(15,40]", Age))

# split by age
split <- factor(plot_df$Age, levels = c("[0,15]", "(15,40]"))
col_fun1 <- list("High-grade glioma/astrocytoma (WHO grade III/IV)" = "lightseagreen",
                 "Diffuse Midline Glioma" = "darkgreen",
                 "Astrocytoma;Oligoastrocytoma" = "mediumorchid2",
                 "Astrocytoma" = "brown2", 
                 "Glioblastoma" = "orange",
                 "Pleomorphic xanthoastrocytoma" = "#005082",
                 "Male" = "navy",
                 "Female" = "deeppink4",
                 "[0,15]" = "gold",
                 "(15,40]" = "purple",
                 "Initial CNS Tumor" = "#cee397",
                 "Progressive" = "#827397",
                 "Recurrence" = "#363062",
                 "Second Malignancy" = "#005082",
                 "Treatment naive" = "lightgray",
                 "Post-treatment" = "gray50",
                 "Post-mortem" = "black",
                 "Cortical" = "magenta",
                 "Other/Multiple locations/NOS" = "pink",
                 "Midline" = "purple",
                 "Cerebellar" = "navy")

# generate plot
circos.clear()
pdf(file = file.path(output_dir, "hope_clinical_data_availability_diagnosis_age_two_groups.pdf"), width = 10, height = 10)
circos.par(start.degree = 30, gap.degree = 1, points.overflow.warning = FALSE)
circos.heatmap(plot_df, 
               split = split, 
               col = unlist(col_fun1), 
               track.height = 0.4, 
               bg.border = "gray50", bg.lwd = 2, cell.lwd = 4,
               show.sector.labels = F, cell.border = "white")
# add border colors to sectors
for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = 1)
  circos.rect(CELL_META$cell.xlim[1], CELL_META$cell.ylim[1],
              CELL_META$cell.xlim[2], CELL_META$cell.ylim[2],
              col = NA, border = unlist(col_fun1)[sn])
}
lgd_diagnosis = Legend(title = "Diagnosis", 
                       at = c("Pleomorphic xanthoastrocytoma",
                              "Diffuse Midline Glioma",
                              "Astrocytoma;Oligoastrocytoma",
                              "High-grade glioma/astrocytoma (WHO grade III/IV)",
                              "Astrocytoma", 
                              "Glioblastoma"), 
                       legend_gp = gpar(fill = c("lightseagreen", "darkgreen", "mediumorchid2", "brown2", "orange", "#005082")))
lgd_gender = Legend(title = "Sex", 
                    at = c("Male", "Female"), 
                    legend_gp = gpar(fill = c("navy", "deeppink4")))
lgd_age = Legend(title = "Age", 
                 at = c("[0,15]", "(15,40]"), 
                 legend_gp = gpar(fill = c("gold", "purple")))
lgd_dtype = Legend(title = "Diagnosis Type",
                   at = c("Initial CNS Tumor", "Progressive", "Recurrence", "Second Malignancy") ,
                   legend_gp = gpar(fill = c("#cee397", "#827397", "#363062", "#005082")))
lgd_annot = Legend(title = "Annotation",
                   at = c("Treatment naive", "Post-treatment", "Post-mortem"),
                   legend_gp = gpar(fill = c("lightgray", "gray50", "black")))
lgd_tumor_location = Legend(title = "Tumor Location",
                            at = c("Cortical", "Other/Multiple locations/NOS", "Midline", "Cerebellar"),
                            legend_gp = gpar(fill = c("magenta", "pink", "purple", "navy")))
h = dev.size()[2]
circle_size = unit(1, "snpc")
lgd_list = packLegend(lgd_age, lgd_gender, max_height = unit(0.9*h, "inch"), direction = "horizontal")
lgd_list1 = packLegend(lgd_diagnosis, max_height = unit(0.9*h, "inch"), direction = "horizontal")
lgd_list2 = packLegend(lgd_dtype, lgd_annot, max_height = unit(0.9*h, "inch"), direction = "horizontal")
lgd_list3 = packLegend(lgd_tumor_location, max_height = unit(0.9*h, "inch"), direction = "horizontal")
draw(lgd_list, x = unit(125, "mm"), y = unit(170, "mm")) 
draw(lgd_list1, x = unit(130, "mm"), y = unit(140, "mm")) 
draw(lgd_list2, x = unit(123, "mm"), y = unit(110, "mm")) 
draw(lgd_list3, x = unit(114, "mm"), y = unit(85, "mm")) 
circos.clear()
dev.off()

# 2) Three age groups

plot_df <- annot

# split by age
split <- factor(plot_df$Age, levels = c("[0,15]", "(15,26]", "(26,40]"))
col_fun1 <- list("High-grade glioma/astrocytoma (WHO grade III/IV)" = "lightseagreen",
                 "Diffuse Midline Glioma" = "darkgreen",
                 "Astrocytoma;Oligoastrocytoma" = "mediumorchid2",
                 "Astrocytoma" = "brown2", 
                 "Glioblastoma" = "orange",
                 "Pleomorphic xanthoastrocytoma" = "#005082",
                 "Male" = "navy",
                 "Female" = "deeppink4",
                 "[0,15]" = "gold",
                 "(15,26]" = "purple",
                 "(26,40]" = "darkgreen",
                 "Initial CNS Tumor" = "#cee397",
                 "Progressive" = "#827397",
                 "Recurrence" = "#363062",
                 "Second Malignancy" = "#005082",
                 "Treatment naive" = "lightgray",
                 "Post-treatment" = "gray50",
                 "Post-mortem" = "black",
                 "Cortical" = "magenta",
                 "Other/Multiple locations/NOS" = "pink",
                 "Midline" = "purple",
                 "Cerebellar" = "navy")

# generate plot
circos.clear()
pdf(file = file.path(output_dir, "hope_clinical_data_availability_diagnosis_age_three_groups.pdf"), width = 10, height = 10)
circos.par(start.degree = 30, gap.degree = 1, points.overflow.warning = FALSE)
circos.heatmap(plot_df, 
               split = split, 
               col = unlist(col_fun1), 
               track.height = 0.4, 
               bg.border = "gray50", bg.lwd = 2, cell.lwd = 4,
               show.sector.labels = F, cell.border = "white")
# add border colors to sectors
for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = 1)
  circos.rect(CELL_META$cell.xlim[1], CELL_META$cell.ylim[1],
              CELL_META$cell.xlim[2], CELL_META$cell.ylim[2],
              col = NA, border = unlist(col_fun1)[sn])
}
lgd_diagnosis = Legend(title = "Diagnosis", 
                       at = c("Pleomorphic xanthoastrocytoma",
                              "Diffuse Midline Glioma",
                              "Astrocytoma;Oligoastrocytoma",
                              "High-grade glioma/astrocytoma (WHO grade III/IV)",
                              "Astrocytoma", 
                              "Glioblastoma"), 
                       legend_gp = gpar(fill = c("lightseagreen", "darkgreen", "mediumorchid2", "brown2", "orange", "#005082")))
lgd_gender = Legend(title = "Sex", 
                    at = c("Male", "Female"), 
                    legend_gp = gpar(fill = c("navy", "deeppink4")))
lgd_age = Legend(title = "Age", 
                 at = c("[0,15]", "(15,26]", "(26,40]"), 
                 legend_gp = gpar(fill = c("gold", "purple", "darkgreen")))
lgd_dtype = Legend(title = "Diagnosis Type",
                   at = c("Initial CNS Tumor", "Progressive", "Recurrence", "Second Malignancy") ,
                   legend_gp = gpar(fill = c("#cee397", "#827397", "#363062", "#005082")))
lgd_annot = Legend(title = "Annotation",
                   at = c("Treatment naive", "Post-treatment", "Post-mortem"),
                   legend_gp = gpar(fill = c("lightgray", "gray50", "black")))
lgd_tumor_location = Legend(title = "Tumor Location",
                            at = c("Cortical", "Other/Multiple locations/NOS", "Midline", "Cerebellar"),
                            legend_gp = gpar(fill = c("magenta", "pink", "purple", "navy")))
h = dev.size()[2]
circle_size = unit(1, "snpc")
lgd_list = packLegend(lgd_age, lgd_gender, max_height = unit(0.9*h, "inch"), direction = "horizontal")
lgd_list1 = packLegend(lgd_diagnosis, max_height = unit(0.9*h, "inch"), direction = "horizontal")
lgd_list2 = packLegend(lgd_dtype, lgd_annot, max_height = unit(0.9*h, "inch"), direction = "horizontal")
lgd_list3 = packLegend(lgd_tumor_location, max_height = unit(0.9*h, "inch"), direction = "horizontal")
draw(lgd_list, x = unit(125, "mm"), y = unit(170, "mm")) 
draw(lgd_list1, x = unit(130, "mm"), y = unit(140, "mm")) 
draw(lgd_list2, x = unit(123, "mm"), y = unit(110, "mm")) 
draw(lgd_list3, x = unit(114, "mm"), y = unit(85, "mm")) 
circos.clear()
dev.off()

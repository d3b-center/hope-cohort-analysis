# Function: data availability heatmap with diagnosis as top level annotation (with Age as continuous variable)

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
                "Age_at_Initial_Diagnosis" = "HOPE_Age_at_Initial_Diagnosis",
                "Gender" = "HARMONY_Gender",
                "Diagnosis" = "HOPE_diagnosis",
                "Diagnosis Type" = "HOPE_diagnosis_type",
                "Annotation" = "HOPE_sample_annotation",
                "Tumor Location" = "HOPE_Tumor.Location.condensed") %>%
  dplyr::mutate(Age_at_Initial_Diagnosis = Age_at_Initial_Diagnosis/365) %>%
  dplyr::select(sample_id, Age, Age_at_Initial_Diagnosis, Gender, Diagnosis, `Diagnosis Type`, Annotation, `Tumor Location`) %>%
  unique() %>%
  column_to_rownames('sample_id') %>%
  dplyr::arrange(Age, Age_at_Initial_Diagnosis, Gender, Diagnosis, `Diagnosis Type`, Annotation, `Tumor Location`)

# merge age into 2 groups
plot_df <- annot %>%
  mutate(Age = as.character(Age)) %>%
  mutate(Age = ifelse(Age %in% c("(15,26]", "(26,40]"), "(15,40]", Age))

# variable for splitting the circular heatmap
split <- factor(plot_df$Age, levels = c("[0,15]", "(15,40]"))

# colour mapping for factor variables
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

# color mapping for numeric age variable
col_fun_age <- colorRamp2(breaks = c(min(plot_df$Age_at_Initial_Diagnosis, na.rm = T), median(plot_df$Age_at_Initial_Diagnosis, na.rm = T), max(plot_df$Age_at_Initial_Diagnosis, na.rm = T)),
                          colors = c("lightskyblue1","skyblue","dodgerblue4"))

# generate plot 
circos.clear()
pdf(file = file.path(output_dir, "hope_clinical_data_availability_diagnosis_age_continuous.pdf"), width = 10, height = 10)
circos.par(start.degree = 30, gap.degree = 1, points.overflow.warning = FALSE)
circos.heatmap(mat = plot_df %>% dplyr::select(Age), 
               split = split,
               col = list("[0,15]" = "white", "(15,40]" = "white"), 
               track.height = 0.001, 
               cell.border = "white",
               bg.lwd = 2, cell.lwd = 2)
circos.heatmap(mat = plot_df %>% dplyr::select(Age_at_Initial_Diagnosis), 
               col = col_fun_age, 
               track.height = 0.06, 
               cell.border = "white",
               bg.lwd = 2, cell.lwd = 2)
circos.heatmap(plot_df %>% 
                 dplyr::select(-c(Age, Age_at_Initial_Diagnosis)), 
               col = unlist(col_fun1),
               track.height = 0.3, 
               cell.border = "white",
               bg.lwd = 2, cell.lwd = 2)

# legends
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
lgd_age = Legend(title = "Age (years)", col_fun = col_fun_age, at = c(0, 20, 39.9), labels = c(0, 20, 39.9))
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

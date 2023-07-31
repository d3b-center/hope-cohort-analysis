# data availability heatmap (with Age as continuous variable)
suppressPackageStartupMessages({
  library(reshape2)
  library(tidyverse) 
  library(ComplexHeatmap)
  library(circlize) 
})

# output directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
output_dir <- file.path(root_dir, "results", "data_plots")
dir.create(output_dir, recursive = T, showWarnings = F)

# updated clinical data from Mateusz
annot <- readr::read_tsv(file.path(data_dir, "hopeonly_clinical_table_011823.tsv"))
annot$diagnosis_type[annot$diagnosis_type == "Recurrent"] = "Recurrence"
annot$diagnosis_type[annot$diagnosis_type == "recurrent"] = "Recurrence"
annot$diagnosis_type[annot$diagnosis_type == "Recurrent, residual"] = "Recurrence"
annot$diagnosis_type[annot$diagnosis_type == "Primary"] = "Initial CNS Tumor"

# get Age from Nicole's file
age_info <- readxl::read_xlsx(file.path(data_dir, "clini_m_030722-for_Komal.xlsx"))
annot <- annot %>%
  inner_join(age_info %>% dplyr::select(age.class, id), by = c("Sample_ID" = "id")) %>%
  dplyr::rename("Age" = "age.class")
annot$Age <- factor(annot$Age, levels = c("[0,15]", "(15,40]"))

# select columns of interest
annot <- annot %>%
  dplyr::select(Sample_ID, Age, Age_at_Initial_Diagnosis, Gender, diagnosis, diagnosis_type, sample_annotation, Tumor.Location.condensed) %>%
  dplyr::mutate(Age_at_Initial_Diagnosis = Age_at_Initial_Diagnosis/365) %>%
  dplyr::rename("Diagnosis" = "diagnosis",
                "Diagnosis Type" = "diagnosis_type",
                "Annotation" = "sample_annotation",
                "Tumor Location" = "Tumor.Location.condensed") %>%
  column_to_rownames('Sample_ID') %>%
  dplyr::arrange(Age, Age_at_Initial_Diagnosis, Gender, Diagnosis, `Diagnosis Type`, Annotation, `Tumor Location`)

# variable for splitting the circular heatmap
split <- factor(annot$Age, levels = c("[0,15]", "(15,40]"))

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

# color mapping for numeric variable
col_fun_age <- colorRamp2(breaks = c(min(annot$Age_at_Initial_Diagnosis), median(annot$Age_at_Initial_Diagnosis), max(annot$Age_at_Initial_Diagnosis)),
                          colors = c("lightskyblue1","skyblue","dodgerblue4"))

pdf(file = file.path(output_dir, "hope_cohort_data_availability_clinical_v2.pdf"), width = 10, height = 10)
circos.clear()
circos.par(start.degree = 30, gap.degree = 1, points.overflow.warning = FALSE)
circos.heatmap(mat = annot %>% dplyr::select(Age), 
               split = split,
               col = list("[0,15]" = "white", "(15,40]" = "white"), 
               track.height = 0.001, 
               cell.border = "white",
               bg.lwd = 2, cell.lwd = 2)
circos.heatmap(mat = annot %>% dplyr::select(Age_at_Initial_Diagnosis), 
               col = col_fun_age, 
               track.height = 0.06, 
               cell.border = "white",
               bg.lwd = 2, cell.lwd = 2)
circos.heatmap(annot %>% 
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

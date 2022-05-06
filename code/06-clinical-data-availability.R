# data availability heatmap

suppressPackageStartupMessages({
  library(reshape2)
  library(tidyverse) 
  library(ComplexHeatmap)
  library(circlize) 
})

# only clinical data
annot <- readxl::read_xlsx('data/clini_m_030722-for_Komal.xlsx')
annot <- annot %>%
  filter(Sample_ID != "7316-4065") %>%
  dplyr::select(Sample_ID, Diagnosis_demoted, age.class, Gender, Diagnosis.Type_demoted, Sample.annotation, Tumor.Location.condensed3) %>%
  column_to_rownames('Sample_ID')
split <- annot$Diagnosis_demoted
col_fun1 <- list("High-grade glioma/astrocytoma (WHO grade III/IV)"="lightseagreen",
                 "Astrocytoma;Oligoastrocytoma" = "mediumorchid2",
                 "Astrocytoma" = "brown2", 
                 "Glioblastoma" = "orange",
                 # "Low-grade glioma/astrocytoma (WHO grade I/II)" = "blue2",
                 "Male" = "navy",
                 "Female" = "deeppink4",
                 "[0,15]" = "gold",
                 "(15,40]" = "purple",
                 "Progressive" = "#827397", 
                 "Initial CNS Tumor" = "#cee397",
                 "Recurrence" = "#363062",
                 "Second Malignancy" = "#005082",
                 "Treatment naive" = "lightgray",
                 "Post-treatment" = "gray50",
                 "Post-mortem" = "black",
                 "Cortical" = "magenta",
                 "Other/Multiple locations/NOS" = "pink",
                 "Midline" = "purple",
                 "Cerebellar" = "navy")
circos.clear()
pdf(file = "results/hope_cohort_data_availability_clinical.pdf", width = 10, height = 10)
circos.par(start.degree = 30, gap.degree = 1, points.overflow.warning = FALSE)
circos.heatmap(annot, 
               split = split, 
               col = unlist(col_fun1), 
               track.height = 0.4, 
               bg.border = "gray50", bg.lwd = 2,
               show.sector.labels = F, cell.border = "white")
# add border colors to sectors
for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = 1)
  circos.rect(CELL_META$cell.xlim[1], CELL_META$cell.ylim[1],
              CELL_META$cell.xlim[2], CELL_META$cell.ylim[2],
              col = NA, border = unlist(col_fun1)[sn])
}
lgd_diagnosis = Legend(title = "Diagnosis", 
                       at = c("High-grade glioma/astrocytoma (WHO grade III/IV)",
                              # "Low-grade glioma/astrocytoma (WHO grade I/II)",
                              "Astrocytoma;Oligoastrocytoma",
                              "Astrocytoma",
                              "Glioblastoma"), 
                       legend_gp = gpar(fill = c("lightseagreen", "blue2", "mediumorchid2", "brown2", "orange")))
lgd_gender = Legend(title = "Sex", 
                 at = c("Male", "Female"), 
                 legend_gp = gpar(fill = c("navy", "deeppink4")))
lgd_age = Legend(title = "Age", 
                 at = c("[0,15]", "(15,40]"), 
                 legend_gp = gpar(fill = c("gold", "purple")))
lgd_dtype = Legend(title = "Diagnosis_type",
                   at = c("Progressive", "Initial CNS Tumor",
                          "Recurrence", "Second Malignancy"),
                   legend_gp = gpar(fill = c("#827397", "#cee397", "#363062", "#005082")))
lgd_annot = Legend(title = "Annotation",
                   at = c("Treatment naive", "Post-treatment", "Post-mortem"),
                   legend_gp = gpar(fill = c("lightgray", "gray50", "black")))
lgd_tumor_location = Legend(title = "Tumor Location",
                            at = c("Cortical", "Other/Multiple locations/NOS",
                                   "Midline", "Cerebellar"),
                            legend_gp = gpar(fill = c("magenta", "pink", "purple", "navy")))
h = dev.size()[2]
circle_size = unit(1, "snpc")
lgd_list = packLegend(lgd_diagnosis, max_height = unit(0.9*h, "inch"), direction = "horizontal")
lgd_list1 = packLegend(lgd_age, lgd_gender, max_height = unit(0.9*h, "inch"), direction = "horizontal")
lgd_list2 = packLegend(lgd_dtype, lgd_annot, max_height = unit(0.9*h, "inch"), direction = "horizontal")
lgd_list3 = packLegend(lgd_tumor_location, max_height = unit(0.9*h, "inch"), direction = "horizontal")
draw(lgd_list, x = unit(125, "mm"), y = unit(160, "mm")) 
draw(lgd_list1, x = unit(100, "mm"), y = unit(137, "mm")) 
draw(lgd_list2, x = unit(116, "mm"), y = unit(115, "mm")) 
draw(lgd_list3, x = unit(106, "mm"), y = unit(88, "mm")) 
circos.clear()
dev.off()





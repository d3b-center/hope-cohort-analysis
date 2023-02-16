# data availability heatmap

suppressPackageStartupMessages({
  library(reshape2)
  library(tidyverse) 
  library(ComplexHeatmap)
  library(circlize) 
})

# only clinical data
annot <- readxl::read_xlsx('data/clini_m_030722-for_Komal.xlsx')
annot_grade <- readxl::read_xlsx('data/Project Hope - free text diagnosis -GradesKomalCircos.xlsx')

annot <- annot %>%
  inner_join(annot_grade, by = c("Sample_ID" = "SDG-ID")) %>%
  filter(Sample_ID != "7316-4065") %>%
  dplyr::rename("WHO_Grade" = "WHO Grade") %>%
  dplyr::select(Sample_ID, WHO_Grade, age.class, Gender, Diagnosis.Type_demoted, Sample.annotation, Tumor.Location.condensed3) %>%
  column_to_rownames('Sample_ID') %>%
  arrange(WHO_Grade, age.class, Gender, Diagnosis.Type_demoted, Sample.annotation, Tumor.Location.condensed3)
split <- factor(annot$WHO_Grade, levels = c("3","4", "N/A"))
col_fun1 <- list("3" = "#ffac59",
                 "4" = "#ff5959",
                 "N/A" = "#FFFFFF",
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
pdf(file = "results/hope_cohort_data_availability_clinical_with_WHOGrade.pdf", width = 10, height = 10)
circos.par(start.degree = 30, gap.degree = 1, points.overflow.warning = FALSE)
circos.heatmap(annot, 
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
lgd_grade = Legend(title = "WHO_Grade",
                   at = c("3", "4", "N/A"),
                   legend_gp = gpar(fill = c("#ffac59", "#ff5959", "#FFFFFF")))
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
lgd_list = packLegend(lgd_grade, max_height = unit(0.9*h, "inch"), direction = "horizontal")
lgd_list1 = packLegend(lgd_age, lgd_gender, max_height = unit(0.9*h, "inch"), direction = "horizontal")
lgd_list2 = packLegend(lgd_dtype, lgd_annot, max_height = unit(0.9*h, "inch"), direction = "horizontal")
lgd_list3 = packLegend(lgd_tumor_location, max_height = unit(0.9*h, "inch"), direction = "horizontal")
draw(lgd_list, x = unit(96, "mm"), y = unit(156, "mm")) 
draw(lgd_list1, x = unit(104, "mm"), y = unit(137, "mm")) 
draw(lgd_list2, x = unit(120, "mm"), y = unit(115, "mm")) 
draw(lgd_list3, x = unit(110, "mm"), y = unit(88, "mm")) 
circos.clear()
dev.off()

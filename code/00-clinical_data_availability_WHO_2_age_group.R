# data availability heatmap

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
annot = readr::read_tsv(file.path(data_dir, "hopeonly_clinical_table_011823.tsv"))

# get Age from Nicole's file
age_info <- readxl::read_xlsx(file.path(data_dir, "clini_m_030722-for_Komal.xlsx"))
annot <- annot %>%
  left_join(age_info %>% dplyr::select(age.class, id), by = c("Sample_ID" = "id")) %>%
  dplyr::rename("Age" = "age.class")

annot <- annot %>%
  dplyr::select(Sample_ID, WHO.Grade, Age, Gender, diagnosis_type_simple, sample_annotation, Tumor.Location.condensed) %>%
  dplyr::rename("Diagnosis Type" = "diagnosis_type_simple",
                "Annotation" = "sample_annotation",
                "Tumor Location" = "Tumor.Location.condensed",
                "WHO Grade" = "WHO.Grade") %>%
  column_to_rownames('Sample_ID') %>%
  arrange(`WHO Grade`, Age, Gender, `Diagnosis Type`, Annotation, `Tumor Location`)
annot$`WHO Grade`[is.na(annot$`WHO Grade`)] <- "NA"
split <- factor(annot$`WHO Grade`, levels = c("1-2?", "2", "3", "4", "NA"))
col_fun1 <- list("1-2?" = "#fceea7",
                 "2" = "#edd34e",
                 "3" = "#ffac59",
                 "4" = "#ff5959",
                 "NA" = "lightgray",
                 "Male" = "navy",
                 "Female" = "deeppink4",
                 "[0,15]" = "gold",
                 "(15,40]" = "purple",
                 "Primary" = "#cee397",
                 "Recurrent" = "#363062",
                 "Treatment naive" = "lightgray",
                 "Post-treatment" = "gray50",
                 "Post-mortem" = "black",
                 "Cortical" = "magenta",
                 "Other/Multiple locations/NOS" = "pink",
                 "Midline" = "purple",
                 "Cerebellar" = "navy")
circos.clear()
pdf(file = file.path(output_dir, "hope_cohort_data_availability_clinical_with_WHOGrade_2_age_groups.pdf"), width = 10, height = 10)
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
                   at = c("1-2?", "2", "3", "4", "N/A"),
                   legend_gp = gpar(fill = c("#fceea7", "#edd34e", "#ffac59", "#ff5959", "lightgray")))
lgd_gender = Legend(title = "Sex", 
                    at = c("Male", "Female"), 
                    legend_gp = gpar(fill = c("navy", "deeppink4")))
lgd_age = Legend(title = "Age", 
                 at = c("[0,15]", "(15,40]"), 
                 legend_gp = gpar(fill = c("gold", "purple")))
lgd_dtype = Legend(title = "Diagnosis_type",
                   at = c("Primary", "Recurrent"),
                   legend_gp = gpar(fill = c("#cee397", "#363062")))
lgd_annot = Legend(title = "Annotation",
                   at = c("Treatment naive", "Post-treatment", "Post-mortem"),
                   legend_gp = gpar(fill = c("lightgray", "gray50", "black")))
lgd_tumor_location = Legend(title = "Tumor Location",
                            at = c("Cortical", "Other/Multiple locations/NOS", "Midline", "Cerebellar"),
                            legend_gp = gpar(fill = c("magenta", "pink", "purple", "navy")))
h = dev.size()[2]
circle_size = unit(1, "snpc")
lgd_list = packLegend(lgd_grade, lgd_age, lgd_gender, max_height = unit(0.9*h, "inch"), direction = "horizontal")
lgd_list2 = packLegend(lgd_dtype, lgd_annot, max_height = unit(0.9*h, "inch"), direction = "horizontal")
lgd_list3 = packLegend(lgd_tumor_location, max_height = unit(0.9*h, "inch"), direction = "horizontal")
draw(lgd_list, x = unit(120, "mm"), y = unit(150, "mm")) 
draw(lgd_list2, x = unit(120, "mm"), y = unit(120, "mm")) 
draw(lgd_list3, x = unit(115, "mm"), y = unit(95, "mm")) 
circos.clear()
dev.off()

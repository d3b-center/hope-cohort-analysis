# Function: data availability heatmap with two age groups

# load libraries
suppressPackageStartupMessages({
  library(reshape2)
  library(tidyverse) 
  library(ComplexHeatmap)
  library(circlize) 
})

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analyses_dir <- file.path(root_dir, "analyses", "data-availability")
input_dir <- file.path(analyses_dir, "input")
output_dir <- file.path(analyses_dir, "results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# read histologies
annot <- read_tsv(file.path(data_dir, "Hope-GBM-histologies.tsv"))
annot <- annot %>%
  filter(!is.na(HOPE_diagnosis))

# fix HOPE_diagnosis_type
annot <- annot %>%
  mutate(HOPE_diagnosis_type = case_when(
    HOPE_diagnosis_type %in% c("Recurrent", "recurrent", "Recurrent, residual")  ~ "Recurrence",
    HOPE_diagnosis_type == "Primary" ~ "Initial CNS Tumor",
    TRUE ~ as.character(HOPE_diagnosis_type)))

# fix WHO Grade
annot <- annot %>%
  mutate(HARMONY_WHO.Grade = ifelse(HARMONY_WHO.Grade == "1-2?", NA, HARMONY_WHO.Grade))

# add WHO Grade to Diagnosis 
diagnosis_who_grade_map <- annot %>% 
  filter(!is.na(HOPE_diagnosis)) %>%
  mutate(Diagnosis = gsub(" [(].*", "", HOPE_diagnosis)) %>%
  group_by(Diagnosis, HOPE_diagnosis) %>% 
  summarise(HARMONY_WHO.Grade = toString(as.roman(sort(unique(na.omit(HARMONY_WHO.Grade)))))) %>%
  mutate(HARMONY_WHO.Grade = gsub(", ", "/", HARMONY_WHO.Grade)) %>%
  mutate(HARMONY_WHO.Grade = paste0("(WHO grade ", HARMONY_WHO.Grade, ")")) %>%
  mutate(Diagnosis = paste(Diagnosis, HARMONY_WHO.Grade)) %>%
  dplyr::select(HOPE_diagnosis, Diagnosis)

# add new Diagnosis to annotation 
annot <- annot %>%
  left_join(diagnosis_who_grade_map)

# subset to columns of interest
annot <- annot %>%
  dplyr::rename("Age" = "HARMONY_age_class_derived",
                "Gender" = "HARMONY_Gender",
                "Diagnosis Type" = "HOPE_diagnosis_type",
                "Annotation" = "HOPE_sample_annotation",
                "Tumor Location" = "HOPE_Tumor.Location.condensed") %>%
  dplyr::select(sample_id, Age, Gender, Diagnosis, `Diagnosis Type`, Annotation, `Tumor Location`) %>%
  unique() %>%
  column_to_rownames('sample_id') %>%
  dplyr::arrange(Age, Gender, Diagnosis, `Diagnosis Type`, Annotation, `Tumor Location`)

# plot with two age groups
# merge age into 2 groups
plot_df <- annot %>%
  mutate(Age = as.character(Age)) %>%
  mutate(Age = ifelse(Age %in% c("(15,26]", "(26,40]"), "(15,40]", Age))

# split by age
split <- factor(plot_df$Age, levels = c("[0,15]", "(15,40]"))

# color codes
col_diagnosis <- c("High-grade glioma/astrocytoma (WHO grade III/IV)" = "lightseagreen",
                   "Diffuse Midline Glioma (WHO grade III/IV)" = "darkgreen",
                   "Astrocytoma;Oligoastrocytoma (WHO grade III)" = "mediumorchid2",
                   "Astrocytoma (WHO grade III/IV)" = "#5fff57", 
                   "Glioblastoma (WHO grade IV)" = "#f268d6",
                   "Pleomorphic xanthoastrocytoma (WHO grade II/III)" = "#005082")
col_gender <- c("Male" = "#0707CF",
                "Female" = "#CC0303")
col_age <- c("[0,15]" = "#C7E9C0",
             "(15,40]" = "#238B45")
col_dtype <- c("Initial CNS Tumor" = "#cee397",
               "Progressive" = "#827397",
               "Recurrence" = "#363062",
               "Second Malignancy" = "#005082")
col_annot <- c("Treatment naive" = "lightgray",
               "Post-treatment" = "gray50",
               "Post-mortem" = "black")
col_tumor_loc <- c("Cortical" = "#D4806C",
                   "Other/Multiple locations/NOS" = "#7C8F97",
                   "Midline" = "#344C68",
                   "Cerebellar" = "#94004C")
col_fun1 = list(c(col_diagnosis, col_gender, col_age, col_dtype, col_annot, col_tumor_loc))[[1]] %>% 
  unlist() %>% as.list()

# generate plot
circos.clear()
pdf(file = file.path(output_dir, "hope_clinical_data_availability_age_two_groups.pdf"), width = 10, height = 10)
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
                       at = names(col_diagnosis), 
                       legend_gp = gpar(fill = col_diagnosis))
lgd_gender = Legend(title = "Sex", 
                    at = names(col_gender), 
                    legend_gp = gpar(fill = col_gender))
lgd_age = Legend(title = "Age", 
                 at = names(col_age), 
                 legend_gp = gpar(fill = col_age))
lgd_dtype = Legend(title = "Diagnosis Type",
                   at = names(col_dtype),
                   legend_gp = gpar(fill = col_dtype))
lgd_annot = Legend(title = "Annotation",
                   at = names(col_annot),
                   legend_gp = gpar(fill = col_annot))
lgd_tumor_location = Legend(title = "Tumor Location",
                            at = names(col_tumor_loc),
                            legend_gp = gpar(fill = col_tumor_loc))
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

# Function: data availability heatmap with continous age variable

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
input_dir <- file.path(analyses_dir, "input")
output_dir <- file.path(analyses_dir, "results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data", "v1")
analyses_dir <- file.path(root_dir, "analyses", "data-availability")
output_dir <- file.path(analyses_dir, "results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# read histologies
annot <- read_tsv(file.path(data_dir, "Hope-GBM-histologies.tsv"))
annot <- annot %>%
  filter(!is.na(HOPE_diagnosis))

# filter to HOPE annotation binary matrix
binary_matrix <- read_tsv(file.path(input_dir, "compare_HOPE_v2plot.annotation.txt"))
binary_matrix <- binary_matrix %>%
  filter(Remove == 0)
annot <- annot %>%
  filter(sample_id %in% binary_matrix$sample_id)

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
                "Age_at_Initial_Diagnosis" = "HOPE_Age_at_Initial_Diagnosis",
                "Gender" = "HARMONY_Gender",
                "Diagnosis Type" = "HOPE_diagnosis_type",
                "Annotation" = "HOPE_sample_annotation",
                "Tumor Location" = "HOPE_Tumor.Location.condensed") %>%
  dplyr::mutate(Age_at_Initial_Diagnosis = Age_at_Initial_Diagnosis/365) %>%
  dplyr::select(sample_id, Age, Age_at_Initial_Diagnosis, Gender, Diagnosis, `Diagnosis Type`, Annotation, `Tumor Location`) %>%
  unique() %>%
  column_to_rownames('sample_id') %>%
  dplyr::arrange(Age, Age_at_Initial_Diagnosis, Gender, Diagnosis, `Diagnosis Type`, Annotation, `Tumor Location`)

# merge age into 2 groups (for separation)
plot_df <- annot %>%
  mutate(Age = as.character(Age)) %>%
  mutate(Age = ifelse(Age %in% c("(15,26]", "(26,40]"), "(15,40]", Age))

# variable for splitting the circular heatmap
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

# color mapping for numeric age variable
col_fun_age <- colorRamp2(breaks = c(min(plot_df$Age_at_Initial_Diagnosis, na.rm = T), median(plot_df$Age_at_Initial_Diagnosis, na.rm = T), max(plot_df$Age_at_Initial_Diagnosis, na.rm = T)),
                          colors = c("lightskyblue1","skyblue","dodgerblue4"))

# generate plot 
circos.clear()
pdf(file = file.path(output_dir, "hope_clinical_data_availability_age_continuous.pdf"), width = 10, height = 10)
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
                       at = names(col_diagnosis), 
                       legend_gp = gpar(fill = col_diagnosis))
lgd_gender = Legend(title = "Sex", 
                    at = names(col_gender), 
                    legend_gp = gpar(fill = col_gender))
lgd_age = Legend(title = "Age (years)", col_fun = col_fun_age, at = c(0, 20, 39.9), labels = c(0, 20, 39.9))
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

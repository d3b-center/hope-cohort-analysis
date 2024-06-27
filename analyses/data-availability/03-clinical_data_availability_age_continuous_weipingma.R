# Function: data availability heatmap with continous age variable
# This version has been updated by Weiping Ma

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
  filter(!is.na(HOPE_diagnosis)) %>%
  filter(sample_id != "7316-3625")

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
# col_fun_age <- colorRamp2(breaks = c(min(plot_df$Age_at_Initial_Diagnosis, na.rm = T), median(plot_df$Age_at_Initial_Diagnosis, na.rm = T), max(plot_df$Age_at_Initial_Diagnosis, na.rm = T)),
#                           colors = c("lightskyblue1","skyblue","dodgerblue4"))

col_fun_age <- colorRamp2(breaks = summary(plot_df$Age_at_Initial_Diagnosis),
                          colors = c("#146CF6", "#188AF0", "#00B7D8", "#00D4B0", "#00E54B", "#00F800"))

###################################### Added by Weiping Ma ###################################### 
### horizontal group name text order names

### add 6 empty slots to the data matrix
k = 6
plot_df2 = rbind(matrix(NA,ncol = ncol(plot_df),nrow = k,dimnames = list(1:k,colnames(plot_df))),
                 plot_df)
plot_df2$Age[is.na(plot_df2$Age)] = '[0]'
col_fun2 = col_fun1
col_fun2$`[0]` = 'white'
# col_fun_age

split2 <- factor(plot_df2$Age, levels = c('[0]',"[0,15]", "(15,40]"))

plot_df2$`Tumor Location`[plot_df2$`Tumor Location` == "Other/Multiple locations/NOS"] = 'Other/Multiple/NOS'

names(col_tumor_loc)[2] = 'Other/Multiple/NOS'

col_tumor_loc

# generate plot 
circos.clear()
pdf(file = file.path(output_dir, "hope_clinical_data_availability_age_continuous_weipingma.pdf"), width = 10, height = 10)
circos.par(start.degree = 90+(360-5*3-89-k)/(92+k)*k/2+0.5+k/2-1, gap.degree = 1, gap.after = c(5), points.overflow.warning = FALSE)
circos.heatmap(mat = plot_df2 %>% dplyr::select(Age), 
               split = split2,
               col = list("[0]" = "white","[0,15]" = "white", "(15,40]" = "white"), 
               track.height = 0.001, 
               cell.border = "white",
               bg.lwd = 2, cell.lwd = 2)

circos.heatmap(plot_df2 %>% 
                 dplyr::select(c('Diagnosis Type','Tumor Location','Annotation','Diagnosis','Gender')), 
               col = unlist(col_fun2),
               track.height = 0.3, 
               cell.border = "white",
               bg.lwd = 2, cell.lwd = 2)

circos.heatmap(mat = plot_df2 %>% dplyr::select(Age_at_Initial_Diagnosis), 
               col = col_fun_age, 
               track.height = 0.06, 
               cell.border = "white",
               bg.lwd = 2, cell.lwd = 2)

highlight.sector("[0]", col = "white")

cex.name = 1.2
# facing.name = 'downward'
facing.name = "bending.outside"
circos.trackText(x = k/2,y = 0.5,
                 labels = 'Age',
                 cex = cex.name, sectors = '[0]',track = 3,facing = facing.name,niceFacing = TRUE)

circos.trackText(x = rep(k/2,5),y = 4:0+0.5,
                 labels = c('Diagnosis Type','Tumor Location','Treatment','Histology','Sex'),
                 cex = cex.name, sectors = rep('[0]',5),track = 2,facing = facing.name,niceFacing = TRUE)

# legends
lgd_diagnosis = Legend(title = "Histology", 
                       at = names(col_diagnosis), 
                       legend_gp = gpar(fill = col_diagnosis))
lgd_gender = Legend(title = "Sex", 
                    at = names(col_gender), 
                    legend_gp = gpar(fill = col_gender))
lgd_age = Legend(title = "Age (years)", col_fun = col_fun_age, at = c(0, 15, 39.9), labels = c(0, 15, 39.9))
lgd_dtype = Legend(title = "Diagnosis Type",
                   at = names(col_dtype),
                   legend_gp = gpar(fill = col_dtype))
lgd_annot = Legend(title = "Treatment",
                   at = names(col_annot),
                   legend_gp = gpar(fill = col_annot))
lgd_tumor_location = Legend(title = "Tumor Location",
                            at = names(col_tumor_loc),
                            legend_gp = gpar(fill = col_tumor_loc))

gap.temp = lgd_diagnosis@grob$vp$width-lgd_tumor_location@grob$vp$width-lgd_dtype@grob$vp$width

h0 = lgd_age@grob$vp$height
h1 = lgd_diagnosis@grob$vp$height
h2 = lgd_tumor_location@grob$vp$height
h3 = lgd_annot@grob$vp$height

h = dev.size()[2]
circle_size = unit(1, "snpc")
lgd_list = packLegend(lgd_age, lgd_gender, max_height = unit(0.9*h, "inch"), direction = "horizontal",
                      column_gap = unit(0.5, "inch"))
lgd_list1 = packLegend(lgd_diagnosis, max_height = unit(0.9*h, "inch"), direction = "horizontal")
lgd_list2 = packLegend(lgd_tumor_location, lgd_dtype, max_height = unit(0.9*h, "inch"), 
                       direction = "horizontal",column_gap = gap.temp)
lgd_list3 = packLegend(lgd_annot, max_height = unit(0.9*h, "inch"), direction = "horizontal")


lgd.x = unit(h/2, "inch")
lgd.gap = unit(5, "mm")
lgd.y3 = (unit(254, "mm")-(lgd.gap*3+h0+h1+h2))/2

draw(lgd_list, x = lgd.x, y = lgd.y3+lgd.gap*3+h0/2+h1+h2+h3/2) 
draw(lgd_list1, x = lgd.x, y = lgd.y3+lgd.gap*2+h1/2+h2+h3/2) 
draw(lgd_list2, x = lgd.x, y = lgd.y3+lgd.gap*1+h2/2+h3/2) 
draw(lgd_list3, x = lgd.x, y = lgd.y3) 

cex.n = 0.8
# x.n = unit(100,'mm')
x.n = -0.05
text(x = x.n, y = 0.39, 'AYA', cex = cex.n)
text(x = x.n, y = 0.35, '(n = 34)', cex = cex.n)

text(x = x.n, y = 0.30, 'PED', cex = cex.n)
text(x = x.n, y = 0.26, '(n = 59)', cex = cex.n)

# text(x = 0, 0.93, 'PED', cex = cex.n)

circos.clear()
dev.off()


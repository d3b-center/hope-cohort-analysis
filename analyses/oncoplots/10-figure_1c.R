# Function: Script to generate Figure 1C
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(tidyverse)
  library(circlize)
})

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analyses_dir <- file.path(root_dir, "analyses", "oncoplots")
input_dir <- file.path(analyses_dir, "results")
output_dir <- file.path(analyses_dir, "results", "cascade_plots_sinai_data")
dir.create(output_dir, recursive = T, showWarnings = F)

# read matrix
dat <- read_tsv(file.path(analyses_dir, "input", "harmonized_mutation_data_090224.txt"))
dat <- dat %>%
  column_to_rownames("geneID")
dat[dat == 0] <- ""
dat[dat != ""] <- "Mutation"
keep = apply(
  dat,
  1,
  FUN = function(x)
    length(unique(x)) == 2
)
dat = dat[keep, ]

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w, h, gp = gpar(fill = "#ffffff", col = "#595959"))
  },
  Mutation = function(x, y, w, h) {
    grid.rect(x,
              y,
              w - unit(0.3, "mm"),
              h - unit(0.3, "mm"),
              gp = gpar(fill = "black", col = NA))
  }
)
col = c("Mutation" = "black")

# read annotation
annot_info <- read.delim(
  file.path(input_dir, "annotation_add_tumor_only.txt"),
  header = TRUE,
  check.names = TRUE
)
annot_info <- annot_info %>%
  filter(Sample %in% colnames(dat), Sequencing_Experiment != "RNA-Seq") %>%
  dplyr::select(-c(CNS_region)) %>%
  remove_rownames() %>%
  column_to_rownames('Sample') %>%
  as.data.frame()
samples_to_use <- intersect(rownames(annot_info), colnames(dat))
dat <- dat[, samples_to_use]
annot_info <- annot_info[samples_to_use, ]

# convert to two age groups
annot_info <- annot_info %>%
  mutate(Age = as.character(Age)) %>%
  mutate(Age = ifelse(Age %in% c("(15,26]", "(26,40]"), "(15,40]", Age))
annot_info$Age <- factor(annot_info$Age, levels = c("[0,15]", "(15,40]"))

# color code for TMB
col_fun_tmb = colorRamp2(c(0, max(annot_info$TMB, na.rm = T)), c("white", "magenta3"))

# remove Molecular_Subtype
annot_info$Molecular_Subtype <- NULL

# replace NA values with "Flag"
annot_info$Cancer_Group[is.na(annot_info$Cancer_Group)] <- "Flag"

# remove underscore from column names
colnames(annot_info) <- gsub("_", " ", colnames(annot_info))

# oncoprint annotation color code
ha = HeatmapAnnotation(
  df = annot_info %>% dplyr::select(-c(`Sequencing Experiment`)),
  col = list(
    TMB = col_fun_tmb,
    Diagnosis = c(
      "High-grade glioma/astrocytoma (WHO grade III/IV)" = "lightseagreen",
      "Diffuse Midline Glioma (WHO grade III/IV)" = "darkgreen",
      "Astrocytoma;Oligoastrocytoma (WHO grade III)" = "mediumorchid2",
      "Astrocytoma (WHO grade III/IV)" = "#5fff57",
      "Glioblastoma (WHO grade IV)" = "#f268d6",
      "Pleomorphic xanthoastrocytoma (WHO grade III)" = "#005082"
    ),
    `Molecular Subtype` = c(
      "DMG, H3 K28" = "#053061",
      "DHG, H3 G35, TP53" = "#A6761D",
      "HGG, H3 wildtype" = "#4393c3",
      "HGG, H3 wildtype, TP53" = "darkgreen",
      "DMG, H3 K28, TP53" = "#BC80BD",
      "HGG, IDH, TP53" = "#FFFF99",
      "IHG, NTRK-altered, TP53"  = "#E7298A",
      "IHG, NTRK-altered" = "#f4a582",
      "IHG, ROS1-altered" = "#d6604d",
      "IHG, ALK-altered" = "#E31A1C",
      "PXA" = "#67001f",
      "HGG, IDH" = "#B3DE69",
      "Flag" = "#f1f1f1"
    ),
    `Diagnosis Type` = c(
      "Initial CNS Tumor" = "#cee397",
      "Progressive" = "#827397",
      "Recurrence" = "#363062",
      "Second Malignancy" = "#005082"
    ),
    `Tumor Location` = c(
      "Cortical" = "#D4806C",
      "Other/Multiple locations/NOS" = "#7C8F97",
      "Midline" = "#344C68",
      "Cerebellar" = "#94004C"
    ),
    `CNS region` = c(
      "Posterior fossa" = "#D4806C",
      "Other" = "#7C8F97",
      "Midline" = "#344C68",
      "Hemispheric" = "#94004C",
      "Mixed" = "darkgreen"
    ),
    `Cancer Group` = c(
      "(DMG) Diffuse Midline Glioma" = "#053061",
      "(DHG) Diffuse Hemispheric Glioma" = "#A6761D",
      "(HGG) High Grade Glioma (not otherwise specified)" = "#4393c3",
      "(IHG) Infantile Hemispheric Glioma" = "#E7298A",
      "(PXA) Pleomorphic Xanthoastrocytoma" = "#67001f",
      "Flag" = "#f1f1f1"
    ),
    Sex = c("Male" = "#0707CF", "Female" = "#CC0303"),
    Age = c("[0,15]" = "#C7E9C0", "(15,40]" = "#238B45")
  ),
  annotation_name_gp = gpar(fontsize = 10),
  gp = gpar(col = "#595959"),
  simple_anno_size = unit(4, "mm"),
  annotation_name_side = "left"
)

# order by Age
ht = oncoPrint(
  dat,
  get_type = function(x)
    strsplit(x, ";")[[1]],
  alter_fun = alter_fun,
  column_split = annot_info$Age,
  col = col,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 10),
  row_names_gp = gpar(fontsize = 10),
  column_names_side = "top",
  top_annotation = ha,
  right_annotation = NULL,
  row_names_side = "left",
  pct_side = "right",
  remove_empty_rows = TRUE,
  heatmap_legend_param = list(
    title = "Alteration",
    nrow = 9,
    title_position = "topleft",
    direction = "horizontal",
    at = c("Mutation"),
    labels = c("SNV")
  )
)
pdf(
  file = file.path(output_dir, "figure_1c.pdf"),
  width = 16,
  height = 8
)
draw(
  ht,
  merge_legend = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()

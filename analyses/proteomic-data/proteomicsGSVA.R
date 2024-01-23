
#load libraries
library(tidyverse)
library(readr)
library(tibble)
library(msigdbr)
library(GSVA)
#


#Magrittr pipe
`%>%` <- dplyr::`%>%`
#


#set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
#data directory
data_dir <- file.path(root_dir, "data")
#directory containing input files
inputfile_dir <- "./output"
#directory to store output
output_dir <- "./GSEAoutput"
#if the output directory does not exist, create it
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
#


#read histology file to data frame
histologyfile <- file.path(data_dir, "Hope-GBM-histologies-base.tsv")
HISTdata <- read.table(histologyfile, sep="\t", quote="\"", header=TRUE)
#


#WCP data input file
quant_data_file <- file.path(inputfile_dir, "hope-protein-imputed-prot-expression.tsv")
if (!file.exists(quant_data_file)){
  stop("input quantification file does not exist")
}
#


#read expression data to data frame
quant_df <- read.table(quant_data_file, sep="\t", quote="\"", header=TRUE, check.names=FALSE)
#get only quant values and convert to matrix
quant_data<-quant_df[,-c(1,2,3)]
#change column names from kids first IDs (BS_....) to sample IDs (7316-...)
#use histology file to get corresponding IDs
quant_data_colnames<-colnames(quant_data)
for (i in 1:length(quant_data_colnames)){
  currcolname<-quant_data_colnames[i]
  currsamplename<-HISTdata$sample_id[HISTdata$Kids_First_Biospecimen_ID==currcolname]
  #if there is no sample ID keep the original column name
  if ( length(currsamplename)==0 ){
    currsamplename<-currcolname
  }
  #change relevant data frame column name
  colnames(quant_data)[which(names(quant_data) == currcolname)] <- currsamplename
}
#convert to matrix for use in GSEA function
quant_data_matrix <- as.matrix(quant_data)
quant_data_matrix <- apply(quant_data_matrix, 2 ,as.numeric)
#add gene symbols as matrix row names
row.names(quant_data_matrix)<-quant_df$ApprovedGeneSymbol
#


#specify annotation category: KEGG, HALLMARK, BIOCARTA
annotation_category<-"KEGG"
#human geneset genes from `migsdbr` package
#KEGG and BIOCARTA have proteasome pathway annotation
#HALLMARK does not have proteasome pathway annotation
if ( annotation_category == "KEGG" ){
  human_genesets  <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
}else if ( annotation_category == "HALLMARK" ){
  human_genesets  <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
}else if ( annotation_category == "BIOCARTA" ){
  human_genesets  <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:BIOCARTA")
}else {
  stop("choose HALLMARK, KEGG, OR BIOCARTA AS ANNOTATION CATEGORY")
}
#
#reformat to geneset name and gene symbol tibble and list
human_genesets_twocols <- human_genesets %>% dplyr::select(gs_name, human_gene_symbol)
human_genesets_genelist <- base::split(human_genesets_twocols$human_gene_symbol, list(human_genesets_twocols$gs_name))
#


#calculate the GSVA scores
gsea_scores <- GSVA::gsva(quant_data_matrix,
                               human_genesets_genelist,
                               method = "gsva",
                               min.sz=1, max.sz=1500, #arguments from K. Rathi
                               parallel.sz = 8, #for the bigger dataset, this ensures this won't crash due to memory problems
                               mx.diff = TRUE)  #setting this argument to TRUE computes Gaussian-distributed scores (bimodal score distribution if FALSE)
#


#reformat
gsea_scores_df <- as.data.frame(gsea_scores) %>%
  rownames_to_column(var = "geneset_name")
#


#output file
output_file <- paste("Hope_proteome_imputed_data_GSVA_", annotation_category, ".tsv", sep="")
gsva_output_file <- file.path(output_dir, output_file)
#


#write output to file
write.table(gsea_scores_df, gsva_output_file, sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
#
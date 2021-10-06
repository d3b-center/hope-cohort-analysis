library(edgeR)

ss_diffexpr <- function(expr, norm_method = c("uq", "tmm"), housekeeping_genes = NULL){
  
  # create groups
  sample_subject <- grep('sample_of_interest', colnames(expr))
  ref_subjects <- grep('sample_of_interest', colnames(expr), invert = TRUE)
  expr <- expr[,c(ref_subjects, sample_subject)]
  group <- factor(c(rep(1, times = length(ref_subjects)), rep(2, times = length(sample_subject))))
  
  # build DGEList
  y <- DGEList(counts = expr, group = group)
  
  # filter out low expression genes
  keep = filterByExpr(y)
  y <- y[keep,,keep.lib.sizes=FALSE]
  
  # normalization
  if(norm_method == "uq"){
    # Perform UQ-pgq2 on y$counts and put normalized data back into DGEList object
    y$counts <- uq_pgQ2(X = y$counts)
  } else if(norm_method == "tmm"){
    # TMM normalization
    y <- calcNormFactors(object = y)
  }
  
  if(!is.null(housekeeping_genes)){
    # copy normalized DGEList object (new object y1) and replace group variable with only one treatment group
    y1 = y
    y1$samples$group <- 1
    
    # housekeeping genes index
    housekeeping_genes <- which(rownames(y1$counts) %in% housekeeping_genes)
    
    # estimate common dispersion from a set of housekeeping genes (linked below).
    y0 <- estimateDisp(y1[housekeeping_genes,], trend="none", tagwise=FALSE)
    y$common.dispersion <- y0$common.dispersion
  } else{
    # set bcv to 0.5
    y$common.dispersion <- 0.5
  }
  
  ### Make design matrix 
  design <- model.matrix(~group)
  fit <- glmFit(y = y, design = design)
  lrt <- glmLRT(fit)
  
  # get table of differentially expressed genes and perform FDR adjustment on p-value. 
  # filter to FDR < 0.05 to get significant genes
  top_tags <- topTags(lrt, adjust.method = "BH", n = Inf)
  top_tags <- top_tags$table %>%
    filter(FDR < 0.05)
  top_tags$diff_expr <- ifelse(top_tags$logFC > 0, 'up', 'down')
  
  # add tpm value to top_tags
  top_tags <- top_tags %>%
    rownames_to_column('genes') %>%
    dplyr::select(genes, logFC, diff_expr) %>%
    column_to_rownames('genes')
  
  return(top_tags) 
}

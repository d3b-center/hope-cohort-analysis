# compute correlations and pvalues for two variables in the clinical file
dat <- readxl::read_xlsx('data/Sequenced RNA additional selectionMK_ with extraction data-2.xlsx', skip = 1, sheet = 2)
dat <- dat[,1:22]
pairs_corr <- list(c("Concentration...19", "sequenicg passed/failed"),
                   c("RIN", "sequenicg passed/failed"),
                   c("260/230", "sequenicg passed/failed"),
                   c("Contamination %", "sequenicg passed/failed"),
                   c("RQS", "sequenicg passed/failed"))

output <- data.frame()
for(i in 1:length(pairs_corr)){
  cols <- pairs_corr[[i]]
  tmp <- dat %>%
    dplyr::select(cols)
  colnames(tmp) <- c("num_var", "cat_var")
  tmp$num_var <- as.numeric(tmp$num_var)
  ks_out <- kruskal.test(num_var ~ cat_var, data = tmp)
  kruskal_wallis_pval <- ks_out$p.value
  df <- data.frame(num_var = cols[1], cat_var = cols[2], kruskal_wallis_pval)
  output <- rbind(output, df)
}

num_pairs_corr <- list(c("260/230", "Concentration...19"),
                   c("Concentration...19", "Contamination %"),
                   c("Reads Aligned in Pairs",	"Contamination %"),
                   c("RQS",	"Contamination %"),
                   c("260/230", "Reads Aligned in Pairs"),
                   c("Concentration...19", "Reads Aligned in Pairs"),
                   c("RIN", "Reads Aligned in Pairs"),
                   c("Concentration...19",	"RQS"),
                   c("RIN", "Concentration...19"))		

num_vars_output <- data.frame()
for(i in 1:length(num_pairs_corr)){
  cols <- num_pairs_corr[[i]]
  tmp <- dat %>%
    dplyr::select(cols)
  colnames(tmp) <- c("num_var1", "num_var2")
  tmp$num_var1 <- as.numeric(tmp$num_var1)
  tmp$num_var2 <- as.numeric(tmp$num_var2)
  corr_out <- cor.test(x = tmp$num_var1, y = tmp$num_var2)
  df <- data.frame(num_var = cols[1], cat_var = cols[2], corr_out$p.value, corr_out$estimate)
  num_vars_output <- rbind(num_vars_output, df)
}

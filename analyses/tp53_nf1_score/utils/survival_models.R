# Functions for conducting survival analyses
#
# C. Savonen and SJ Spielman for ALSF - CCDL
#
# 2019, 2022
#

# Attach this package
library(survminer)

# Magrittr pipe
`%>%` <- dplyr::`%>%`

efs_palette <- c("#E69F01", "#56B4E8", "#009E74", "#F0E441",
                 "#0072B3", "#D55E01", "#CC79A6", "#999998")

survival_analysis <- function(metadata,
                              ind_var,
                              test = "kap.meier",
                              ind_data = NULL,
                              metadata_sample_col = "Kids_First_Biospecimen_ID_DNA",
                              ind_data_sample_col = NULL,
                              days_col = "OS_days",
                              status_col = "OS_status") {
  # Given the overall survival information, and an independent variable will
  # run survival analysis and return a list with 1) the model fit object and
  # 2) a summary table 3) the data used.
  #
  # Args:
  #   metadata: a data.frame that contains columns OS_status and OS_days to use
  #             for the survival model. This also assumes "LIVING" and "DECEASED"
  #             are the two statuses. THis will be converted to a numeric variable
  #             for use with `survival` R functions. Samples with NAs are dropped.
  #   ind_var: a character string noting the name of the independent variable to
  #            test as a predictor for survival.
  #   test: a character string noting which test to use. Supported choices:
  #         "kap.meier", "log.rank", "cox.reg"
  #   ind_data: If the data for the independent variable needed for the test is not
  #         in metadata, add it here. This assumes it will be a data.frame with
  #         either a "Tumor_Sample_Barcode", "Kids_First_Biospecimen_ID", or
  #         "sample" column which to inner_join with metadata by.
  #   ind_data_sample_col : A character string that states the name of the column in
  #                `ind_data` that contains the sample IDs that match to the
  #                metadata.
  #   metadata_sample_col : A character string that states the name of the column in
  #                `metadata` that contains the sample IDs that match to the
  #                ind_data. Default is "Kids_First_Biospecimen_ID".
  #   os_days_col : A character string that states the name of the column in
  #                `metadata` that contains the overall survival days information.
  #                 Default is "OS_days".
  #   os_status_col : A character string that states the name of the column in
  #                `metadata` that contains the overall survival status information
  #                 The data in this column will be converted to numeric if it
  #                 is not already.  Default is "OS_status".
  #
  # Returns:
  # A list with three objects: 1) the original model fit object,
  #                            2) the summary table
  #                            3) the original data.frame
  
  ####################### Check the options supplied ###########################
  # List the tests
  supported_tests <- c("kap.meier", "log.rank", "cox.reg")
  
  # Check that it is a supported test
  if (!(test %in% supported_tests)) {
    stop(paste(
      test, "is not a supported test. Please specify one of the following:",
      paste(supported_tests, collapse = ", ")
    ))
  }
  
  # Don't continue if there's no variable specified
  if (is.null(ind_var)) {
    stop("No variable has been supplied to test with using the `ind_var` argument.
         Stopping.")
  }
  
  # List the columns we need
  needed_cols <- c(status_col, days_col)
  
  # Get logical vector indicating which are in metadata
  found_cols <- (needed_cols %in% colnames(metadata))
  
  # If not all the columns are found, stop
  if (!all(found_cols)) {
    stop(cat(
      "The following column names specified for overall survival information: \n",
      paste(needed_cols[which(found_cols)], collapse = "\n"),
      "\n ...were not found in the specified metadata data.frame.",
      "Check your `status_col` and `days_col` arguments."
    ))
  }
  
  ######################## Set up the survival variables #######################
  # Pull out this data
  status <- metadata %>%
    dplyr::pull({{status_col}})

  # Code status variable as 0 (no event) or 1 (event)
  metadata <- metadata %>%
    dplyr::mutate(
      !!status_col := case_when(
        status %in% c("LIVING", "NO EVENT") ~ FALSE,
        status %in% c("DECEASED", "EVENT") ~ TRUE
      ))
      
      
    #   !!status_col := case_when(
    #     status_col == "OS_status" ~ as.numeric(
    #       factor(status, levels = c("LIVING", "DECEASED"))),
    #     status_col == "EFS_status" ~ as.numeric(
    #       factor(status, levels = c("NO EVENT", "EVENT")))
    # ))
  
  ############################ Set up the ind data #############################
  # If other ind_data has been supplied, attempt to join it to metadata
  if (!is.null(ind_data)) {
    
    # Check that a sample column was specified
    if (is.null(ind_data_sample_col)) {
      stop("A `ind_data` data.frame was specified but the column with the sample information
           was not specified using `ind_data_sample_col`.")
    }
    
    # List the columns we need
    found_cols <- c(
      ind_data_sample_col %in% colnames(ind_data),
      ind_var %in% colnames(ind_data),
      metadata_sample_col %in% colnames(metadata)
    )
    
    # If not all the columns are found, stop
    if (!all(found_cols)) {
      stop(cat(
        "The following column names specified: \n",
        paste(needed_cols[which(found_cols)], collapse = "\n"),
        "\n ...were not found in the specified data.frames.",
        "Check your `metadata_sample_col`, `ind_data_sample_col`, and `ind_var` arguments."
      ))
    }
    
    # Set up the enquosures
    ind_data_sample_col <- enquo(ind_data_sample_col)
    metadata_sample_col <- enquo(metadata_sample_col)
    
    # Set up the sample by thing
    sample_by <- set_names(quo_name(ind_data_sample_col), quo_name(metadata_sample_col))
    
    # Join this ind_data to the metadata
    metadata <- dplyr::inner_join(metadata, ind_data,
                                  by = sample_by
    )
  } else {
    # Set up
    metadata_sample_col <- enquo(metadata_sample_col)
    
    # Look for the independent variable columns in the metadata
    found_cols <- (ind_var %in% colnames(metadata))
    
    # If not all the columns are found, stop
    if (!all(found_cols)) {
      stop(cat(
        "The following column names specified for the independent variable(s): \n",
        ind_var[which(found_cols)],
        "Check your `ind_var` argument."
      ))
    }
  }
  
  ############################ Set up model and data ###########################
  # Extract independent variables
  ind_var_df <- metadata %>%
    dplyr::select(!!metadata_sample_col, {{days_col}}, !!status_col, ind_var)
  # dplyr::mutate(status_col = case_when(
  #    status_col %in% c("LIVING", "NO EVENT") ~ FALSE,
  #    status_col %in% c("DECEASED", "EVENT") ~ TRUE
  #  ))
  
  # For the model need a plus sign for separating multiple independent variables
  ind_var <- paste0(ind_var, collapse = "+")
  
  # Piece together a model
  model <- paste0(
    "survival::Surv(",
    days_col, ", ", status_col, ") ~ ",
    ind_var
  )
  
  ################################# Do the test! ###############################
  # Print out what the model is
  message(paste("Testing model:", model, "with", test))
  
  # Run the appropriate test
  if (test == "kap.meier") {
    # Make the model
    fit <- survival::survfit(
      formula(model),
      data = ind_var_df
    )
  }
  if (test == "log.rank") {
    # Make the model
    fit <- survival::survdiff(
      formula(model),
      data = ind_var_df
    )
    # Obtain p value for Chi-Squared stat
    fit$p.value <- pchisq(fit$chisq, df = 1, lower = FALSE)
  }
  if (test == "cox.reg") {
    # Make the model
    fit <- survival::coxph(
      formula(model),
      data = ind_var_df
    )
  }
  # Tidy up the model object with broom
  table <- broom::tidy(fit)
  
  # Restore the model in this slot so the ggsurvplot function can find it
  fit$call$formula <- formula(model)
  
  # Return both the fit object and the table
  return(list(model = fit, table = table, original_data = ind_var_df))
}




fit_save_model <- function(df,
                           terms,
                           output_file,
                           model_type = "univariate",
                           test = "cox",
                           years_col = "OS_years",
                           status_col = "OS_status") {
  # df: A data frame that contains columns:
  #   - `OS_status` with character values "LIVING" and "DECEASED", or numeric values 0 and 1 respectively
  #   - `OS_years` giving survival time in years
  # terms: A string providing RHS of model equation. For univariate models, this should be a column in `df`. For multivariate models, this should contain the full RHS
  # output_file: Path to RDS file where exported model should be saved
  # model_type: String "univariate" or "multivariate", where default is "univariate"
  # test: String indicating test, generally either "cox" or "log.rank"
  
  # Check model_type
  if (!(model_type %in% c("univariate", "multivariate"))) {
    stop("Must specify `model_type` as 'univariate' or 'multivariate'")
  }
  
  
  # Fit the model
  if (model_type == "univariate") {
    fit <- survival_analysis(
      df,
      ind_var = terms,
      test    = "cox.reg",
      metadata_sample_col = "Kids_First_Biospecimen_ID",
      os_days_col = years_col # we want to use years, not days
    )
      } else if (model_type == "multivariate") {
        
        if (years_col == "EFS_years"){
          
          # Recode OS_Status (for univariate, `survival_analysis()` does the recoding)
          df <- df %>%
            mutate(EFS_status = ifelse(EFS_status == "NO EVENT", 0, 1))
          
          # Fit model
          fitted_multi <- survival::coxph(
            formula(
              paste0("survival::Surv(time = EFS_years, event = EFS_status) ~ ", terms)
            ),
            data = df
          )
          
        } else if (years_col == "OS_years"){
          
          # Recode OS_Status (for univariate, `survival_analysis()` does the recoding)
          df <- df %>%
            mutate(OS_status = ifelse(OS_status == "LIVING", 0, 1))
          
          # Fit model
          fitted_multi <- survival::coxph(
            formula(
              paste0("survival::Surv(time = OS_years, event = OS_status) ~ ", terms)
            ),
            data = df
          )
          
        }

    # Set up list object to match parts of `survival_analysis()` output we need
    fit <- list(
      model = fitted_multi,
      table = broom::tidy(fitted_multi)
    )
  }
  
  # Save model RDS
  readr::write_rds(fit$model,
                   # remove spaces in case (common for cancer groups)
                   str_replace_all(output_file, " ", "-"))
  
  # Return the model fit table to be printed out
  fit$table
}






plotKM <- function(model,
                   variable,
                   combined = FALSE, 
                   title) {
  
  if (combined == FALSE){
    
    if ("OS_days" %in% names(model$original_data)){
      
      event_type <- "OS"
      
      diff_obj <- survdiff(survival::Surv(OS_days, OS_status) ~ cpgPLP_status,  
                              model$original_data)
      diff_pvalue <- 1 - pchisq(diff_obj$chisq, length(diff_obj$n) - 1)
      diff_pvalue_formatted <- format(
        signif(diff_pvalue, 2),
        scientific = FALSE)
      
    }
    
    if ("EFS_days" %in% names(model$original_data)){
      
      event_type <- "EFS"
      
      diff_obj <- survdiff(survival::Surv(EFS_days, EFS_status) ~ cpgPLP_status,  
                               model$original_data)
      diff_pvalue <- 1 - pchisq(diff_obj$chisq, length(diff_obj$n) - 1)
      diff_pvalue_formatted <- format(
        signif(diff_pvalue, 2),
        scientific = FALSE)
      
    }
    
    levels <- model$original_data %>%
      pull(variable) %>%
      unique()
    levels <- levels[order(levels)]
    
    colors <- colorblindr::palette_OkabeIto[1:length(levels)]
    lines <- c(rep("solid", length(levels)), 
              rep("dashed", length(levels)))
    labels <- glue::glue("{event_type}:{levels}")
    
      
    km_plot <- survminer::ggsurvplot(fit = model$model, 
                                     data = model$original_data,
                                     palette = colors,
                                     linetype = lines,
                                     risk.table = TRUE,
                                     xlim = c(0, 4000),
                                     break.time.by = 500,
                                     ggtheme = ggpubr::theme_pubr(),
                                     title = title,
                                     legend = "right",
                                     xlab = "Time (days)",
                                     legend.title = '',
                                     legend.labs = labels,
                                     risk.table.y.text.col = TRUE,
                                     risk.table.y.text = TRUE) 
    
    km_plot_graph <- km_plot$plot + 
      ggplot2::annotate("text", 
                        200, 0.15, 
                        label = paste0(event_type, " P = ", diff_pvalue_formatted)) +
      theme(legend.text = element_text(size = 16, color = "black", face = "bold")) +
      cowplot::background_grid()
    
    km_plot_table <- km_plot$table
    
    # Re-combine plot and table
    km_final <- km_plot_graph/km_plot_table +
      plot_layout(widths = c(2, 1), 
                  heights = c(2, 0.75))
    
    return(km_final)
    
  } 
  
  if (combined == TRUE){
    
    fit <- list(OS = model[[1]]$model,
                EFS = model[[2]]$model)
    
    data_os <- model[[1]]$original_data
    
    data_efs <- model[[2]]$original_data
    
    variable_os <- data_os %>%
      pull(variable)
    
    levels_os <- unique(variable_os[!is.na(data_os$OS_days)][order(variable_os[!is.na(data_os$OS_days)])])
    levels_os <- levels_os[!is.na(levels_os)]
    
    os_palette <- colorblindr::palette_OkabeIto[1:length(levels_os)]
    
    variable_efs <- data_efs %>%
      pull(variable)
    
    levels_efs <- unique(variable_efs[!is.na(data_efs$EFS_days)][order(variable_efs[!is.na(data_efs$EFS_days)])])
    levels_efs <- levels_efs[!is.na(levels_efs)]
    
    colors <- c(os_palette[1:length(levels_os)],
                efs_palette[1:length(levels_efs)])
    
    lines <- c(rep("solid", length(levels_os)), 
               rep("dashed", length(levels_efs)))
    
    labels <- c(glue::glue("OS:{levels_os}"),
                glue::glue("EFS:{levels_efs}"))
    
    diff_os_obj <- survdiff(survival::Surv(OS_days, OS_status) ~ variable_os,  
                            data_os)
    diff_os_pvalue <- 1 - pchisq(diff_os_obj$chisq, length(diff_os_obj$n) - 1)
    diff_os_pvalue_formatted <- format(
      signif(diff_os_pvalue, 2),
      scientific = FALSE)
    
    diff_efs_obj <- survdiff(survival::Surv(EFS_days, EFS_status) ~ variable_efs,  
                             data_efs)
    diff_efs_pvalue <- 1 - pchisq(diff_efs_obj$chisq, length(diff_efs_obj$n) - 1)
    diff_efs_pvalue_formatted <- format(
      signif(diff_efs_pvalue, 2),
      scientific = FALSE)
    
    km_plot <- survminer::ggsurvplot(fit = fit, 
                                     data = data_efs,
                                     palette = colors,
                                     linetype = lines,
                                     combine = TRUE,
                                     risk.table = TRUE,
                                     xlim = c(0, 4000),
                                     break.time.by = 500,
                                     ggtheme = ggpubr::theme_pubr(),
                                     title = title,
                                     legend = "right",
                                     xlab = "Time (days)",
                                     legend.title = '',
                                     legend.labs = labels,
                                     risk.table.y.text.col = TRUE,
                                     risk.table.y.text = TRUE) 
    
    km_plot_graph <- km_plot$plot + 
      ggplot2::annotate("text", 
                        300, 0.15, 
                        label = paste0("OS P = ", diff_os_pvalue_formatted)) +
      ggplot2::annotate("text", 
                        300, 0.10, 
                        label = paste0("EFS P = ", diff_efs_pvalue_formatted)) +
      theme(legend.key.size = unit(1, 'cm')) +
      cowplot::background_grid()
    
    km_plot_table <- km_plot$table
    
    # Re-combine plot and table
    km_final <- km_plot_graph/km_plot_table +
      plot_layout(widths = c(2, 1), 
                  heights = c(2, 0.75))
    
    return(km_final)
    
  }
}





plotForest <- function(model) {
  
  # Determine if OS or EFS model 
  
  event_type <- ifelse(grepl("OS", model$formula[2]), 
                       "OS",
                       "EFS")
  
  # Set up ordering and labels for y-axis
  term_order <- rev(paste0(unlist(lapply(names(model$xlevels), function(x) rep(x, length(model$xlevels[[x]])))),
                           as.vector(unlist(model$xlevels))))
  
  term_labels <- rev(as.vector(unlist(model$xlevels)))
  
  numeric_terms <- names(model$coefficients)[!names(model$coefficients) %in% term_order]
  
  term_order <- c(numeric_terms, term_order)
  
  term_labels <- c(numeric_terms, term_labels)
  
  survival_n <- broom::glance(model) %>%
    select(n, nevent)
  
  # Convert survival model result to data frame, and exponentiate estimates/CIs to get HRs
  survival_df <- broom::tidy(model) %>%
    
    # Add references
    add_row(term = term_order[!term_order %in% broom::tidy(model)$term], 
            estimate = 0) %>%
    mutate(
      conf.low = exp(estimate-std.error),
      conf.high = exp(estimate+std.error),
      estimate = exp(estimate),
      # significance indicator column for filling points.
      # Note T/F these are strings for type compatibility with "REF"
      significant = case_when(p.value <= 0.05 ~ "TRUE", 
                              p.value > 0.05 ~ "FALSE", 
                              is.na(p.value) ~ "REF"),
      # y-axis factor re-labeling
      term = factor(term, 
                    levels = term_order,
                    labels = term_labels)
    ) %>%
    filter(estimate > 1e-5 & estimate < 5000)
  
  forest_plot <- ggplot(survival_df) +
    aes(x = estimate, y = term, fill = significant
    ) + 
    # add CI first so line doesn't cover open point
    geom_errorbarh(
      aes(xmin = conf.low,xmax = conf.high,
      ), height = 0.15, linewidth = 0.65) + 
    geom_point(size = 3.5, shape = 23) +
    # Point fill based on sigificance
    scale_fill_manual(
      values = c("FALSE" = "white", 
                 "TRUE" = "black",
                 "REF" = "gray"),
      guide = FALSE # turn off legend
    ) + 
    # Vertical guiding line at 1
    geom_vline(xintercept = 1, linetype = 3
    ) +
    labs(x = "Hazard Ratio ± 95% CI", y = "",
         subtitle = glue::glue("{event_type}: N = {survival_n$n} with {survival_n$nevent} events")
    ) + 
    # log-scale the x-axis
    scale_x_log10() +
    ggpubr::theme_pubr() + 
    theme(
      plot.subtitle = element_text(face = "bold")
    ) +
    # grid makes it easier to follow lines
    cowplot::background_grid()
  
  # Accompanying panel with sample sizes, P-values, etc.
  
  # prepare data for panel
  # note this warning is OK and EXPECTED because there is no CI for the reference group: 
  #    Removed 2 rows containing missing values (geom_text). 
  survival_df_spread <- survival_df %>%
    mutate(
      # Clean pvalues into labels. 
      p_string = if_else(
        p.value >= 0.001, 
        paste0("P = ",round(p.value,3)),
        "P < 0.001"
      ),
      # round to 2 digits and create single string with "hr (low-high)"
      conf.low = signif(conf.low, 2),
      conf.high = signif(conf.high, 2),
      estimate = signif(estimate, 2),
      hr_ci = glue::glue("{estimate} ({conf.low} - {conf.high})")
    ) %>%
    select(term, hr_ci, p_string) %>%
    # this throws a warning but it's ok
    # format tibble for plotting
    gather(hr_ci:p_string, key = "name", value = "value") %>%
    #remove values for reference
    mutate(value = value)
  
  labels_panel <- ggplot(survival_df_spread) +
    aes(x = name, y = term, label = value) + 
    geom_text(hjust = 0, size = 3) +
    labs(
      # hack!
      subtitle = paste0("               ",
                        "HR (95% CI)        P-value")
    ) +
    ggpubr::theme_pubr() + 
    # remove axes.
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      # -26 is as low as we can go before plot starts to get coverd
      plot.margin = margin(6, 0, 36, -25, unit = "pt"),
      plot.subtitle = element_text(face = "bold")
    ) 
  
  forest_panels <- cowplot::plot_grid(forest_plot, labels_panel, nrow = 1, rel_widths = c(1,0.5), scale = 0.95)
  
  print(forest_panels)
}

  
  
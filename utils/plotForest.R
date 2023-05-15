# adapted from https://github.com/d3b-center/pbta-splicing/blob/main/analyses/survival/util/survival_models.R
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

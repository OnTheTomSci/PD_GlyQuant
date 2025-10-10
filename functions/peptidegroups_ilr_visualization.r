# ILR Coordinate Visualization Functions
# Functions for visualizing ILR analysis results across multiple proteins and coordinates

#' Create heatmap of ILR coordinate significance
#'
#' @param ilr_results List of ILR results from analyze_multiple_proteins_ilr()
#' @param value_type What to display: "p_value", "effect", or "both"
#' @param p_threshold P-value threshold for significance (default: 0.05)
#' @return ggplot2 heatmap object
#' @export
plot_ilr_heatmap <- function(ilr_results, 
                             value_type = "both",
                             p_threshold = 0.05) {
  
  # Compile all results into one data frame
  all_tests <- bind_rows(lapply(names(ilr_results), function(protein_id) {
    result <- ilr_results[[protein_id]]
    if (!is.null(result$statistical_tests) && nrow(result$statistical_tests) > 0) {
      result$statistical_tests %>%
        mutate(protein = protein_id)
    } else {
      NULL
    }
  }))
  
  if (nrow(all_tests) == 0) {
    stop("No statistical test results found in ilr_results")
  }
  
  # Prepare data for heatmap
  heatmap_data <- all_tests %>%
    mutate(
      # Format for display
      p_display = sprintf("%.3f", p_value_adj),
      effect_display = sprintf("%.2f", difference),
      significant = p_value_adj < p_threshold,
      # Create combined label
      cell_label = case_when(
        value_type == "p_value" ~ p_display,
        value_type == "effect" ~ effect_display,
        value_type == "both" ~ paste0(effect_display, "\n(p=", p_display, ")"),
        TRUE ~ ""
      ),
      # Extract numeric part of ILR coordinate for proper ordering
      ilr_num = as.numeric(stringr::str_extract(ilr_coordinate, "\\d+"))
    ) %>%
    # Order ILR coordinates numerically
    mutate(ilr_coordinate = factor(ilr_coordinate, 
                                   levels = unique(ilr_coordinate[order(ilr_num)])))
  
  # Create heatmap
  if (value_type == "p_value") {
    # Color by p-value
    plot <- ggplot2::ggplot(heatmap_data, 
                            ggplot2::aes(x = ilr_coordinate, y = protein, fill = -log10(p_value_adj))) +
      ggplot2::geom_tile(color = "white", size = 0.5) +
      ggplot2::geom_text(ggplot2::aes(label = ifelse(significant, "*", "")), 
                         size = 6, color = "black") +
      ggplot2::scale_fill_gradient2(
        low = "lightgray", mid = "yellow", high = "red",
        midpoint = -log10(0.05),
        name = "-log10(p-adj)",
        limits = c(0, max(-log10(heatmap_data$p_value_adj), 2))
      )
  } else if (value_type == "effect") {
    # Color by effect size
    max_abs_effect <- max(abs(heatmap_data$difference))
    plot <- ggplot2::ggplot(heatmap_data, 
                            ggplot2::aes(x = ilr_coordinate, y = protein, fill = difference)) +
      ggplot2::geom_tile(color = "white", size = 0.5) +
      ggplot2::geom_text(ggplot2::aes(label = ifelse(significant, "*", "")), 
                         size = 6, color = "black") +
      ggplot2::scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B",
        midpoint = 0,
        name = "Effect Size\n(MECFS - Healthy)",
        limits = c(-max_abs_effect, max_abs_effect)
      )
  } else {
    # Combined: color by effect, asterisk for significance
    max_abs_effect <- max(abs(heatmap_data$difference))
    plot <- ggplot2::ggplot(heatmap_data, 
                            ggplot2::aes(x = ilr_coordinate, y = protein, fill = difference)) +
      ggplot2::geom_tile(color = "white", size = 0.5) +
      ggplot2::geom_point(data = filter(heatmap_data, significant),
                          shape = 8, size = 4, color = "black") +
      ggplot2::scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B",
        midpoint = 0,
        name = "Effect Size\n(MECFS - Healthy)",
        limits = c(-max_abs_effect, max_abs_effect)
      )
  }
  
  plot <- plot +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right"
    ) +
    ggplot2::labs(
      title = "ILR Coordinate Analysis Results",
      subtitle = "* indicates p-adj < 0.05",
      x = "ILR Coordinate",
      y = "Protein"
    )
  
  return(plot)
}


#' Create forest plot of ILR coordinate effects
#'
#' @param ilr_results List of ILR results
#' @param show_only_significant Show only significant results (default: FALSE)
#' @param p_threshold P-value threshold (default: 0.05)
#' @return ggplot2 forest plot object
#' @export
plot_ilr_forest <- function(ilr_results,
                            show_only_significant = FALSE,
                            p_threshold = 0.05) {
  
  # Compile all results
  all_tests <- bind_rows(lapply(names(ilr_results), function(protein_id) {
    result <- ilr_results[[protein_id]]
    if (!is.null(result$statistical_tests) && nrow(result$statistical_tests) > 0) {
      result$statistical_tests %>%
        mutate(
          protein = protein_id,
          protein_coord = paste(protein, ilr_coordinate, sep = " - ")
        )
    } else {
      NULL
    }
  }))
  
  if (nrow(all_tests) == 0) {
    stop("No statistical test results found")
  }
  
  # Add significance
  all_tests <- all_tests %>%
    mutate(
      significant = p_value_adj < p_threshold,
      sig_label = ifelse(significant, "Significant", "Not Significant")
    )
  
  # Filter if requested
  if (show_only_significant) {
    all_tests <- all_tests %>% filter(significant)
    
    if (nrow(all_tests) == 0) {
      warning("No significant results to plot")
      return(NULL)
    }
  }
  
  # Sort by effect size
  all_tests <- all_tests %>%
    arrange(desc(abs(difference))) %>%
    mutate(protein_coord = factor(protein_coord, levels = protein_coord))
  
  # Create forest plot
  plot <- ggplot2::ggplot(all_tests, 
                          ggplot2::aes(x = difference, y = protein_coord, color = sig_label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_lower, xmax = ci_upper), 
                            height = 0.3, alpha = 0.7) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(values = c("Significant" = "#E49CB1", 
                                           "Not Significant" = "gray60")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8),
      legend.position = "top",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    ) +
    ggplot2::labs(
      title = "ILR Coordinate Effect Sizes with 95% Confidence Intervals",
      x = "Effect Size (MECFS - Healthy)",
      y = "Protein - ILR Coordinate",
      color = NULL
    )
  
  return(plot)
}


#' Create dot plot showing significant coordinates
#'
#' @param ilr_results List of ILR results
#' @param p_threshold P-value threshold (default: 0.05)
#' @return ggplot2 dot plot object
#' @export
plot_ilr_dotplot <- function(ilr_results, p_threshold = 0.05) {
  
  # Compile results
  all_tests <- bind_rows(lapply(names(ilr_results), function(protein_id) {
    result <- ilr_results[[protein_id]]
    if (!is.null(result$statistical_tests) && nrow(result$statistical_tests) > 0) {
      result$statistical_tests %>%
        mutate(protein = protein_id)
    } else {
      NULL
    }
  }))
  
  if (nrow(all_tests) == 0) {
    stop("No results to plot")
  }
  
  # Add significance and direction
  all_tests <- all_tests %>%
    mutate(
      significant = p_value_adj < p_threshold,
      direction = case_when(
        !significant ~ "Not Significant",
        difference > 0 ~ "Increased in MECFS",
        difference < 0 ~ "Decreased in MECFS"
      ),
      abs_difference = abs(difference),
      neg_log_p = -log10(p_value_adj),
      # Extract numeric part of ILR coordinate for proper ordering
      ilr_num = as.numeric(stringr::str_extract(ilr_coordinate, "\\d+"))
    ) %>%
    # Order ILR coordinates numerically
    mutate(ilr_coordinate = factor(ilr_coordinate, 
                                   levels = unique(ilr_coordinate[order(ilr_num)])))
  
  # Create dot plot
  plot <- ggplot2::ggplot(all_tests, 
                          ggplot2::aes(x = ilr_coordinate, y = protein)) +
    ggplot2::geom_point(ggplot2::aes(size = abs_difference, 
                                     color = direction,
                                     alpha = significant)) +
    ggplot2::scale_size_continuous(
      name = "Abs(Effect Size)",
      range = c(2, 10)
    ) +
    ggplot2::scale_color_manual(
      values = c("Increased in MECFS" = "#E49CB1",
                 "Decreased in MECFS" = "#9DD4CC",
                 "Not Significant" = "gray80"),
      name = "Direction"
    ) +
    ggplot2::scale_alpha_manual(
      values = c("TRUE" = 1, "FALSE" = 0.3),
      guide = "none"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    ) +
    ggplot2::labs(
      title = "ILR Coordinate Analysis - Dot Plot",
      subtitle = paste("Significant results (p <", p_threshold, ") shown with full opacity"),
      x = "ILR Coordinate",
      y = "Protein"
    )
  
  return(plot)
}


#' Create multi-panel boxplots for significant coordinates
#'
#' @param ilr_results List of ILR results
#' @param p_threshold P-value threshold (default: 0.05)
#' @param max_plots Maximum number of plots to create (default: 20)
#' @return List of ggplot2 objects
#' @export
plot_ilr_significant_boxplots <- function(ilr_results, 
                                          p_threshold = 0.05,
                                          max_plots = 20) {
  
  # Find all significant protein-coordinate combinations
  significant_combos <- list()
  
  for (protein_id in names(ilr_results)) {
    result <- ilr_results[[protein_id]]
    
    if (!is.null(result$statistical_tests)) {
      sig_coords <- result$statistical_tests %>%
        filter(p_value_adj < p_threshold)
      
      if (nrow(sig_coords) > 0) {
        for (i in 1:nrow(sig_coords)) {
          combo_name <- paste(protein_id, sig_coords$ilr_coordinate[i], sep = "_")
          significant_combos[[combo_name]] <- list(
            protein = protein_id,
            coordinate = sig_coords$ilr_coordinate[i],
            p_value = sig_coords$p_value_adj[i],
            effect = sig_coords$difference[i]
          )
        }
      }
    }
  }
  
  if (length(significant_combos) == 0) {
    cat("No significant results found to plot\n")
    return(NULL)
  }
  
  # Sort by p-value and limit
  significant_combos <- significant_combos[order(sapply(significant_combos, function(x) x$p_value))]
  significant_combos <- head(significant_combos, max_plots)
  
  cat(sprintf("Creating boxplots for %d significant protein-coordinate combinations\n", 
              length(significant_combos)))
  
  # Create plots
  plots <- list()
  
  for (combo_name in names(significant_combos)) {
    combo <- significant_combos[[combo_name]]
    result <- ilr_results[[combo$protein]]
    
    # Create boxplot
    plot <- ggplot2::ggplot(result$ilr_transformed,
                            ggplot2::aes(x = group, 
                                         y = .data[[combo$coordinate]], 
                                         fill = group)) +
      ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = 1) +
      ggplot2::geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
      ggplot2::scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 11, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9)
      ) +
      ggplot2::labs(
        title = sprintf("%s - %s", combo$protein, combo$coordinate),
        subtitle = sprintf("p-adj = %.3e, effect = %.2f", combo$p_value, combo$effect),
        x = "Group",
        y = "ILR Value"
      )
    
    plots[[combo_name]] <- plot
  }
  
  return(plots)
}


#' Create comprehensive summary plot combining multiple views
#'
#' @param ilr_results List of ILR results
#' @param output_file Optional filename to save combined plot
#' @param p_threshold P-value threshold (default: 0.05)
#' @return Combined ggplot object (using patchwork)
#' @export
plot_ilr_summary <- function(ilr_results, 
                             output_file = NULL,
                             p_threshold = 0.05) {
  
  # Requires patchwork for combining plots
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required. Install with: install.packages('patchwork')")
  }
  
  # Create individual plots
  cat("Creating heatmap...\n")
  heatmap <- plot_ilr_heatmap(ilr_results, value_type = "effect", p_threshold = p_threshold)
  
  cat("Creating forest plot...\n")
  forest <- plot_ilr_forest(ilr_results, show_only_significant = TRUE, p_threshold = p_threshold)
  
  cat("Creating dot plot...\n")
  dotplot <- plot_ilr_dotplot(ilr_results, p_threshold = p_threshold)
  
  # Combine plots
  if (!is.null(forest) && nrow(forest$data) > 0) {
    combined_plot <- (heatmap / forest) | dotplot
    combined_plot <- combined_plot + 
      patchwork::plot_annotation(
        title = "ILR Compositional Analysis - Summary",
        theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5))
      )
  } else {
    combined_plot <- heatmap / dotplot
    combined_plot <- combined_plot + 
      patchwork::plot_annotation(
        title = "ILR Compositional Analysis - Summary",
        subtitle = "No significant results for forest plot",
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
          plot.subtitle = ggplot2::element_text(hjust = 0.5)
        )
      )
  }
  
  # Save if requested
  if (!is.null(output_file)) {
    cat(sprintf("Saving combined plot to %s...\n", output_file))
    ggplot2::ggsave(output_file, combined_plot, width = 16, height = 12, dpi = 300)
  }
  
  return(combined_plot)
}


#' Create interactive table of ILR results
#'
#' @param ilr_results List of ILR results
#' @param output_file Optional filename to save table
#' @return Data frame with all results
#' @export
create_ilr_results_table <- function(ilr_results, output_file = NULL) {
  
  # Compile all results
  results_table <- bind_rows(lapply(names(ilr_results), function(protein_id) {
    result <- ilr_results[[protein_id]]
    
    if (!is.null(result$statistical_tests) && nrow(result$statistical_tests) > 0) {
      result$statistical_tests %>%
        mutate(
          protein = protein_id,
          n_compositions = result$transformation_info$n_compositions,
          n_samples = result$transformation_info$n_samples
        ) %>%
        select(protein, ilr_coordinate, n_compositions, n_samples,
               mean_healthy, mean_mecfs, difference, 
               t_statistic, p_value, p_value_adj, significant, ci_lower, ci_upper)
    } else {
      NULL
    }
  }))
  
  if (!is.null(output_file)) {
    write.csv(results_table, output_file, row.names = FALSE)
    cat(sprintf("Results table saved to %s\n", output_file))
  }
  
  return(results_table)
}



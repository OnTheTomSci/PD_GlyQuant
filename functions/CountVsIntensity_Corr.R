

#' Correlate relative abundances calculated from PSM counts vs intensities
#' @param psm_relative_abundance_df The relative abundance dataframe from PSM counts
#' @param intensity_relative_abundance_df The relative abundance dataframe from intensities
#' @param feature_col The column name containing the feature categories
#' @param output_file The output CSV file path for correlation results
#' @param plot_file The output plot file path for correlation plots
#' @return A dataframe with correlation results
correlate_psm_intensity_abundances <- function(psm_relative_abundance_df, 
                                              intensity_relative_abundance_df, 
                                              feature_col, 
                                              output_file, 
                                              plot_file) {
  # Load required libraries
  library(ggplot2)
  library(viridis)
  library(dplyr)
  
  # Validate that the required columns exist in both dataframes
  required_cols <- c("sample", "disease_status", feature_col, "relative_percentage")
  missing_cols_psm <- setdiff(required_cols, names(psm_relative_abundance_df))
  missing_cols_intensity <- setdiff(required_cols, names(intensity_relative_abundance_df))
  
  if (length(missing_cols_psm) > 0) {
    stop(paste("Missing required columns in PSM data:", paste(missing_cols_psm, collapse = ", ")))
  }
  if (length(missing_cols_intensity) > 0) {
    stop(paste("Missing required columns in intensity data:", paste(missing_cols_intensity, collapse = ", ")))
  }
  
  # Get unique feature categories
  unique_features <- unique(psm_relative_abundance_df[[feature_col]])
  
  # Initialize results dataframe
  results <- data.frame(
    feature = character(),
    pearson_correlation = numeric(),
    spearman_correlation = numeric(),
    pearson_p_value = numeric(),
    spearman_p_value = numeric(),
    n_samples = numeric(),
    mean_psm_abundance = numeric(),
    mean_intensity_abundance = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Create a list to store plots
  plot_list <- list()
  
  # For each feature category, calculate correlation
  for (feature_value in unique_features) {
    # Filter data for this feature
    psm_data <- psm_relative_abundance_df %>% 
      dplyr::filter(!!rlang::sym(feature_col) == feature_value) %>%
      dplyr::select(sample, disease_status, relative_percentage) %>%
      dplyr::rename(psm_relative_percentage = relative_percentage)
    
    intensity_data <- intensity_relative_abundance_df %>% 
      dplyr::filter(!!rlang::sym(feature_col) == feature_value) %>%
      dplyr::select(sample, disease_status, relative_percentage) %>%
      dplyr::rename(intensity_relative_percentage = relative_percentage)
    
    # Join the data
    combined_data <- psm_data %>%
      dplyr::inner_join(intensity_data, by = c("sample", "disease_status"))
    
    # Skip if there are not enough samples
    if (nrow(combined_data) < 3) {
      warning(paste("Not enough samples for feature", feature_value, 
                   "(n =", nrow(combined_data), ")"))
      next
    }
    
    # Calculate correlations
    pearson_result <- cor.test(combined_data$psm_relative_percentage, 
                              combined_data$intensity_relative_percentage, 
                              method = "pearson")
    spearman_result <- cor.test(combined_data$psm_relative_percentage, 
                               combined_data$intensity_relative_percentage, 
                               method = "spearman")
    
    # Calculate means
    mean_psm <- mean(combined_data$psm_relative_percentage, na.rm = TRUE)
    mean_intensity <- mean(combined_data$intensity_relative_percentage, na.rm = TRUE)
    
    # Add results to dataframe
    results <- rbind(results, data.frame(
      feature = feature_value,
      pearson_correlation = pearson_result$estimate,
      spearman_correlation = spearman_result$estimate,
      pearson_p_value = pearson_result$p.value,
      spearman_p_value = spearman_result$p.value,
      n_samples = nrow(combined_data),
      mean_psm_abundance = mean_psm,
      mean_intensity_abundance = mean_intensity,
      stringsAsFactors = FALSE
    ))
    
    # Create correlation plot
    p <- ggplot(combined_data, 
                aes(x = psm_relative_percentage, 
                    y = intensity_relative_percentage, 
                    color = disease_status)) +
      geom_point(alpha = 0.7, size = 2) +
      geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
      scale_color_viridis(discrete = TRUE, option = "D") +
      theme_minimal() +
      theme(
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "top",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "cm")
      ) +
      labs(title = paste("PSM vs Intensity Correlation:", feature_value),
           x = "PSM Relative Abundance (%)",
           y = "Intensity Relative Abundance (%)",
           color = "Disease Status") +
      annotate("text", 
               x = max(combined_data$psm_relative_percentage, na.rm = TRUE) * 0.1,
               y = max(combined_data$intensity_relative_percentage, na.rm = TRUE) * 0.9,
               label = paste("Pearson r =", round(pearson_result$estimate, 3),
                           "\nSpearman ρ =", round(spearman_result$estimate, 3)),
               hjust = 0, vjust = 1, size = 4)
    
    # Store plot
    plot_list[[feature_value]] <- p
  }
  
  # Arrange results by absolute Pearson correlation
  results <- results %>%
    dplyr::arrange(desc(abs(pearson_correlation)))
  
  # Save results
  write.csv(results, 
            file = output_file, 
            row.names = FALSE)
  
  # Create combined plot if there are multiple features
  if (length(plot_list) > 1) {
    # Combine plots using patchwork
    library(patchwork)
    
    # Create a grid of plots
    combined_plot <- wrap_plots(plot_list, ncol = 2)
    
    # Save combined plot
    ggsave(plot_file, 
           combined_plot, 
           width = 16, 
           height = ceiling(length(plot_list) / 2) * 8, 
           dpi = 300,
           bg = "white")
  } else if (length(plot_list) == 1) {
    # Save single plot
    ggsave(plot_file, 
           plot_list[[1]], 
           width = 10, 
           height = 8, 
           dpi = 300,
           bg = "white")
  }
  
  # Print summary
  cat("\nCorrelation Summary for", feature_col, ":\n")
  cat("Total features analyzed:", nrow(results), "\n")
  cat("Features with significant Pearson correlation (p < 0.05):", 
      sum(results$pearson_p_value < 0.05), "\n")
  cat("Features with significant Spearman correlation (p < 0.05):", 
      sum(results$spearman_p_value < 0.05), "\n")
  
  # Print top 5 correlations
  cat("\nTop 5 Features by Pearson Correlation:\n")
  top_5 <- results %>% head(5)
  for (i in 1:nrow(top_5)) {
    cat(i, ". ", top_5$feature[i], 
        " (Pearson r = ", round(top_5$pearson_correlation[i], 3),
        ", p = ", format(top_5$pearson_p_value[i], scientific = TRUE, digits = 2),
        ")\n", sep = "")
  }
  
  return(list(results = results, plots = plot_list))
}

# Example usage (you'll need to create the intensity-based relative abundance data first)
# correlate_psm_intensity_abundances(
#   psm_relative_abundance_df = composition_relative_abundance,
#   intensity_relative_abundance_df = intensity_composition_relative_abundance,
#   feature_col = "glycan_composition",
#   output_file = "output_data/correlation_psm_intensity_composition.csv",
#   plot_file = "output_data/correlation_plots/psm_intensity_composition_correlation.png"
# )

# Print summary
cat("\nCorrelation Analysis Summary:\n")
cat("Function created to correlate PSM-based vs intensity-based relative abundances.\n")
cat("Use this function with your intensity-based relative abundance data.\n")
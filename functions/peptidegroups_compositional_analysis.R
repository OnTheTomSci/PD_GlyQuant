# Peptide Groups Compositional Analysis Functions
# This module contains functions for ILR-based compositional analysis of glycan data

library(tidyverse)
library(compositions)
library(stats)

#' Calculate ILR transformation for a single protein's glycan compositions
#'
#' @param data_long Long format data with abundance, sample, group, protein_accessions, glycan_composition columns
#' @param protein_id Protein accession to analyze
#' @return List containing ILR results and metadata
calculate_protein_glycan_ilr <- function(data_long, protein_id) {
  
  # Filter data for the specific protein
  protein_data <- data_long %>%
    filter(protein_accessions == protein_id) %>%
    filter(!is.na(glycan_composition), abundance > 0)
  
  if (nrow(protein_data) == 0) {
    stop(paste("No data found for protein:", protein_id))
  }
  
  # Check if we have multiple compositions
  compositions <- unique(protein_data$glycan_composition)
  if (length(compositions) < 2) {
    stop(paste("Need at least 2 glycan compositions for protein:", protein_id, 
               "Found only:", length(compositions)))
  }
  
  # Create wide format data (samples × compositions)
  # First aggregate abundances by sample and composition (in case of multiple peptides)
  aggregated_data <- protein_data %>%
    select(sample, group, glycan_composition, abundance) %>%
    group_by(sample, group, glycan_composition) %>%
    summarise(abundance = sum(abundance, na.rm = TRUE), .groups = 'drop')
  
  wide_data <- aggregated_data %>%
    spread(key = glycan_composition, value = abundance, fill = 0)
  
  # Extract composition matrix and metadata
  composition_matrix <- as.matrix(wide_data[, 3:ncol(wide_data)])
  rownames(composition_matrix) <- wide_data$sample
  group_info <- wide_data$group
  
  # Calculate relative abundances (proportions)
  relative_abundances <- compositions::acomp(composition_matrix)
  
  # Calculate ILR transformation
  # For compositional data with n parts, we get n-1 ILR coordinates
  ilr_coords <- compositions::ilr(relative_abundances)
  
  # Create ILR coordinate names
  n_compositions <- ncol(composition_matrix)
  ilr_names <- paste0("ILR_", 1:(n_compositions-1))
  colnames(ilr_coords) <- ilr_names
  
  # Combine with metadata
  ilr_transformed <- data.frame(
    sample = wide_data$sample,
    group = group_info,
    ilr_coords
  )
  
  # Create relative abundance data frame for reference
  relative_abundance_df <- data.frame(
    sample = wide_data$sample,
    group = group_info,
    relative_abundances
  ) %>%
    pivot_longer(cols = -c(sample, group), 
                 names_to = "glycan_composition", 
                 values_to = "relative_abundance")
  
  # Summary information
  composition_summary <- data.frame(
    composition = compositions,
    mean_relative_abundance = colMeans(relative_abundances),
    min_relative_abundance = apply(relative_abundances, 2, min),
    max_relative_abundance = apply(relative_abundances, 2, max)
  ) %>%
    arrange(desc(mean_relative_abundance))
  
  transformation_info <- list(
    protein_id = protein_id,
    n_compositions = n_compositions,
    n_samples = nrow(wide_data),
    n_groups = length(unique(group_info)),
    compositions = compositions
  )
  
  return(list(
    ilr_transformed = ilr_transformed,
    relative_abundance = relative_abundance_df,
    composition_summary = composition_summary,
    transformation_info = transformation_info,
    raw_compositions = composition_matrix
  ))
}

#' Test ILR coordinate differences between groups
#'
#' @param ilr_result Result from calculate_protein_glycan_ilr()
#' @param p_adjust_method Method for p-value adjustment (default: "BH")
#' @return Data frame with statistical test results
test_ilr_differences <- function(ilr_result, p_adjust_method = "BH") {
  
  ilr_data <- ilr_result$ilr_transformed
  group_info <- ilr_data$group
  
  # Get ILR coordinate names
  ilr_cols <- colnames(ilr_data)[grepl("^ILR_", colnames(ilr_data))]
  
  if (length(ilr_cols) == 0) {
    stop("No ILR coordinates found in the data")
  }
  
  # Test each ILR coordinate
  test_results <- data.frame()
  
  for (coord in ilr_cols) {
    # Extract ILR coordinate values
    healthy_values <- ilr_data[group_info == "Healthy", coord]
    mecfs_values <- ilr_data[group_info == "MECFS", coord]
    
    # Skip if insufficient data
    if (length(healthy_values) < 2 || length(mecfs_values) < 2) {
      next
    }
    
    # Perform t-test
    tryCatch({
      t_test <- t.test(mecfs_values, healthy_values)
      
      # Calculate effect size (Cohen's d)
      pooled_sd <- sqrt(((length(healthy_values) - 1) * var(healthy_values) + 
                        (length(mecfs_values) - 1) * var(mecfs_values)) / 
                       (length(healthy_values) + length(mecfs_values) - 2))
      cohens_d <- (mean(mecfs_values) - mean(healthy_values)) / pooled_sd
      
      # Calculate confidence intervals for the difference
      se_diff <- sqrt(var(healthy_values)/length(healthy_values) + var(mecfs_values)/length(mecfs_values))
      ci_lower <- (mean(mecfs_values) - mean(healthy_values)) - qt(0.975, t_test$parameter) * se_diff
      ci_upper <- (mean(mecfs_values) - mean(healthy_values)) + qt(0.975, t_test$parameter) * se_diff
      
      # Store results
      test_results <- rbind(test_results, data.frame(
        ilr_coordinate = coord,
        healthy_mean = mean(healthy_values),
        mecfs_mean = mean(mecfs_values),
        difference = mean(mecfs_values) - mean(healthy_values),
        healthy_sd = sd(healthy_values),
        mecfs_sd = sd(mecfs_values),
        cohens_d = cohens_d,
        t_statistic = t_test$statistic,
        p_value = t_test$p.value,
        df = t_test$parameter,
        n_healthy = length(healthy_values),
        n_mecfs = length(mecfs_values),
        ci_lower = ci_lower,
        ci_upper = ci_upper
      ))
    }, error = function(e) {
      # If t-test fails, record NA values
      test_results <<- rbind(test_results, data.frame(
        ilr_coordinate = coord,
        healthy_mean = mean(healthy_values, na.rm = TRUE),
        mecfs_mean = mean(mecfs_values, na.rm = TRUE),
        difference = mean(mecfs_values, na.rm = TRUE) - mean(healthy_values, na.rm = TRUE),
        healthy_sd = sd(healthy_values, na.rm = TRUE),
        mecfs_sd = sd(mecfs_values, na.rm = TRUE),
        cohens_d = NA,
        t_statistic = NA,
        p_value = NA,
        df = NA,
        n_healthy = length(healthy_values),
        n_mecfs = length(mecfs_values),
        ci_lower = NA,
        ci_upper = NA
      ))
    })
  }
  
  # Adjust p-values
  if (nrow(test_results) > 0 && any(!is.na(test_results$p_value))) {
    test_results$p_value_adj <- p.adjust(test_results$p_value, method = p_adjust_method)
    test_results$significant <- test_results$p_value_adj < 0.05
  } else {
    test_results$p_value_adj <- NA
    test_results$significant <- FALSE
  }
  
  return(test_results)
}

#' Create boxplot for a specific ILR coordinate
#'
#' @param ilr_result Result from calculate_protein_glycan_ilr()
#' @param ilr_coordinate Name of the ILR coordinate to plot
#' @return ggplot2 boxplot object
plot_ilr_boxplot <- function(ilr_result, ilr_coordinate) {
  
  ilr_data <- ilr_result$ilr_transformed
  protein_id <- ilr_result$transformation_info$protein_id
  
  # Check if coordinate exists
  if (!ilr_coordinate %in% colnames(ilr_data)) {
    stop(paste("ILR coordinate", ilr_coordinate, "not found in data"))
  }
  
  # Create plot
  plot <- ggplot(ilr_data, aes_string(x = "group", y = ilr_coordinate, fill = "group")) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.8, size = 2) +
    scale_fill_manual(values = c("Healthy" = "#2E8B57", "MECFS" = "#DC143C")) +
    labs(
      title = paste("ILR Coordinate:", ilr_coordinate),
      subtitle = paste("Protein:", protein_id),
      x = "Group",
      y = paste("ILR Coordinate Value:", ilr_coordinate)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "none",
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12, face = "bold")
    )
  
  return(plot)
}

#' Analyze multiple proteins using ILR transformation
#'
#' @param data_long Long format data with abundance, sample, group, protein_accessions, glycan_composition columns
#' @param protein_ids Vector of protein accessions to analyze
#' @param output_dir Directory to save output files (optional)
#' @return List containing results for each protein
analyze_multiple_proteins_ilr <- function(data_long, protein_ids, output_dir = NULL) {
  
  cat("Starting ILR analysis for", length(protein_ids), "proteins...\n")
  
  results <- list()
  successful_analyses <- 0
  failed_analyses <- 0
  
  for (i in seq_along(protein_ids)) {
    protein_id <- protein_ids[i]
    
    cat(sprintf("Processing %d/%d: %s\n", i, length(protein_ids), protein_id))
    
    tryCatch({
      # Calculate ILR transformation
      ilr_result <- calculate_protein_glycan_ilr(data_long, protein_id)
      
      # Perform statistical tests
      statistical_tests <- test_ilr_differences(ilr_result)
      
      # Add statistical tests to result
      ilr_result$statistical_tests <- statistical_tests
      
      # Store result
      results[[protein_id]] <- ilr_result
      
      # Save individual files if output directory specified
      if (!is.null(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        
        # Save ILR transformed data
        write.csv(ilr_result$ilr_transformed, 
                  file.path(output_dir, paste0(protein_id, "_ilr_data.csv")), 
                  row.names = FALSE)
        
        # Save relative abundance data
        write.csv(ilr_result$relative_abundance, 
                  file.path(output_dir, paste0(protein_id, "_relative_abundance.csv")), 
                  row.names = FALSE)
        
        # Save statistical tests
        write.csv(statistical_tests, 
                  file.path(output_dir, paste0(protein_id, "_statistical_tests.csv")), 
                  row.names = FALSE)
        
        # Save composition summary
        write.csv(ilr_result$composition_summary, 
                  file.path(output_dir, paste0(protein_id, "_composition_summary.csv")), 
                  row.names = FALSE)
        
        # Create and save plots for significant coordinates
        if (nrow(statistical_tests) > 0 && any(statistical_tests$significant, na.rm = TRUE)) {
          sig_coords <- statistical_tests$ilr_coordinate[statistical_tests$significant]
          
          for (coord in sig_coords) {
            plot <- plot_ilr_boxplot(ilr_result, coord)
            ggsave(file.path(output_dir, paste0(protein_id, "_", coord, "_plot.png")),
                   plot, width = 8, height = 6, dpi = 300)
          }
        }
      }
      
      successful_analyses <- successful_analyses + 1
      
    }, error = function(e) {
      cat(sprintf("  Error analyzing %s: %s\n", protein_id, e$message))
      failed_analyses <- failed_analyses + 1
    })
  }
  
  cat(sprintf("\nAnalysis complete: %d successful, %d failed\n", 
              successful_analyses, failed_analyses))
  
  if (!is.null(output_dir)) {
    cat(sprintf("Results saved to: %s\n", output_dir))
  }
  
  return(results)
}

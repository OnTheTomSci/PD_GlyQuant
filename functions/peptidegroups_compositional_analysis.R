# Peptide Groups Compositional Analysis Functions
# This module contains functions for CLR and ILR-based compositional analysis of glycan data

library(tidyverse)
library(compositions)
library(stats)

#' Calculate CLR transformation for a multiple protein's glycan compositions
#' 
#' @param data_long Long format data with abundance, sample, group, protein_accessions, glycan_composition columns
#' @return List containing CLR results and metadata
calculate_multiple_proteins_glycan_clr <- function(data_long) {
  
  cat("Starting CLR transformation for glycan compositions across proteins...\n")
  
  # Filter out missing data
  clean_data <- data_long %>%
    filter(!is.na(glycan_composition), 
           !is.na(protein_accessions),
           abundance > 0)
  
  if (nrow(clean_data) == 0) {
    stop("No valid data found after filtering")
  }
  
  # Aggregate abundances by sample, protein, and glycan composition
  # This sums abundances if there are multiple peptides with the same composition
  aggregated_data <- clean_data %>%
    select(sample, group, protein_accessions, glycan_composition, abundance) %>%
    group_by(sample, group, protein_accessions, glycan_composition) %>%
    summarise(abundance = sum(abundance, na.rm = TRUE), .groups = 'drop')
  
  # Calculate relative abundance (proportions) within each protein and sample
  # This normalizes so that for each protein-sample combination, compositions sum to 1
  relative_abundance_data <- aggregated_data %>%
    group_by(sample, protein_accessions) %>%
    mutate(
      protein_total = sum(abundance, na.rm = TRUE),
      relative_abundance = abundance / protein_total
    ) %>%
    ungroup()
  
  # For CLR transformation, we need to convert to wide format
  # Create a composite key for unique combinations
  wide_data <- relative_abundance_data %>%
    mutate(sample_protein = paste(sample, protein_accessions, sep = "_")) %>%
    select(sample_protein, sample, group, protein_accessions, glycan_composition, relative_abundance) %>%
    spread(key = glycan_composition, value = relative_abundance, fill = 0)
  
  # Extract the composition matrix
  metadata_cols <- c("sample_protein", "sample", "group", "protein_accessions")
  composition_matrix <- as.matrix(wide_data[, !colnames(wide_data) %in% metadata_cols])
  rownames(composition_matrix) <- wide_data$sample_protein
  
  # Apply CLR transformation
  # CLR = log(x_i) - mean(log(x))
  # Using the compositions package for proper handling
  relative_abundances_acomp <- compositions::acomp(composition_matrix)
  clr_transformed_matrix <- compositions::clr(relative_abundances_acomp)
  
  # Convert CLR results back to data frame with metadata
  clr_transformed_df <- data.frame(
    sample = wide_data$sample,
    group = wide_data$group,
    protein_accessions = wide_data$protein_accessions,
    clr_transformed_matrix
  )
  
  # Convert to long format for easier analysis
  clr_transformed_long <- clr_transformed_df %>%
    pivot_longer(
      cols = -c(sample, group, protein_accessions),
      names_to = "glycan_composition",
      values_to = "clr_value"
    )
  
  # Create summary statistics
  clr_summary <- clr_transformed_long %>%
    group_by(protein_accessions, glycan_composition) %>%
    summarise(
      mean_clr = mean(clr_value, na.rm = TRUE),
      sd_clr = sd(clr_value, na.rm = TRUE),
      median_clr = median(clr_value, na.rm = TRUE),
      min_clr = min(clr_value, na.rm = TRUE),
      max_clr = max(clr_value, na.rm = TRUE),
      n_samples = n(),
      .groups = 'drop'
    ) %>%
    arrange(protein_accessions, desc(mean_clr))
  
  # Summary by protein
  protein_summary <- aggregated_data %>%
    group_by(protein_accessions) %>%
    summarise(
      n_compositions = n_distinct(glycan_composition),
      n_samples = n_distinct(sample),
      total_abundance = sum(abundance, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(n_compositions))
  
  # Metadata
  transformation_info <- list(
    n_proteins = n_distinct(clean_data$protein_accessions),
    n_samples = n_distinct(clean_data$sample),
    n_compositions = n_distinct(clean_data$glycan_composition),
    n_groups = n_distinct(clean_data$group),
    total_observations = nrow(clr_transformed_long)
  )
  
  cat(sprintf("CLR transformation complete:\n"))
  cat(sprintf("  - Proteins analyzed: %d\n", transformation_info$n_proteins))
  cat(sprintf("  - Samples: %d\n", transformation_info$n_samples))
  cat(sprintf("  - Unique glycan compositions: %d\n", transformation_info$n_compositions))
  cat(sprintf("  - Total CLR values: %d\n", transformation_info$total_observations))
  
  return(list(
    clr_transformed = clr_transformed_df,
    clr_transformed_long = clr_transformed_long,
    relative_abundance = relative_abundance_data,
    clr_summary = clr_summary,
    protein_summary = protein_summary,
    transformation_info = transformation_info,
    raw_matrix = composition_matrix,
    clr_matrix = clr_transformed_matrix
  ))
}

#' Test CLR protein composition differences between groups
#'
#' @param clr_result Result from calculate_multiple_proteins_glycan_clr()
#' @param group1_name Name of the first group (default: "Healthy")
#' @param group2_name Name of the second group (default: "MECFS")
#' @param p_adjust_method Method for p-value adjustment (default: "BH")
#' @return Data frame with statistical test results for each protein-glycan combination
test_clr_protein_composition_differences <- function(clr_result, 
                                                     group1_name = "Healthy", 
                                                     group2_name = "MECFS",
                                                     p_adjust_method = "BH") {
  
  cat("Testing CLR differences between", group1_name, "and", group2_name, "...\n")
  
  # Extract the long format CLR data
  clr_data <- clr_result$clr_transformed_long
  
  # Get unique protein-glycan combinations
  protein_glycan_combos <- clr_data %>%
    select(protein_accessions, glycan_composition) %>%
    distinct() %>%
    arrange(protein_accessions, glycan_composition)
  
  cat(sprintf("Testing %d protein-glycan combinations...\n", nrow(protein_glycan_combos)))
  
  # Initialize results data frame
  test_results <- data.frame()
  
  # Test each protein-glycan combination
  for (i in 1:nrow(protein_glycan_combos)) {
    protein_id <- protein_glycan_combos$protein_accessions[i]
    glycan_comp <- protein_glycan_combos$glycan_composition[i]
    
    # Filter data for this specific protein-glycan combination
    combo_data <- clr_data %>%
      filter(protein_accessions == protein_id,
             glycan_composition == glycan_comp)
    
    # Extract values for each group
    group1_values <- combo_data %>%
      filter(group == group1_name) %>%
      pull(clr_value)
    
    group2_values <- combo_data %>%
      filter(group == group2_name) %>%
      pull(clr_value)
    
    # Skip if insufficient data
    if (length(group1_values) < 2 || length(group2_values) < 2) {
      next
    }
    
    # Perform t-test and calculate statistics
    tryCatch({
      t_test <- t.test(group2_values, group1_values)
      
      # Calculate effect size (Cohen's d)
      pooled_sd <- sqrt(((length(group1_values) - 1) * var(group1_values) + 
                        (length(group2_values) - 1) * var(group2_values)) / 
                       (length(group1_values) + length(group2_values) - 2))
      
      cohens_d <- if(pooled_sd > 0) {
        (mean(group2_values) - mean(group1_values)) / pooled_sd
      } else {
        NA
      }
      
      # Calculate fold change in CLR space
      clr_difference <- mean(group2_values) - mean(group1_values)
      
      # Calculate confidence intervals for the difference
      se_diff <- sqrt(var(group1_values)/length(group1_values) + 
                     var(group2_values)/length(group2_values))
      ci_lower <- clr_difference - qt(0.975, t_test$parameter) * se_diff
      ci_upper <- clr_difference + qt(0.975, t_test$parameter) * se_diff
      
      # Store results
      test_results <- rbind(test_results, data.frame(
        protein_accessions = protein_id,
        glycan_composition = glycan_comp,
        group1_name = group1_name,
        group2_name = group2_name,
        group1_mean_clr = mean(group1_values),
        group2_mean_clr = mean(group2_values),
        clr_difference = clr_difference,
        group1_sd = sd(group1_values),
        group2_sd = sd(group2_values),
        cohens_d = cohens_d,
        effect_size_interpretation = case_when(
          is.na(cohens_d) ~ "NA",
          abs(cohens_d) < 0.2 ~ "negligible",
          abs(cohens_d) < 0.5 ~ "small",
          abs(cohens_d) < 0.8 ~ "medium",
          TRUE ~ "large"
        ),
        t_statistic = t_test$statistic,
        p_value = t_test$p.value,
        df = t_test$parameter,
        n_group1 = length(group1_values),
        n_group2 = length(group2_values),
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        stringsAsFactors = FALSE
      ))
    }, error = function(e) {
      # If t-test fails, record NA values
      test_results <<- rbind(test_results, data.frame(
        protein_accessions = protein_id,
        glycan_composition = glycan_comp,
        group1_name = group1_name,
        group2_name = group2_name,
        group1_mean_clr = mean(group1_values, na.rm = TRUE),
        group2_mean_clr = mean(group2_values, na.rm = TRUE),
        clr_difference = mean(group2_values, na.rm = TRUE) - mean(group1_values, na.rm = TRUE),
        group1_sd = sd(group1_values, na.rm = TRUE),
        group2_sd = sd(group2_values, na.rm = TRUE),
        cohens_d = NA,
        effect_size_interpretation = "NA",
        t_statistic = NA,
        p_value = NA,
        df = NA,
        n_group1 = length(group1_values),
        n_group2 = length(group2_values),
        ci_lower = NA,
        ci_upper = NA,
        stringsAsFactors = FALSE
      ))
    })
  }
  
  # Adjust p-values for multiple testing
  if (nrow(test_results) > 0 && any(!is.na(test_results$p_value))) {
    test_results$p_value_adj <- p.adjust(test_results$p_value, method = p_adjust_method)
    test_results$significant <- test_results$p_value_adj < 0.05
    test_results$significance_level <- case_when(
      is.na(test_results$p_value_adj) ~ "NS",
      test_results$p_value_adj < 0.001 ~ "***",
      test_results$p_value_adj < 0.01 ~ "**",
      test_results$p_value_adj < 0.05 ~ "*",
      TRUE ~ "NS"
    )
  } else {
    test_results$p_value_adj <- NA
    test_results$significant <- FALSE
    test_results$significance_level <- "NS"
  }
  
  # Sort by p-value (most significant first)
  test_results <- test_results %>%
    arrange(p_value_adj, protein_accessions)
  
  # Print summary
  n_significant <- sum(test_results$significant, na.rm = TRUE)
  n_total <- nrow(test_results)
  
  cat(sprintf("\nStatistical testing complete:\n"))
  cat(sprintf("  - Total tests: %d\n", n_total))
  cat(sprintf("  - Significant (adj p < 0.05): %d (%.1f%%)\n", 
              n_significant, 100 * n_significant / n_total))
  
  if (n_significant > 0) {
    cat(sprintf("  - Significant by level:\n"))
    cat(sprintf("      p < 0.001 (***): %d\n", sum(test_results$significance_level == "***", na.rm = TRUE)))
    cat(sprintf("      p < 0.01 (**): %d\n", sum(test_results$significance_level == "**", na.rm = TRUE)))
    cat(sprintf("      p < 0.05 (*): %d\n", sum(test_results$significance_level == "*", na.rm = TRUE)))
    
    # Count by effect size among significant results
    sig_results <- test_results %>% filter(significant)
    cat(sprintf("  - Effect sizes (among significant):\n"))
    cat(sprintf("      Large: %d\n", sum(sig_results$effect_size_interpretation == "large", na.rm = TRUE)))
    cat(sprintf("      Medium: %d\n", sum(sig_results$effect_size_interpretation == "medium", na.rm = TRUE)))
    cat(sprintf("      Small: %d\n", sum(sig_results$effect_size_interpretation == "small", na.rm = TRUE)))
  }
  
  return(test_results)
}

#' Create dot plot showing CLR differences across proteins and glycan compositions
#'
#' @param clr_stats_result Result from test_clr_protein_composition_differences()
#' @param p_threshold P-value threshold for significance (default: 0.05)
#' @param top_n_proteins Optional: show only top N proteins by number of significant glycans (default: NULL shows all)
#' @param min_samples_per_group Minimum number of samples required per group (default: NULL for no filtering)
#' @return ggplot2 dot plot object
plot_clr_dotplot <- function(clr_stats_result, 
                             p_threshold = 0.05,
                             top_n_proteins = NULL,
                             min_samples_per_group = NULL) {
  
  library(ggplot2)
  
  if (nrow(clr_stats_result) == 0) {
    stop("No results to plot")
  }
  
  # Prepare data
  plot_data <- clr_stats_result %>%
    mutate(
      significant = p_value_adj < p_threshold,
      direction = case_when(
        !significant ~ "Not Significant",
        clr_difference > 0 ~ paste("Increased in", group2_name),
        clr_difference < 0 ~ paste("Decreased in", group2_name)
      ),
      abs_clr_diff = abs(clr_difference),
      neg_log_p = -log10(p_value_adj)
    )
  
  # Filter by minimum samples per group if specified
  if (!is.null(min_samples_per_group)) {
    # First, identify protein-glycan combinations with enough samples in both groups
    valid_protein_glycans <- plot_data %>%
      filter(n_group1 >= min_samples_per_group & n_group2 >= min_samples_per_group) %>%
      select(protein_accessions, glycan_composition) %>%
      distinct()
    
    # Count how many valid glycans each protein has
    proteins_with_enough_glycans <- valid_protein_glycans %>%
      group_by(protein_accessions) %>%
      summarise(n_valid_glycans = n(), .groups = 'drop') %>%
      filter(n_valid_glycans >= 3) %>%  # Require at least 3 glycans
      pull(protein_accessions)
    
    if (length(proteins_with_enough_glycans) == 0) {
      stop(sprintf("No proteins found with at least 3 glycans (each with >= %d samples per group)", 
                   min_samples_per_group))
    }
    
    # Filter to only keep valid protein-glycan combinations
    plot_data <- plot_data %>%
      semi_join(valid_protein_glycans, by = c("protein_accessions", "glycan_composition")) %>%
      filter(protein_accessions %in% proteins_with_enough_glycans)
    
    cat(sprintf("Filtered to %d proteins with >= 3 glycans (each glycan found in >= %d samples per group)\n", 
                length(proteins_with_enough_glycans), min_samples_per_group))
  }
  
  # Optionally filter to top proteins
  if (!is.null(top_n_proteins)) {
    top_proteins <- plot_data %>%
      filter(significant) %>%
      group_by(protein_accessions) %>%
      summarise(n_sig = n(), .groups = 'drop') %>%
      arrange(desc(n_sig)) %>%
      head(top_n_proteins) %>%
      pull(protein_accessions)
    
    plot_data <- plot_data %>%
      filter(protein_accessions %in% top_proteins)
    
    if (nrow(plot_data) == 0) {
      stop("No significant results in top proteins")
    }
  }
  
  # Order proteins by number of significant glycans
  protein_order <- plot_data %>%
    group_by(protein_accessions) %>%
    summarise(
      n_significant = sum(significant),
      mean_abs_diff = mean(abs_clr_diff, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(n_significant), desc(mean_abs_diff)) %>%
    pull(protein_accessions)
  
  # Order glycans by name
  glycan_order <- sort(unique(plot_data$glycan_composition))
  
  plot_data <- plot_data %>%
    mutate(
      protein_accessions = factor(protein_accessions, levels = protein_order),
      glycan_composition = factor(glycan_composition, levels = glycan_order)
    )
  
  # Create dot plot (FLIPPED: protein on x-axis, glycan on y-axis)
  plot <- ggplot(plot_data, aes(x = protein_accessions, y = glycan_composition)) +
    geom_point(aes(size = abs_clr_diff, 
                   color = direction,
                   alpha = significant)) +
    scale_size_continuous(
      name = "Abs(CLR Difference)",
      range = c(1, 8),
      breaks = c(0.5, 1, 2, 3)
    ) +
    scale_color_manual(
      values = c(
        "Increased in MECFS" = "#E49CB1",
        "Decreased in MECFS" = "#9DD4CC",
        "Increased in Healthy" = "#E49CB1",
        "Decreased in Healthy" = "#9DD4CC",
        "Not Significant" = "gray80"
      ),
      name = "Direction"
    ) +
    scale_alpha_manual(
      values = c("TRUE" = 1, "FALSE" = 0.3),
      guide = "none"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 7),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "CLR Compositional Analysis - Dot Plot",
      subtitle = paste("Significant results (p <", p_threshold, ") shown with full opacity"),
      x = "Protein",
      y = "Glycan Composition"
    )
  
  # Add count of significant results to subtitle
  n_sig <- sum(plot_data$significant, na.rm = TRUE)
  n_total <- nrow(plot_data)
  
  subtitle_text <- sprintf("Significant: %d/%d (%.1f%%) | p-adj < %.3f",
                          n_sig, n_total, 100*n_sig/n_total, p_threshold)
  
  if (!is.null(min_samples_per_group)) {
    subtitle_text <- paste0(subtitle_text, 
                           sprintf(" | Proteins with ≥3 glycans (each in ≥%d samples/group)", 
                                   min_samples_per_group))
  }
  
  plot <- plot +
    labs(subtitle = subtitle_text)
  
  return(plot)
}

#' Create heatmap of CLR differences
#'
#' @param clr_stats_result Result from test_clr_protein_composition_differences()
#' @param value_type What to display: "p_value", "effect", or "clr_difference" (default)
#' @param p_threshold P-value threshold for significance (default: 0.05)
#' @param top_n_proteins Optional: show only top N proteins (default: NULL shows all)
#' @return ggplot2 heatmap object
plot_clr_heatmap <- function(clr_stats_result,
                             value_type = "clr_difference",
                             p_threshold = 0.05,
                             top_n_proteins = NULL) {
  
  library(ggplot2)
  
  if (nrow(clr_stats_result) == 0) {
    stop("No results to plot")
  }
  
  # Prepare data
  heatmap_data <- clr_stats_result %>%
    mutate(significant = p_value_adj < p_threshold)
  
  # Optionally filter to top proteins
  if (!is.null(top_n_proteins)) {
    top_proteins <- heatmap_data %>%
      filter(significant) %>%
      group_by(protein_accessions) %>%
      summarise(n_sig = n(), .groups = 'drop') %>%
      arrange(desc(n_sig)) %>%
      head(top_n_proteins) %>%
      pull(protein_accessions)
    
    heatmap_data <- heatmap_data %>%
      filter(protein_accessions %in% top_proteins)
  }
  
  # Order proteins by significance
  protein_order <- heatmap_data %>%
    group_by(protein_accessions) %>%
    summarise(
      n_significant = sum(significant),
      mean_abs_diff = mean(abs(clr_difference), na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(n_significant), desc(mean_abs_diff)) %>%
    pull(protein_accessions)
  
  heatmap_data <- heatmap_data %>%
    mutate(protein_accessions = factor(protein_accessions, levels = protein_order))
  
  # Create heatmap based on value type
  if (value_type == "p_value") {
    # Color by p-value
    plot <- ggplot(heatmap_data, 
                   aes(x = glycan_composition, y = protein_accessions, fill = -log10(p_value_adj))) +
      geom_tile(color = "white", size = 0.5) +
      geom_text(aes(label = ifelse(significant, "*", "")), 
                size = 4, color = "black") +
      scale_fill_gradient2(
        low = "lightgray", mid = "yellow", high = "red",
        midpoint = -log10(0.05),
        name = "-log10(p-adj)",
        limits = c(0, max(-log10(heatmap_data$p_value_adj), 2))
      )
  } else if (value_type == "effect" || value_type == "clr_difference") {
    # Color by CLR difference
    max_abs_diff <- max(abs(heatmap_data$clr_difference), na.rm = TRUE)
    plot <- ggplot(heatmap_data, 
                   aes(x = glycan_composition, y = protein_accessions, fill = clr_difference)) +
      geom_tile(color = "white", size = 0.5) +
      geom_point(data = filter(heatmap_data, significant),
                 shape = 8, size = 3, color = "black") +
      scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B",
        midpoint = 0,
        name = paste0("CLR Difference\n(", unique(heatmap_data$group2_name)[1], " - ", 
                     unique(heatmap_data$group1_name)[1], ")"),
        limits = c(-max_abs_diff, max_abs_diff)
      )
  } else {
    stop("value_type must be 'p_value', 'effect', or 'clr_difference'")
  }
  
  plot <- plot +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "right"
    ) +
    labs(
      title = "CLR Compositional Analysis - Heatmap",
      subtitle = "* indicates p-adj < 0.05",
      x = "Glycan Composition",
      y = "Protein"
    )
  
  return(plot)
}

#' Create forest plot of CLR differences
#'
#' @param clr_stats_result Result from test_clr_protein_composition_differences()
#' @param show_only_significant Show only significant results (default: TRUE)
#' @param p_threshold P-value threshold (default: 0.05)
#' @param top_n Optional: show only top N results (default: 50)
#' @return ggplot2 forest plot object
plot_clr_forest <- function(clr_stats_result,
                            show_only_significant = TRUE,
                            p_threshold = 0.05,
                            top_n = 50) {
  
  library(ggplot2)
  
  if (nrow(clr_stats_result) == 0) {
    stop("No results to plot")
  }
  
  # Prepare data
  forest_data <- clr_stats_result %>%
    mutate(
      significant = p_value_adj < p_threshold,
      sig_label = ifelse(significant, "Significant", "Not Significant"),
      protein_glycan = paste0(protein_accessions, "\n", glycan_composition)
    )
  
  # Filter if requested
  if (show_only_significant) {
    forest_data <- forest_data %>% filter(significant)
    
    if (nrow(forest_data) == 0) {
      warning("No significant results to plot")
      return(NULL)
    }
  }
  
  # Limit to top N by p-value
  if (!is.null(top_n) && nrow(forest_data) > top_n) {
    forest_data <- forest_data %>%
      arrange(p_value_adj) %>%
      head(top_n)
  }
  
  # Sort by effect size
  forest_data <- forest_data %>%
    arrange(desc(abs(clr_difference))) %>%
    mutate(protein_glycan = factor(protein_glycan, levels = protein_glycan))
  
  # Create forest plot
  plot <- ggplot(forest_data, 
                 aes(x = clr_difference, y = protein_glycan, color = sig_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), 
                   height = 0.3, alpha = 0.7) +
    geom_point(size = 3) +
    scale_color_manual(values = c("Significant" = "#E49CB1", 
                                  "Not Significant" = "gray60")) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 7),
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    ) +
    labs(
      title = "CLR Differences with 95% Confidence Intervals",
      subtitle = sprintf("Top %d results (by p-value)", nrow(forest_data)),
      x = paste0("CLR Difference (", unique(forest_data$group2_name)[1], " - ", 
                unique(forest_data$group1_name)[1], ")"),
      y = "Protein - Glycan Composition",
      color = NULL
    )
  
  return(plot)
}

#' Create comprehensive summary plot of CLR analysis
#'
#' @param clr_stats_result Result from test_clr_protein_composition_differences()
#' @param output_file Optional filename to save combined plot
#' @param p_threshold P-value threshold (default: 0.05)
#' @param top_n_proteins Show top N proteins in heatmap/dotplot (default: 30)
#' @param min_samples_per_group Minimum samples per group for filtering (default: NULL)
#' @return Combined ggplot object (using patchwork)
plot_clr_summary <- function(clr_stats_result, 
                            output_file = NULL,
                            p_threshold = 0.05,
                            top_n_proteins = 30,
                            min_samples_per_group = NULL) {
  
  # Requires patchwork for combining plots
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required. Install with: install.packages('patchwork')")
  }
  
  library(patchwork)
  
  # Check if there are any significant results
  n_significant <- sum(clr_stats_result$significant, na.rm = TRUE)
  
  # Create individual plots
  cat("Creating heatmap...\n")
  heatmap <- plot_clr_heatmap(clr_stats_result, 
                              value_type = "clr_difference", 
                              p_threshold = p_threshold,
                              top_n_proteins = if(n_significant > 0) top_n_proteins else NULL)
  
  cat("Creating forest plot...\n")
  forest <- plot_clr_forest(clr_stats_result, 
                            show_only_significant = TRUE, 
                            p_threshold = p_threshold,
                            top_n = 30)
  
  cat("Creating dot plot...\n")
  dotplot <- plot_clr_dotplot(clr_stats_result, 
                              p_threshold = p_threshold,
                              top_n_proteins = if(n_significant > 0) top_n_proteins else NULL,
                              min_samples_per_group = min_samples_per_group)
  
  # Combine plots
  if (!is.null(forest)) {
    combined_plot <- (heatmap / forest) | dotplot
    combined_plot <- combined_plot + 
      plot_annotation(
        title = "CLR Compositional Analysis - Summary",
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
      )
  } else {
    combined_plot <- heatmap / dotplot
    combined_plot <- combined_plot + 
      plot_annotation(
        title = "CLR Compositional Analysis - Summary",
        subtitle = "No significant results for forest plot",
        theme = theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)
        )
      )
  }
  
  # Save if requested
  if (!is.null(output_file)) {
    cat(sprintf("Saving combined plot to %s...\n", output_file))
    ggsave(output_file, combined_plot, width = 18, height = 12, dpi = 300)
  }
  
  return(combined_plot)
}

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


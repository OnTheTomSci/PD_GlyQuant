# Peptide Groups Glycan Analysis Functions
# This module contains functions for analyzing glycan compositions and pseudo-glycomics

library(tidyverse)
library(ggplot2)
library(patchwork)

#' Analyze glycan compositions and create visualizations
#' 
#' @param data Long format data with abundance, sample, glycan_composition columns
#' @param output_dir Directory to save output files
#' @param figures_dir Directory to save figures
#' @return List containing analysis results and plots
analyze_glycan_compositions <- function(data, output_dir = "output_data/peptidegroups_intensity", 
                                       figures_dir = "figures") {
  
  # Create output directories if they don't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Calculate relative abundance by glycan composition and sample
  relative_abundance <- data %>%
    group_by(sample, glycan_composition) %>%
    summarise(
      total_abundance = sum(abundance, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    group_by(sample) %>%
    mutate(
      sample_total = sum(total_abundance, na.rm = TRUE),
      relative_abundance = (total_abundance / sample_total) * 100
    ) %>%
    ungroup() %>%
    # Add group classification based on sample names
    mutate(
      group = case_when(
        str_starts(sample, "hc") ~ "Healthy",
        str_starts(sample, "m") ~ "MECFS",
        TRUE ~ "Unknown"
      )
    )
  
  # Calculate mean and standard error for each glycan composition
  glycan_summary <- relative_abundance %>%
    group_by(glycan_composition) %>%
    summarise(
      mean_relative_abundance = mean(relative_abundance, na.rm = TRUE),
      se_relative_abundance = sd(relative_abundance, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
  
  # Save the glycan summary
  write_csv(glycan_summary, file.path(output_dir, "glycan_composition_peptidegroups_intensity_summary.csv"))
  
  # Create bar plot with error bars
  glycan_plot <- plot_glycan_composition_barplot(glycan_summary, 
                                                title = "Relative Abundance of Glycan Compositions")
  
  # Save the plot
  ggsave(file.path(figures_dir, "glycan_composition_relative_abundance.png"), 
         glycan_plot, width = 12, height = 8, dpi = 300)
  
  cat("Glycan composition plot saved to:", file.path(figures_dir, "glycan_composition_relative_abundance.png"), "\n")
  
  # Calculate the sum of top 5 glycans for each sample
  top_5_by_sample <- relative_abundance %>%
    group_by(sample) %>%
    arrange(desc(relative_abundance), .by_group = TRUE) %>%
    slice_head(n = 5) %>%
    summarise(
      top_5_sum = sum(relative_abundance),
      top_5_glycans = paste(glycan_composition, collapse = ", ")
    )
  
  # Calculate the overall top 5 glycan compositions across all samples
  overall_top_5 <- relative_abundance %>%
    group_by(glycan_composition) %>%
    summarise(
      mean_abundance = mean(relative_abundance, na.rm = TRUE),
      sd_abundance = sd(relative_abundance, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(mean_abundance)) %>%
    slice_head(n = 5)
  
  # Print the overall top 5 glycan compositions
  cat("\nOverall Top 5 Most Abundant Glycan Compositions:\n")
  for(i in 1:nrow(overall_top_5)) {
    cat(sprintf("%d. %s (Mean: %.2f%% ± %.2f%% SD)\n",
                i,
                overall_top_5$glycan_composition[i],
                overall_top_5$mean_abundance[i],
                overall_top_5$sd_abundance[i]))
  }
  
  # Perform t-test and F-test for each glycan composition between disease groups
  glycan_ttests <- relative_abundance %>%
    group_by(glycan_composition) %>%
    summarise(
      # F-test for variance
      f_stat = tryCatch({
        var.test(relative_abundance[group == "MECFS"], 
                 relative_abundance[group == "Healthy"])$statistic
      }, error = function(e) NA_real_),
      f_pvalue = tryCatch({
        var.test(relative_abundance[group == "MECFS"],
                 relative_abundance[group == "Healthy"])$p.value
      }, error = function(e) NA_real_),
      # t-test 
      t_stat = tryCatch({
        t.test(relative_abundance ~ group)$statistic
      }, error = function(e) NA_real_),
      p_value = tryCatch({
        t.test(relative_abundance ~ group)$p.value
      }, error = function(e) NA_real_),
      mean_mecfs = mean(relative_abundance[group == "MECFS"], na.rm = TRUE),
      mean_healthy = mean(relative_abundance[group == "Healthy"], na.rm = TRUE), 
      sd_mecfs = sd(relative_abundance[group == "MECFS"], na.rm = TRUE),
      sd_healthy = sd(relative_abundance[group == "Healthy"], na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    # Add BH correction
    mutate(
      p_value_adj = p.adjust(p_value, method = "BH"),
      significant = p_value_adj < 0.05,
      fold_change = mean_mecfs / mean_healthy,
      variance_equal = f_pvalue >= 0.05
    ) %>%
    arrange(p_value_adj)
  
  # Save t-test results
  write.csv(glycan_ttests, file.path(output_dir, "glycan_composition_ttests.csv"), row.names = TRUE)
  
  # Print significant results
  cat("\nSignificant differences in glycan compositions (p < 0.05):\n")
  significant_glycans <- glycan_ttests %>% filter(significant)
  
  if(nrow(significant_glycans) > 0) {
    for(i in 1:nrow(significant_glycans)) {
      cat(sprintf("\n%d. %s\n", i, significant_glycans$glycan_composition[i]))
      cat(sprintf("   MECFS: %.2f%% ± %.2f%% SD\n", 
                  significant_glycans$mean_mecfs[i],
                  significant_glycans$sd_mecfs[i]))
      cat(sprintf("   Healthy: %.2f%% ± %.2f%% SD\n", 
                  significant_glycans$mean_healthy[i],
                  significant_glycans$sd_healthy[i]))
      cat(sprintf("   p-value: %.3e\n", significant_glycans$p_value[i]))
      cat(sprintf("   Fold change: %.2f\n", significant_glycans$fold_change[i]))
    }
  } else {
    cat("No significant differences found between MECFS and Healthy groups.\n")
  }
  
  # Find glycan compositions with fold change > 2 or < 0.5 (1/2)
  high_fold_changes <- glycan_ttests %>%
    filter(fold_change > 2 | fold_change < (1/2)) %>%
    arrange(desc(fold_change))
  
  # Print summary
  cat("\nGlycan compositions with >2-fold change between groups:\n")
  cat("Total count:", nrow(high_fold_changes), "\n")
  
  if(nrow(high_fold_changes) > 0) {
    for(i in 1:nrow(high_fold_changes)) {
      cat(sprintf("\n%d. %s\n", i, high_fold_changes$glycan_composition[i]))
      cat(sprintf("   Fold change: %.2f\n", high_fold_changes$fold_change[i]))
      cat(sprintf("   MECFS: %.2f%% ± %.2f%% SD\n", 
                  high_fold_changes$mean_mecfs[i],
                  high_fold_changes$sd_mecfs[i]))
      cat(sprintf("   Healthy: %.2f%% ± %.2f%% SD\n", 
                  high_fold_changes$mean_healthy[i],
                  high_fold_changes$sd_healthy[i]))
      cat(sprintf("   Adjusted p-value: %.3e\n", high_fold_changes$p_value_adj[i]))
    }
  } else {
    cat("No glycan compositions had >2-fold change between groups.\n")
  }
  
  # Create a bar plot of the top 5 glycan compositions
  top_5_plot <- ggplot(overall_top_5, 
                       aes(x = reorder(glycan_composition, mean_abundance), 
                           y = mean_abundance)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
    geom_errorbar(aes(ymin = mean_abundance - sd_abundance,
                      ymax = mean_abundance + sd_abundance),
                  width = 0.2) +
    labs(title = "Top 5 Most Abundant Glycan Compositions",
         x = "Glycan Composition",
         y = "Mean Relative Abundance (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))
  
  # Save the plot
  ggsave(file.path(figures_dir, "top_5_glycan_compositions.png"), top_5_plot,
         width = 8, height = 6, dpi = 300)
  
  cat("\nTop 5 glycan compositions plot saved to:", file.path(figures_dir, "top_5_glycan_compositions.png"), "\n")
  
  # Calculate relative abundance by glycan class and sample
  relative_abundance_class <- data %>%
    group_by(sample, glycan_class) %>%
    summarise(
      total_abundance = sum(abundance, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    group_by(sample) %>%
    mutate(
      sample_total = sum(total_abundance, na.rm = TRUE),
      relative_abundance = (total_abundance / sample_total) * 100
    ) %>%
    ungroup()
  
  # Create boxplot for glycan classes
  glycan_class_plot <- ggplot(relative_abundance_class, aes(x = glycan_class, y = relative_abundance)) +
    geom_boxplot(fill = "lightblue", alpha = 0.7, outlier.shape = 1) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
    labs(title = "Relative Abundance by Glycan Class",
         x = "Glycan Class",
         y = "Relative Abundance (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))
  
  # Save the boxplot
  ggsave(file.path(figures_dir, "glycan_class_relative_abundance_boxplot.png"), glycan_class_plot, 
         width = 10, height = 8, dpi = 300)
  
  cat("Glycan class boxplot saved to:", file.path(figures_dir, "glycan_class_relative_abundance_boxplot.png"), "\n")
  
  return(list(
    relative_abundance = relative_abundance,
    glycan_summary = glycan_summary,
    top_5_by_sample = top_5_by_sample,
    overall_top_5 = overall_top_5,
    glycan_ttests = glycan_ttests,
    significant_glycans = significant_glycans,
    high_fold_changes = high_fold_changes,
    relative_abundance_class = relative_abundance_class,
    plots = list(
      glycan_plot = glycan_plot,
      top_5_plot = top_5_plot,
      glycan_class_plot = glycan_class_plot
    )
  ))
}

#' Create a barplot for glycan compositions
#' 
#' @param glycan_summary Data frame with glycan composition summary statistics
#' @param title Plot title
#' @return ggplot object
plot_glycan_composition_barplot <- function(glycan_summary, title = "Relative Abundance of Glycan Compositions") {
  
  ggplot(glycan_summary, aes(x = reorder(glycan_composition, -mean_relative_abundance), 
                             y = mean_relative_abundance)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
    geom_errorbar(aes(ymin = mean_relative_abundance - se_relative_abundance, 
                      ymax = mean_relative_abundance + se_relative_abundance), 
                  width = 0.2, position = position_dodge(0.9)) +
    labs(title = title,
         x = "Glycan Composition",
         y = "Relative Abundance (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          plot.title = element_text(hjust = 0.5))
}

#' Analyze pseudo-glycomics composition
#' 
#' @param data Long format data with abundance, sample, glycan_composition, group columns
#' @param output_dir Directory to save output files
#' @return List containing pseudo-glycomics analysis results
analyze_pseudo_glycomics <- function(data, output_dir = "output_data/peptidegroups_intensity") {
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Add group classification based on sample names if not already present
  if (!"group" %in% colnames(data)) {
    data <- data %>%
      mutate(group = case_when(
        str_starts(sample, "hc") ~ "Healthy",
        str_starts(sample, "m") ~ "MECFS",
        TRUE ~ "Unknown"
      ))
  }
  
  pseudo_glycomics <- data %>%
    group_by(group) %>%
    # First calculate total abundance per group
    mutate(total_abundance = sum(abundance, na.rm = TRUE)) %>%
    # Then group by both group and glycan composition
    group_by(glycan_composition, group) %>%
    summarise(
      glycan_abundance = sum(abundance, na.rm = TRUE),
      total_group_abundance = first(total_abundance),
      relative_percentage = (glycan_abundance / total_group_abundance) * 100,
      .groups = 'drop'
    ) %>%
    arrange(group, desc(relative_percentage))
  
  # Print summary
  cat("\nRelative Glycan Composition Percentages by Group:\n")
  print(pseudo_glycomics)
  
  # Save pseudo-glycomics data
  write.csv(pseudo_glycomics, file.path(output_dir, "pseudo_glycomics_composition.csv"), row.names = FALSE)
  
  cat("\nPseudo-glycomics composition data saved to:", file.path(output_dir, "pseudo_glycomics_composition.csv"), "\n")
  
  # Calculate differences between groups
  pseudo_glycomics_wide <- pseudo_glycomics %>%
    select(glycan_composition, group, relative_percentage) %>%
    pivot_wider(names_from = group, values_from = relative_percentage, names_prefix = "percentage_") %>%
    mutate(
      difference_mecfs_healthy = percentage_MECFS - percentage_Healthy,
      fold_change = percentage_MECFS / percentage_Healthy,
      log2_fold_change = log2(fold_change)
    ) %>%
    arrange(desc(abs(difference_mecfs_healthy)))
  
  # Save differences
  write.csv(pseudo_glycomics_wide, file.path(output_dir, "pseudo_glycomics_group_differences.csv"), row.names = FALSE)
  
  cat("Group differences saved to:", file.path(output_dir, "pseudo_glycomics_group_differences.csv"), "\n")
  
  # Print top differences
  cat("\nTop 10 glycan compositions with largest differences between MECFS and Healthy:\n")
  top_differences <- pseudo_glycomics_wide %>%
    head(10)
  
  for(i in 1:nrow(top_differences)) {
    cat(sprintf("\n%d. %s\n", i, top_differences$glycan_composition[i]))
    cat(sprintf("   MECFS: %.2f%%\n", top_differences$percentage_MECFS[i]))
    cat(sprintf("   Healthy: %.2f%%\n", top_differences$percentage_Healthy[i]))
    cat(sprintf("   Difference: %.2f%%\n", top_differences$difference_mecfs_healthy[i]))
    cat(sprintf("   Fold change: %.2f\n", top_differences$fold_change[i]))
  }
  
  return(list(
    pseudo_glycomics = pseudo_glycomics,
    pseudo_glycomics_wide = pseudo_glycomics_wide,
    top_differences = top_differences
  ))
}

#' Create comprehensive glycan analysis plots
#' 
#' @param data Long format data with abundance, sample, glycan_composition columns
#' @param figures_dir Directory to save figures
#' @return List of ggplot objects
create_glycan_analysis_plots <- function(data, figures_dir = "figures") {
  
  # Create figures directory if it doesn't exist
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  
  plots <- list()
  
  # Create histogram of peptide scores
  score_hist <- ggplot(data, aes(x = log2(pep_2d_by_search_engine_a2_pmi_byonic))) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, color = "black") +
    labs(title = "Distribution of log2 Peptide Pep2d Scores",
         x = "log2 Peptide Pep2d Score", 
         y = "Count") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Save the histogram
  ggsave(file.path(figures_dir, "peptidegroups_intensity/peptide_log2_pep2d_histogram.png"), 
         score_hist, width = 8, height = 6, dpi = 300)
  
  cat("Peptide scores histogram saved to:", file.path(figures_dir, "peptidegroups_intensity/peptide_log2_pep2d_histogram.png"), "\n")
  
  # Create boxplot of pep2d scores by sample
  pep2d_sample_plot <- ggplot(data, 
                              aes(x = sample, y = log2(pep_2d_by_search_engine_a2_pmi_byonic))) +
    geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.shape = 1) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
    labs(title = "Peptide log2 Pep2d Score Distribution by Sample",
         x = "Sample", 
         y = "log2 Pep2d Score") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
  
  # Save the pep2d score boxplot
  ggsave(file.path(figures_dir, "peptidegroups_intensity/log2_pep2d_score_by_sample_boxplot.png"),
         pep2d_sample_plot, width = 12, height = 8, dpi = 300)
  
  cat("Pep2d scores by sample boxplot saved to:", file.path(figures_dir, "peptidegroups_intensity/log2_pep2d_score_by_sample_boxplot.png"), "\n")
  
  # Create boxplot of abundances by sample
  sample_abundance_plot <- ggplot(data, aes(x = sample, y = log2(abundance))) +
    geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.shape = 1) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
    labs(title = "Abundance log2 Distribution by Sample",
         x = "Sample",
         y = "Abundance log2") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
  
  # Save the sample abundance boxplot
  ggsave(file.path(figures_dir, "peptidegroups_intensity/sample_abundance_boxplot.png"), 
         sample_abundance_plot, width = 12, height = 8, dpi = 300)
  
  cat("Sample abundance boxplot saved to:", file.path(figures_dir, "peptidegroups_intensity/sample_abundance_boxplot.png"), "\n")
  
  plots$score_histogram <- score_hist
  plots$pep2d_sample_plot <- pep2d_sample_plot
  plots$sample_abundance_plot <- sample_abundance_plot
  
  return(plots)
}

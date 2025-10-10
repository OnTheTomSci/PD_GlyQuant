# Peptide Groups Modification Statistics Functions
# This module contains functions for analyzing fucosylation and sialylation patterns

library(tidyverse)
library(ggplot2)
library(effectsize)

#' Analyze sample-level fucosylation with statistical tests
#' 
#' @param data Long format data with abundance, sample, group, contains_Fuc columns
#' @param output_dir Directory to save output files
#' @param figures_dir Directory to save figures
#' @return List containing fucosylation analysis results
analyze_sample_fucosylation <- function(data, output_dir = "output_data/peptidegroups_intensity", 
                                       figures_dir = "figures/peptidegroups_intensity") {
  
  # Create output directories if they don't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Calculate fucosylation percentages by sample
  fucosylation_by_sample <- data %>%
    # Group by sample and group to get total abundance per sample
    group_by(sample, group) %>%
    mutate(total_sample_abundance = sum(abundance, na.rm = TRUE)) %>%
    # Calculate Fuc percentage per sample using existing contains_Fuc column
    summarise(
      fuc_abundance = sum(abundance[contains_Fuc == TRUE], na.rm = TRUE),
      total_abundance = first(total_sample_abundance),
      fuc_percentage = (fuc_abundance / total_abundance) * 100,
      .groups = 'drop'
    )
  
  # Perform t-test
  t_test_result <- t.test(fuc_percentage ~ group, data = fucosylation_by_sample)
  
  # Calculate F-test for variance (homogeneity of variance test)
  healthy_var <- var(fucosylation_by_sample$fuc_percentage[fucosylation_by_sample$group == "Healthy"])
  mecfs_var <- var(fucosylation_by_sample$fuc_percentage[fucosylation_by_sample$group == "MECFS"])
  
  # F-test for variance (larger variance in numerator)
  if (healthy_var >= mecfs_var) {
    f_ratio <- healthy_var / mecfs_var
    df1 <- sum(fucosylation_by_sample$group == "Healthy") - 1
    df2 <- sum(fucosylation_by_sample$group == "MECFS") - 1
  } else {
    f_ratio <- mecfs_var / healthy_var
    df1 <- sum(fucosylation_by_sample$group == "MECFS") - 1
    df2 <- sum(fucosylation_by_sample$group == "Healthy") - 1
  }
  
  f_p_value <- 2 * (1 - pf(f_ratio, df1, df2))  # Two-tailed test
  
  # Calculate Cohen's d effect size
  cohens_d_result <- cohens_d(fuc_percentage ~ group, data = fucosylation_by_sample)
  cohens_d_value <- cohens_d_result$Cohens_d
  
  # Create summary statistics
  fucosylation_summary <- fucosylation_by_sample %>%
    group_by(group) %>%
    summarise(
      n = n(),
      mean_fuc_percentage = mean(fuc_percentage, na.rm = TRUE),
      sd_fuc_percentage = sd(fuc_percentage, na.rm = TRUE),
      se_fuc_percentage = sd_fuc_percentage / sqrt(n),
      ci_lower = mean_fuc_percentage - (1.96 * se_fuc_percentage),
      ci_upper = mean_fuc_percentage + (1.96 * se_fuc_percentage),
      median_fuc = median(fuc_percentage),
      min_fuc = min(fuc_percentage),
      max_fuc = max(fuc_percentage),
      .groups = 'drop'
    )
  
  # Create plots
  plots <- create_modification_boxplot(fucosylation_by_sample, fucosylation_summary,
                                      t_test_result, f_ratio, f_p_value, cohens_d_value,
                                      modification_type = "Fucosylation",
                                      y_var = "fuc_percentage",
                                      output_path = file.path(figures_dir, "fucosylation_comparison.png"))
  
  # Save detailed summary
  write.csv(fucosylation_summary, file.path(output_dir, "fucosylation_detailed_summary_with_se.csv"), row.names = FALSE)
  
  # Print summary
  cat("\nFucosylation Analysis Summary:\n")
  cat("T-test p-value:", t_test_result$p.value, "\n")
  cat("F-test p-value:", f_p_value, "\n")
  cat("Cohen's d:", cohens_d_value, "\n")
  
  # Interpret Cohen's d
  cohens_d_interpretation <- case_when(
    abs(cohens_d_value) < 0.2 ~ "negligible",
    abs(cohens_d_value) < 0.5 ~ "small",
    abs(cohens_d_value) < 0.8 ~ "medium",
    TRUE ~ "large"
  )
  
  cat("Effect size interpretation:", cohens_d_interpretation, "effect\n")
  
  return(list(
    fucosylation_by_sample = fucosylation_by_sample,
    fucosylation_summary = fucosylation_summary,
    t_test_result = t_test_result,
    f_test = list(f_ratio = f_ratio, p_value = f_p_value, df1 = df1, df2 = df2),
    cohens_d_result = cohens_d_result,
    cohens_d_interpretation = cohens_d_interpretation,
    plots = plots
  ))
}

#' Analyze sample-level sialylation with statistical tests
#' 
#' @param data Long format data with abundance, sample, group, contains_NeuAc columns
#' @param output_dir Directory to save output files
#' @param figures_dir Directory to save figures
#' @return List containing sialylation analysis results
analyze_sample_sialylation <- function(data, output_dir = "output_data/peptidegroups_intensity", 
                                      figures_dir = "figures/peptidegroups_intensity") {
  
  # Create output directories if they don't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Calculate NeuAc percentages by sample
  neuac_by_sample <- data %>%
    # Group by sample and group to get total abundance per sample
    group_by(sample, group) %>%
    mutate(total_sample_abundance = sum(abundance, na.rm = TRUE)) %>%
    # Calculate NeuAc percentage per sample using existing contains_NeuAc column
    summarise(
      neuac_abundance = sum(abundance[contains_NeuAc == TRUE], na.rm = TRUE),
      total_abundance = first(total_sample_abundance),
      neuac_percentage = (neuac_abundance / total_abundance) * 100,
      .groups = 'drop'
    )
  
  # Perform t-test
  t_test_result <- t.test(neuac_percentage ~ group, data = neuac_by_sample)
  
  # Calculate F-test for variance (homogeneity of variance test)
  healthy_var <- var(neuac_by_sample$neuac_percentage[neuac_by_sample$group == "Healthy"])
  mecfs_var <- var(neuac_by_sample$neuac_percentage[neuac_by_sample$group == "MECFS"])
  
  # F-test for variance (larger variance in numerator)
  if (healthy_var >= mecfs_var) {
    f_ratio <- healthy_var / mecfs_var
    df1 <- sum(neuac_by_sample$group == "Healthy") - 1
    df2 <- sum(neuac_by_sample$group == "MECFS") - 1
  } else {
    f_ratio <- mecfs_var / healthy_var
    df1 <- sum(neuac_by_sample$group == "MECFS") - 1
    df2 <- sum(neuac_by_sample$group == "Healthy") - 1
  }
  
  f_p_value <- 2 * (1 - pf(f_ratio, df1, df2))  # Two-tailed test
  
  # Calculate Cohen's d effect size
  cohens_d_result <- cohens_d(neuac_percentage ~ group, data = neuac_by_sample)
  cohens_d_value <- cohens_d_result$Cohens_d
  
  # Create summary statistics
  neuac_summary <- neuac_by_sample %>%
    group_by(group) %>%
    summarise(
      n = n(),
      mean_neuac_percentage = mean(neuac_percentage, na.rm = TRUE),
      sd_neuac_percentage = sd(neuac_percentage, na.rm = TRUE),
      se_neuac_percentage = sd_neuac_percentage / sqrt(n),
      ci_lower = mean_neuac_percentage - (1.96 * se_neuac_percentage),
      ci_upper = mean_neuac_percentage + (1.96 * se_neuac_percentage),
      median_neuac = median(neuac_percentage),
      min_neuac = min(neuac_percentage),
      max_neuac = max(neuac_percentage),
      .groups = 'drop'
    )
  
  # Create plots
  plots <- create_modification_boxplot(neuac_by_sample, neuac_summary,
                                      t_test_result, f_ratio, f_p_value, cohens_d_value,
                                      modification_type = "Sialylation (NeuAc)",
                                      y_var = "neuac_percentage",
                                      output_path = file.path(figures_dir, "sialylation_comparison.png"))
  
  # Save detailed summary
  write.csv(neuac_summary, file.path(output_dir, "sialylation_detailed_summary_with_se.csv"), row.names = FALSE)
  
  # Print summary
  cat("\nSialylation Analysis Summary:\n")
  cat("T-test p-value:", t_test_result$p.value, "\n")
  cat("F-test p-value:", f_p_value, "\n")
  cat("Cohen's d:", cohens_d_value, "\n")
  
  # Interpret Cohen's d
  cohens_d_interpretation <- case_when(
    abs(cohens_d_value) < 0.2 ~ "negligible",
    abs(cohens_d_value) < 0.5 ~ "small",
    abs(cohens_d_value) < 0.8 ~ "medium",
    TRUE ~ "large"
  )
  
  cat("Effect size interpretation:", cohens_d_interpretation, "effect\n")
  
  return(list(
    neuac_by_sample = neuac_by_sample,
    neuac_summary = neuac_summary,
    t_test_result = t_test_result,
    f_test = list(f_ratio = f_ratio, p_value = f_p_value, df1 = df1, df2 = df2),
    cohens_d_result = cohens_d_result,
    cohens_d_interpretation = cohens_d_interpretation,
    plots = plots
  ))
}

#' Create reusable boxplot with statistical annotations
#' 
#' @param data Data frame with group and percentage columns
#' @param summary_data Summary statistics data frame
#' @param t_test_result Result from t.test()
#' @param f_ratio F-ratio from F-test
#' @param f_p_value P-value from F-test
#' @param cohens_d_value Cohen's d value
#' @param modification_type Type of modification (e.g., "Fucosylation", "Sialylation")
#' @param y_var Name of the y-variable column
#' @param output_path Path to save the plot
#' @return List of ggplot objects
create_modification_boxplot <- function(data, summary_data, t_test_result, f_ratio, f_p_value, 
                                       cohens_d_value, modification_type, y_var, output_path) {
  
  # Create labels for the plot
  t_test_label <- sprintf("t-test: p = %.3g", t_test_result$p.value)
  f_test_label <- sprintf("F-var = %.3f, p = %.3g", f_ratio, f_p_value)
  cohens_d_label <- sprintf("Cohen's d = %.3f", cohens_d_value)
  
  # Get the y values for positioning
  y_values <- data[[y_var]]
  max_y <- max(y_values, na.rm = TRUE)
  
  # Create boxplot with individual points
  boxplot <- ggplot(data, aes_string(x = "group", y = y_var, fill = "group")) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
    scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    labs(
      title = paste(modification_type, "Percentage by Group"),
      x = "Group",
      y = paste(modification_type, "Percentage (%)")
    ) +
    # Add statistical annotations
    annotate("text", 
             x = 1.5, 
             y = max_y * 1.15,
             label = t_test_label,
             size = 3.5) +
    annotate("text", 
             x = 1.5, 
             y = max_y * 1.05,
             label = f_test_label,
             size = 3.5) +
    annotate("text", 
             x = 1.5, 
             y = max_y * 0.95,
             label = cohens_d_label,
             size = 3.5)
  
  # Create bar plot with error bars (mean ± SE)
  mean_col <- paste0("mean_", str_remove(y_var, "_percentage"), "_percentage")
  se_col <- paste0("se_", str_remove(y_var, "_percentage"), "_percentage")
  
  barplot_se <- ggplot(summary_data, aes_string(x = "group", y = mean_col, fill = "group")) +
    geom_col(alpha = 0.7, width = 0.6) +
    geom_errorbar(aes_string(ymin = paste0(mean_col, " - ", se_col),
                             ymax = paste0(mean_col, " + ", se_col)),
                  width = 0.2, size = 1, color = "black") +
    scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    labs(
      title = paste(modification_type, "Percentage by Group (Mean ± SE)"),
      x = "Group",
      y = paste(modification_type, "Percentage (%)")
    ) +
    # Add sample size annotations
    geom_text(aes_string(label = "paste('n =', n)"), 
              vjust = -0.5, size = 3.5) +
    # Add mean value annotations
    geom_text(aes_string(label = paste0("sprintf('Mean = %.2f%%', ", mean_col, ")")), 
              vjust = -1.5, size = 3.5) +
    # Add SE annotations
    geom_text(aes_string(label = paste0("sprintf('SE = %.2f%%', ", se_col, ")")), 
              vjust = -2.5, size = 3.5)
  
  # Create bar plot with 95% confidence intervals
  barplot_ci <- ggplot(summary_data, aes_string(x = "group", y = mean_col, fill = "group")) +
    geom_col(alpha = 0.7, width = 0.6) +
    geom_errorbar(aes_string(ymin = "ci_lower", ymax = "ci_upper"),
                  width = 0.2, size = 1, color = "black") +
    scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    labs(
      title = paste(modification_type, "Percentage by Group (Mean ± 95% CI)"),
      x = "Group",
      y = paste(modification_type, "Percentage (%)")
    ) +
    # Add sample size annotations
    geom_text(aes_string(label = "paste('n =', n)"), 
              vjust = -0.5, size = 3.5) +
    # Add mean value annotations
    geom_text(aes_string(label = paste0("sprintf('Mean = %.2f%%', ", mean_col, ")")), 
              vjust = -1.5, size = 3.5) +
    # Add CI annotations
    geom_text(aes(label = sprintf("95%% CI = [%.2f, %.2f]", ci_lower, ci_upper)), 
              vjust = -2.5, size = 3)
  
  # Save plots
  ggsave(str_replace(output_path, "\\.png$", ".png"), boxplot, width = 8, height = 6, dpi = 300)
  ggsave(str_replace(output_path, "\\.png$", "_mean_se.png"), barplot_se, width = 8, height = 6, dpi = 300)
  ggsave(str_replace(output_path, "\\.png$", "_mean_ci.png"), barplot_ci, width = 8, height = 6, dpi = 300)
  
  cat("Plots saved:\n")
  cat("-", output_path, "\n")
  cat("-", str_replace(output_path, "\\.png$", "_mean_se.png"), "\n")
  cat("-", str_replace(output_path, "\\.png$", "_mean_ci.png"), "\n")
  
  return(list(
    boxplot = boxplot,
    barplot_se = barplot_se,
    barplot_ci = barplot_ci
  ))
}

#' Analyze glycan class relative abundance by sample
#' 
#' @param data Long format data with abundance, sample, group, glycan_class columns
#' @param output_dir Directory to save output files
#' @param figures_dir Directory to save figures
#' @return List containing glycan class analysis results
analyze_glycan_class_by_sample <- function(data, output_dir = "output_data/peptidegroups_intensity", 
                                          figures_dir = "figures") {
  
  # Create output directories if they don't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Calculate relative abundance of glycan classes by sample
  glycan_class_by_sample <- data %>%
    # Group by sample, group, and glycan_class to get abundances
    group_by(sample, group, glycan_class) %>%
    summarise(
      group_abundance = sum(abundance, na.rm = TRUE),
      .groups = 'keep'
    ) %>%
    # Calculate total abundance per sample for percentage
    group_by(sample, group) %>%
    mutate(
      total_sample_abundance = sum(group_abundance, na.rm = TRUE),
      class_percentage = (group_abundance / total_sample_abundance) * 100
    ) %>%
    ungroup()
  
  # Create boxplot for glycan classes
  glycan_class_plot <- ggplot(glycan_class_by_sample, 
                              aes(x = glycan_class, y = class_percentage, fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 1, position = position_dodge(0.8)) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.5, position = position_dodge(0.8)) +
    scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
    labs(title = "Relative Abundance of Glycan Classes by Group",
         x = "Glycan Class",
         y = "Relative Abundance (%)",
         fill = "Group") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
  
  # Save the boxplot
  ggsave(file.path(figures_dir, "glycan_class_by_group_boxplot.png"), glycan_class_plot, 
         width = 12, height = 8, dpi = 300)
  
  cat("Glycan class boxplot saved to:", file.path(figures_dir, "glycan_class_by_group_boxplot.png"), "\n")
  
  # Calculate summary statistics by glycan class and group
  glycan_class_summary <- glycan_class_by_sample %>%
    group_by(glycan_class, group) %>%
    summarise(
      n = n(),
      mean_percentage = mean(class_percentage, na.rm = TRUE),
      sd_percentage = sd(class_percentage, na.rm = TRUE),
      se_percentage = sd_percentage / sqrt(n),
      median_percentage = median(class_percentage, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Save summary statistics
  write.csv(glycan_class_summary, file.path(output_dir, "glycan_class_by_group_summary.csv"), row.names = FALSE)
  
  return(list(
    glycan_class_by_sample = glycan_class_by_sample,
    glycan_class_summary = glycan_class_summary,
    plot = glycan_class_plot
  ))
}

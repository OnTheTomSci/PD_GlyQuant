# Peptide Groups Modification Statistics Functions
# This module contains functions for analyzing fucosylation and sialylation patterns

library(tidyverse)
library(ggplot2)
library(effectsize)
library(effsize)


#' Analyze sample-level fucosylation with statistical tests
#' 
#' @param data Long format data with abundance, sample, group, contains_Fuc columns
#' @param output_dir Directory to save output files
#' @param figures_dir Directory to save figures
#' @return List containing fucosylation analysis results
analyze_sample_fucosylation <- function(data, output_dir = "output_data/peptidegroups_intensity/puesdo_glycomics", 
                                       figures_dir = "figures/peptidegroups_intensity/puesdo_glycomics") {
  
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
  
  # Interpret Cohen's d
  cohens_d_interpretation <- case_when(
    abs(cohens_d_value) < 0.2 ~ "negligible",
    abs(cohens_d_value) < 0.5 ~ "small",
    abs(cohens_d_value) < 0.8 ~ "medium",
    TRUE ~ "large"
  )
  
  # Create comprehensive statistical test results
  # Handle potential missing or differently structured cohens_d_result fields
  cohens_d_p_value <- if("p_value" %in% names(cohens_d_result)) cohens_d_result$p_value else NA
  cohens_d_ci_low <- if("CI_low" %in% names(cohens_d_result)) cohens_d_result$CI_low else NA
  cohens_d_ci_high <- if("CI_high" %in% names(cohens_d_result)) cohens_d_result$CI_high else NA
  
  statistical_results <- data.frame(
    test_type = c("t_test", "f_test", "cohens_d"),
    statistic_name = c("t_statistic", "f_ratio", "cohens_d"),
    statistic_value = c(t_test_result$statistic, f_ratio, cohens_d_value),
    p_value = c(t_test_result$p.value, f_p_value, cohens_d_p_value),
    df1 = c(t_test_result$parameter, df1, NA),
    df2 = c(NA, df2, NA),
    confidence_interval_lower = c(t_test_result$conf.int[1], NA, cohens_d_ci_low),
    confidence_interval_upper = c(t_test_result$conf.int[2], NA, cohens_d_ci_high),
    effect_size_interpretation = c(NA, NA, cohens_d_interpretation),
    description = c("Independent samples t-test comparing fucosylation percentages between groups",
                   "F-test for homogeneity of variance",
                   "Cohen's d effect size measure")
  )
  
  # Save detailed summary
  write.csv(fucosylation_summary, file.path(output_dir, "fucosylation_detailed_summary_with_se.csv"), row.names = FALSE)
  
  # Save comprehensive statistical test results
  write.csv(statistical_results, file.path(output_dir, "fucosylation_statistical_tests_results.csv"), row.names = FALSE)
  
  # Save raw sample-level data
  write.csv(fucosylation_by_sample, file.path(output_dir, "fucosylation_by_sample_data.csv"), row.names = FALSE)
  
  # Print summary
  cat("\nFucosylation Analysis Summary:\n")
  cat("T-test p-value:", t_test_result$p.value, "\n")
  cat("F-test p-value:", f_p_value, "\n")
  cat("Cohen's d:", cohens_d_value, "\n")
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
  
  # Interpret Cohen's d
  cohens_d_interpretation <- case_when(
    abs(cohens_d_value) < 0.2 ~ "negligible",
    abs(cohens_d_value) < 0.5 ~ "small",
    abs(cohens_d_value) < 0.8 ~ "medium",
    TRUE ~ "large"
  )
  
  # Create comprehensive statistical test results
  # Handle potential missing or differently structured cohens_d_result fields
  cohens_d_p_value <- if("p_value" %in% names(cohens_d_result)) cohens_d_result$p_value else NA
  cohens_d_ci_low <- if("CI_low" %in% names(cohens_d_result)) cohens_d_result$CI_low else NA
  cohens_d_ci_high <- if("CI_high" %in% names(cohens_d_result)) cohens_d_result$CI_high else NA
  
  statistical_results <- data.frame(
    test_type = c("t_test", "f_test", "cohens_d"),
    statistic_name = c("t_statistic", "f_ratio", "cohens_d"),
    statistic_value = c(t_test_result$statistic, f_ratio, cohens_d_value),
    p_value = c(t_test_result$p.value, f_p_value, cohens_d_p_value),
    df1 = c(t_test_result$parameter, df1, NA),
    df2 = c(NA, df2, NA),
    confidence_interval_lower = c(t_test_result$conf.int[1], NA, cohens_d_ci_low),
    confidence_interval_upper = c(t_test_result$conf.int[2], NA, cohens_d_ci_high),
    effect_size_interpretation = c(NA, NA, cohens_d_interpretation),
    description = c("Independent samples t-test comparing sialylation percentages between groups",
                   "F-test for homogeneity of variance",
                   "Cohen's d effect size measure")
  )
  
  # Save detailed summary
  write.csv(neuac_summary, file.path(output_dir, "sialylation_detailed_summary_with_se.csv"), row.names = FALSE)
  
  # Save comprehensive statistical test results
  write.csv(statistical_results, file.path(output_dir, "sialylation_statistical_tests_results.csv"), row.names = FALSE)
  
  # Save raw sample-level data
  write.csv(neuac_by_sample, file.path(output_dir, "sialylation_by_sample_data.csv"), row.names = FALSE)
  
  # Print summary
  cat("\nSialylation Analysis Summary:\n")
  cat("T-test p-value:", t_test_result$p.value, "\n")
  cat("F-test p-value:", f_p_value, "\n")
  cat("Cohen's d:", cohens_d_value, "\n")
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

#' Analyze glycan class statistical differences between groups
#' 
#' @param data Long format data with abundance, sample, group, glycan_class columns
#' @param output_dir Directory to save output files
#' @param figures_dir Directory to save figures
#' @return List containing glycan class analysis results
analyze_glycan_class_statistics <- function(data, output_dir = "output_data/peptidegroups_intensity", 
                                           figures_dir = "figures/peptidegroups_intensity") {
  
  # Create output directories if they don't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Calculate glycan class percentages by sample
  glycan_class_by_sample <- data %>%
    # Group by sample and group to get total abundance per sample
    group_by(sample, group) %>%
    mutate(total_sample_abundance = sum(abundance, na.rm = TRUE)) %>%
    # Group by glycan class and calculate percentages
    group_by(sample, group, glycan_class) %>%
    summarise(
      class_abundance = sum(abundance, na.rm = TRUE),
      total_abundance = first(total_sample_abundance),
      class_percentage = (class_abundance / total_abundance) * 100,
      .groups = 'drop'
    )
  
  # Perform statistical tests for each glycan class
  glycan_class_stats <- glycan_class_by_sample %>%
    group_by(glycan_class) %>%
    summarise(
      # Sample sizes
      n_healthy = sum(group == "Healthy"),
      n_mecfs = sum(group == "MECFS"),
      
      # Descriptive statistics
      mean_healthy = mean(class_percentage[group == "Healthy"], na.rm = TRUE),
      mean_mecfs = mean(class_percentage[group == "MECFS"], na.rm = TRUE),
      sd_healthy = sd(class_percentage[group == "Healthy"], na.rm = TRUE),
      sd_mecfs = sd(class_percentage[group == "MECFS"], na.rm = TRUE),
      se_healthy = sd_healthy / sqrt(n_healthy),
      se_mecfs = sd_mecfs / sqrt(n_mecfs),
      
      # Statistical tests
      t_statistic = tryCatch({
        t.test(class_percentage ~ group)$statistic
      }, error = function(e) NA_real_),
      t_p_value = tryCatch({
        t.test(class_percentage ~ group)$p.value
      }, error = function(e) NA_real_),
      
      # F-test for variance
      f_statistic = tryCatch({
        if(sd_healthy > 0 && sd_mecfs > 0) {
          var.test(class_percentage ~ group)$statistic
        } else NA_real_
      }, error = function(e) NA_real_),
      f_p_value = tryCatch({
        if(sd_healthy > 0 && sd_mecfs > 0) {
          var.test(class_percentage ~ group)$p.value
        } else NA_real_
      }, error = function(e) NA_real_),
      
      # Cohen's d effect size
      cohens_d = tryCatch({
        if(sd_healthy > 0 && sd_mecfs > 0) {
          cohens_d_result <- cohens_d(class_percentage ~ group)
          cohens_d_result$Cohens_d
        } else NA_real_
      }, error = function(e) NA_real_),
      
      .groups = 'drop'
    ) %>%
    mutate(
      # Calculate additional statistics
      fold_change = mean_mecfs / mean_healthy,
      log2_fold_change = log2(fold_change),
      mean_difference = mean_mecfs - mean_healthy,
      
      # Effect size interpretation
      effect_size_interpretation = case_when(
        is.na(cohens_d) ~ "Cannot calculate",
        abs(cohens_d) < 0.2 ~ "negligible",
        abs(cohens_d) < 0.5 ~ "small",
        abs(cohens_d) < 0.8 ~ "medium",
        TRUE ~ "large"
      ),
      
      # Significance flags
      significant_t_test = t_p_value < 0.05 & !is.na(t_p_value),
      significant_f_test = f_p_value < 0.05 & !is.na(f_p_value),
      high_fold_change = abs(log2_fold_change) > 1,  # >2-fold or <0.5-fold
      
      # Confidence intervals (95%)
      ci_lower_healthy = mean_healthy - (1.96 * se_healthy),
      ci_upper_healthy = mean_healthy + (1.96 * se_healthy),
      ci_lower_mecfs = mean_mecfs - (1.96 * se_mecfs),
      ci_upper_mecfs = mean_mecfs + (1.96 * se_mecfs)
    ) %>%
    arrange(t_p_value)
  
  # Apply multiple testing correction
  glycan_class_stats$t_p_value_adj <- p.adjust(glycan_class_stats$t_p_value, method = "BH")
  glycan_class_stats$significant_t_test_adj <- glycan_class_stats$t_p_value_adj < 0.05
  
  # Create summary statistics by group
  group_summary <- glycan_class_by_sample %>%
    group_by(group, glycan_class) %>%
    summarise(
      n = n(),
      mean_percentage = mean(class_percentage, na.rm = TRUE),
      sd_percentage = sd(class_percentage, na.rm = TRUE),
      se_percentage = sd_percentage / sqrt(n),
      median_percentage = median(class_percentage),
      min_percentage = min(class_percentage),
      max_percentage = max(class_percentage),
      .groups = 'drop'
    ) %>%
    pivot_wider(names_from = group, 
                values_from = c(n, mean_percentage, sd_percentage, se_percentage, 
                               median_percentage, min_percentage, max_percentage),
                names_sep = "_")
  
  # Create comprehensive results table
  comprehensive_results <- glycan_class_stats %>%
    left_join(group_summary, by = "glycan_class") %>%
    select(
      glycan_class,
      n_healthy, n_mecfs,
      mean_healthy, mean_mecfs,
      sd_healthy, sd_mecfs,
      se_healthy, se_mecfs,
      median_percentage_Healthy, median_percentage_MECFS,
      min_percentage_Healthy, min_percentage_MECFS,
      max_percentage_Healthy, max_percentage_MECFS,
      ci_lower_healthy, ci_upper_healthy,
      ci_lower_mecfs, ci_upper_mecfs,
      mean_difference, fold_change, log2_fold_change,
      t_statistic, t_p_value,
      f_statistic, f_p_value,
      cohens_d, effect_size_interpretation,
      significant_t_test, significant_t_test_adj, significant_f_test, high_fold_change
    )
  
  # Save all results
  write.csv(glycan_class_by_sample, file.path(output_dir, "glycan_class_by_sample_data.csv"), row.names = FALSE)
  write.csv(glycan_class_stats, file.path(output_dir, "glycan_class_statistical_tests.csv"), row.names = FALSE)
  write.csv(comprehensive_results, file.path(output_dir, "glycan_class_comprehensive_results.csv"), row.names = FALSE)
  
  # Create plots
  plots <- create_glycan_class_plots(glycan_class_by_sample, glycan_class_stats, figures_dir)
  
  # Print summary
  cat("\nGlycan Class Analysis Summary:\n")
  cat("Total glycan classes analyzed:", nrow(glycan_class_stats), "\n")
  cat("Significant differences (p < 0.05):", sum(glycan_class_stats$significant_t_test, na.rm = TRUE), "\n")
  cat("Significant after multiple testing correction:", sum(glycan_class_stats$significant_t_test_adj, na.rm = TRUE), "\n")
  cat("High fold changes (>2-fold or <0.5-fold):", sum(glycan_class_stats$high_fold_change, na.rm = TRUE), "\n")
  
  # Print significant results
  significant_classes <- glycan_class_stats %>%
    filter(significant_t_test_adj) %>%
    arrange(t_p_value_adj)
  
  if(nrow(significant_classes) > 0) {
    cat("\nSignificantly different glycan classes (after multiple testing correction):\n")
    for(i in 1:nrow(significant_classes)) {
      cat(sprintf("\n%d. %s\n", i, significant_classes$glycan_class[i]))
      cat(sprintf("   Healthy: %.2f%% ± %.2f%% SD\n", 
                  significant_classes$mean_healthy[i],
                  significant_classes$sd_healthy[i]))
      cat(sprintf("   MECFS: %.2f%% ± %.2f%% SD\n", 
                  significant_classes$mean_mecfs[i],
                  significant_classes$sd_mecfs[i]))
      cat(sprintf("   Adjusted p-value: %.3e\n", significant_classes$t_p_value_adj[i]))
      cat(sprintf("   Fold change: %.2f\n", significant_classes$fold_change[i]))
      cat(sprintf("   Effect size: %s\n", significant_classes$effect_size_interpretation[i]))
    }
  }
  
  return(list(
    glycan_class_by_sample = glycan_class_by_sample,
    glycan_class_stats = glycan_class_stats,
    comprehensive_results = comprehensive_results,
    significant_classes = significant_classes,
    plots = plots
  ))
}

#' Create plots for glycan class analysis
#' 
#' @param glycan_class_by_sample Sample-level data
#' @param glycan_class_stats Statistical results
#' @param figures_dir Directory to save figures
#' @return List of ggplot objects
create_glycan_class_plots <- function(glycan_class_by_sample, glycan_class_stats, figures_dir) {
  
  # Create boxplot
  boxplot <- ggplot(glycan_class_by_sample, aes(x = glycan_class, y = class_percentage, fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 1) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
    scale_fill_manual(values = c("Healthy" = "#2E8B57", "MECFS" = "#DC143C")) +
    labs(title = "Glycan Class Distribution by Group",
         x = "Glycan Class",
         y = "Relative Abundance (%)",
         fill = "Group") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))
  
  # Save boxplot
  ggsave(file.path(figures_dir, "glycan_class_comparison.png"), boxplot, 
         width = 12, height = 8, dpi = 300)
  
 
  return(list(
    boxplot = boxplot,
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
#' @param data Long format data with columns: abundance, sample, group, glycan_class
#' @param output_dir Directory to save output files (CSVs)
#' @param figures_dir Directory to save figures
#' @return List with per-sample table, summary table, between-group test table, and ggplot object
#' @details
#' - Between-group tests are computed per glycan_class.
#' - Welch's t-test is used by default (robust to unequal variances).
#' - F-test (var.test) reports variance ratio (s^2_Healthy / s^2_MECFS).
#' - Effect size is Hedges' g (bias-corrected Cohen's d).
#' - If a class has <2 observations in either group, stats are returned as NA.
analyze_glycan_class_by_sample <- function(data,
                                           output_dir = "output_data/peptidegroups_intensity",
                                           figures_dir = "figures") {

  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("purrr", quietly = TRUE)
  requireNamespace("broom", quietly = TRUE)
  if (!requireNamespace("effsize", quietly = TRUE)) {
    stop("Package 'effsize' is required for Hedges' g. Please install.packages('effsize').")
  }

  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(broom)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  glycan_class_by_sample <- data %>%
    group_by(sample, group, glycan_class) %>%
    summarise(
      group_abundance = sum(abundance, na.rm = TRUE),
      .groups = 'keep'
    ) %>%
    group_by(sample, group) %>%
    mutate(
      total_sample_abundance = sum(group_abundance, na.rm = TRUE),
      class_percentage = (group_abundance / total_sample_abundance) * 100
    ) %>%
    ungroup()

  glycan_class_plot <- ggplot(glycan_class_by_sample,
                              aes(x = glycan_class, y = class_percentage, fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 1, position = position_dodge(0.8)) +
    geom_jitter(alpha = 0.6, size = 1.5,
                position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2)) +
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

  ggsave(file.path(figures_dir, "glycan_class_by_group_boxplot.png"), glycan_class_plot,
         width = 12, height = 8, dpi = 300)

  message("Glycan class boxplot saved to: ",
          file.path(figures_dir, "glycan_class_by_group_boxplot.png"))

  glycan_class_summary <- glycan_class_by_sample %>%
    group_by(glycan_class, group) %>%
    summarise(
      n = dplyr::n(),
      mean_percentage = mean(class_percentage, na.rm = TRUE),
      sd_percentage   = sd(class_percentage, na.rm = TRUE),
      se_percentage   = sd_percentage / sqrt(n),
      median_percentage = median(class_percentage, na.rm = TRUE),
      .groups = 'drop'
    )

  write.csv(glycan_class_summary,
            file.path(output_dir, "glycan_class_by_group_summary.csv"),
            row.names = FALSE)

  compute_tests <- function(df) {
    df <- df %>% dplyr::filter(group %in% c("Healthy", "MECFS"))

    n_H <- sum(df$group == "Healthy")
    n_M <- sum(df$group == "MECFS")

    m_H <- mean(df$class_percentage[df$group == "Healthy"], na.rm = TRUE)
    m_M <- mean(df$class_percentage[df$group == "MECFS"], na.rm = TRUE)
    sd_H <- sd(df$class_percentage[df$group == "Healthy"], na.rm = TRUE)
    sd_M <- sd(df$class_percentage[df$group == "MECFS"], na.rm = TRUE)

    out <- list(
      n_Healthy = n_H, n_MECFS = n_M,
      mean_Healthy = m_H, mean_MECFS = m_M,
      sd_Healthy = sd_H, sd_MECFS = sd_M,
      t_stat = NA_real_, t_df = NA_real_, t_p_value = NA_real_,
      t_ci_low = NA_real_, t_ci_high = NA_real_,
      hedges_g = NA_real_,
      F_stat = NA_real_, F_num_df = NA_real_, F_den_df = NA_real_, F_p_value = NA_real_
    )

    if (n_H >= 2 && n_M >= 2 && all(is.finite(df$class_percentage))) {
      t_res <- tryCatch(
        t.test(class_percentage ~ group, data = df, var.equal = FALSE, conf.level = 0.95),
        error = function(e) NULL
      )
      if (!is.null(t_res)) {
        out$t_stat    <- unname(t_res$statistic)
        out$t_df      <- unname(t_res$parameter)
        out$t_p_value <- unname(t_res$p.value)
        out$t_ci_low  <- unname(t_res$conf.int[1])
        out$t_ci_high <- unname(t_res$conf.int[2])
      }

      if (is.finite(sd_H) && is.finite(sd_M) && sd_H > 0 && sd_M > 0) {
        f_res <- tryCatch(
          var.test(class_percentage ~ group, data = df),
          error = function(e) NULL
        )
        if (!is.null(f_res)) {
          out$F_stat    <- unname(f_res$statistic)    # s^2_H / s^2_M
          out$F_num_df  <- unname(f_res$parameter[1]) # df1
          out$F_den_df  <- unname(f_res$parameter[2]) # df2
          out$F_p_value <- unname(f_res$p.value)
        }
      }

      g_res <- tryCatch(
        effsize::cohen.d(class_percentage ~ group, data = df,
                         hedges.correction = TRUE, na.rm = TRUE),
        error = function(e) NULL
      )
      if (!is.null(g_res) && length(g_res$estimate) == 1) {
        out$hedges_g <- unname(g_res$estimate)
      }
    }

    tibble::as_tibble(out)
  }

  test_results <- glycan_class_by_sample %>%
    group_by(glycan_class) %>%
    group_modify(~ compute_tests(.x)) %>%
    ungroup() %>%
    mutate(direction = dplyr::case_when(
      is.finite(mean_Healthy) & is.finite(mean_MECFS) & (mean_MECFS > mean_Healthy) ~ "MECFS>Healthy",
      is.finite(mean_Healthy) & is.finite(mean_MECFS) & (mean_MECFS < mean_Healthy) ~ "Healthy>MECFS",
      TRUE ~ NA_character_
    ))

  out_tests_path <- file.path(output_dir, "glycan_class_between_group_tests.csv")
  write.csv(test_results, out_tests_path, row.names = FALSE)
  message("Between-group test results written to: ", out_tests_path)

  list(
    glycan_class_by_sample = glycan_class_by_sample,
    glycan_class_summary   = glycan_class_summary,
    glycan_class_tests     = test_results,
    plot = glycan_class_plot
  )
}

analyze_sample_fucosylation(data = glyco_peptide_groups_long, output_dir = "output_data/peptidegroups_intensity/puesdo_glycomics", figures_dir = "figures/peptidegroups_intensity/puesdo_glycomics")
analyze_sample_sialylation(data = glyco_peptide_groups_long, output_dir = "output_data/peptidegroups_intensity/puesdo_glycomics", figures_dir = "figures/peptidegroups_intensity/puesdo_glycomics")
analyze_glycan_class_by_sample(data = glyco_peptide_groups_long, output_dir = "output_data/peptidegroups_intensity/puesdo_glycomics", figures_dir = "figures/peptidegroups_intensity/puesdo_glycomics")
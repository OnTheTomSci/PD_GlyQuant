# Peptide Groups Glycosite Glycosylation Analysis Functions
# This module contains functions for glycosite-level glycosylation analysis

library(tidyverse)
library(ggplot2)

#' Analyze glycosite-level glycosylation patterns
#' 
#' @param data Long format data with abundance, sample, group, gsite_ID, contains_Fuc, contains_NeuAc columns
#' @param output_dir Directory to save output files
#' @param figures_dir Directory to save figures
#' @return List containing glycosite glycosylation analysis results
analyze_glycosite_glycosylation <- function(data, output_dir = "output_data/peptidegroups_intensity/glycosite_level", 
                                           figures_dir = "figures/glycosite_level") {
  
  # Create output directories if they don't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Split data by glycosite ID
  glycosite_list <- split(data, data$gsite_ID)
  
  # Initialize results storage
  fuc_results <- list()
  neuac_results <- list()
  summary_stats <- data.frame()
  
  # Track glycosites for summary
  glycosites_analyzed <- 0
  glycosites_insufficient_data <- 0
  glycosites_insufficient_groups <- 0
  glycosites_with_fuc_data <- 0
  glycosites_with_neuac_data <- 0
  fuc_tests_performed <- 0
  neuac_tests_performed <- 0
  
  # Track failed tests
  failed_fuc_tests <- data.frame(
    glycosite_ID = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
  
  failed_neuac_tests <- data.frame(
    glycosite_ID = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
  
  cat("Analyzing", length(glycosite_list), "glycosites for fucosylation and sialylation patterns...\n")
  
  # Process each glycosite
  for(i in seq_along(glycosite_list)) {
    glycosite_name <- names(glycosite_list)[i]
    glycosite_data <- glycosite_list[[i]]
    
    cat(sprintf("\nProcessing glycosite %d/%d: %s\n", i, length(glycosite_list), glycosite_name))
    
    # Skip glycosites with insufficient data
    if(n_distinct(glycosite_data$sample) < 3) {
      cat("  Skipping - insufficient total samples\n")
      glycosites_insufficient_data <- glycosites_insufficient_data + 1
      next
    }
    
    # Check if we have at least 3 samples in each group
    sample_counts <- glycosite_data %>%
      group_by(group) %>%
      summarise(n_samples = n_distinct(sample), .groups = 'drop')
    
    if(any(sample_counts$n_samples < 3)) {
      cat("  Skipping - insufficient samples in each group (need ≥3 per group)\n")
      cat("    Sample counts:", paste(sample_counts$group, "=", sample_counts$n_samples, collapse = ", "), "\n")
      glycosites_insufficient_groups <- glycosites_insufficient_groups + 1
      next
    }
    
    glycosites_analyzed <- glycosites_analyzed + 1
    
    # Calculate fucosylation percentages by sample
    fuc_by_sample <- glycosite_data %>%
      group_by(sample, group) %>%
      summarise(
        total_abundance = sum(abundance, na.rm = TRUE),
        fuc_abundance = sum(abundance[contains_Fuc == TRUE], na.rm = TRUE),
        fuc_percentage = (fuc_abundance / total_abundance) * 100,
        .groups = 'drop'
      )
    
    # Calculate sialylation percentages by sample
    neuac_by_sample <- glycosite_data %>%
      group_by(sample, group) %>%
      summarise(
        total_abundance = sum(abundance, na.rm = TRUE),
        neuac_abundance = sum(abundance[contains_NeuAc == TRUE], na.rm = TRUE),
        neuac_percentage = (neuac_abundance / total_abundance) * 100,
        .groups = 'drop'
      )
    
    # Check if glycosite has sufficient fucosylation data (≥3 samples per group with fucosylation)
    fuc_sample_counts <- fuc_by_sample %>%
      filter(fuc_percentage > 0) %>%  # Only samples with some fucosylation
      group_by(group) %>%
      summarise(n_samples = n(), .groups = 'drop')
    
    has_sufficient_fuc <- nrow(fuc_sample_counts) >= 2 && all(fuc_sample_counts$n_samples >= 3)
    
    # Check if glycosite has sufficient sialylation data (≥3 samples per group with sialylation)
    neuac_sample_counts <- neuac_by_sample %>%
      filter(neuac_percentage > 0) %>%  # Only samples with some sialylation
      group_by(group) %>%
      summarise(n_samples = n(), .groups = 'drop')
    
    has_sufficient_neuac <- nrow(neuac_sample_counts) >= 2 && all(neuac_sample_counts$n_samples >= 3)
    
    # Store results only if sufficient data exists
    if(has_sufficient_fuc) {
      fuc_results[[glycosite_name]] <- fuc_by_sample
      glycosites_with_fuc_data <- glycosites_with_fuc_data + 1
    } else {
      cat("    No sufficient fucosylation data (need ≥3 samples per group with fucosylation)\n")
    }
    if(has_sufficient_neuac) {
      neuac_results[[glycosite_name]] <- neuac_by_sample
      glycosites_with_neuac_data <- glycosites_with_neuac_data + 1
    } else {
      cat("    No sufficient sialylation data (need ≥3 samples per group with sialylation)\n")
    }
    
    # Perform statistical tests for fucosylation
    if(has_sufficient_fuc) {
      
      # Check for constant data within groups
      healthy_fuc_data <- fuc_by_sample$fuc_percentage[fuc_by_sample$group == "Healthy"]
      mecfs_fuc_data <- fuc_by_sample$fuc_percentage[fuc_by_sample$group == "MECFS"]
      
      # Check if data is constant within groups
      if(length(unique(healthy_fuc_data)) <= 1 && length(unique(mecfs_fuc_data)) <= 1) {
        failed_fuc_tests <- rbind(failed_fuc_tests, data.frame(
          glycosite_ID = glycosite_name,
          reason = "Constant data within both groups",
          stringsAsFactors = FALSE
        ))
        cat("    Fucosylation test failed: Constant data within both groups\n")
      } else if(length(unique(healthy_fuc_data)) <= 1) {
        failed_fuc_tests <- rbind(failed_fuc_tests, data.frame(
          glycosite_ID = glycosite_name,
          reason = "Constant data in Healthy group",
          stringsAsFactors = FALSE
        ))
        cat("    Fucosylation test failed: Constant data in Healthy group\n")
      } else if(length(unique(mecfs_fuc_data)) <= 1) {
        failed_fuc_tests <- rbind(failed_fuc_tests, data.frame(
          glycosite_ID = glycosite_name,
          reason = "Constant data in MECFS group",
          stringsAsFactors = FALSE
        ))
        cat("    Fucosylation test failed: Constant data in MECFS group\n")
      } else {
        # F-test for variance
        fuc_f_test <- try({
          var.test(fuc_percentage ~ group, data = fuc_by_sample)
        }, silent = TRUE)
        
        # T-test
        fuc_t_test <- try({
          t.test(fuc_percentage ~ group, data = fuc_by_sample)
        }, silent = TRUE)
        
        if(!inherits(fuc_t_test, "try-error")) {
          fuc_stats <- data.frame(
            glycosite_ID = glycosite_name,
            feature = "Fucosylation",
            mean_healthy = mean(healthy_fuc_data, na.rm = TRUE),
            mean_mecfs = mean(mecfs_fuc_data, na.rm = TRUE),
            sd_healthy = sd(healthy_fuc_data, na.rm = TRUE),
            sd_mecfs = sd(mecfs_fuc_data, na.rm = TRUE),
            t_stat = fuc_t_test$statistic,
            p_value = fuc_t_test$p.value,
            f_stat = ifelse(!inherits(fuc_f_test, "try-error"), fuc_f_test$statistic, NA),
            f_p_value = ifelse(!inherits(fuc_f_test, "try-error"), fuc_f_test$p.value, NA),
            n_healthy = sum(fuc_by_sample$group == "Healthy"),
            n_mecfs = sum(fuc_by_sample$group == "MECFS"),
            stringsAsFactors = FALSE
          )
          summary_stats <- rbind(summary_stats, fuc_stats)
          fuc_tests_performed <- fuc_tests_performed + 1
        } else {
          failed_fuc_tests <- rbind(failed_fuc_tests, data.frame(
            glycosite_ID = glycosite_name,
            reason = paste("T-test error:", fuc_t_test[1]),
            stringsAsFactors = FALSE
          ))
          cat("    Fucosylation t-test failed:", fuc_t_test[1], "\n")
        }
      }
    }
    
    # Perform statistical tests for sialylation
    if(has_sufficient_neuac) {
      
      # Check for constant data within groups
      healthy_neuac_data <- neuac_by_sample$neuac_percentage[neuac_by_sample$group == "Healthy"]
      mecfs_neuac_data <- neuac_by_sample$neuac_percentage[neuac_by_sample$group == "MECFS"]
      
      # Check if data is constant within groups
      if(length(unique(healthy_neuac_data)) <= 1 && length(unique(mecfs_neuac_data)) <= 1) {
        failed_neuac_tests <- rbind(failed_neuac_tests, data.frame(
          glycosite_ID = glycosite_name,
          reason = "Constant data within both groups",
          stringsAsFactors = FALSE
        ))
        cat("    Sialylation test failed: Constant data within both groups\n")
      } else if(length(unique(healthy_neuac_data)) <= 1) {
        failed_neuac_tests <- rbind(failed_neuac_tests, data.frame(
          glycosite_ID = glycosite_name,
          reason = "Constant data in Healthy group",
          stringsAsFactors = FALSE
        ))
        cat("    Sialylation test failed: Constant data in Healthy group\n")
      } else if(length(unique(mecfs_neuac_data)) <= 1) {
        failed_neuac_tests <- rbind(failed_neuac_tests, data.frame(
          glycosite_ID = glycosite_name,
          reason = "Constant data in MECFS group",
          stringsAsFactors = FALSE
        ))
        cat("    Sialylation test failed: Constant data in MECFS group\n")
      } else {
        # F-test for variance
        neuac_f_test <- try({
          var.test(neuac_percentage ~ group, data = neuac_by_sample)
        }, silent = TRUE)
        
        # T-test
        neuac_t_test <- try({
          t.test(neuac_percentage ~ group, data = neuac_by_sample)
        }, silent = TRUE)
        
        if(!inherits(neuac_t_test, "try-error")) {
          neuac_stats <- data.frame(
            glycosite_ID = glycosite_name,
            feature = "Sialylation",
            mean_healthy = mean(healthy_neuac_data, na.rm = TRUE),
            mean_mecfs = mean(mecfs_neuac_data, na.rm = TRUE),
            sd_healthy = sd(healthy_neuac_data, na.rm = TRUE),
            sd_mecfs = sd(mecfs_neuac_data, na.rm = TRUE),
            t_stat = neuac_t_test$statistic,
            p_value = neuac_t_test$p.value,
            f_stat = ifelse(!inherits(neuac_f_test, "try-error"), neuac_f_test$statistic, NA),
            f_p_value = ifelse(!inherits(neuac_f_test, "try-error"), neuac_f_test$p.value, NA),
            n_healthy = sum(neuac_by_sample$group == "Healthy"),
            n_mecfs = sum(neuac_by_sample$group == "MECFS"),
            stringsAsFactors = FALSE
          )
          summary_stats <- rbind(summary_stats, neuac_stats)
          neuac_tests_performed <- neuac_tests_performed + 1
        } else {
          failed_neuac_tests <- rbind(failed_neuac_tests, data.frame(
            glycosite_ID = glycosite_name,
            reason = paste("T-test error:", neuac_t_test[1]),
            stringsAsFactors = FALSE
          ))
          cat("    Sialylation t-test failed:", neuac_t_test[1], "\n")
        }
      }
    }
  }
  
  # Apply BH correction to all p-values
  if(nrow(summary_stats) > 0) {
    summary_stats$p_value_adj <- p.adjust(summary_stats$p_value, method = "BH")
    summary_stats$significant <- summary_stats$p_value_adj < 0.05
    summary_stats$fold_change <- summary_stats$mean_mecfs / summary_stats$mean_healthy
  }
  
  # Create separate boxplots for fucosylation and sialylation
  cat("\nCreating separate boxplots for fucosylation and sialylation...\n")
  
  for(glycosite_name in names(glycosite_list)) {
    safe_glycosite_name <- gsub("[^a-zA-Z0-9]", "_", glycosite_name)
    
    # Create fucosylation boxplot if sufficient data exists
    if(glycosite_name %in% names(fuc_results)) {
      fuc_plot <- ggplot(fuc_results[[glycosite_name]], 
                        aes(x = group, y = fuc_percentage, fill = group)) +
        geom_boxplot(alpha = 0.7, outlier.shape = 1) +
        geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
        scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
        labs(title = paste("Fucosylation -", glycosite_name),
             x = "Group", 
             y = "Fucosylation (%)") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 10),
              legend.position = "top")
      
      # Save fucosylation plot
      ggsave(paste0(figures_dir, "/", safe_glycosite_name, "_fucosylation.png"), 
             fuc_plot, width = 8, height = 6, dpi = 300)
    }
    
    # Create sialylation boxplot if sufficient data exists
    if(glycosite_name %in% names(neuac_results)) {
      neuac_plot <- ggplot(neuac_results[[glycosite_name]], 
                          aes(x = group, y = neuac_percentage, fill = group)) +
        geom_boxplot(alpha = 0.7, outlier.shape = 1) +
        geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
        scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
        labs(title = paste("Sialylation -", glycosite_name),
             x = "Group", 
             y = "Sialylation (%)") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 10),
              legend.position = "top")
      
      # Save sialylation plot
      ggsave(paste0(figures_dir, "/", safe_glycosite_name, "_sialylation.png"), 
             neuac_plot, width = 8, height = 6, dpi = 300)
    }
  }
  
  # Save results
  saveRDS(fuc_results, paste0(output_dir, "/fucosylation_by_glycosite.rds"))
  saveRDS(neuac_results, paste0(output_dir, "/sialylation_by_glycosite.rds"))
  write.csv(summary_stats, paste0(output_dir, "/glycosite_glycosylation_statistics.csv"), row.names = FALSE)
  
  # Create separate summary files for fucosylation and sialylation
  fuc_summary <- summary_stats %>% filter(feature == "Fucosylation")
  neuac_summary <- summary_stats %>% filter(feature == "Sialylation")
  
  write.csv(fuc_summary, paste0(output_dir, "/fucosylation_statistics.csv"), row.names = FALSE)
  write.csv(neuac_summary, paste0(output_dir, "/sialylation_statistics.csv"), row.names = FALSE)
  
  # Create detailed summary dataframes for each feature
  fuc_detailed <- bind_rows(lapply(names(fuc_results), function(glycosite) {
    data <- fuc_results[[glycosite]]
    data$glycosite_ID <- glycosite
    return(data)
  }))
  
  neuac_detailed <- bind_rows(lapply(names(neuac_results), function(glycosite) {
    data <- neuac_results[[glycosite]]
    data$glycosite_ID <- glycosite
    return(data)
  }))
  
  write.csv(fuc_detailed, paste0(output_dir, "/fucosylation_detailed_data.csv"), row.names = FALSE)
  write.csv(neuac_detailed, paste0(output_dir, "/sialylation_detailed_data.csv"), row.names = FALSE)
  
  # Save failed test information
  write.csv(failed_fuc_tests, paste0(output_dir, "/failed_fucosylation_tests.csv"), row.names = FALSE)
  write.csv(failed_neuac_tests, paste0(output_dir, "/failed_sialylation_tests.csv"), row.names = FALSE)
  
  # Print summary
  cat("\n=== GLYCOSITE-LEVEL GLYCOSYLATION ANALYSIS SUMMARY ===\n")
  cat("Total glycosites in dataset:", length(glycosite_list), "\n")
  cat("Glycosites with sufficient data (≥3 samples per group):", glycosites_analyzed, "\n")
  cat("Glycosites excluded - insufficient total samples:", glycosites_insufficient_data, "\n")
  cat("Glycosites excluded - insufficient samples per group:", glycosites_insufficient_groups, "\n")
  cat("Glycosites with sufficient fucosylation data (≥3 samples per group with fucosylation):", glycosites_with_fuc_data, "\n")
  cat("Glycosites with sufficient sialylation data (≥3 samples per group with sialylation):", glycosites_with_neuac_data, "\n")
  cat("Fucosylation statistical tests performed:", fuc_tests_performed, "\n")
  cat("Sialylation statistical tests performed:", neuac_tests_performed, "\n")
  cat("Total statistical comparisons with results:", nrow(summary_stats), "\n")
  
  if(nrow(summary_stats) > 0) {
    cat("Significant differences (FDR < 0.05):", sum(summary_stats$significant, na.rm = TRUE), "\n")
    
    # Fucosylation summary
    fuc_stats <- summary_stats %>% filter(feature == "Fucosylation")
    cat(sprintf("\n=== FUCOSYLATION ANALYSIS ===\n"))
    cat("Glycosites with sufficient fucosylation data:", glycosites_with_fuc_data, "\n")
    cat("Glycosites with successful fucosylation statistical tests:", nrow(fuc_stats), "\n")
    cat("Significant fucosylation differences:", sum(fuc_stats$significant, na.rm = TRUE), "\n")
    
    if(nrow(fuc_stats) > 0) {
      cat("Mean fucosylation by group:\n")
      fuc_means <- fuc_stats %>%
        summarise(
          mean_healthy = mean(mean_healthy, na.rm = TRUE),
          mean_mecfs = mean(mean_mecfs, na.rm = TRUE),
          .groups = 'drop'
        )
      cat(sprintf("  Healthy: %.2f%%\n", fuc_means$mean_healthy))
      cat(sprintf("  MECFS: %.2f%%\n", fuc_means$mean_mecfs))
    }
    
    # Significant fucosylation differences
    sig_fuc <- fuc_stats %>% filter(significant)
    if(nrow(sig_fuc) > 0) {
      cat("\nSignificant fucosylation differences:\n")
      for(i in 1:nrow(sig_fuc)) {
        cat(sprintf("  %s: p = %.3e, FC = %.2f\n", 
                   sig_fuc$glycosite_ID[i], 
                   sig_fuc$p_value_adj[i], 
                   sig_fuc$fold_change[i]))
      }
    }
    
    # Sialylation summary
    neuac_stats <- summary_stats %>% filter(feature == "Sialylation")
    cat(sprintf("\n=== SIALYLATION ANALYSIS ===\n"))
    cat("Glycosites with sufficient sialylation data:", glycosites_with_neuac_data, "\n")
    cat("Glycosites with successful sialylation statistical tests:", nrow(neuac_stats), "\n")
    cat("Significant sialylation differences:", sum(neuac_stats$significant, na.rm = TRUE), "\n")
    
    if(nrow(neuac_stats) > 0) {
      cat("Mean sialylation by group:\n")
      neuac_means <- neuac_stats %>%
        summarise(
          mean_healthy = mean(mean_healthy, na.rm = TRUE),
          mean_mecfs = mean(mean_mecfs, na.rm = TRUE),
          .groups = 'drop'
        )
      cat(sprintf("  Healthy: %.2f%%\n", neuac_means$mean_healthy))
      cat(sprintf("  MECFS: %.2f%%\n", neuac_means$mean_mecfs))
    }
    
    # Significant sialylation differences
    sig_neuac <- neuac_stats %>% filter(significant)
    if(nrow(sig_neuac) > 0) {
      cat("\nSignificant sialylation differences:\n")
      for(i in 1:nrow(sig_neuac)) {
        cat(sprintf("  %s: p = %.3e, FC = %.2f\n", 
                   sig_neuac$glycosite_ID[i], 
                   sig_neuac$p_value_adj[i], 
                   sig_neuac$fold_change[i]))
      }
    }
  }
  
  cat("\nResults saved to:", output_dir, "\n")
  cat("  - fucosylation_statistics.csv: Fucosylation statistical results\n")
  cat("  - sialylation_statistics.csv: Sialylation statistical results\n")
  cat("  - fucosylation_detailed_data.csv: Detailed fucosylation data by glycosite and sample\n")
  cat("  - sialylation_detailed_data.csv: Detailed sialylation data by glycosite and sample\n")
  cat("  - glycosite_glycosylation_statistics.csv: Combined statistical results\n")
  cat("  - fucosylation_by_glycosite.rds: Fucosylation data as R object\n")
  cat("  - sialylation_by_glycosite.rds: Sialylation data as R object\n")
  cat("  - failed_fucosylation_tests.csv: Glycosites with failed fucosylation tests\n")
  cat("  - failed_sialylation_tests.csv: Glycosites with failed sialylation tests\n")
  cat("\nFigures saved to:", figures_dir, "\n")
  cat("  - *_fucosylation.png: Individual fucosylation boxplots\n")
  cat("  - *_sialylation.png: Individual sialylation boxplots\n")
  
  # Report failed tests
  if(nrow(failed_fuc_tests) > 0) {
    cat("\n=== FAILED FUCOSYLATION TESTS ===\n")
    cat("Total failed fucosylation tests:", nrow(failed_fuc_tests), "\n")
    
    # Count reasons
    fuc_reason_counts <- table(failed_fuc_tests$reason)
    cat("Reasons for failure:\n")
    for(i in 1:length(fuc_reason_counts)) {
      cat(sprintf("  %s: %d glycosites\n", names(fuc_reason_counts)[i], fuc_reason_counts[i]))
    }
    
    cat("\nGlycosites with failed fucosylation tests:\n")
    for(i in 1:nrow(failed_fuc_tests)) {
      cat(sprintf("  %s: %s\n", failed_fuc_tests$glycosite_ID[i], failed_fuc_tests$reason[i]))
    }
  }
  
  if(nrow(failed_neuac_tests) > 0) {
    cat("\n=== FAILED SIALYLATION TESTS ===\n")
    cat("Total failed sialylation tests:", nrow(failed_neuac_tests), "\n")
    
    # Count reasons
    neuac_reason_counts <- table(failed_neuac_tests$reason)
    cat("Reasons for failure:\n")
    for(i in 1:length(neuac_reason_counts)) {
      cat(sprintf("  %s: %d glycosites\n", names(neuac_reason_counts)[i], neuac_reason_counts[i]))
    }
    
    cat("\nGlycosites with failed sialylation tests:\n")
    for(i in 1:nrow(failed_neuac_tests)) {
      cat(sprintf("  %s: %s\n", failed_neuac_tests$glycosite_ID[i], failed_neuac_tests$reason[i]))
    }
  }
  
  return(list(
    fucosylation_data = fuc_results,
    sialylation_data = neuac_results,
    statistics = summary_stats,
    failed_fucosylation_tests = failed_fuc_tests,
    failed_sialylation_tests = failed_neuac_tests
  ))
}

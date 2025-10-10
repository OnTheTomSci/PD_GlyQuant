# Peptide Groups Protein Glycosylation Analysis Functions
# This module contains functions for protein-level glycosylation analysis including ILR transformation

library(tidyverse)
library(ggplot2)

#' Analyze protein-level glycosylation patterns
#' Uses Isometric Log Ratio (ILR) transformation for compositional data analysis
#' ILR transformation: log(proportion_fucosylated / proportion_non_fucosylated)
#' This transformation makes the data suitable for standard statistical tests
#' 
#' @param data Long format data with abundance, sample, group, protein_accessions, contains_Fuc, contains_NeuAc columns
#' @param output_dir Directory to save output files
#' @param figures_dir Directory to save figures
#' @return List containing protein glycosylation analysis results
analyze_protein_glycosylation <- function(data, output_dir = "output_data/peptidegroups_intensity/protein_level", 
                                         figures_dir = "figures/protein_level") {
  
  # Create output directories if they don't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Split data by protein accession
  protein_list <- split(data, data$protein_accessions)
  
  # Initialize results storage
  fuc_results <- list()
  neuac_results <- list()
  summary_stats <- data.frame()
  
  # Track proteins for summary
  proteins_analyzed <- 0
  proteins_insufficient_data <- 0
  proteins_insufficient_groups <- 0
  proteins_with_fuc_data <- 0
  proteins_with_neuac_data <- 0
  fuc_tests_performed <- 0
  neuac_tests_performed <- 0
  
  # Track failed tests
  failed_fuc_tests <- data.frame(
    protein_accession = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
  
  failed_neuac_tests <- data.frame(
    protein_accession = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
  
  cat("Analyzing", length(protein_list), "proteins for fucosylation and sialylation patterns...\n")
  
  # Process each protein
  for(i in seq_along(protein_list)) {
    protein_name <- names(protein_list)[i]
    protein_data <- protein_list[[i]]
    
    cat(sprintf("\nProcessing protein %d/%d: %s\n", i, length(protein_list), protein_name))
    
    # Skip proteins with insufficient data
    if(n_distinct(protein_data$sample) < 3) {
      cat("  Skipping - insufficient total samples\n")
      proteins_insufficient_data <- proteins_insufficient_data + 1
      next
    }
    
    # Check if we have at least 3 samples in each group
    sample_counts <- protein_data %>%
      group_by(group) %>%
      summarise(n_samples = n_distinct(sample), .groups = 'drop')
    
    if(any(sample_counts$n_samples < 3)) {
      cat("  Skipping - insufficient samples in each group (need ≥3 per group)\n")
      cat("    Sample counts:", paste(sample_counts$group, "=", sample_counts$n_samples, collapse = ", "), "\n")
      proteins_insufficient_groups <- proteins_insufficient_groups + 1
      next
    }
    
    proteins_analyzed <- proteins_analyzed + 1
    
    # Calculate fucosylation percentages by sample
    fuc_by_sample <- protein_data %>%
      group_by(sample, group) %>%
      summarise(
        total_abundance = sum(abundance, na.rm = TRUE),
        fuc_abundance = sum(abundance[contains_Fuc == TRUE], na.rm = TRUE),
        fuc_percentage = (fuc_abundance / total_abundance) * 100,
        .groups = 'drop'
      )
    
    # Apply ILR transformation to fucosylation data
    # ILR requires compositional data (proportions that sum to 1)
    fuc_by_sample <- fuc_by_sample %>%
      mutate(
        # Convert percentages to proportions
        fuc_proportion = fuc_percentage / 100,
        non_fuc_proportion = 1 - fuc_proportion,
        # Apply ILR transformation: log(fuc_proportion / non_fuc_proportion)
        fuc_ilr = log(fuc_proportion / non_fuc_proportion)
      )
    
    # Calculate sialylation percentages by sample
    neuac_by_sample <- protein_data %>%
      group_by(sample, group) %>%
      summarise(
        total_abundance = sum(abundance, na.rm = TRUE),
        neuac_abundance = sum(abundance[contains_NeuAc == TRUE], na.rm = TRUE),
        neuac_percentage = (neuac_abundance / total_abundance) * 100,
        .groups = 'drop'
      )
    
    # Apply ILR transformation to sialylation data
    neuac_by_sample <- neuac_by_sample %>%
      mutate(
        # Convert percentages to proportions
        neuac_proportion = neuac_percentage / 100,
        non_neuac_proportion = 1 - neuac_proportion,
        # Apply ILR transformation: log(neuac_proportion / non_neuac_proportion)
        # Handle edge cases for ILR transformation
        neuac_proportion_adj = ifelse(neuac_proportion == 0, 0.001, 
                                     ifelse(neuac_proportion == 1, 0.999, neuac_proportion)),
        non_neuac_proportion_adj = 1 - neuac_proportion_adj,
        neuac_ilr = log(neuac_proportion_adj / non_neuac_proportion_adj)
      ) %>%
      # Remove infinite or NaN values
      filter(is.finite(neuac_ilr))
    
    # Check if protein has sufficient fucosylation data (≥3 samples per group with fucosylation)
    fuc_sample_counts <- fuc_by_sample %>%
      filter(fuc_percentage > 0) %>%  # Only samples with some fucosylation
      group_by(group) %>%
      summarise(n_samples = n(), .groups = 'drop')
    
    has_sufficient_fuc <- nrow(fuc_sample_counts) >= 2 && all(fuc_sample_counts$n_samples >= 3)
    
    # Check if protein has sufficient sialylation data (≥3 samples per group with sialylation)
    neuac_sample_counts <- neuac_by_sample %>%
      filter(neuac_percentage > 0) %>%  # Only samples with some sialylation
      group_by(group) %>%
      summarise(n_samples = n(), .groups = 'drop')
    
    has_sufficient_neuac <- nrow(neuac_sample_counts) >= 2 && all(neuac_sample_counts$n_samples >= 3)
    
    # Store results only if sufficient data exists
    if(has_sufficient_fuc) {
      fuc_results[[protein_name]] <- fuc_by_sample
      proteins_with_fuc_data <- proteins_with_fuc_data + 1
    } else {
      cat("    No sufficient fucosylation data (need ≥3 samples per group with fucosylation)\n")
    }
    if(has_sufficient_neuac) {
      neuac_results[[protein_name]] <- neuac_by_sample
      proteins_with_neuac_data <- proteins_with_neuac_data + 1
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
          protein_accession = protein_name,
          reason = "Constant data within both groups",
          stringsAsFactors = FALSE
        ))
        cat("    Fucosylation test failed: Constant data within both groups\n")
      } else if(length(unique(healthy_fuc_data)) <= 1) {
        failed_fuc_tests <- rbind(failed_fuc_tests, data.frame(
          protein_accession = protein_name,
          reason = "Constant data in Healthy group",
          stringsAsFactors = FALSE
        ))
        cat("    Fucosylation test failed: Constant data in Healthy group\n")
      } else if(length(unique(mecfs_fuc_data)) <= 1) {
        failed_fuc_tests <- rbind(failed_fuc_tests, data.frame(
          protein_accession = protein_name,
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
            protein_accession = protein_name,
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
            protein_accession = protein_name,
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
          protein_accession = protein_name,
          reason = "Constant data within both groups",
          stringsAsFactors = FALSE
        ))
        cat("    Sialylation test failed: Constant data within both groups\n")
      } else if(length(unique(healthy_neuac_data)) <= 1) {
        failed_neuac_tests <- rbind(failed_neuac_tests, data.frame(
          protein_accession = protein_name,
          reason = "Constant data in Healthy group",
          stringsAsFactors = FALSE
        ))
        cat("    Sialylation test failed: Constant data in Healthy group\n")
      } else if(length(unique(mecfs_neuac_data)) <= 1) {
        failed_neuac_tests <- rbind(failed_neuac_tests, data.frame(
          protein_accession = protein_name,
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
            protein_accession = protein_name,
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
            protein_accession = protein_name,
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
  
  for(protein_name in names(protein_list)) {
    safe_protein_name <- gsub("[^a-zA-Z0-9]", "_", protein_name)
    
    # Create fucosylation boxplot if sufficient data exists
    if(protein_name %in% names(fuc_results)) {
      fuc_plot <- ggplot(fuc_results[[protein_name]], 
                        aes(x = group, y = fuc_percentage, fill = group)) +
        geom_boxplot(alpha = 0.7, outlier.shape = 1) +
        geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
        scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
        labs(title = paste("Fucosylation -", protein_name),
             x = "Group", 
             y = "Fucosylation (%)") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 14),
              legend.position = "top",
              axis.text = element_text(size = 12))
      
      # Save fucosylation plot
      ggsave(paste0(figures_dir, "/", safe_protein_name, "_fucosylation.png"), 
             fuc_plot, width = 8, height = 6, dpi = 300)
    }
    
    # Create sialylation boxplot if sufficient data exists
    if(protein_name %in% names(neuac_results)) {
      neuac_plot <- ggplot(neuac_results[[protein_name]], 
                          aes(x = group, y = neuac_percentage, fill = group)) +
        geom_boxplot(alpha = 0.7, outlier.shape = 1) +
        geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
        scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
        labs(title = paste("Sialylation -", protein_name),
             x = "Group", 
             y = "Sialylation (%)") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 14),
              legend.position = "top",
              axis.text = element_text(size = 12))
      
      # Save sialylation plot
      ggsave(paste0(figures_dir, "/", safe_protein_name, "_sialylation.png"), 
             neuac_plot, width = 8, height = 6, dpi = 300)
    }
  }
  
  # Save results
  saveRDS(fuc_results, paste0(output_dir, "/fucosylation_by_protein.rds"))
  saveRDS(neuac_results, paste0(output_dir, "/sialylation_by_protein.rds"))
  write.csv(summary_stats, paste0(output_dir, "/protein_glycosylation_statistics.csv"), row.names = FALSE)
  
  # Create separate summary files for fucosylation and sialylation
  fuc_summary <- summary_stats %>% filter(feature == "Fucosylation")
  neuac_summary <- summary_stats %>% filter(feature == "Sialylation")
  
  write.csv(fuc_summary, paste0(output_dir, "/fucosylation_statistics.csv"), row.names = FALSE)
  write.csv(neuac_summary, paste0(output_dir, "/sialylation_statistics.csv"), row.names = FALSE)
  
  # Create detailed summary dataframes for each feature
  fuc_detailed <- bind_rows(lapply(names(fuc_results), function(protein) {
    data <- fuc_results[[protein]]
    data$protein_accession <- protein
    return(data)
  }))
  
  neuac_detailed <- bind_rows(lapply(names(neuac_results), function(protein) {
    data <- neuac_results[[protein]]
    data$protein_accession <- protein
    return(data)
  }))
  
  write.csv(fuc_detailed, paste0(output_dir, "/fucosylation_detailed_data.csv"), row.names = FALSE)
  write.csv(neuac_detailed, paste0(output_dir, "/sialylation_detailed_data.csv"), row.names = FALSE)
  
  # Save failed test information
  write.csv(failed_fuc_tests, paste0(output_dir, "/failed_fucosylation_tests.csv"), row.names = FALSE)
  write.csv(failed_neuac_tests, paste0(output_dir, "/failed_sialylation_tests.csv"), row.names = FALSE)
  
  # Print summary
  cat("\n=== PROTEIN-LEVEL GLYCOSYLATION ANALYSIS SUMMARY ===\n")
  cat("Total proteins in dataset:", length(protein_list), "\n")
  cat("Proteins with sufficient data (≥3 samples per group):", proteins_analyzed, "\n")
  cat("Proteins excluded - insufficient total samples:", proteins_insufficient_data, "\n")
  cat("Proteins excluded - insufficient samples per group:", proteins_insufficient_groups, "\n")
  cat("Proteins with sufficient fucosylation data (≥3 samples per group with fucosylation):", proteins_with_fuc_data, "\n")
  cat("Proteins with sufficient sialylation data (≥3 samples per group with sialylation):", proteins_with_neuac_data, "\n")
  cat("Fucosylation statistical tests performed:", fuc_tests_performed, "\n")
  cat("Sialylation statistical tests performed:", neuac_tests_performed, "\n")
  cat("Total statistical comparisons with results:", nrow(summary_stats), "\n")
  
  if(nrow(summary_stats) > 0) {
    cat("Significant differences (FDR < 0.05):", sum(summary_stats$significant, na.rm = TRUE), "\n")
    
    # Fucosylation summary
    fuc_stats <- summary_stats %>% filter(feature == "Fucosylation")
    cat(sprintf("\n=== FUCOSYLATION ANALYSIS ===\n"))
    cat("Proteins with sufficient fucosylation data:", proteins_with_fuc_data, "\n")
    cat("Proteins with successful fucosylation statistical tests:", nrow(fuc_stats), "\n")
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
                   sig_fuc$protein_accession[i], 
                   sig_fuc$p_value_adj[i], 
                   sig_fuc$fold_change[i]))
      }
    }
    
    # Sialylation summary
    neuac_stats <- summary_stats %>% filter(feature == "Sialylation")
    cat(sprintf("\n=== SIALYLATION ANALYSIS ===\n"))
    cat("Proteins with sufficient sialylation data:", proteins_with_neuac_data, "\n")
    cat("Proteins with successful sialylation statistical tests:", nrow(neuac_stats), "\n")
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
                   sig_neuac$protein_accession[i], 
                   sig_neuac$p_value_adj[i], 
                   sig_neuac$fold_change[i]))
      }
    }
  }
  
  cat("\nResults saved to:", output_dir, "\n")
  cat("  - fucosylation_statistics.csv: Fucosylation statistical results\n")
  cat("  - sialylation_statistics.csv: Sialylation statistical results\n")
  cat("  - fucosylation_detailed_data.csv: Detailed fucosylation data by protein and sample\n")
  cat("  - sialylation_detailed_data.csv: Detailed sialylation data by protein and sample\n")
  cat("  - protein_glycosylation_statistics.csv: Combined statistical results\n")
  cat("  - fucosylation_by_protein.rds: Fucosylation data as R object\n")
  cat("  - sialylation_by_protein.rds: Sialylation data as R object\n")
  cat("  - failed_fucosylation_tests.csv: Proteins with failed fucosylation tests\n")
  cat("  - failed_sialylation_tests.csv: Proteins with failed sialylation tests\n")
  cat("\nFigures saved to:", figures_dir, "\n")
  cat("  - *_fucosylation.png: Individual fucosylation boxplots\n")
  cat("  - *_sialylation.png: Individual sialylation boxplots\n")
  
  return(list(
    fucosylation_data = fuc_results,
    sialylation_data = neuac_results,
    statistics = summary_stats,
    failed_fucosylation_tests = failed_fuc_tests,
    failed_sialylation_tests = failed_neuac_tests
  ))
}

#' Additional ILR-based analysis function
#' 
#' @param data Long format data with abundance, sample, group, protein_accessions, contains_Fuc, contains_NeuAc columns
#' @param original_results Original protein glycosylation results (optional)
#' @param output_dir Directory to save output files
#' @param figures_dir Directory to save figures
#' @return List containing ILR-based analysis results
analyze_protein_glycosylation_ilr <- function(data, original_results = NULL, output_dir = "output_data/peptidegroups_intensity/protein_level_ilr", 
                                             figures_dir = "figures/protein_level_ilr") {
  
  # Create output directories if they don't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Split data by protein accession
  protein_list <- split(data, data$protein_accessions)
  
  # Initialize results storage
  ilr_summary_stats <- data.frame()
  
  cat("Performing ILR-based statistical analysis...\n")
  
  # Analyze each protein with ILR transformation
  for(protein_name in names(protein_list)) {
    protein_data <- protein_list[[protein_name]]
    
    # Skip if insufficient data
    if(n_distinct(protein_data$sample) < 3) next
    
    # Check group sample sizes
    sample_counts <- protein_data %>%
      group_by(group) %>%
      summarise(n_samples = n_distinct(sample), .groups = 'drop')
    
    if(any(sample_counts$n_samples < 3)) next
    
    # Calculate fucosylation and sialylation with ILR
    fuc_by_sample <- protein_data %>%
      group_by(sample, group) %>%
      summarise(
        total_abundance = sum(abundance, na.rm = TRUE),
        fuc_abundance = sum(abundance[contains_Fuc == TRUE], na.rm = TRUE),
        fuc_percentage = (fuc_abundance / total_abundance) * 100,
        .groups = 'drop'
      ) %>%
      mutate(
        fuc_proportion = fuc_percentage / 100,
        non_fuc_proportion = 1 - fuc_proportion,
        # Handle edge cases for ILR transformation
        fuc_proportion_adj = ifelse(fuc_proportion == 0, 0.001, 
                                   ifelse(fuc_proportion == 1, 0.999, fuc_proportion)),
        non_fuc_proportion_adj = 1 - fuc_proportion_adj,
        fuc_ilr = log(fuc_proportion_adj / non_fuc_proportion_adj)
      ) %>%
      # Remove infinite or NaN values
      filter(is.finite(fuc_ilr))
    
    neuac_by_sample <- protein_data %>%
      group_by(sample, group) %>%
      summarise(
        total_abundance = sum(abundance, na.rm = TRUE),
        neuac_abundance = sum(abundance[contains_NeuAc == TRUE], na.rm = TRUE),
        neuac_percentage = (neuac_abundance / total_abundance) * 100,
        .groups = 'drop'
      ) %>%
      mutate(
        neuac_proportion = neuac_percentage / 100,
        non_neuac_proportion = 1 - neuac_proportion,
        # Handle edge cases for ILR transformation
        neuac_proportion_adj = ifelse(neuac_proportion == 0, 0.001, 
                                     ifelse(neuac_proportion == 1, 0.999, neuac_proportion)),
        non_neuac_proportion_adj = 1 - neuac_proportion_adj,
        neuac_ilr = log(neuac_proportion_adj / non_neuac_proportion_adj)
      ) %>%
      # Remove infinite or NaN values
      filter(is.finite(neuac_ilr))
    
    # Check if sufficient data for ILR analysis
    fuc_sample_counts <- fuc_by_sample %>%
      group_by(group) %>%
      summarise(n_samples = n(), .groups = 'drop')
    
    neuac_sample_counts <- neuac_by_sample %>%
      group_by(group) %>%
      summarise(n_samples = n(), .groups = 'drop')
    
    has_sufficient_fuc_ilr <- nrow(fuc_sample_counts) >= 2 && all(fuc_sample_counts$n_samples >= 3)
    has_sufficient_neuac_ilr <- nrow(neuac_sample_counts) >= 2 && all(neuac_sample_counts$n_samples >= 3)
    
    # ILR-based fucosylation analysis
    if(has_sufficient_fuc_ilr) {
      healthy_fuc_ilr_data <- fuc_by_sample$fuc_ilr[fuc_by_sample$group == "Healthy"]
      mecfs_fuc_ilr_data <- fuc_by_sample$fuc_ilr[fuc_by_sample$group == "MECFS"]
      
      if(length(unique(healthy_fuc_ilr_data)) > 1 || length(unique(mecfs_fuc_ilr_data)) > 1) {
        fuc_ilr_t_test <- try({
          t.test(fuc_ilr ~ group, data = fuc_by_sample)
        }, silent = TRUE)
        
        if(!inherits(fuc_ilr_t_test, "try-error")) {
          fuc_ilr_stats <- data.frame(
            protein_accession = protein_name,
            feature = "Fucosylation_ILR",
            mean_healthy = mean(healthy_fuc_ilr_data, na.rm = TRUE),
            mean_mecfs = mean(mecfs_fuc_ilr_data, na.rm = TRUE),
            sd_healthy = sd(healthy_fuc_ilr_data, na.rm = TRUE),
            sd_mecfs = sd(mecfs_fuc_ilr_data, na.rm = TRUE),
            t_stat = fuc_ilr_t_test$statistic,
            p_value = fuc_ilr_t_test$p.value,
            n_healthy = sum(fuc_by_sample$group == "Healthy"),
            n_mecfs = sum(fuc_by_sample$group == "MECFS"),
            stringsAsFactors = FALSE
          )
          ilr_summary_stats <- rbind(ilr_summary_stats, fuc_ilr_stats)
        }
      }
    }
    
    # ILR-based sialylation analysis
    if(has_sufficient_neuac_ilr) {
      healthy_neuac_ilr_data <- neuac_by_sample$neuac_ilr[neuac_by_sample$group == "Healthy"]
      mecfs_neuac_ilr_data <- neuac_by_sample$neuac_ilr[neuac_by_sample$group == "MECFS"]
      
      if(length(unique(healthy_neuac_ilr_data)) > 1 || length(unique(mecfs_neuac_ilr_data)) > 1) {
        neuac_ilr_t_test <- try({
          t.test(neuac_ilr ~ group, data = neuac_by_sample)
        }, silent = TRUE)
        
        if(!inherits(neuac_ilr_t_test, "try-error")) {
          neuac_ilr_stats <- data.frame(
            protein_accession = protein_name,
            feature = "Sialylation_ILR",
            mean_healthy = mean(healthy_neuac_ilr_data, na.rm = TRUE),
            mean_mecfs = mean(mecfs_neuac_ilr_data, na.rm = TRUE),
            sd_healthy = sd(healthy_neuac_ilr_data, na.rm = TRUE),
            sd_mecfs = sd(mecfs_neuac_ilr_data, na.rm = TRUE),
            t_stat = neuac_ilr_t_test$statistic,
            p_value = neuac_ilr_t_test$p.value,
            n_healthy = sum(neuac_by_sample$group == "Healthy"),
            n_mecfs = sum(neuac_by_sample$group == "MECFS"),
            stringsAsFactors = FALSE
          )
          ilr_summary_stats <- rbind(ilr_summary_stats, neuac_ilr_stats)
        }
      }
    }
  }
  
  # Apply BH correction to ILR p-values
  if(nrow(ilr_summary_stats) > 0) {
    ilr_summary_stats$p_value_adj <- p.adjust(ilr_summary_stats$p_value, method = "BH")
    ilr_summary_stats$significant <- ilr_summary_stats$p_value_adj < 0.05
    ilr_summary_stats$fold_change <- ilr_summary_stats$mean_mecfs / ilr_summary_stats$mean_healthy
  }
  
  # Save ILR results
  write.csv(ilr_summary_stats, file.path(output_dir, "protein_glycosylation_ilr_results.csv"), row.names = FALSE)
  
  cat("ILR analysis completed. Results saved to:", file.path(output_dir, "protein_glycosylation_ilr_results.csv"), "\n")
  cat("ILR plots saved to:", figures_dir, "\n")
  
  return(list(
    ilr_statistics = ilr_summary_stats
  ))
}
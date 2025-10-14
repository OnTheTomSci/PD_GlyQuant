# Example CLR Analysis Script
# This script demonstrates how to use the CLR transformation and statistical testing functions

# Load required libraries
library(tidyverse)
library(compositions)

# Source required functions
source("functions/peptidegroups_preprocessing.R")
source("functions/peptidegroups_compositional_analysis.R")

# ====== Step 1: Load and preprocess data ======
cat("Step 1: Loading and preprocessing peptide groups data...\n")
data <- load_and_preprocess_peptidegroups(
  study_info_path = "input_data/10S_MECFS_GPEPS_250125_StudyInformation.txt",
  peptide_groups_path = "input_data/10S_MECFS_GPEPS_250125_PeptideGroups.txt",
  glycan_class_map_path = "input_data/glycan_class_map.csv",
  fdr_threshold = 0.001
)

# ====== Step 2: Calculate CLR transformation ======
cat("\nStep 2: Calculating CLR transformation...\n")
clr_results <- calculate_multiple_proteins_glycan_clr(data$glyco_peptide_groups_long)

# ====== Step 3: Perform statistical testing ======
cat("\nStep 3: Testing CLR differences between groups...\n")
clr_stats <- test_clr_protein_composition_differences(
  clr_result = clr_results,
  group1_name = "Healthy",
  group2_name = "MECFS",
  p_adjust_method = "BH"
)

# ====== Step 4: Save results ======
cat("\nStep 4: Saving results...\n")

# Create output directory
dir.create("output_data/clr_analysis", recursive = TRUE, showWarnings = FALSE)

# Save CLR transformed data (long format)
write.csv(clr_results$clr_transformed_long, 
          "output_data/clr_analysis/clr_transformed_long.csv", 
          row.names = FALSE)

# Save CLR transformed data (wide format)
write.csv(clr_results$clr_transformed, 
          "output_data/clr_analysis/clr_transformed_wide.csv", 
          row.names = FALSE)

# Save relative abundance data
write.csv(clr_results$relative_abundance, 
          "output_data/clr_analysis/relative_abundance.csv", 
          row.names = FALSE)

# Save CLR summary statistics
write.csv(clr_results$clr_summary, 
          "output_data/clr_analysis/clr_summary.csv", 
          row.names = FALSE)

# Save protein summary
write.csv(clr_results$protein_summary, 
          "output_data/clr_analysis/protein_summary.csv", 
          row.names = FALSE)

# Save statistical test results
write.csv(clr_stats, 
          "output_data/clr_analysis/clr_statistical_tests.csv", 
          row.names = FALSE)

# Save only significant results
significant_results <- clr_stats %>%
  filter(significant == TRUE) %>%
  arrange(p_value_adj)

write.csv(significant_results, 
          "output_data/clr_analysis/clr_significant_results.csv", 
          row.names = FALSE)

# ====== Step 5: Display summary ======
cat("\n====== Analysis Summary ======\n")
cat(sprintf("Total proteins analyzed: %d\n", clr_results$transformation_info$n_proteins))
cat(sprintf("Total glycan compositions: %d\n", clr_results$transformation_info$n_compositions))
cat(sprintf("Total samples: %d\n", clr_results$transformation_info$n_samples))
cat(sprintf("\nStatistical tests performed: %d\n", nrow(clr_stats)))
cat(sprintf("Significant differences found: %d (%.1f%%)\n", 
            nrow(significant_results), 
            100 * nrow(significant_results) / nrow(clr_stats)))

if (nrow(significant_results) > 0) {
  cat("\nTop 10 most significant protein-glycan combinations:\n")
  print(significant_results %>%
          select(protein_accessions, glycan_composition, clr_difference, 
                 cohens_d, effect_size_interpretation, p_value_adj, significance_level) %>%
          head(10))
}

cat("\n====== Files saved to output_data/clr_analysis/ ======\n")
cat("- clr_transformed_long.csv: Long format CLR values\n")
cat("- clr_transformed_wide.csv: Wide format CLR values\n")
cat("- relative_abundance.csv: Relative abundances before CLR\n")
cat("- clr_summary.csv: Summary statistics for each protein-glycan\n")
cat("- protein_summary.csv: Summary by protein\n")
cat("- clr_statistical_tests.csv: All statistical test results\n")
cat("- clr_significant_results.csv: Only significant results\n")

cat("\nAnalysis complete!\n")


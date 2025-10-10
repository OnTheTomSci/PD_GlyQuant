# Example: ILR Compositional Analysis of Glycan Data
# This script demonstrates how to use the new ILR transformation functions

library(tidyverse)

# Load required modules
source("functions/peptidegroups_preprocessing.R")
source("functions/peptidegroups_compositional_analysis.R")

cat("================================================================================\n")
cat("  ILR COMPOSITIONAL ANALYSIS - EXAMPLE\n")
cat("================================================================================\n\n")

# ============================================================================
# STEP 1: Load and Prepare Data
# ============================================================================

cat("Step 1: Loading data...\n")
data <- load_and_preprocess_peptidegroups("input_data")
glyco_long <- data$glyco_peptide_groups_long

cat("  Data loaded successfully!\n")
cat(sprintf("  Total rows: %d\n", nrow(glyco_long)))
cat(sprintf("  Total proteins: %d\n\n", n_distinct(glyco_long$protein_accessions)))

# ============================================================================
# STEP 2: Identify Candidate Proteins
# ============================================================================

cat("Step 2: Finding proteins with multiple glycan compositions...\n")

protein_summary <- glyco_long %>%
  group_by(protein_accessions, gene_name) %>%
  summarise(
    n_compositions = n_distinct(glycan_composition),
    n_samples = n_distinct(sample),
    total_abundance = sum(abundance, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(n_compositions >= 3) %>%  # Need at least 3 for meaningful analysis
  arrange(desc(n_compositions))

cat(sprintf("  Found %d proteins with ≥3 glycan compositions\n\n", nrow(protein_summary)))
cat("  All proteins sorted by number of compositions:\n")
print(protein_summary)

# ============================================================================
# STEP 3: Single Protein Example - Detailed Analysis
# ============================================================================

cat("\n================================================================================\n")
cat("Step 3: Detailed analysis of one protein\n")
cat("================================================================================\n\n")

# Select the protein with most compositions
example_protein <- protein_summary$protein_accessions[1]
example_gene <- protein_summary$gene_name[1]

cat(sprintf("Analyzing: %s (%s)\n\n", example_protein, example_gene))

# Calculate ILR transformation
result <- calculate_protein_glycan_ilr(
  data_long = glyco_long,
  protein_id = example_protein
)

# Display composition summary
cat("Glycan Composition Summary:\n")
print(result$composition_summary)

# Display ILR transformed data (first 10 rows)
cat("\n\nILR-Transformed Data (first 10 samples):\n")
print(head(result$ilr_transformed, 10))

# Perform statistical tests
cat("\n\nStatistical Tests on ILR Coordinates:\n")
stats <- test_ilr_differences(result)
print(stats)

# Create and display plot
cat("\n\nCreating visualization...\n")
plot <- plot_ilr_boxplot(result, "ILR_1")
ggsave("example_ilr_plot.png", plot, width = 8, height = 6, dpi = 300)
cat("  Plot saved to: example_ilr_plot.png\n")

# ============================================================================
# STEP 4: Batch Analysis of Top Proteins
# ============================================================================

cat("\n================================================================================\n")
cat(sprintf("Step 4: Batch analysis of ALL %d proteins\n", nrow(protein_summary))) 
cat("================================================================================\n\n")

# Get all proteins that passed the filter (≥3 compositions)
all_proteins <- protein_summary$protein_accessions

cat(sprintf("Analyzing %d proteins with ≥3 glycan compositions...\n", length(all_proteins)))
cat("(This may take a few minutes)\n\n")

# Perform batch analysis
batch_results <- analyze_multiple_proteins_ilr(
  data_long = glyco_long,
  protein_ids = all_proteins,
  output_dir = "output_data/example_ilr_batch"
)

cat("\n\n=== BATCH ANALYSIS COMPLETE ===\n")
cat(sprintf("✓ Successfully analyzed %d proteins\n", length(batch_results)))
cat(sprintf("✓ Results saved to: output_data/example_ilr_batch/\n"))
cat(sprintf("✓ Created %d result files per protein (ilr_tests, relative_abundance, etc.)\n", 3))

# ============================================================================
# STEP 5: Summarize Findings
# ============================================================================

cat("\n================================================================================\n")
cat("Step 5: Summary of significant findings\n")
cat("================================================================================\n\n")

# Compile summary of significant results
summary_table <- data.frame()

for (protein_id in names(batch_results)) {
  result <- batch_results[[protein_id]]
  
  if (!is.null(result$statistical_tests)) {
    sig_tests <- result$statistical_tests %>%
      filter(significant == TRUE)
    
    if (nrow(sig_tests) > 0) {
      summary_table <- rbind(summary_table, data.frame(
        Protein = protein_id,
        N_Compositions = result$transformation_info$n_compositions,
        N_Samples = result$transformation_info$n_samples,
        N_Significant_Coords = nrow(sig_tests),
        Min_P_Adj = min(sig_tests$p_value_adj),
        Max_Abs_Difference = max(abs(sig_tests$difference))
      ))
    }
  }
}

if (nrow(summary_table) > 0) {
  cat("Proteins with significant compositional differences:\n\n")
  print(summary_table)
  
  write.csv(summary_table, "example_ilr_summary.csv", row.names = FALSE)
  cat("\n\nSummary saved to: example_ilr_summary.csv\n")
} else {
  cat("No significant compositional differences detected in analyzed proteins.\n")
}

# ============================================================================
# STEP 6: Detailed Look at Most Significant Protein
# ============================================================================

if (nrow(summary_table) > 0) {
  cat("\n================================================================================\n")
  cat("Step 6: Detailed examination of most significant protein\n")
  cat("================================================================================\n\n")
  
  # Get protein with lowest p-value
  most_sig_protein <- summary_table$Protein[which.min(summary_table$Min_P_Adj)]
  most_sig_result <- batch_results[[most_sig_protein]]
  
  cat(sprintf("Most significant protein: %s\n\n", most_sig_protein))
  
  # Show composition breakdown by group
  cat("Relative abundance by group:\n")
  group_comparison <- most_sig_result$relative_abundance %>%
    group_by(glycan_composition, group) %>%
    summarise(
      mean_rel_abund = mean(relative_abundance) * 100,
      sd_rel_abund = sd(relative_abundance) * 100,
      .groups = 'drop'
    ) %>%
    arrange(desc(mean_rel_abund))
  
  print(group_comparison)
  
  # Show significant ILR coordinates
  cat("\n\nSignificant ILR coordinates:\n")
  sig_coords <- most_sig_result$statistical_tests %>%
    filter(significant == TRUE) %>%
    arrange(p_value_adj)
  
  print(sig_coords)
  
  # Create plots for all significant coordinates
  cat("\n\nCreating plots for significant coordinates...\n")
  for (i in 1:nrow(sig_coords)) {
    coord <- sig_coords$ilr_coordinate[i]
    plot <- plot_ilr_boxplot(most_sig_result, coord)
    
    filename <- sprintf("example_%s_%s.png", most_sig_protein, coord)
    ggsave(filename, plot, width = 8, height = 6, dpi = 300)
    cat(sprintf("  Saved: %s\n", filename))
  }
}

# ============================================================================
# STEP 7: Create Volcano Plot
# ============================================================================

cat("\n================================================================================\n")
cat("Step 7: Creating volcano plot of ILR results\n") 
cat("================================================================================\n\n")

# Prepare data for volcano plot
volcano_data <- summary_table %>%
  mutate(
    Label = ifelse(Min_P_Adj < 0.05 & abs(Max_Abs_Difference) > 1, Protein, ""),
    Significant = Min_P_Adj < 0.05 & abs(Max_Abs_Difference) > 1
  )

# Create volcano plot
volcano_plot <- EnhancedVolcano(
  volcano_data,
  lab = volcano_data$Label,
  x = 'Max_Abs_Difference',
  y = 'Min_P_Adj',
  title = 'Volcano Plot of ILR Analysis Results',
  subtitle = NULL,
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 3,
  labSize = 4,
  col = c('grey', 'grey', 'grey', '#E49CB1'),
  colAlpha = 0.5,
  legendPosition = 'right',
  legendLabSize = 12,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  xlab = "Maximum Effect Size",
  ylab = "-log10(Minimum Adjusted p-value)"
)

# Save the plot
ggsave("example_ilr_volcano.png", volcano_plot, width = 120, height = 80, dpi = 300, units = "mm")
cat("Created volcano plot: example_ilr_volcano.png\n\n")


# ============================================================================
# COMPLETE
# ============================================================================

cat("\n================================================================================\n")
cat("  ANALYSIS COMPLETE\n")
cat("================================================================================\n\n")

cat("Generated files:\n")
cat("  - example_ilr_plot.png              (Single protein visualization)\n")
cat("  - example_ilr_summary.csv           (Summary table)\n")
cat("  - output_data/example_ilr_batch/    (Batch analysis results)\n")
cat("  - example_<protein>_ILR_*.png       (Significant coordinate plots)\n")
cat("\n")

# ============================================================================
# List all analyzed proteins
# ============================================================================

cat("================================================================================\n")
cat(sprintf("ALL %d PROTEINS ANALYZED:\n", length(batch_results)))
cat("================================================================================\n\n")

# Create detailed list with composition counts
analyzed_proteins <- protein_summary %>%
  mutate(
    analyzed = protein_accessions %in% names(batch_results),
    has_sig_results = sapply(protein_accessions, function(p) {
      if (p %in% names(batch_results)) {
        result <- batch_results[[p]]
        if (!is.null(result$statistical_tests)) {
          return(any(result$statistical_tests$significant, na.rm = TRUE))
        }
      }
      return(FALSE)
    })
  )

# Print summary
cat(sprintf("Total proteins with ≥3 compositions: %d\n", nrow(analyzed_proteins)))
cat(sprintf("Successfully analyzed: %d\n", sum(analyzed_proteins$analyzed)))
cat(sprintf("With significant results: %d\n\n", sum(analyzed_proteins$has_sig_results)))

# Print list
cat("Protein List (showing composition count and significance):\n")
cat("------------------------------------------------------------\n")
for (i in 1:nrow(analyzed_proteins)) {
  sig_marker <- ifelse(analyzed_proteins$has_sig_results[i], "***", "   ")
  cat(sprintf("%2d. %s %s | %s | %d compositions\n",
              i,
              sig_marker,
              analyzed_proteins$protein_accessions[i],
              ifelse(is.na(analyzed_proteins$gene_name[i]), "Unknown", analyzed_proteins$gene_name[i]),
              analyzed_proteins$n_compositions[i]))
}

cat("\n*** = Has significant ILR coordinate differences\n\n")

# Save batch_results for later use
cat("Saving batch_results for visualization scripts...\n")
saveRDS(batch_results, "batch_results.rds")
cat("  ✓ Saved: batch_results.rds\n\n")

cat("Next steps:\n")
cat("  1. Review the summary table (example_ilr_summary.csv)\n")
cat("  2. Check proteins marked with *** for significant changes\n")
cat("  3. Examine individual protein plots in output_data/example_ilr_batch/\n")
cat("  4. Run visualize_all_proteins.R to create comprehensive visualizations\n")
cat("     → source('visualize_all_proteins.R')\n")
cat("\n")


# CLR Analysis Visualization Script
# Demonstrates different ways to visualize CLR statistical test results

library(tidyverse)
library(ggplot2)
library(patchwork)

# Source required functions
source("functions/peptidegroups_preprocessing.R")
source("functions/peptidegroups_compositional_analysis.R")

cat("================================================================================\n")
cat("  CLR COMPOSITIONAL ANALYSIS - VISUALIZATION WORKFLOW\n")
cat("================================================================================\n\n")

# ============================================================================
# STEP 1: Load data and run CLR analysis
# ============================================================================

cat("Step 1: Loading data and performing CLR analysis...\n")

# Load preprocessed data
data <- load_and_preprocess_peptidegroups(
  study_info_path = "input_data/10S_MECFS_GPEPS_250125_StudyInformation.txt",
  peptide_groups_path = "input_data/10S_MECFS_GPEPS_250125_PeptideGroups.txt",
  glycan_class_map_path = "input_data/glycan_class_map.csv",
  fdr_threshold = 0.001
)

# Calculate CLR transformation
cat("\nCalculating CLR transformation...\n")
clr_results <- calculate_multiple_proteins_glycan_clr(data$glyco_peptide_groups_long)

# Perform statistical testing
cat("\nTesting CLR differences between groups...\n")
clr_stats <- test_clr_protein_composition_differences(
  clr_result = clr_results,
  group1_name = "Healthy",
  group2_name = "MECFS",
  p_adjust_method = "BH"
)

cat(sprintf("\nAnalysis complete: %d protein-glycan combinations tested\n\n", nrow(clr_stats)))

# ============================================================================
# STEP 2: Create dot plot (primary visualization)
# ============================================================================

cat("================================================================================\n")
cat("VISUALIZATION 1: DOT PLOT\n")
cat("Shows all proteins and glycan compositions\n")
cat("================================================================================\n\n")

# Full dot plot with proteins that have at least 3 samples per group
dotplot_all <- plot_clr_dotplot(
  clr_stats_result = clr_stats,
  p_threshold = 0.05,
  top_n_proteins = NULL,  # Show all proteins
  min_samples_per_group = 3  # Require at least 3 samples per group
)

ggsave("clr_dotplot_all.png", dotplot_all, 
       width = 16, height = 12, dpi = 300)
cat("✓ Saved: clr_dotplot_all.png\n")
cat("  - Dot size = Magnitude of CLR difference\n")
cat("  - Dot color = Direction (pink = increased, teal = decreased)\n")
cat("  - Opacity = Significance (solid = significant)\n")
cat("  - Filtered to proteins with ≥3 glycans (each in ≥3 samples per group)\n\n")

# Dot plot with top 30 proteins (only if there are significant results)
n_significant <- sum(clr_stats$significant, na.rm = TRUE)
if (n_significant > 0) {
  dotplot_top30 <- plot_clr_dotplot(
    clr_stats_result = clr_stats,
    p_threshold = 0.05,
    top_n_proteins = 30  # Show top 30 proteins
  )
  
  ggsave("clr_dotplot_top30.png", dotplot_top30, 
         width = 14, height = 10, dpi = 300)
  cat("✓ Saved: clr_dotplot_top30.png\n")
  cat("  - Shows top 30 proteins by number of significant glycans\n\n")
} else {
  cat("⚠ Skipped: clr_dotplot_top30.png (no significant results)\n\n")
}

# ============================================================================
# STEP 3: Create heatmap
# ============================================================================

cat("================================================================================\n")
cat("VISUALIZATION 2: HEATMAP\n")
cat("Shows CLR differences as a heatmap\n")
cat("================================================================================\n\n")

# Heatmap colored by CLR difference
heatmap_clr <- plot_clr_heatmap(
  clr_stats_result = clr_stats,
  value_type = "clr_difference",
  p_threshold = 0.05,
  top_n_proteins = 30
)

ggsave("clr_heatmap_difference.png", heatmap_clr,
       width = 12, height = 10, dpi = 300)
cat("✓ Saved: clr_heatmap_difference.png\n")
cat("  - Color = CLR difference (blue = decreased, red = increased)\n")
cat("  - Star (*) = Significant (p-adj < 0.05)\n\n")

# Heatmap colored by p-value
heatmap_pval <- plot_clr_heatmap(
  clr_stats_result = clr_stats,
  value_type = "p_value",
  p_threshold = 0.05,
  top_n_proteins = 30
)

ggsave("clr_heatmap_pvalue.png", heatmap_pval,
       width = 12, height = 10, dpi = 300)
cat("✓ Saved: clr_heatmap_pvalue.png\n")
cat("  - Color = -log10(p-value) (darker = more significant)\n\n")

# ============================================================================
# STEP 4: Create forest plot
# ============================================================================

cat("================================================================================\n")
cat("VISUALIZATION 3: FOREST PLOT\n")
cat("Shows CLR differences with 95% confidence intervals\n")
cat("================================================================================\n\n")

# Forest plot with significant results only
forest_sig <- plot_clr_forest(
  clr_stats_result = clr_stats,
  show_only_significant = TRUE,
  p_threshold = 0.05,
  top_n = 50  # Top 50 most significant
)

if (!is.null(forest_sig)) {
  ggsave("clr_forest_significant.png", forest_sig,
         width = 10, height = 14, dpi = 300)
  cat("✓ Saved: clr_forest_significant.png\n")
  cat("  - Dots = CLR difference\n")
  cat("  - Lines = 95% confidence intervals\n")
  cat("  - Vertical line at 0 = no change\n\n")
} else {
  cat("  No significant results to plot\n\n")
}

# Forest plot with all results (top 50 by p-value)
forest_all <- plot_clr_forest(
  clr_stats_result = clr_stats,
  show_only_significant = FALSE,
  p_threshold = 0.05,
  top_n = 50
)

if (!is.null(forest_all)) {
  ggsave("clr_forest_all.png", forest_all,
         width = 10, height = 14, dpi = 300)
  cat("✓ Saved: clr_forest_all.png\n\n")
}

# ============================================================================
# STEP 5: Create comprehensive summary plot
# ============================================================================

cat("================================================================================\n")
cat("VISUALIZATION 4: COMPREHENSIVE SUMMARY\n")
cat("Combines heatmap, forest plot, and dot plot\n")
cat("================================================================================\n\n")

summary_plot <- plot_clr_summary(
  clr_stats_result = clr_stats,
  output_file = "clr_comprehensive_summary.png",
  p_threshold = 0.05,
  top_n_proteins = 30,
  min_samples_per_group = 3
)

cat("✓ Saved: clr_comprehensive_summary.png\n")
cat("  - Combined view of all visualization types\n\n")

# ============================================================================
# STEP 6: Save results tables
# ============================================================================

cat("================================================================================\n")
cat("SAVING RESULTS TABLES\n")
cat("================================================================================\n\n")

# Save all results
write.csv(clr_stats, "clr_all_results_table.csv", row.names = FALSE)
cat("✓ Saved: clr_all_results_table.csv\n")

# Save significant results only
significant_results <- clr_stats %>%
  filter(significant == TRUE) %>%
  arrange(p_value_adj)

write.csv(significant_results, "clr_significant_results_table.csv", row.names = FALSE)
cat("✓ Saved: clr_significant_results_table.csv\n")

# Create summary by protein
protein_summary <- clr_stats %>%
  group_by(protein_accessions) %>%
  summarise(
    n_glycans_tested = n(),
    n_significant = sum(significant, na.rm = TRUE),
    pct_significant = 100 * n_significant / n_glycans_tested,
    mean_abs_clr_diff = mean(abs(clr_difference), na.rm = TRUE),
    max_abs_clr_diff = max(abs(clr_difference), na.rm = TRUE),
    min_p_value_adj = min(p_value_adj, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(n_significant), min_p_value_adj)

write.csv(protein_summary, "clr_protein_summary.csv", row.names = FALSE)
cat("✓ Saved: clr_protein_summary.csv\n")

# ============================================================================
# STEP 7: Print summary statistics
# ============================================================================

cat("\n================================================================================\n")
cat("SUMMARY STATISTICS\n")
cat("================================================================================\n\n")

cat(sprintf("Total protein-glycan combinations tested: %d\n", nrow(clr_stats)))
cat(sprintf("Significant differences (p-adj < 0.05): %d (%.1f%%)\n",
            nrow(significant_results),
            100 * nrow(significant_results) / nrow(clr_stats)))

if (nrow(significant_results) > 0) {
  cat("\nSignificance levels:\n")
  cat(sprintf("  p < 0.001 (***): %d\n", sum(clr_stats$significance_level == "***", na.rm = TRUE)))
  cat(sprintf("  p < 0.01 (**): %d\n", sum(clr_stats$significance_level == "**", na.rm = TRUE)))
  cat(sprintf("  p < 0.05 (*): %d\n", sum(clr_stats$significance_level == "*", na.rm = TRUE)))
  
  cat("\nEffect sizes (among significant):\n")
  cat(sprintf("  Large: %d\n", sum(significant_results$effect_size_interpretation == "large", na.rm = TRUE)))
  cat(sprintf("  Medium: %d\n", sum(significant_results$effect_size_interpretation == "medium", na.rm = TRUE)))
  cat(sprintf("  Small: %d\n", sum(significant_results$effect_size_interpretation == "small", na.rm = TRUE)))
  
  cat("\nTop 10 most significant protein-glycan combinations:\n")
  print(significant_results %>%
          select(protein_accessions, glycan_composition, clr_difference, 
                 cohens_d, p_value_adj, significance_level) %>%
          head(10))
  
  cat("\nProteins with most significant glycan differences:\n")
  print(protein_summary %>%
          filter(n_significant > 0) %>%
          select(protein_accessions, n_glycans_tested, n_significant, 
                 pct_significant, min_p_value_adj) %>%
          head(10))
}

cat("\n================================================================================\n")
cat("VISUALIZATION COMPLETE!\n")
cat("================================================================================\n\n")

cat("Generated files:\n")
cat("  Dot plots:\n")
cat("    - clr_dotplot_all.png\n")
cat("    - clr_dotplot_top30.png\n")
cat("  Heatmaps:\n")
cat("    - clr_heatmap_difference.png\n")
cat("    - clr_heatmap_pvalue.png\n")
cat("  Forest plots:\n")
cat("    - clr_forest_significant.png\n")
cat("    - clr_forest_all.png\n")
cat("  Summary:\n")
cat("    - clr_comprehensive_summary.png\n")
cat("  Tables:\n")
cat("    - clr_all_results_table.csv\n")
cat("    - clr_significant_results_table.csv\n")
cat("    - clr_protein_summary.csv\n\n")


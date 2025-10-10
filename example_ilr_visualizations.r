# Example: Visualizing ILR Analysis Results
# Demonstrates different ways to plot all ILR coordinates together

library(tidyverse)
library(patchwork)  # For combining plots

# Load required modules
source("functions/peptidegroups_preprocessing.R")
source("functions/peptidegroups_compositional_analysis.R")
source("functions/peptidegroups_ilr_visualization.R")

cat("================================================================================\n")
cat("  ILR RESULTS VISUALIZATION EXAMPLES\n")
cat("================================================================================\n\n")

# ============================================================================
# STEP 1: Load existing results or run analysis
# ============================================================================

cat("Step 1: Loading or generating ILR results...\n\n")

# Option A: If you already ran the analysis, load the saved results
# (The analyze_multiple_proteins_ilr function returns results in memory)
# For this example, we'll run a quick analysis

data <- load_and_preprocess_peptidegroups("input_data")
glyco_long <- data$glyco_peptide_groups_long

# Find proteins with multiple compositions
protein_candidates <- glyco_long %>%
  group_by(protein_accessions) %>%
  summarise(n_compositions = n_distinct(glycan_composition), .groups = 'drop') %>%
  filter(n_compositions >= 3) %>%
  head(10)

# Run ILR analysis (or use your existing results)
cat("Running ILR analysis on sample proteins...\n")
ilr_results <- analyze_multiple_proteins_ilr(
  glyco_long,
  protein_candidates$protein_accessions,
  output_dir = NULL  # Don't save files for this example
)

cat(sprintf("Analyzed %d proteins\n\n", length(ilr_results)))

# ============================================================================
# STEP 2: HEATMAP - Best for overview of all coordinates
# ============================================================================

cat("================================================================================\n")
cat("VISUALIZATION 1: HEATMAP\n")
cat("Shows all proteins × ILR coordinates with effect sizes\n")
cat("================================================================================\n\n")

# Create heatmap showing effect sizes
heatmap_effect <- plot_ilr_heatmap(
  ilr_results,
  value_type = "effect",  # Color by effect size
  p_threshold = 0.05
)

# Save
ggsave("ilr_heatmap_effect.png", heatmap_effect, 
       width = 10, height = 8, dpi = 300)
cat("✓ Saved: ilr_heatmap_effect.png\n")
cat("  - Rows = Proteins\n")
cat("  - Columns = ILR Coordinates\n")
cat("  - Color = Effect size (blue = decreased, red = increased)\n")
cat("  - Star (*) = Significant (p < 0.05)\n\n")

# Alternative: Heatmap showing p-values
heatmap_pval <- plot_ilr_heatmap(
  ilr_results,
  value_type = "p_value",  # Color by p-value
  p_threshold = 0.05
)

ggsave("ilr_heatmap_pvalue.png", heatmap_pval,
       width = 10, height = 8, dpi = 300)
cat("✓ Saved: ilr_heatmap_pvalue.png\n")
cat("  - Color = -log10(p-value) (darker = more significant)\n\n")

# ============================================================================
# STEP 3: FOREST PLOT - Best for showing effect sizes with uncertainty
# ============================================================================

cat("================================================================================\n")
cat("VISUALIZATION 2: FOREST PLOT\n")
cat("Shows effect sizes with 95% confidence intervals\n")
cat("================================================================================\n\n")

# Forest plot with all results
forest_all <- plot_ilr_forest(
  ilr_results,
  show_only_significant = FALSE
)

ggsave("ilr_forest_all.png", forest_all,
       width = 8, height = 12, dpi = 300)
cat("✓ Saved: ilr_forest_all.png\n")
cat("  - Dots = Effect size\n")
cat("  - Lines = 95% confidence intervals\n")
cat("  - Vertical line at 0 = no change\n\n")

# Forest plot with only significant results
forest_sig <- plot_ilr_forest(
  ilr_results,
  show_only_significant = TRUE,
  p_threshold = 0.05
)

if (!is.null(forest_sig)) {
  ggsave("ilr_forest_significant.png", forest_sig,
         width = 8, height = 6, dpi = 300)
  cat("✓ Saved: ilr_forest_significant.png\n")
  cat("  - Shows only significant results\n\n")
} else {
  cat("  No significant results to plot\n\n")
}

# ============================================================================
# STEP 4: DOT PLOT - Best for showing patterns across proteins
# ============================================================================

cat("================================================================================\n")
cat("VISUALIZATION 3: DOT PLOT\n")
cat("Shows direction and magnitude of changes\n")
cat("================================================================================\n\n")

dotplot <- plot_ilr_dotplot(
  ilr_results,
  p_threshold = 0.05
)

ggsave("ilr_dotplot.png", dotplot,
       width = 200, height = 240, dpi = 300, units = "mm")
cat("✓ Saved: ilr_dotplot.png\n")
cat("  - Dot size = Magnitude of effect\n")
cat("  - Dot color = Direction (pink = increased, teal = decreased)\n")
cat("  - Opacity = Significance (solid = significant)\n\n")

# ============================================================================
# STEP 5: BOXPLOTS - Best for detailed inspection of significant results
# ============================================================================

cat("================================================================================\n")
cat("VISUALIZATION 4: BOXPLOTS FOR SIGNIFICANT COORDINATES\n")
cat("Individual plots showing actual data distribution\n")
cat("================================================================================\n\n")

boxplots <- plot_ilr_significant_boxplots(
  ilr_results,
  p_threshold = 0.05,
  max_plots = 9  # Create up to 9 plots
)

if (!is.null(boxplots) && length(boxplots) > 0) {
  # Combine into grid using patchwork
  if (length(boxplots) <= 9) {
    ncol <- 3
    combined_boxplots <- wrap_plots(boxplots, ncol = ncol)
    
    ggsave("ilr_boxplots_grid.png", combined_boxplots,
           width = 12, height = 4 * ceiling(length(boxplots) / ncol), dpi = 300)
    cat(sprintf("✓ Saved: ilr_boxplots_grid.png (%d plots)\n", length(boxplots)))
    cat("  - Shows distribution of ILR values by group\n")
    cat("  - Only significant protein-coordinate combinations\n\n")
  }
  
  # Also save individual plots
  dir.create("ilr_individual_boxplots", showWarnings = FALSE)
  for (name in names(boxplots)) {
    ggsave(file.path("ilr_individual_boxplots", paste0(name, ".png")),
           boxplots[[name]], width = 6, height = 5, dpi = 300)
  }
  cat(sprintf("✓ Saved: %d individual boxplots to ilr_individual_boxplots/\n\n", 
              length(boxplots)))
} else {
  cat("  No significant results to plot\n\n")
}

# ============================================================================
# STEP 6: COMPREHENSIVE SUMMARY - Combines multiple views
# ============================================================================

cat("================================================================================\n")
cat("VISUALIZATION 5: COMPREHENSIVE SUMMARY PLOT\n")
cat("Combines heatmap, forest plot, and dot plot\n")
cat("================================================================================\n\n")

summary_plot <- plot_ilr_summary(
  ilr_results,
  output_file = "ilr_comprehensive_summary.png",
  p_threshold = 0.05
)

cat("✓ Saved: ilr_comprehensive_summary.png\n")
cat("  - Multi-panel figure with all views\n")
cat("  - Best for presentations or reports\n\n")

# ============================================================================
# STEP 7: RESULTS TABLE - For detailed inspection
# ============================================================================

cat("================================================================================\n")
cat("EXPORT: COMPLETE RESULTS TABLE\n")
cat("================================================================================\n\n")

results_table <- create_ilr_results_table(
  ilr_results,
  output_file = "ilr_complete_results_table.csv"
)

cat("✓ Saved: ilr_complete_results_table.csv\n")
cat(sprintf("  - %d rows (all protein-coordinate combinations)\n", nrow(results_table)))
cat("  - Sortable and filterable in Excel\n\n")

# Print summary
sig_count <- sum(results_table$significant, na.rm = TRUE)
cat(sprintf("Summary: %d/%d coordinate tests were significant (%.1f%%)\n", 
            sig_count, nrow(results_table), 
            100 * sig_count / nrow(results_table)))

# ============================================================================
# RECOMMENDATIONS BY USE CASE
# ============================================================================

cat("\n================================================================================\n")
cat("RECOMMENDATIONS FOR YOUR USE CASE\n")
cat("================================================================================\n\n")

cat("FOR PRESENTATIONS:\n")
cat("  → Use: ilr_comprehensive_summary.png\n")
cat("  → Shows everything in one figure\n\n")

cat("FOR MANUSCRIPTS:\n")
cat("  → Main figure: ilr_heatmap_effect.png\n")
cat("  → Supplement: ilr_forest_significant.png\n")
cat("  → Supplement: ilr_boxplots_grid.png\n\n")

cat("FOR EXPLORATION:\n")
cat("  → Start with: ilr_dotplot.png (see patterns)\n")
cat("  → Then: ilr_heatmap_effect.png (identify specific coordinates)\n")
cat("  → Finally: Individual boxplots (detailed inspection)\n\n")

cat("FOR IDENTIFYING CANDIDATES:\n")
cat("  → ilr_forest_significant.png\n")
cat("  → Sorted by effect size, easy to spot top hits\n\n")

cat("FOR DETAILED ANALYSIS:\n")
cat("  → ilr_complete_results_table.csv\n")
cat("  → Filter and sort in Excel/R for specific proteins\n\n")

# ============================================================================
# INTERPRETATION GUIDE
# ============================================================================

cat("================================================================================\n")
cat("INTERPRETATION GUIDE\n")
cat("================================================================================\n\n")

cat("WHAT EACH PLOT TELLS YOU:\n\n")

cat("1. HEATMAP:\n")
cat("   • Quick overview of which protein-coordinate combinations changed\n")
cat("   • Color shows direction (blue/red) and magnitude\n")
cat("   • Stars show statistical significance\n")
cat("   • Patterns: Clusters of proteins with similar changes\n\n")

cat("2. FOREST PLOT:\n")
cat("   • Effect sizes with uncertainty (confidence intervals)\n")
cat("   • Wide intervals = less certain estimates\n")
cat("   • Intervals crossing zero = not significant\n")
cat("   • Sorted by effect size for easy ranking\n\n")

cat("3. DOT PLOT:\n")
cat("   • Size = strength of effect\n")
cat("   • Color = direction of change\n")
cat("   • Opacity = significance level\n")
cat("   • Good for spotting coordinate-specific patterns\n\n")

cat("4. BOXPLOTS:\n")
cat("   • Actual distribution of ILR values\n")
cat("   • Shows group separation visually\n")
cat("   • Individual points reveal outliers\n")
cat("   • Best for detailed inspection of top hits\n\n")

# ============================================================================
# COMPLETE
# ============================================================================

cat("================================================================================\n")
cat("VISUALIZATION COMPLETE\n")
cat("================================================================================\n\n")

cat("Generated files:\n")
cat(" ilr_heatmap_effect.png\n")
cat(" ilr_heatmap_pvalue.png\n")
cat(" ilr_forest_all.png\n")
cat("  ilr_forest_significant.png\n")
cat(" ilr_dotplot.png\n")
cat(" ilr_boxplots_grid.png\n")
cat("  ilr_comprehensive_summary.png\n")
cat("  ilr_individual_boxplots/ (folder)\n")
cat("  ilr_complete_results_table.csv\n\n")

cat("Next steps:\n")
cat("  1. Open ilr_comprehensive_summary.png for overview\n")
cat("  2. Check ilr_complete_results_table.csv for details\n")
cat("  3. Examine individual boxplots for significant hits\n")
cat("  4. Use findings to select proteins for validation\n\n")


# Quick Visualization of All 53 Protein ILR Results
# Works with existing CSV files in output_data/example_ilr_batch/

library(tidyverse)
library(patchwork)

cat("================================================================================\n")
cat("  VISUALIZING ALL 53 PROTEIN ILR RESULTS\n")
cat("  (Using saved CSV files)\n")
cat("================================================================================\n\n")

# ============================================================================
# STEP 1: Load all results from CSV files
# ============================================================================

cat("Step 1: Loading results from CSV files...\n")

# Find all ilr_tests.csv files
test_files <- list.files("output_data/example_ilr_batch", 
                         pattern = "_ilr_tests.csv$", 
                         full.names = TRUE)

cat(sprintf("  Found %d protein result files\n\n", length(test_files)))

# Load all test results
all_results <- bind_rows(lapply(test_files, function(file) {
  protein_id <- gsub("_ilr_tests.csv", "", basename(file))
  
  tests <- read.csv(file, stringsAsFactors = FALSE)
  
  if (nrow(tests) > 0) {
    tests$protein <- protein_id
    return(tests)
  } else {
    return(NULL)
  }
}))

cat(sprintf("  Loaded %d total ILR coordinate tests\n", nrow(all_results)))
cat(sprintf("  From %d proteins\n\n", n_distinct(all_results$protein)))

# ============================================================================
# STEP 2: Create Heatmap - Effect Sizes
# ============================================================================

cat("Creating Heatmap of Effect Sizes...\n")

# Prepare heatmap data
heatmap_data <- all_results %>%
  mutate(
    significant = p_value_adj < 0.05,
    sig_marker = ifelse(significant, "*", ""),
    # Extract numeric part for proper ordering
    ilr_num = as.numeric(str_extract(ilr_coordinate, "\\d+"))
  ) %>%
  # Order ILR coordinates numerically (ILR_1, ILR_2, ... not ILR_1, ILR_10, ILR_2)
  mutate(ilr_coordinate = factor(ilr_coordinate, 
                                 levels = unique(ilr_coordinate[order(ilr_num)])))

# Create heatmap
heatmap <- ggplot(heatmap_data, 
                  aes(x = ilr_coordinate, y = protein, fill = difference)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = sig_marker), size = 5, color = "black") +
  scale_fill_gradient2(
    low = "#2166AC", 
    mid = "white", 
    high = "#B2182B",
    midpoint = 0,
    name = "Effect Size\n(MECFS - Healthy)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = "ILR Coordinate Analysis - All 53 Proteins",
    subtitle = "* indicates p-adj < 0.05",
    x = "ILR Coordinate",
    y = "Protein"
  )

ggsave("heatmap_all_53_proteins.png", heatmap, 
       width = 12, height = max(10, n_distinct(all_results$protein) * 0.25), 
       dpi = 300)

cat("  ✓ Saved: heatmap_all_53_proteins.png\n\n")

# ============================================================================
# STEP 3: Create Forest Plot - Significant Results Only
# ============================================================================

cat("Creating Forest Plot (Significant Results)...\n")

sig_results <- all_results %>%
  filter(p_value_adj < 0.05) %>%
  mutate(protein_coord = paste(protein, ilr_coordinate, sep = " - ")) %>%
  arrange(desc(abs(difference)))

if (nrow(sig_results) > 0) {
  sig_results$protein_coord <- factor(sig_results$protein_coord, 
                                      levels = sig_results$protein_coord)
  
  forest <- ggplot(sig_results, 
                   aes(x = difference, y = protein_coord)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), 
                   height = 0.3, alpha = 0.7, color = "#E49CB1") +
    geom_point(size = 3, color = "#E49CB1") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
      title = "Significant ILR Coordinates - Effect Sizes",
      subtitle = sprintf("%d significant results (p-adj < 0.05)", nrow(sig_results)),
      x = "Effect Size (MECFS - Healthy)",
      y = "Protein - ILR Coordinate"
    )
  
  ggsave("forest_plot_significant_all.png", forest, 
         width = 10, height = max(6, nrow(sig_results) * 0.25), 
         dpi = 300)
  
  cat(sprintf("  ✓ Saved: forest_plot_significant_all.png\n"))
  cat(sprintf("    - Shows %d significant results\n\n", nrow(sig_results)))
} else {
  cat("  ! No significant results found\n\n")
}

# ============================================================================
# STEP 4: Create Dot Plot
# ============================================================================

cat("Creating Dot Plot...\n")

plot_data <- all_results %>%
  mutate(
    significant = p_value_adj < 0.05,
    direction = case_when(
      !significant ~ "Not Significant",
      difference > 0 ~ "Increased in MECFS",
      difference < 0 ~ "Decreased in MECFS"
    ),
    abs_difference = abs(difference),
    # Extract numeric part for proper ordering
    ilr_num = as.numeric(str_extract(ilr_coordinate, "\\d+"))
  ) %>%
  # Order ILR coordinates numerically
  mutate(ilr_coordinate = factor(ilr_coordinate, 
                                 levels = unique(ilr_coordinate[order(ilr_num)])))

dotplot <- ggplot(plot_data, 
                  aes(x = ilr_coordinate, y = protein)) +
  geom_point(aes(size = abs_difference, 
                 color = direction,
                 alpha = significant)) +
  scale_size_continuous(
    name = "Abs(Effect Size)",
    range = c(2, 10)
  ) +
  scale_color_manual(
    values = c("Increased in MECFS" = "#E49CB1",
               "Decreased in MECFS" = "#9DD4CC",
               "Not Significant" = "gray80"),
    name = "Direction"
  ) +
  scale_alpha_manual(
    values = c("TRUE" = 1, "FALSE" = 0.3),
    guide = "none"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "ILR Coordinate Analysis - Dot Plot (All 53 Proteins)",
    subtitle = "Larger dots = larger effect, Solid = significant",
    x = "ILR Coordinate",
    y = "Protein"
  )

ggsave("dotplot_all_53_proteins.png", dotplot, 
       width = 200, height = max(240, n_distinct(plot_data$protein) * 0.25), 
       dpi = 300, units = "mm")

cat("  ✓ Saved: dotplot_all_53_proteins.png\n\n")

# ============================================================================
# STEP 5: Summary Statistics
# ============================================================================

cat("================================================================================\n")
cat("  SUMMARY STATISTICS\n")
cat("================================================================================\n\n")

n_proteins <- n_distinct(all_results$protein)
n_total_tests <- nrow(all_results)
n_significant <- sum(all_results$p_value_adj < 0.05, na.rm = TRUE)

cat(sprintf("Total proteins: %d\n", n_proteins))
cat(sprintf("Total ILR coordinate tests: %d\n", n_total_tests))
cat(sprintf("Significant tests (p-adj < 0.05): %d (%.1f%%)\n\n", 
            n_significant, 100 * n_significant / n_total_tests))

# Proteins with significant results
proteins_with_sig <- all_results %>%
  filter(p_value_adj < 0.05) %>%
  group_by(protein) %>%
  summarise(
    n_sig_coords = n(),
    min_p = min(p_value_adj),
    max_effect = max(abs(difference)),
    mean_effect = mean(abs(difference)),
    .groups = 'drop'
  ) %>%
  arrange(min_p)

cat(sprintf("Proteins with significant results: %d/%d (%.1f%%)\n\n", 
            nrow(proteins_with_sig), n_proteins,
            100 * nrow(proteins_with_sig) / n_proteins))

cat("Top 10 proteins by significance:\n")
cat("────────────────────────────────────────────────────────────────────\n")
top_10 <- head(proteins_with_sig, 10)
for (i in 1:nrow(top_10)) {
  cat(sprintf("%2d. %s | %d coords | p=%.2e | max effect=%.2f\n",
              i,
              top_10$protein[i],
              top_10$n_sig_coords[i],
              top_10$min_p[i],
              top_10$max_effect[i]))
}
cat("\n")

# Export summary
write.csv(proteins_with_sig, "proteins_with_significant_ilr.csv", row.names = FALSE)
cat("  ✓ Saved: proteins_with_significant_ilr.csv\n\n")

# Export complete table
write.csv(all_results, "all_ilr_results_table.csv", row.names = FALSE)
cat("  ✓ Saved: all_ilr_results_table.csv\n")
cat(sprintf("    - %d rows with all protein-coordinate tests\n\n", nrow(all_results)))

# ============================================================================
# COMPLETE
# ============================================================================

cat("================================================================================\n")
cat("  VISUALIZATION COMPLETE\n")
cat("================================================================================\n\n")

cat("Generated visualizations:\n")
cat("  📊 heatmap_all_53_proteins.png          - Overview of all results\n")
cat("  📊 forest_plot_significant_all.png      - Significant results ranked\n")
cat("  📊 dotplot_all_53_proteins.png          - Pattern visualization\n")
cat("  📋 all_ilr_results_table.csv            - Complete results table\n")
cat("  📋 proteins_with_significant_ilr.csv    - Summary of significant proteins\n\n")

cat("Recommended viewing order:\n")
cat("  1. heatmap_all_53_proteins.png    - See all results at once\n")
cat("  2. forest_plot_significant_all.png - Check significant hits\n")
cat("  3. dotplot_all_53_proteins.png     - Understand patterns\n")
cat("  4. Check CSV files for detailed numbers\n\n")

cat("Next steps:\n")
cat("  • Open the PNG files to view all protein results\n")
cat("  • Review proteins_with_significant_ilr.csv for top candidates\n")
cat("  • Examine detailed plots in output_data/example_ilr_batch/\n\n")


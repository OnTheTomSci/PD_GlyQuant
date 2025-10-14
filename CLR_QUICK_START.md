# CLR Analysis - Quick Start Guide

## Overview
This guide provides a quick reference for performing CLR (Centered Log-Ratio) compositional analysis on glycan data.

## Files Created
1. **Functions**: `functions/peptidegroups_compositional_analysis.R` (updated with new functions)
2. **Example Analysis**: `example_clr_analysis.R`
3. **Visualization Script**: `CLR_visualisations.R`
4. **Documentation**: `CLR_ANALYSIS_README.md`

---

## Quick Usage

### Option 1: Run Complete Analysis (Recommended)

```r
# Run the complete CLR analysis workflow
source("CLR_visualisations.R")
```

This will:
- Load and preprocess your data
- Calculate CLR transformations
- Perform statistical testing
- Generate all visualizations
- Save all results to CSV files

### Option 2: Step-by-Step Analysis

```r
library(tidyverse)
library(compositions)

# Load functions
source("functions/peptidegroups_preprocessing.R")
source("functions/peptidegroups_compositional_analysis.R")

# 1. Load data
data <- load_and_preprocess_peptidegroups()

# 2. Calculate CLR transformation
clr_results <- calculate_multiple_proteins_glycan_clr(
  data$glyco_peptide_groups_long
)

# 3. Test differences between groups
clr_stats <- test_clr_protein_composition_differences(
  clr_result = clr_results,
  group1_name = "Healthy",
  group2_name = "MECFS"
)

# 4. Create visualizations
dotplot <- plot_clr_dotplot(clr_stats, top_n_proteins = 30)
heatmap <- plot_clr_heatmap(clr_stats, top_n_proteins = 30)
forest <- plot_clr_forest(clr_stats, top_n = 50)
summary <- plot_clr_summary(clr_stats, output_file = "clr_summary.png")

# 5. Save results
write.csv(clr_stats, "clr_results.csv", row.names = FALSE)
```

---

## Available Functions

### Analysis Functions

1. **`calculate_multiple_proteins_glycan_clr(data_long)`**
   - Calculates relative abundances and CLR transformation
   - Input: Long format data with columns: `abundance`, `sample`, `group`, `protein_accessions`, `glycan_composition`
   - Output: List with CLR-transformed data and summaries

2. **`test_clr_protein_composition_differences(clr_result, group1_name, group2_name)`**
   - Tests CLR differences between two groups
   - Performs t-tests with multiple testing correction
   - Calculates effect sizes (Cohen's d)
   - Output: Data frame with statistical test results

### Visualization Functions

3. **`plot_clr_dotplot(clr_stats_result, p_threshold, top_n_proteins)`**
   - **Best for**: Overview of all protein-glycan combinations
   - Shows magnitude and direction of changes

4. **`plot_clr_heatmap(clr_stats_result, value_type, p_threshold, top_n_proteins)`**
   - **Best for**: Pattern recognition across proteins
   - Can show CLR differences or p-values

5. **`plot_clr_forest(clr_stats_result, show_only_significant, p_threshold, top_n)`**
   - **Best for**: Detailed view of effect sizes with uncertainty
   - Shows confidence intervals

6. **`plot_clr_summary(clr_stats_result, output_file, p_threshold, top_n_proteins)`**
   - **Best for**: Comprehensive overview
   - Combines all three visualization types

---

## Output Files

### From `CLR_visualisations.R`:

**Visualizations:**
- `clr_dotplot_all.png` - All proteins
- `clr_dotplot_top30.png` - Top 30 proteins
- `clr_heatmap_difference.png` - Heatmap of CLR differences
- `clr_heatmap_pvalue.png` - Heatmap of p-values
- `clr_forest_significant.png` - Forest plot of significant results
- `clr_forest_all.png` - Forest plot of all results
- `clr_comprehensive_summary.png` - Combined view

**Data Tables:**
- `clr_all_results_table.csv` - All statistical test results
- `clr_significant_results_table.csv` - Significant results only
- `clr_protein_summary.csv` - Summary by protein

### From `example_clr_analysis.R`:

**In `output_data/clr_analysis/`:**
- `clr_transformed_long.csv` - CLR values (long format)
- `clr_transformed_wide.csv` - CLR values (wide format)
- `relative_abundance.csv` - Relative abundances before CLR
- `clr_summary.csv` - Summary statistics per protein-glycan
- `clr_statistical_tests.csv` - All statistical tests
- `clr_significant_results.csv` - Significant results only

---

## Interpreting Results

### CLR Values
- **Positive CLR**: Composition more abundant than average
- **Negative CLR**: Composition less abundant than average
- **CLR = 0**: Composition at geometric mean

### Statistical Significance
- **\*\*\***: p < 0.001 (highly significant)
- **\*\***: p < 0.01 (very significant)
- **\***: p < 0.05 (significant)

### Effect Sizes (Cohen's d)
- **Large**: |d| ≥ 0.8
- **Medium**: 0.5 ≤ |d| < 0.8
- **Small**: 0.2 ≤ |d| < 0.5
- **Negligible**: |d| < 0.2

---

## Common Workflows

### 1. Exploratory Analysis
```r
# Run complete analysis
source("CLR_visualisations.R")

# Review clr_comprehensive_summary.png
# Check clr_protein_summary.csv for proteins of interest
```

### 2. Focus on Specific Proteins
```r
# Filter results for specific proteins
protein_of_interest <- "PROTEIN_ID"

results_filtered <- clr_stats %>%
  filter(protein_accessions == protein_of_interest)

# Create focused visualization
dotplot_focused <- plot_clr_dotplot(
  results_filtered,
  top_n_proteins = NULL
)
```

### 3. Export for Publication
```r
# Create high-resolution figures
summary_plot <- plot_clr_summary(
  clr_stats,
  output_file = "Figure_CLR_Analysis.png",
  top_n_proteins = 30
)

# Export publication-ready table
publication_table <- clr_stats %>%
  filter(significant == TRUE) %>%
  select(protein_accessions, glycan_composition, 
         clr_difference, cohens_d, p_value_adj) %>%
  arrange(p_value_adj)

write.csv(publication_table, "Table_CLR_Significant.csv", row.names = FALSE)
```

---

## Troubleshooting

**Issue**: "No significant results to plot"
- **Solution**: Lower p_threshold or check if there are truly no differences

**Issue**: Plot is too crowded
- **Solution**: Use `top_n_proteins` parameter to show fewer proteins

**Issue**: Error about missing packages
- **Solution**: Install required packages:
  ```r
  install.packages(c("tidyverse", "compositions", "ggplot2", "patchwork"))
  ```

---

## Next Steps

1. **For detailed methodology**: See `CLR_ANALYSIS_README.md`
2. **For comparisons**: Compare with ILR results from `example_ilr_analysis.R`
3. **For further analysis**: Use the CLR values in downstream multivariate analyses

---

**Questions?** Refer to the full documentation in `CLR_ANALYSIS_README.md`


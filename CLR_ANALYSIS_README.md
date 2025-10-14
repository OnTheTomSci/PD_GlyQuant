# CLR Analysis Functions - User Guide

## Overview

Two new functions have been added to analyze glycan compositions using Centered Log-Ratio (CLR) transformation:

1. `calculate_multiple_proteins_glycan_clr()` - Performs CLR transformation on glycan compositions
2. `test_clr_protein_composition_differences()` - Statistically tests differences between groups

## Function 1: calculate_multiple_proteins_glycan_clr()

### Purpose
Calculates relative abundances and CLR transformation for glycan compositions across multiple proteins.

### Input
- `data_long`: Long format data with columns:
  - `abundance`: Glycopeptide abundance values
  - `sample`: Sample identifiers
  - `group`: Group labels (e.g., "Healthy", "MECFS")
  - `protein_accessions`: Protein identifiers
  - `glycan_composition`: Glycan composition strings

### Process
1. **Filters** out missing/invalid data
2. **Aggregates** abundances by sample, protein, and glycan composition
3. **Calculates relative abundances**: For each protein in each sample, converts glycan abundances to proportions (sum = 1)
4. **Applies CLR transformation**: `CLR = log(x_i) - mean(log(x))` using the `compositions` package
5. **Generates summaries**: Creates statistics for each protein-glycan combination

### Output
A list containing:
- `clr_transformed`: Wide format data frame (samples × glycan compositions)
- `clr_transformed_long`: Long format CLR values (easy for plotting/analysis)
- `relative_abundance`: Relative abundances before transformation
- `clr_summary`: Mean, SD, median, min, max CLR values per protein-glycan
- `protein_summary`: Number of compositions and samples per protein
- `transformation_info`: Metadata (counts of proteins, samples, compositions)
- `raw_matrix`: Original composition matrix
- `clr_matrix`: CLR-transformed matrix

### Example Usage
```r
source("functions/peptidegroups_preprocessing.R")
source("functions/peptidegroups_compositional_analysis.R")

# Load data
data <- load_and_preprocess_peptidegroups()

# Apply CLR transformation
clr_results <- calculate_multiple_proteins_glycan_clr(data$glyco_peptide_groups_long)

# Access results
clr_table <- clr_results$clr_transformed_long
summary_stats <- clr_results$clr_summary
```

---

## Function 2: test_clr_protein_composition_differences()

### Purpose
Statistically tests whether CLR-transformed glycan compositions differ between two groups.

### Input Parameters
- `clr_result`: Output from `calculate_multiple_proteins_glycan_clr()`
- `group1_name`: Name of first group (default: "Healthy")
- `group2_name`: Name of second group (default: "MECFS")
- `p_adjust_method`: Multiple testing correction method (default: "BH" for Benjamini-Hochberg)

### Statistical Tests Performed
For each protein-glycan combination:
1. **Two-sample t-test**: Tests if mean CLR values differ between groups
2. **Cohen's d effect size**: Standardized difference (small: 0.2-0.5, medium: 0.5-0.8, large: >0.8)
3. **95% Confidence intervals**: For the CLR difference
4. **Multiple testing correction**: Adjusts p-values using specified method

### Output
A data frame with one row per protein-glycan combination containing:
- `protein_accessions`: Protein identifier
- `glycan_composition`: Glycan composition
- `group1_name`, `group2_name`: Group identifiers
- `group1_mean_clr`, `group2_mean_clr`: Mean CLR values per group
- `clr_difference`: Group2 mean - Group1 mean
- `group1_sd`, `group2_sd`: Standard deviations
- `cohens_d`: Effect size
- `effect_size_interpretation`: "negligible", "small", "medium", or "large"
- `t_statistic`: t-test statistic
- `p_value`: Raw p-value
- `p_value_adj`: Adjusted p-value (for multiple testing)
- `significant`: TRUE if p_value_adj < 0.05
- `significance_level`: "***" (p<0.001), "**" (p<0.01), "*" (p<0.05), or "NS"
- `df`: Degrees of freedom
- `n_group1`, `n_group2`: Sample sizes
- `ci_lower`, `ci_upper`: 95% confidence interval bounds

### Example Usage
```r
# Test differences between groups
clr_stats <- test_clr_protein_composition_differences(
  clr_result = clr_results,
  group1_name = "Healthy",
  group2_name = "MECFS",
  p_adjust_method = "BH"
)

# Filter significant results
significant <- clr_stats %>%
  filter(significant == TRUE) %>%
  arrange(p_value_adj)

# View top results
head(significant)

# Filter by effect size
large_effects <- significant %>%
  filter(effect_size_interpretation == "large")
```

---

## Complete Workflow Example

See `example_clr_analysis.R` for a complete workflow that:
1. Loads and preprocesses data
2. Calculates CLR transformation
3. Performs statistical testing
4. Saves all results to CSV files
5. Prints summary statistics

To run:
```r
source("example_clr_analysis.R")
```

### Output Files
All results are saved to `output_data/clr_analysis/`:
- `clr_transformed_long.csv`: Long format CLR values
- `clr_transformed_wide.csv`: Wide format CLR values
- `relative_abundance.csv`: Relative abundances (before CLR)
- `clr_summary.csv`: Summary statistics
- `protein_summary.csv`: Per-protein summary
- `clr_statistical_tests.csv`: All test results
- `clr_significant_results.csv`: Significant results only

---

## Interpretation Guide

### CLR Values
- **CLR values** represent log-ratios relative to the geometric mean
- **Positive CLR**: Composition is more abundant than average
- **Negative CLR**: Composition is less abundant than average
- **CLR = 0**: Composition is at the geometric mean

### CLR Differences
- **Positive difference**: Group2 has higher relative abundance than Group1
- **Negative difference**: Group2 has lower relative abundance than Group1

### Effect Sizes (Cohen's d)
- **Negligible**: |d| < 0.2
- **Small**: 0.2 ≤ |d| < 0.5
- **Medium**: 0.5 ≤ |d| < 0.8
- **Large**: |d| ≥ 0.8

### Statistical Significance
- **\*\*\***: p < 0.001 (highly significant)
- **\*\***: p < 0.01 (very significant)
- **\***: p < 0.05 (significant)
- **NS**: Not significant

---

## Why Use CLR?

### Advantages
1. **Compositional data**: Properly handles data that sum to a constant (relative abundances)
2. **Removes closure**: Avoids spurious correlations from compositional constraint
3. **Interpretable**: Differences in CLR space represent fold-changes in ratios
4. **Multivariate**: Each glycan composition gets its own CLR value (unlike ILR which reduces dimensionality)

### When to Use CLR vs ILR
- **CLR**: When you want to test each glycan composition individually
- **ILR**: When you want to reduce dimensionality and analyze overall compositional shifts

Both approaches are valid for compositional data analysis!

---

## Visualization Functions

Four visualization functions are available to display CLR analysis results:

### 1. plot_clr_dotplot()

Creates a dot plot showing CLR differences across proteins and glycan compositions.

**Parameters:**
- `clr_stats_result`: Result from `test_clr_protein_composition_differences()`
- `p_threshold`: P-value threshold for significance (default: 0.05)
- `top_n_proteins`: Show only top N proteins (default: NULL shows all)

**Features:**
- Dot size = Magnitude of CLR difference
- Dot color = Direction (pink = increased, teal = decreased)
- Opacity = Significance (solid = significant, faded = not significant)
- Proteins ordered by number of significant glycans

**Example:**
```r
dotplot <- plot_clr_dotplot(clr_stats, p_threshold = 0.05, top_n_proteins = 30)
ggsave("clr_dotplot.png", dotplot, width = 14, height = 10, dpi = 300)
```

### 2. plot_clr_heatmap()

Creates a heatmap of CLR differences or p-values.

**Parameters:**
- `clr_stats_result`: Result from statistical testing
- `value_type`: "clr_difference", "effect", or "p_value" (default: "clr_difference")
- `p_threshold`: P-value threshold (default: 0.05)
- `top_n_proteins`: Show only top N proteins (default: NULL)

**Features:**
- Color scale shows CLR differences (blue = decreased, red = increased) or p-values
- Star (*) marks significant results
- Proteins ordered by number of significant changes

**Example:**
```r
# Heatmap colored by CLR difference
heatmap <- plot_clr_heatmap(clr_stats, value_type = "clr_difference", top_n_proteins = 30)

# Heatmap colored by p-value
heatmap_p <- plot_clr_heatmap(clr_stats, value_type = "p_value", top_n_proteins = 30)
```

### 3. plot_clr_forest()

Creates a forest plot showing CLR differences with 95% confidence intervals.

**Parameters:**
- `clr_stats_result`: Result from statistical testing
- `show_only_significant`: Show only significant results (default: TRUE)
- `p_threshold`: P-value threshold (default: 0.05)
- `top_n`: Show only top N results (default: 50)

**Features:**
- Dots = CLR difference (point estimate)
- Error bars = 95% confidence intervals
- Vertical line at 0 = no change
- Results sorted by effect size

**Example:**
```r
forest <- plot_clr_forest(clr_stats, show_only_significant = TRUE, top_n = 50)
ggsave("clr_forest.png", forest, width = 10, height = 14, dpi = 300)
```

### 4. plot_clr_summary()

Creates a comprehensive summary combining heatmap, forest plot, and dot plot.

**Parameters:**
- `clr_stats_result`: Result from statistical testing
- `output_file`: Optional filename to save (default: NULL)
- `p_threshold`: P-value threshold (default: 0.05)
- `top_n_proteins`: Show top N proteins in heatmap/dotplot (default: 30)

**Features:**
- Combined view using patchwork package
- Automatically saves if output_file is specified
- Shows all three visualization types in one figure

**Example:**
```r
summary_plot <- plot_clr_summary(
  clr_stats, 
  output_file = "clr_comprehensive_summary.png",
  p_threshold = 0.05,
  top_n_proteins = 30
)
```

### Complete Visualization Workflow

See `CLR_visualisations.R` for a complete workflow that:
1. Loads and preprocesses data
2. Calculates CLR transformation
3. Performs statistical testing
4. Creates all visualization types
5. Saves results tables
6. Prints summary statistics

To run the complete workflow:
```r
source("CLR_visualisations.R")
```

This will generate:
- **Dot plots**: `clr_dotplot_all.png`, `clr_dotplot_top30.png`
- **Heatmaps**: `clr_heatmap_difference.png`, `clr_heatmap_pvalue.png`
- **Forest plots**: `clr_forest_significant.png`, `clr_forest_all.png`
- **Summary**: `clr_comprehensive_summary.png`
- **Tables**: `clr_all_results_table.csv`, `clr_significant_results_table.csv`, `clr_protein_summary.csv`

---

## References

- Aitchison, J. (1982). The statistical analysis of compositional data. Journal of the Royal Statistical Society: Series B (Methodological), 44(2), 139-160.
- Gloor, G. B., et al. (2017). Microbiome datasets are compositional: and this is not optional. Frontiers in microbiology, 8, 2224.


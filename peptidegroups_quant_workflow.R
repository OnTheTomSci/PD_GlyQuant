# Peptide Groups Quantification Workflow
# Main orchestration script for modularized peptide groups analysis
# This script loads all function modules and executes the complete analysis pipeline

# Load required libraries
library(tidyverse)
library(here)
library(janitor)
library(stringr)
library(readr)
library(dplyr)
library(ggplot2)
library(broom)
library(viridis)
library(patchwork)
library(EnhancedVolcano)
library(effectsize)

# Set working directory to the project root
setwd("/Users/thomasreilly/Desktop/PD_GlyQuant")

# Source all function modules
source("functions/peptidegroups_preprocessing.R")
source("functions/peptidegroups_volcano.R")
source("functions/peptidegroups_glycan_analysis.R")
source("functions/peptidegroups_modification_stats.R")
source("functions/peptidegroups_protein_glycosylation.R")
source("functions/peptidegroups_glycosite_glycosylation.R")

cat("=== PEPTIDE GROUPS QUANTIFICATION WORKFLOW ===\n")
cat("Starting modularized analysis pipeline...\n\n")

# =============================================================================
# 1. DATA LOADING AND PREPROCESSING
# =============================================================================
cat("Step 1: Data Loading and Preprocessing\n")
cat("=====================================\n")

# Load and preprocess data
preprocessed_data <- load_and_preprocess_peptidegroups(
  study_info_path = "input_data/10S_MECFS_GPEPS_250125_StudyInformation.txt",
  peptide_groups_path = "input_data/10S_MECFS_GPEPS_250125_PeptideGroups.txt",
  glycan_class_map_path = "input_data/glycan_class_map.csv",
  fdr_threshold = 0.001
)

# Extract processed data
glyco_peptide_groups <- preprocessed_data$glyco_peptide_groups
glyco_peptide_groups_long <- preprocessed_data$glyco_peptide_groups_long
StudyInformation <- preprocessed_data$StudyInformation
glycan_class_map <- preprocessed_data$glycan_class_map
sample_metadata <- preprocessed_data$sample_metadata

# Add group classification if not already present
if (!"group" %in% colnames(glyco_peptide_groups_long)) {
  glyco_peptide_groups_long <- glyco_peptide_groups_long %>%
    mutate(group = case_when(
      str_starts(sample, "hc") ~ "Healthy",
      str_starts(sample, "m") ~ "MECFS",
      TRUE ~ "Unknown"
    ))
}

cat("✓ Data loading and preprocessing completed\n\n")

# =============================================================================
# 2. VOLCANO PLOT ANALYSES
# =============================================================================
cat("Step 2: Volcano Plot Analyses\n")
cat("=============================\n")

# Create volcano plots at different levels
gene_volcano_results <- create_gene_volcano(
  data = glyco_peptide_groups_long,
  output_path = "figures/peptidegroups_intensity/gene_volcano_plot.png",
  group1_name = "Healthy",
  group2_name = "MECFS",
  p_cutoff = 0.05,
  fc_cutoff = 1
)

protein_glycan_volcano_results <- create_protein_glycan_volcano(
  data = glyco_peptide_groups_long,
  output_path = "figures/peptidegroups_intensity/protein_glycan_volcano_plot.png",
  group1_name = "Healthy",
  group2_name = "MECFS",
  p_cutoff = 0.05,
  fc_cutoff = 1
)

protein_glycosite_glycan_volcano_results <- create_protein_glycosite_glycan_volcano(
  data = glyco_peptide_groups_long,
  output_path = "figures/peptidegroups_intensity/protein_glycosite_glycan_volcano_plot.png",
  group1_name = "Healthy",
  group2_name = "MECFS",
  p_cutoff = 0.05,
  fc_cutoff = 1
)

cat("✓ Volcano plot analyses completed\n\n")

# =============================================================================
# 3. GLYCAN COMPOSITION ANALYSIS
# =============================================================================
cat("Step 3: Glycan Composition Analysis\n")
cat("===================================\n")

# Analyze glycan compositions
glycan_analysis_results <- analyze_glycan_compositions(
  data = glyco_peptide_groups_long,
  output_dir = "output_data/peptidegroups_intensity",
  figures_dir = "figures"
)

# Analyze pseudo-glycomics
pseudo_glycomics_results <- analyze_pseudo_glycomics(
  data = glyco_peptide_groups_long,
  output_dir = "output_data/peptidegroups_intensity"
)

# Create additional glycan analysis plots
glycan_plots <- create_glycan_analysis_plots(
  data = glyco_peptide_groups_long,
  figures_dir = "figures"
)

cat("✓ Glycan composition analysis completed\n\n")

# =============================================================================
# 4. MODIFICATION STATISTICS (SAMPLE-LEVEL)
# =============================================================================
cat("Step 4: Sample-Level Modification Statistics\n")
cat("===========================================\n")

# Analyze sample-level fucosylation
fucosylation_results <- analyze_sample_fucosylation(
  data = glyco_peptide_groups_long,
  output_dir = "output_data/peptidegroups_intensity",
  figures_dir = "figures/peptidegroups_intensity"
)

# Analyze sample-level sialylation
sialylation_results <- analyze_sample_sialylation(
  data = glyco_peptide_groups_long,
  output_dir = "output_data/peptidegroups_intensity",
  figures_dir = "figures/peptidegroups_intensity"
)

# Analyze glycan class by sample
glycan_class_results <- analyze_glycan_class_by_sample(
  data = glyco_peptide_groups_long,
  output_dir = "output_data/peptidegroups_intensity",
  figures_dir = "figures"
)

cat("✓ Sample-level modification statistics completed\n\n")

# =============================================================================
# 5. PROTEIN-LEVEL GLYCOSYLATION ANALYSIS
# =============================================================================
cat("Step 5: Protein-Level Glycosylation Analysis\n")
cat("===========================================\n")

# Analyze protein-level glycosylation
protein_glycosylation_results <- analyze_protein_glycosylation(
  data = glyco_peptide_groups_long,
  output_dir = "output_data/peptidegroups_intensity/protein_level",
  figures_dir = "figures/protein_level"
)

# Run ILR-based analysis
protein_glycosylation_ilr_results <- analyze_protein_glycosylation_ilr(
  data = glyco_peptide_groups_long,
  original_results = protein_glycosylation_results,
  output_dir = "output_data/peptidegroups_intensity/protein_level_ilr",
  figures_dir = "figures/protein_level_ilr"
)

cat("✓ Protein-level glycosylation analysis completed\n\n")

# =============================================================================
# 6. GLYCOSITE-LEVEL GLYCOSYLATION ANALYSIS
# =============================================================================
cat("Step 6: Glycosite-Level Glycosylation Analysis\n")
cat("=============================================\n")

# Analyze glycosite-level glycosylation
glycosite_glycosylation_results <- analyze_glycosite_glycosylation(
  data = glyco_peptide_groups_long,
  output_dir = "output_data/peptidegroups_intensity/glycosite_level",
  figures_dir = "figures/glycosite_level"
)

cat("✓ Glycosite-level glycosylation analysis completed\n\n")

# =============================================================================
# 7. GENE-LEVEL DIFFERENTIAL ABUNDANCE ANALYSIS
# =============================================================================
cat("Step 7: Gene-Level Differential Abundance Analysis\n")
cat("=================================================\n")

# Summarize glycopeptide data to gene level
gene_level_data <- glyco_peptide_groups_long %>%
  group_by(gene_name, sample) %>%
  summarise(
    total_abundance = sum(abundance, na.rm = TRUE),
    .groups = 'drop'
  )

# Map samples to groups
gene_stats <- gene_level_data %>%
  mutate(
    group = case_when(
      str_detect(sample, "^hc") ~ "Healthy",
      str_detect(sample, "^m")  ~ "MECFS",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(gene_name, group) %>%
  summarise(
    mean_abundance = mean(total_abundance, na.rm = TRUE),
    sd_abundance   = sd(total_abundance, na.rm = TRUE),
    n              = n(),
    .groups = 'drop'
  ) %>%
  pivot_wider(
    names_from = group,
    values_from = c(mean_abundance, sd_abundance, n),
    names_glue = "{.value}_{group}"
  )

# Differential statistics
gene_stats <- gene_stats %>%
  mutate(
    log2FC = log2(mean_abundance_MECFS / mean_abundance_Healthy),
    se_diff = sqrt((sd_abundance_Healthy^2 / n_Healthy) +
                   (sd_abundance_MECFS^2 / n_MECFS)),
    df = n_Healthy + n_MECFS - 2,
    t_stat = (mean_abundance_MECFS - mean_abundance_Healthy) / se_diff,
    p_value = 2 * pt(-abs(t_stat), df),
    adj_p_value = p.adjust(p_value, method = "BH")
  )

# Create gene-level volcano plot
gene_level_volcano <- EnhancedVolcano(
  gene_stats,
  lab = gene_stats$gene_name,
  x = "log2FC",
  y = "adj_p_value",
  xlab = "log2(Fold Change) MECFS / Healthy",
  ylab = "-log10(Adjusted P-value)",
  pCutoff = 0.08,
  FCcutoff = 0.5,      # |log2FC| > 0.5
  pointSize = 2.5,
  labSize = 3.5,
  colAlpha = 0.6,
  title = "Gene-Level Glycopeptide Changes",
  subtitle = "MECFS vs Healthy",
  caption = paste("Total genes:", nrow(gene_stats)),
  legendLabels = c("NS", "Log2FC", "Adj.P", "Adj.P & Log2FC"),
  legendPosition = "top",
  drawConnectors = TRUE,
  widthConnectors = 0.5
)

# Save gene-level volcano plot
ggsave("figures/gene_level_volcano_plot_enhanced.png", 
       gene_level_volcano, width = 10, height = 8, dpi = 300)

# Summary
cat("\nGene-level differential abundance analysis:\n")
cat("Total genes analyzed:", nrow(gene_stats), "\n")
cat("Significantly different genes:", sum(gene_stats$adj_p_value < 0.05, na.rm = TRUE), "\n")
cat("Upregulated in MECFS:", sum(gene_stats$adj_p_value < 0.05 & gene_stats$log2FC > 0, na.rm = TRUE), "\n")
cat("Downregulated in MECFS:", sum(gene_stats$adj_p_value < 0.05 & gene_stats$log2FC < 0, na.rm = TRUE), "\n")

cat("✓ Gene-level differential abundance analysis completed\n\n")

# =============================================================================
# 8. SUMMARY AND OUTPUT
# =============================================================================
cat("Step 8: Analysis Summary\n")
cat("=======================\n")

# Create comprehensive summary
analysis_summary <- list(
  preprocessing = list(
    total_glycopeptides = sample_metadata$total_samples,
    unique_genes = length(sample_metadata$common_genes),
    unique_proteins = length(sample_metadata$common_proteins),
    unique_glycosites = length(sample_metadata$common_gsites)
  ),
  volcano_analyses = list(
    gene_level = gene_volcano_results,
    protein_glycan_level = protein_glycan_volcano_results,
    protein_glycosite_glycan_level = protein_glycosite_glycan_volcano_results
  ),
  glycan_analysis = glycan_analysis_results,
  pseudo_glycomics = pseudo_glycomics_results,
  sample_level_modifications = list(
    fucosylation = fucosylation_results,
    sialylation = sialylation_results,
    glycan_class = glycan_class_results
  ),
  protein_level_analysis = protein_glycosylation_results,
  protein_level_ilr_analysis = protein_glycosylation_ilr_results,
  glycosite_level_analysis = glycosite_glycosylation_results,
  gene_level_analysis = list(
    statistics = gene_stats,
    volcano_plot = gene_level_volcano
  )
)

# Save comprehensive summary
saveRDS(analysis_summary, "output_data/peptidegroups_intensity/comprehensive_analysis_summary.rds")

cat("✓ Comprehensive analysis summary saved\n")

# Print final summary
cat("\n=== FINAL ANALYSIS SUMMARY ===\n")
cat("All analyses completed successfully!\n")
cat("Results saved to:\n")
cat("  - output_data/peptidegroups_intensity/\n")
cat("  - figures/\n")
cat("  - input_data/Corr/ (for correlation analyses)\n")
cat("\nKey outputs:\n")
cat("  - Volcano plots at gene, protein-glycan, and protein-glycosite-glycan levels\n")
cat("  - Glycan composition and pseudo-glycomics analysis\n")
cat("  - Sample-level fucosylation and sialylation statistics\n")
cat("  - Protein-level glycosylation analysis with ILR transformation\n")
cat("  - Glycosite-level glycosylation analysis\n")
cat("  - Gene-level differential abundance analysis\n")
cat("\nAnalysis pipeline completed successfully!\n")
cat("=========================================\n")
# Main analysis script
library(tidyverse)
library(patchwork)

# Source all function files
source("functions/data_ingest.R")
source("functions/count_stats.R")
source("functions/intensity_stats.R")
source("functions/ebayes.R")
source("functions/clustering_models.R")
source("functions/plot_funtions.R")
source("functions/pesudo_glycomics.R")

# Main analysis workflow
main_analysis <- function(glycoPSMs) {
  # Get unique proteins
  unique_proteins <- unique(glycoPSMs$protein_accessions)
  
  # Analyze all proteins
  protein_analyses <- list()
  for (protein in unique_proteins) {
    result <- analyze_protein(protein, glycoPSMs)
    if (!is.null(result)) {
      protein_analyses[[protein]] <- result
    }
  }
  
  # Process significant results
  significant_results <- process_significant_results(protein_analyses)
  
  # Save combined plots
  if (!is.null(significant_results$plots)) {
    save_combined_plots(significant_results$plots, 
                       "significant_proteins_fucosylation_boxplots")
  }
  
  # Save statistics
  if (nrow(significant_results$stats) > 0) {
    write.csv(significant_results$stats,
              "output_data/significant_protein_fucosylation_stats.csv",
              row.names = FALSE)
  }
  
  # Print diagnostic information
  print_diagnostic_info(protein_analyses)
  
  return(protein_analyses)
}

# Run the analysis
results <- main_analysis(glycoPSMs)

# Print summary of all analyses
cat("\n=== GLYCAN ANALYSIS SUMMARY ===\n")
cat("All glycan analysis functions have been executed.\n")
cat("Results have been saved to the output_data directory.\n")
cat("Check the following files for detailed results:\n")
cat("- PSM counts by composition and features\n")
cat("- Relative abundance calculations\n")
cat("- Chi-square goodness of fit tests\n")
cat("- Mann-Whitney U tests\n")
cat("- Ranked glycan compositions and features\n")
cat("===============================\n")
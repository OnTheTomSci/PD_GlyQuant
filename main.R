# Main analysis script
library(tidyverse)
library(patchwork)

# Source all function files
source("R/functions/glycan_analysis_functions.R")
source("R/functions/plot_functions.R")
source("R/functions/statistics_functions.R")

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
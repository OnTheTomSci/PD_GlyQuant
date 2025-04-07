library(dplyr)

#' Count PSMs by glycan composition grouped by sample
#' @param data The glycoPSMs dataframe
#' @return A dataframe with PSM counts by glycan composition and sample
count_psms_by_composition <- function(data) {
  psm_counts <- data %>%
    group_by(sample, disease_status, glycan_composition) %>%
    summarise(
      psm_count = n(),
      .groups = 'drop'
    ) %>%
    arrange(sample, desc(psm_count))
  
  return(psm_counts)
}

# Calculate PSM counts
psm_counts_by_composition <- count_psms_by_composition(glycoPSMs)

# Save results
write.csv(psm_counts_by_composition, 
          file = "output_data/psm_counts_by_composition.csv", 
          row.names = FALSE)

# Print summary
cat("\nPSM Count Summary by Glycan Composition:\n")
cat("Total samples:", length(unique(psm_counts_by_composition$sample)), "\n")
cat("Total glycan compositions:", length(unique(psm_counts_by_composition$glycan_composition)), "\n")
cat("Total PSMs:", sum(psm_counts_by_composition$psm_count), "\n")

#' Count PSMs by a specified glycan feature grouped by sample
#' @param data The glycoPSMs dataframe
#' @param feature The column name to group by (e.g., "contains_Fuc", "contains_NeuAc", "glycan_class", "sia_count", "fuc_count")
#' @return A dataframe with PSM counts by the specified feature and sample
count_psms_by_feature <- function(data, feature) {
  # Validate that the feature exists in the data
  if (!feature %in% names(data)) {
    stop(paste("Feature", feature, "not found in the data"))
  }
  
  # Create the grouping expression dynamically
  group_expr <- rlang::sym(feature)
  
  # Count PSMs by the specified feature and sample
  psm_counts <- data %>%
    group_by(sample, disease_status, !!group_expr) %>%
    summarise(
      psm_count = n(),
      .groups = 'drop'
    ) %>%
    arrange(sample, desc(psm_count))
  
  return(psm_counts)
}

# Example usage:
# Count PSMs by fucose presence
fuc_psm_counts <- count_psms_by_feature(glycoPSMs, "contains_Fuc")
write.csv(fuc_psm_counts, 
          file = "output_data/psm_counts_by_fucose.csv", 
          row.names = FALSE)

# Count PSMs by sialic acid presence
sia_psm_counts <- count_psms_by_feature(glycoPSMs, "contains_NeuAc")
write.csv(sia_psm_counts, 
          file = "output_data/psm_counts_by_sialic_acid.csv", 
          row.names = FALSE)

# Count PSMs by glycan class
class_psm_counts <- count_psms_by_feature(glycoPSMs, "glycan_class")
write.csv(class_psm_counts, 
          file = "output_data/psm_counts_by_glycan_class.csv", 
          row.names = FALSE)

# Count PSMs by sialic acid count
sia_count_psm_counts <- count_psms_by_feature(glycoPSMs, "sia_count")
write.csv(sia_count_psm_counts, 
          file = "output_data/psm_counts_by_sia_count.csv", 
          row.names = FALSE)

# Count PSMs by fucose count
fuc_count_psm_counts <- count_psms_by_feature(glycoPSMs, "fuc_count")
write.csv(fuc_count_psm_counts, 
          file = "output_data/psm_counts_by_fuc_count.csv", 
          row.names = FALSE)

# Print summary for each feature
cat("\nPSM Count Summary by Fucose Presence:\n")
cat("Total samples:", length(unique(fuc_psm_counts$sample)), "\n")
cat("Total PSMs:", sum(fuc_psm_counts$psm_count), "\n")

cat("\nPSM Count Summary by Sialic Acid Presence:\n")
cat("Total samples:", length(unique(sia_psm_counts$sample)), "\n")
cat("Total PSMs:", sum(sia_psm_counts$psm_count), "\n")

cat("\nPSM Count Summary by Glycan Class:\n")
cat("Total samples:", length(unique(class_psm_counts$sample)), "\n")
cat("Total glycan classes:", length(unique(class_psm_counts$glycan_class)), "\n")
cat("Total PSMs:", sum(class_psm_counts$psm_count), "\n")

cat("\nPSM Count Summary by Sialic Acid Count:\n")
cat("Total samples:", length(unique(sia_count_psm_counts$sample)), "\n")
cat("Total PSMs:", sum(sia_count_psm_counts$psm_count), "\n")

cat("\nPSM Count Summary by Fucose Count:\n")
cat("Total samples:", length(unique(fuc_count_psm_counts$sample)), "\n")
cat("Total PSMs:", sum(fuc_count_psm_counts$psm_count), "\n")

#' Calculate relative abundance percentages from PSM count data frames
#' @param psm_counts_df The PSM counts dataframe (output from count_psms_by_feature or count_psms_by_composition)
#' @param feature_col The column name containing the feature categories (e.g., "contains_Fuc", "glycan_class")
#' @param output_file The output CSV file path
#' @return A dataframe with relative abundance percentages
calculate_relative_abundance <- function(psm_counts_df, feature_col, output_file) {
  # Validate that the feature column exists in the data
  if (!feature_col %in% names(psm_counts_df)) {
    stop(paste("Feature column", feature_col, "not found in the data"))
  }
  
  # Calculate total PSMs per sample
  sample_totals <- psm_counts_df %>%
    group_by(sample, disease_status) %>%
    summarise(
      total_psms = sum(psm_count),
      .groups = 'drop'
    )
  
  # Calculate relative abundance percentages
  relative_abundance <- psm_counts_df %>%
    # Join with sample totals
    left_join(sample_totals, by = c("sample", "disease_status")) %>%
    # Calculate relative percentage
    mutate(
      relative_percentage = (psm_count / total_psms) * 100
    ) %>%
    # Select and arrange columns
    select(sample, disease_status, !!rlang::sym(feature_col), psm_count, total_psms, relative_percentage) %>%
    arrange(sample, desc(relative_percentage))
  
  # Save results
  write.csv(relative_abundance, 
            file = output_file, 
            row.names = FALSE)
  
  # Print summary
  cat("\nRelative Abundance Summary for", feature_col, ":\n")
  cat("Total samples:", length(unique(relative_abundance$sample)), "\n")
  cat("Total", feature_col, "categories:", length(unique(relative_abundance[[feature_col]])), "\n")
  
  return(relative_abundance)
}

# Calculate relative abundances for each feature
# For glycan composition
composition_relative_abundance <- calculate_relative_abundance(
  psm_counts_by_composition, 
  "glycan_composition", 
  "output_data/relative_abundance_by_composition.csv"
)

# For fucose presence
fuc_relative_abundance <- calculate_relative_abundance(
  fuc_psm_counts, 
  "contains_Fuc", 
  "output_data/relative_abundance_by_fucose.csv"
)

# For sialic acid presence
sia_relative_abundance <- calculate_relative_abundance(
  sia_psm_counts, 
  "contains_NeuAc", 
  "output_data/relative_abundance_by_sialic_acid.csv"
)

# For glycan class
class_relative_abundance <- calculate_relative_abundance(
  class_psm_counts, 
  "glycan_class", 
  "output_data/relative_abundance_by_glycan_class.csv"
)

# For sialic acid count
sia_count_relative_abundance <- calculate_relative_abundance(
  sia_count_psm_counts, 
  "sia_count", 
  "output_data/relative_abundance_by_sia_count.csv"
)

# For fucose count
fuc_count_relative_abundance <- calculate_relative_abundance(
  fuc_count_psm_counts, 
  "fuc_count", 
  "output_data/relative_abundance_by_fuc_count.csv"
)

# Print overall summary
cat("\nOverall Relative Abundance Summary:\n")
cat("All relative abundance calculations completed and saved to output_data directory.\n")


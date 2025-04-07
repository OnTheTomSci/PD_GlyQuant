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

# check for normality of the data 









#' Perform chi-square goodness of fit test on glycan feature PSM count data
#' @param psm_counts_df The PSM counts dataframe (output from count_psms_by_feature or count_psms_by_composition)
#' @param feature_col The column name containing the feature categories (e.g., "contains_Fuc", "glycan_class")
#' @param output_file The output CSV file path for the test results
#' @param expected_distribution Optional named vector of expected proportions (must sum to 1)
#' @return A dataframe with chi-square test results
perform_chi_square_test <- function(psm_counts_df, feature_col, output_file, expected_distribution = NULL) {
  # Validate that the feature column exists in the data
  if (!feature_col %in% names(psm_counts_df)) {
    stop(paste("Feature column", feature_col, "not found in the data"))
  }
  
  # Get unique samples
  unique_samples <- unique(psm_counts_df$sample)
  
  # Initialize results dataframe
  results <- data.frame(
    sample = character(),
    disease_status = character(),
    chi_square_statistic = numeric(),
    p_value = numeric(),
    degrees_of_freedom = numeric(),
    significant = logical(),
    stringsAsFactors = FALSE
  )
  
  # For each sample, perform chi-square test
  for (sample_name in unique_samples) {
    # Filter data for this sample
    sample_data <- psm_counts_df %>% filter(sample == sample_name)
    
    # Get disease status for this sample
    disease_status <- unique(sample_data$disease_status)
    
    # Extract observed counts and feature categories
    observed_counts <- sample_data$psm_count
    feature_categories <- sample_data[[feature_col]]
    
    # If no expected distribution is provided, use equal distribution
    if (is.null(expected_distribution)) {
      # Create equal expected distribution
      expected_distribution <- rep(1/length(observed_counts), length(observed_counts))
      names(expected_distribution) <- as.character(feature_categories)
    } else {
      # Validate that expected distribution has the right categories
      missing_categories <- setdiff(as.character(feature_categories), names(expected_distribution))
      if (length(missing_categories) > 0) {
        warning(paste("Expected distribution missing categories:", 
                     paste(missing_categories, collapse = ", "), 
                     "for sample", sample_name))
        # Add missing categories with zero probability
        for (cat in missing_categories) {
          expected_distribution[cat] <- 0
        }
      }
      
      # Ensure expected distribution sums to 1
      if (abs(sum(expected_distribution) - 1) > 1e-10) {
        warning("Expected distribution does not sum to 1, normalizing")
        expected_distribution <- expected_distribution / sum(expected_distribution)
      }
    }
    
    # Calculate expected counts
    total_count <- sum(observed_counts)
    expected_counts <- expected_distribution[as.character(feature_categories)] * total_count
    
    # Perform chi-square test
    chi_square_result <- chisq.test(observed_counts, p = expected_distribution[as.character(feature_categories)])
    
    # Add results to dataframe
    results <- rbind(results, data.frame(
      sample = sample_name,
      disease_status = disease_status,
      chi_square_statistic = chi_square_result$statistic,
      p_value = chi_square_result$p.value,
      degrees_of_freedom = chi_square_result$parameter,
      significant = chi_square_result$p.value < 0.05,
      stringsAsFactors = FALSE
    ))
  }
  
  # Save results
  write.csv(results, 
            file = output_file, 
            row.names = FALSE)
  
  # Print summary
  cat("\nChi-square Goodness of Fit Test Summary for", feature_col, ":\n")
  cat("Total samples tested:", nrow(results), "\n")
  cat("Samples with significant deviation (p < 0.05):", sum(results$significant), "\n")
  cat("Percentage of samples with significant deviation:", 
      round(sum(results$significant) / nrow(results) * 100, 2), "%\n")
  
  return(results)
}

# Perform chi-square tests for each feature
# For glycan composition (using equal expected distribution)
composition_chi_square <- perform_chi_square_test(
  psm_counts_by_composition, 
  "glycan_composition", 
  "output_data/chi_square_by_composition.csv"
)

# For fucose presence (using expected distribution based on literature or previous studies)
# Example: 60% with fucose, 40% without fucose
fuc_expected <- c("TRUE" = 0.6, "FALSE" = 0.4)
fuc_chi_square <- perform_chi_square_test(
  fuc_psm_counts, 
  "contains_Fuc", 
  "output_data/chi_square_by_fucose.csv",
  fuc_expected
)

# For sialic acid presence (using expected distribution)
# Example: 70% with sialic acid, 30% without sialic acid
sia_expected <- c("TRUE" = 0.7, "FALSE" = 0.3)
sia_chi_square <- perform_chi_square_test(
  sia_psm_counts, 
  "contains_NeuAc", 
  "output_data/chi_square_by_sialic_acid.csv",
  sia_expected
)

# For glycan class (using equal expected distribution)
class_chi_square <- perform_chi_square_test(
  class_psm_counts, 
  "glycan_class", 
  "output_data/chi_square_by_glycan_class.csv"
)

# For sialic acid count (using expected distribution based on literature)
# Example: 40% with 0 sialic acids, 30% with 1, 20% with 2, 10% with 3+
sia_count_expected <- c("Sialic Acid 0" = 0.4, "Sialic Acid 1" = 0.3, 
                        "Sialic Acid 2" = 0.2, "Sialic Acid 3+" = 0.1)
sia_count_chi_square <- perform_chi_square_test(
  sia_count_psm_counts, 
  "sia_count", 
  "output_data/chi_square_by_sia_count.csv",
  sia_count_expected
)

# For fucose count (using expected distribution based on literature)
# Example: 50% with 0 fucose, 30% with 1, 15% with 2, 5% with 3+
fuc_count_expected <- c("Fucose 0" = 0.5, "Fucose 1" = 0.3, 
                        "Fucose 2" = 0.15, "Fucose 3+" = 0.05)
fuc_count_chi_square <- perform_chi_square_test(
  fuc_count_psm_counts, 
  "fuc_count", 
  "output_data/chi_square_by_fuc_count.csv",
  fuc_count_expected
)

# Print overall summary
cat("\nOverall Chi-square Test Summary:\n")
cat("All chi-square goodness of fit tests completed and saved to output_data directory.\n")

#########################################################
# create ordianl dataframes of glycans and glyco feaures 

#' Rank glycan compositions from highest to lowest for each sample
#' @param psm_counts_df The PSM counts dataframe (output from count_psms_by_composition)
#' @param output_file The output CSV file path for the ranked data
#' @return A dataframe with glycan compositions ranked by PSM count for each sample
rank_glycan_compositions <- function(psm_counts_df, output_file) {
  # Validate that the required columns exist
  required_cols <- c("sample", "disease_status", "glycan_composition", "psm_count")
  missing_cols <- setdiff(required_cols, names(psm_counts_df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Rank glycan compositions by PSM count for each sample
  ranked_compositions <- psm_counts_df %>%
    # Group by sample to rank within each sample
    group_by(sample, disease_status) %>%
    # Add rank column (1 = highest PSM count)
    mutate(
      rank = rank(-psm_count, ties.method = "min")
    ) %>%
    # Arrange by sample and rank
    arrange(sample, rank) %>%
    # Ungroup to prepare for output
    ungroup()
  
  # Save results
  write.csv(ranked_compositions, 
            file = output_file, 
            row.names = FALSE)
  
  # Print summary
  cat("\nGlycan Composition Ranking Summary:\n")
  cat("Total samples:", length(unique(ranked_compositions$sample)), "\n")
  cat("Total glycan compositions:", length(unique(ranked_compositions$glycan_composition)), "\n")
  
  # Print top 5 compositions for each sample
  cat("\nTop 5 Glycan Compositions by Sample:\n")
  for (sample_name in unique(ranked_compositions$sample)) {
    cat("\nSample:", sample_name, "\n")
    top_5 <- ranked_compositions %>% 
      filter(sample == sample_name) %>% 
      arrange(rank) %>% 
      head(5)
    
    for (i in 1:nrow(top_5)) {
      cat(i, ". ", top_5$glycan_composition[i], 
          " (PSM count: ", top_5$psm_count[i], 
          ", Rank: ", top_5$rank[i], ")\n", sep = "")
    }
  }
  
  return(ranked_compositions)
}

# Rank glycan compositions
ranked_compositions <- rank_glycan_compositions(
  psm_counts_by_composition, 
  "output_data/ranked_glycan_compositions.csv"
)

# Print overall summary
cat("\nOverall Glycan Composition Ranking Summary:\n")
cat("Glycan compositions ranked by PSM count for each sample and saved to output_data/ranked_glycan_compositions.csv\n")

#' Rank glycan features across samples based on relative abundance
#' @param relative_abundance_df The relative abundance dataframe (output from calculate_relative_abundance)
#' @param feature_col The column name containing the feature categories (e.g., "contains_Fuc", "glycan_class")
#' @param output_file The output CSV file path for the ranked data
#' @return A dataframe with glycan features ranked across samples
rank_features_across_samples <- function(relative_abundance_df, feature_col, output_file) {
  # Validate that the required columns exist
  required_cols <- c("sample", "disease_status", feature_col, "relative_percentage")
  missing_cols <- setdiff(required_cols, names(relative_abundance_df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Get unique feature categories
  unique_features <- unique(relative_abundance_df[[feature_col]])
  
  # Initialize results dataframe
  results <- data.frame(
    sample = character(),
    disease_status = character(),
    feature = character(),
    relative_percentage = numeric(),
    rank = numeric(),
    stringsAsFactors = FALSE
  )
  
  # For each feature category, rank samples by relative percentage
  for (feature_value in unique_features) {
    # Filter data for this feature
    feature_data <- relative_abundance_df %>% 
      filter(!!rlang::sym(feature_col) == feature_value)
    
    # Add rank column (1 = highest relative percentage)
    feature_data <- feature_data %>%
      mutate(
        rank = rank(-relative_percentage, ties.method = "min")
      )
    
    # Add to results
    results <- rbind(results, data.frame(
      sample = feature_data$sample,
      disease_status = feature_data$disease_status,
      feature = feature_value,
      relative_percentage = feature_data$relative_percentage,
      rank = feature_data$rank,
      stringsAsFactors = FALSE
    ))
  }
  
  # Arrange by feature and rank
  results <- results %>%
    arrange(feature, rank)
  
  # Save results
  write.csv(results, 
            file = output_file, 
            row.names = FALSE)
  
  # Print summary
  cat("\nFeature Ranking Across Samples Summary for", feature_col, ":\n")
  cat("Total samples:", length(unique(results$sample)), "\n")
  cat("Total features:", length(unique(results$feature)), "\n")
  
  # Print top 3 samples for each feature
  cat("\nTop 3 Samples for Each Feature:\n")
  for (feature_value in unique_features) {
    cat("\nFeature:", feature_value, "\n")
    top_3 <- results %>% 
      filter(feature == feature_value) %>% 
      arrange(rank) %>% 
      head(3)
    
    for (i in 1:nrow(top_3)) {
      cat(i, ". Sample: ", top_3$sample[i], 
          " (Relative %: ", round(top_3$relative_percentage[i], 2), 
          "%, Rank: ", top_3$rank[i], ")\n", sep = "")
    }
  }
  
  return(results)
}

# Rank features across samples for each glycan feature
# For glycan composition
composition_ranked_across_samples <- rank_features_across_samples(
  composition_relative_abundance, 
  "glycan_composition", 
  "output_data/ranked_composition_across_samples.csv"
)

# For fucose presence
fuc_ranked_across_samples <- rank_features_across_samples(
  fuc_relative_abundance, 
  "contains_Fuc", 
  "output_data/ranked_fucose_across_samples.csv"
)

# For sialic acid presence
sia_ranked_across_samples <- rank_features_across_samples(
  sia_relative_abundance, 
  "contains_NeuAc", 
  "output_data/ranked_sialic_acid_across_samples.csv"
)

# For glycan class
class_ranked_across_samples <- rank_features_across_samples(
  class_relative_abundance, 
  "glycan_class", 
  "output_data/ranked_glycan_class_across_samples.csv"
)

# For sialic acid count
sia_count_ranked_across_samples <- rank_features_across_samples(
  sia_count_relative_abundance, 
  "sia_count", 
  "output_data/ranked_sia_count_across_samples.csv"
)

# For fucose count
fuc_count_ranked_across_samples <- rank_features_across_samples(
  fuc_count_relative_abundance, 
  "fuc_count", 
  "output_data/ranked_fuc_count_across_samples.csv"
)

# Print overall summary
cat("\nOverall Feature Ranking Across Samples Summary:\n")
cat("All feature rankings across samples completed and saved to output_data directory.\n")


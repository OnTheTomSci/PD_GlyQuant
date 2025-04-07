library(dplyr)
#' @param gpeps_dataframe the input glycopsm dataframe with all the glyco features
#' anotated to the data frame
#' @param top_lev_group the top level grouping is t how you subset and grup data to
#' peform glycan feature anlysse at eg. protein or glycosites
#' @param value_col Column name for the measurement values to analyze
#' @param sample_col Column name to use for pivoting to wide format
#' @param group_values Vector of expected group values (default: c("Healthy", "MECFS"))
#' @param min_samples Minimum number of samples required per group (default: 3)
#' @param file_prefix Prefix for output files (default: "Analysis_")

#' @return a matrix of relative aubundances for each top level groupingings and for
#' each sample
#'
#' @export
glyco_matrix <- function(
  gpeps_dataframe,
  top_lev_group,
  glycofeature_group,
  value_col,
  sample_col,
  group_col,
  min_samples = 3,
  file_prefix = "Analysis_") {
  # Convert column names to symbols for dplyr operations
  top_lev_sym <- rlang::sym(top_lev_group)
  glycofeature_sym <- rlang::sym(glycofeature_group)
  value_sym <- rlang::sym(value_col)
  sample_sym <- rlang::sym(sample_col)
  group_sym <- rlang::sym(group_col)

  # Get unique top level grouping values
  unique_top_groups <- unique(gpeps_dataframe[[top_lev_group]])

  # Get unique samples across all data
  unique_samples <- unique(gpeps_dataframe[[sample_col]])

  # Initialize result dataframe to store all results
  result_df <- data.frame()

  # Loop through each unique top level group
  for (current_group in unique_top_groups) {
    # Subset data for the current top level group
    group_data <- gpeps_dataframe %>%
      dplyr::filter(!!top_lev_sym == current_group)

    # Skip if insufficient data
    if (nrow(group_data) == 0) {
      next
    }

    # Get unique glycofeature values for this group
    unique_glycofeatures <- unique(group_data[[glycofeature_group]])

    # Create temporary dataframe to store results for this group
    temp_df <- data.frame(
      top_level_group = rep(current_group, length(unique_glycofeatures)),
      glycofeature = unique_glycofeatures,
      stringsAsFactors = FALSE
    )

    # Add columns for each sample, initialized to NA
    for (sample_id in unique_samples) {
      temp_df[[sample_id]] <- NA
    }

    # Loop through each sample
    for (sample_id in unique_samples) {
      # Subset data for current sample
      sample_data <- group_data %>%
        dplyr::filter(!!sample_sym == sample_id)

      # If no data exists for this sample and top level group, leave as NA
      if (nrow(sample_data) == 0) {
        # Values already initialized to NA, so skip to next sample
        next
      }

      # Calculate total sum for this sample and top level group
      total_sum <- sum(sample_data[[value_col]], na.rm = TRUE)

      # If total sum is zero or NA, leave values as NA and skip to next sample
      if (is.na(total_sum) || total_sum == 0) {
        next
      }

      # Calculate relative abundance for each glycofeature
      for (i in 1:nrow(temp_df)) {
        glyco_feature <- temp_df$glycofeature[i]

        # Sum values for the current glycofeature
        feature_data <- sample_data %>%
          dplyr::filter(!!glycofeature_sym == glyco_feature)

        # If no data for this glycofeature in this sample, leave as NA
        if (nrow(feature_data) == 0) {
          # Value already initialized to NA
          next
        }

        feature_sum <- sum(feature_data[[value_col]], na.rm = TRUE)

        # If feature sum is NA or zero, leave as NA
        if (is.na(feature_sum) || feature_sum == 0) {
          next
        }

        # Calculate relative abundance (percentage)
        relative_abundance <- (feature_sum / total_sum) * 100
        relative_abundance <- format(relative_abundance, scientific = F)
        # Store in result dataframe
        temp_df[i, sample_id] <- relative_abundance
      }
    }

    # Append to main result dataframe
    result_df <- rbind(result_df, temp_df)
  }

  # Convert to matrix format
  # Create row names that combine top_level_group and glycofeature
  rownames <- paste(result_df$top_level_group, result_df$glycofeature, sep = "_")

  # Create the final matrix
  final_matrix <- as.matrix(result_df[, -c(1, 2)]) # Remove the first two columns (top_level_group and glycofeature)
  rownames(final_matrix) <- rownames

  # We can also add group information to the column names if needed
  # This section uses the group_col to add group info to sample names
  if (!is.null(group_col) && group_col %in% colnames(gpeps_dataframe)) {
    # Create a mapping of sample IDs to their group
    sample_to_group <- gpeps_dataframe %>%
      dplyr::select(!!sample_sym, !!group_sym) %>%
      dplyr::distinct()

    # If there are duplicates (a sample appears in multiple groups), take the first occurrence
    sample_to_group <- sample_to_group[!duplicated(sample_to_group[[sample_col]]), ]

    # Create a named vector for easy lookup
    group_lookup <- setNames(
      sample_to_group[[group_col]],
      sample_to_group[[sample_col]]
    )

    # Add group info to column names if available
    new_colnames <- colnames(final_matrix)
    for (i in 1:length(new_colnames)) {
      sample_id <- new_colnames[i]
      if (sample_id %in% names(group_lookup)) {
        group_value <- group_lookup[sample_id]
        new_colnames[i] <- paste(sample_id, group_value, sep = "_")
      }
    }
    colnames(final_matrix) <- new_colnames
  }

  return(final_matrix)
}


protein_gly_class <- glyco_matrix(
  gpeps_dataframe = glycoPSMs,
  top_lev_group = "protein_accessions",
  glycofeature_group = "glycan_class",
  value_col = "intensity",
  sample_col = "sample",
  group_col = "disease_status",
  file_prefix = "protein_gly_class"
)
write.csv(protein_gly_class, file = "output_data/protein_gly_class_RA.csv")

protein_gly_sia <- glyco_matrix(
  gpeps_dataframe = glycoPSMs,
  top_lev_group = "protein_accessions",
  glycofeature_group = "contains_NeuAc",
  value_col = "intensity",
  sample_col = "sample",
  group_col = "disease_status",
  file_prefix = "protein_gly_sia"
)
write.csv(protein_gly_sia, file = "output_data/protein_gly_sia_RA.csv")

protein_gly_fuc <- glyco_matrix(
  gpeps_dataframe = glycoPSMs,
  top_lev_group = "protein_accessions",
  glycofeature_group = "contains_Fuc",
  value_col = "intensity",
  sample_col = "sample",
  group_col = "disease_status",
  file_prefix = "protein_gly_fuc"
)
write.csv(protein_gly_fuc, file = "output_data/protein_gly_fuc_RA.csv")

protein_gly_sia_count <- glyco_matrix(
  gpeps_dataframe = glycoPSMs,
  top_lev_group = "protein_accessions",
  glycofeature_group = "sia_count",
  value_col = "intensity",
  sample_col = "sample",
  group_col = "disease_status",
  file_prefix = "protein_gly_sia_count"
)
write.csv(protein_gly_sia_count, file = "output_data/protein_gly_sia_count_RA.csv")

protein_gly_comp <- glyco_matrix(
  gpeps_dataframe = glycoPSMs,
  top_lev_group = "protein_accessions",
  glycofeature_group = "glycan_composition",
  value_col = "intensity",
  sample_col = "sample",
  group_col = "disease_status",
  file_prefix = "protein_gly_comp"
)
write.csv(protein_gly_comp, file = "output_data/protein_gly_comp_RA.csv")

glycosite_gly_class <- glyco_matrix(
  gpeps_dataframe = glycoPSMs,
  top_lev_group = "gsite_ID",
  glycofeature_group = "glycan_class",
  value_col = "intensity",
  sample_col = "sample",
  group_col = "disease_status",
  file_prefix = "glycosite_gly_class"
)
write.csv(glycosite_gly_class, file = "output_data/glycosite_gly_class_RA.csv")

glycosite_gly_sia <- glyco_matrix(
  gpeps_dataframe = glycoPSMs,
  top_lev_group = "gsite_ID",
  glycofeature_group = "contains_NeuAc",
  value_col = "intensity",
  sample_col = "sample",
  group_col = "disease_status",
  file_prefix = "glycosite_gly_sia"
)
write.csv(glycosite_gly_sia, file = "output_data/glycosite_gly_sia_RA.csv")

glycosite_gly_fuc <- glyco_matrix(
  gpeps_dataframe = glycoPSMs,
  top_lev_group = "gsite_ID",
  glycofeature_group = "contains_Fuc",
  value_col = "intensity",
  sample_col = "sample",
  group_col = "disease_status",
  file_prefix = "glycosite_gly_fuc"
)
write.csv(glycosite_gly_fuc, file = "output_data/glycosite_gly_fuc_RA.csv")

glycosite_gly_sia_count <- glyco_matrix(
  gpeps_dataframe = glycoPSMs,
  top_lev_group = "gsite_ID",
  glycofeature_group = "sia_count",
  value_col = "intensity",
  sample_col = "sample",
  group_col = "disease_status",
  file_prefix = "glycosite_gly_sia_count"
)
write.csv(glycosite_gly_sia_count, file = "output_data/glycosite_gly_sia_count_RA.csv")

glycosite_gly_comp <- glyco_matrix(
  gpeps_dataframe = glycoPSMs,
  top_lev_group = "gsite_ID",
  glycofeature_group = "glycan_composition",
  value_col = "intensity",
  sample_col = "sample",
  group_col = "disease_status",
  file_prefix = "glycosite_gly_comp"
)
write.csv(glycosite_gly_comp, file = "output_data/glycosite_gly_comp_RA.csv")

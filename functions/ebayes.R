## Empirical Bayes Differential Expression Analysis with CLR Transformation
# Load required packages
if (!requireNamespace("limma", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("limma")
}
if (!requireNamespace("compositions", quietly = TRUE)) {
  install.packages("compositions")
}
library(limma)
library(compositions)
library(dplyr)
library(rlang)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Convert glycopeptide long format data to wide matrix with relative abundances
#' @param gpeps_dataframe the input glycopsm dataframe with all the glyco features annotated to the data frame 
#' @param top_lev_group the top level grouping is how you subset and group data to perform glycan feature analysis at eg. protein or glycosites
#' @param glycofeature_group the glycan feature to analyze (e.g., glycan_composition, glycan_class)
#' @param value_col Column name for the measurement values to analyze
#' @param sample_col Column name to use for pivoting to wide format
#' @param group_col Column name containing group information
#' @param min_samples Minimum number of samples required per group (default: 3)
#' @param file_prefix Prefix for output files (default: "Analysis_")
#' @return a matrix of relative abundances for each top level grouping and for each sample 
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
      for (i in seq_len(nrow(temp_df))) {
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
    for (i in seq_along(new_colnames)) {
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

#' Centered Log-Ratio (CLR) Transform Matrix Data using compositions package
#' @param data_matrix A numeric matrix containing the data to be transformed
#' @param pseudo_count A small number to add before log transformation to handle zeros (default: 1)
#' @param na_handling Strategy for handling NA values: 
#'        "keep" (default) - keep NAs as is
#'        "remove" - remove rows with any NAs
#'        "impute_min" - replace NAs with minimum non-NA value in dataset
#'        "impute_zero" - replace NAs with zero before adding pseudo count
#' @param min_non_na Minimum number of non-NA values required per row to keep row (default: 1)
#' @return A matrix of the same dimensions as input with CLR-transformed values
clr_transform_matrix <- function(data_matrix,
                               pseudo_count = 1,
                               na_handling = "keep",
                               min_non_na = 1) {
  
  # Input validation
  if (!is.matrix(data_matrix)) {
    stop("Input must be a matrix")
  }
  
  # Convert matrix to numeric if it isn't already
  data_matrix <- matrix(as.numeric(data_matrix), 
                       nrow = nrow(data_matrix),
                       dimnames = dimnames(data_matrix))
  
  # Handle NA values according to specified strategy
  if (na_handling == "remove") {
    # Count non-NA values per row
    non_na_count <- rowSums(!is.na(data_matrix))
    # Keep only rows with sufficient non-NA values
    data_matrix <- data_matrix[non_na_count >= min_non_na, , drop = FALSE]
    
  } else if (na_handling == "impute_min") {
    # Find minimum non-NA value in entire dataset
    min_value <- min(data_matrix, na.rm = TRUE)
    # Replace NAs with minimum value
    data_matrix[is.na(data_matrix)] <- min_value
    
  } else if (na_handling == "impute_zero") {
    # Replace NAs with 0
    data_matrix[is.na(data_matrix)] <- 0
    
  } else if (na_handling != "keep") {
    warning("Invalid na_handling option. Using 'keep' as default.")
  }
  
  # Add pseudo count to handle zeros
  data_matrix <- data_matrix + pseudo_count
  
  # Initialize matrix for CLR transformed values
  clr_matrix <- matrix(NA, 
                      nrow = nrow(data_matrix), 
                      ncol = ncol(data_matrix),
                      dimnames = dimnames(data_matrix))
  
  # Perform CLR transformation row by row using compositions package
  for (i in seq_len(nrow(data_matrix))) {
    row_data <- data_matrix[i, ]
    
    # Skip rows with all NAs
    if (all(is.na(row_data))) {
      next
    }
    
    # Use compositions package for CLR transformation
    # Convert to acomp (compositional data) and then apply clr
    tryCatch({
      # Create compositional data object
      comp_data <- compositions::acomp(row_data)
      # Apply CLR transformation
      clr_values <- compositions::clr(comp_data)
      clr_matrix[i, ] <- as.numeric(clr_values)
    }, error = function(e) {
      # If compositions package fails, fall back to manual calculation
      warning(paste("CLR transformation failed for row", i, "using fallback method:", e$message))
      geom_mean <- exp(mean(log(row_data), na.rm = TRUE))
      clr_matrix[i, ] <- log(row_data / geom_mean)
    })
  }
  
  return(clr_matrix)
}

# =============================================================================
# MAIN ANALYSIS PIPELINE
# =============================================================================

cat("=== Empirical Bayes Analysis with CLR Transformation ===\n")
cat("Processing glyco_peptide_groups_long data...\n\n")

# Check if glyco_peptide_groups_long exists
if (!exists("glyco_peptide_groups_long")) {
  stop("glyco_peptide_groups_long data not found. Please ensure the data is loaded.")
}

# =============================================================================
# 1. CREATE RELATIVE ABUNDANCE MATRICES
# =============================================================================

cat("Step 1: Creating relative abundance matrices...\n")

# Protein-level glycan composition matrix
cat("Creating protein-level glycan composition matrix...\n")
protein_gly_comp <- glyco_matrix(
  gpeps_dataframe = glyco_peptide_groups_long,
  top_lev_group = "protein_accessions",
  glycofeature_group = "glycan_composition",
  value_col = "abundance",
  sample_col = "sample",
  group_col = "group",
  file_prefix = "protein_gly_comp"
)

# Glycosite-level glycan composition matrix
cat("Creating glycosite-level glycan composition matrix...\n")
glycosite_gly_comp <- glyco_matrix(
  gpeps_dataframe = glyco_peptide_groups_long,
  top_lev_group = "gsite_ID",
  glycofeature_group = "glycan_composition",
  value_col = "abundance",
  sample_col = "sample",
  group_col = "group",
  file_prefix = "glycosite_gly_comp"
)

cat("✓ Relative abundance matrices created\n")
cat("  Protein-level features:", nrow(protein_gly_comp), "\n")
cat("  Glycosite-level features:", nrow(glycosite_gly_comp), "\n\n")

# =============================================================================
# 2. APPLY CLR TRANSFORMATION
# =============================================================================

cat("Step 2: Applying CLR transformation...\n")

# Apply CLR transformation to both matrices
protein_gly_comp_clr <- clr_transform_matrix(protein_gly_comp, pseudo_count = 1)
glycosite_gly_comp_clr <- clr_transform_matrix(glycosite_gly_comp, pseudo_count = 1)

cat("✓ CLR transformation completed\n\n")

# =============================================================================
# 3. EMPIRICAL BAYES ANALYSIS - PROTEIN LEVEL
# =============================================================================

cat("Step 3: Running empirical Bayes analysis - Protein level...\n")

# Extract group information (HC vs M) from column names
# Handle lowercase sample names from glycopeptide_groups_long
groups_protein <- factor(ifelse(grepl("^hc", colnames(protein_gly_comp_clr), ignore.case = TRUE), "HC", "M"))
print("Protein-level group levels:")
print(levels(groups_protein))

# Create design matrix
design_protein <- model.matrix(~0 + groups_protein)
colnames(design_protein) <- levels(groups_protein)
print("Protein-level design matrix column names:")
print(colnames(design_protein))

# Create contrast matrix for HC vs M comparison
contrast.matrix_protein <- makeContrasts(
  M_vs_HC = M - HC,
  levels = design_protein
)

# Fit linear model
fit_protein <- lmFit(protein_gly_comp_clr, design_protein)

# Fit contrasts
fit2_protein <- contrasts.fit(fit_protein, contrast.matrix_protein)
fit2_protein <- eBayes(fit2_protein)

# Get results
results_protein <- topTable(fit2_protein, 
                   coef = "M_vs_HC", 
                   number = Inf,  # Return all results
                   adjust.method = "BH")  # Benjamini-Hochberg correction

# Add more information to results
results_protein$Feature <- rownames(results_protein)
results_protein$Significant <- !is.na(results_protein$adj.P.Val) & results_protein$adj.P.Val < 0.05
results_protein$Direction <- ifelse(!is.na(results_protein$logFC) & results_protein$logFC > 0, "Up in MECFS", "Down in MECFS")

# Print summary
cat("\nProtein-level Differential Analysis Summary:\n")
cat("Total features tested:", nrow(results_protein), "\n")
cat("Features with NA p-values:", sum(is.na(results_protein$adj.P.Val)), "\n")
cat("Significant features (FDR < 0.05):", sum(results_protein$Significant, na.rm = TRUE), "\n")
cat("  Up in MECFS:", sum(results_protein$Significant & results_protein$logFC > 0, na.rm = TRUE), "\n")
cat("  Down in MECFS:", sum(results_protein$Significant & results_protein$logFC < 0, na.rm = TRUE), "\n")

# Save results
write.csv(results_protein, "output_data/protein_glycan_composition_limma_clr_results.csv")
cat("✓ Protein-level results saved to output_data/protein_glycan_composition_limma_clr_results.csv\n\n")

# =============================================================================
# 4. EMPIRICAL BAYES ANALYSIS - GLYCOSITE LEVEL
# =============================================================================

cat("Step 4: Running empirical Bayes analysis - Glycosite level...\n")

# Extract group information (HC vs M) from column names
# Handle lowercase sample names from glycopeptide_groups_long
groups_glycosite <- factor(ifelse(grepl("^hc", colnames(glycosite_gly_comp_clr), ignore.case = TRUE), "HC", "M"))
print("Glycosite-level group levels:")
print(levels(groups_glycosite))

# Create design matrix
design_glycosite <- model.matrix(~0 + groups_glycosite)
colnames(design_glycosite) <- levels(groups_glycosite)
print("Glycosite-level design matrix column names:")
print(colnames(design_glycosite))

# Create contrast matrix for HC vs M comparison
contrast.matrix_glycosite <- makeContrasts(
  M_vs_HC = M - HC,
  levels = design_glycosite
)

# Fit linear model
fit_glycosite <- lmFit(glycosite_gly_comp_clr, design_glycosite)

# Fit contrasts
fit2_glycosite <- contrasts.fit(fit_glycosite, contrast.matrix_glycosite)
fit2_glycosite <- eBayes(fit2_glycosite)

# Get results
results_glycosite <- topTable(fit2_glycosite, 
                   coef = "M_vs_HC", 
                   number = Inf,  # Return all results
                   adjust.method = "BH")  # Benjamini-Hochberg correction

# Add more information to results
results_glycosite$Feature <- rownames(results_glycosite)
results_glycosite$Significant <- !is.na(results_glycosite$adj.P.Val) & results_glycosite$adj.P.Val < 0.05
results_glycosite$Direction <- ifelse(!is.na(results_glycosite$logFC) & results_glycosite$logFC > 0, "Up in MECFS", "Down in MECFS")

# Print summary
cat("\nGlycosite-level Differential Analysis Summary:\n")
cat("Total features tested:", nrow(results_glycosite), "\n")
cat("Features with NA p-values:", sum(is.na(results_glycosite$adj.P.Val)), "\n")
cat("Significant features (FDR < 0.05):", sum(results_glycosite$Significant, na.rm = TRUE), "\n")
cat("  Up in MECFS:", sum(results_glycosite$Significant & results_glycosite$logFC > 0, na.rm = TRUE), "\n")
cat("  Down in MECFS:", sum(results_glycosite$Significant & results_glycosite$logFC < 0, na.rm = TRUE), "\n")

# Save results
write.csv(results_glycosite, "output_data/glycosite_glycan_composition_limma_clr_results.csv")
cat("✓ Glycosite-level results saved to output_data/glycosite_glycan_composition_limma_clr_results.csv\n\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("=== ANALYSIS COMPLETE ===\n")
cat("Pipeline: Raw Abundance → Relative Abundance → CLR Transformation (compositions package) → Empirical Bayes\n")
cat("Results saved:\n")
cat("  - Protein level: output_data/protein_glycan_composition_limma_clr_results.csv\n")
cat("  - Glycosite level: output_data/glycosite_glycan_composition_limma_clr_results.csv\n")
cat("Each result file contains: Feature, logFC, AveExpr, t, P.Value, adj.P.Val, B, Significant, Direction\n")
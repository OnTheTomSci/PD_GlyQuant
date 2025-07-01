# Create and save heatmaps with clustering enabled
# Protein glycan composition heatmap
create_glyco_heatmap(
  protein_gly_comp_log,
  title = "Protein Glycan Composition",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  filename = "figures/protein_glycan_comp_heatmap.png"
)

# Handle NA values before creating the heatmap
glycosite_gly_comp_log_clean <- glycosite_gly_comp_log
glycosite_gly_comp_log_clean[is.na(glycosite_gly_comp_log_clean)] <- 0  # or another appropriate value
glycosite_gly_comp_log_clean[is.infinite(glycosite_gly_comp_log_clean)] <- 0

create_glyco_heatmap(
  glycosite_gly_comp_log_clean,
  title = "Glycosite Glycan Composition",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  filename = "figures/glycosite_glycan_comp_heatmap.png"
)

# Select top glycosites that have complete data
select_top_complete_glycosites <- function(data_matrix, n = 20) {
  # Find rows with complete data (no NAs or Inf)
  complete_rows <- which(apply(data_matrix, 1, function(x) !any(is.na(x) | is.infinite(x))))
  
  if(length(complete_rows) == 0) {
    stop("No rows with complete data found")
  }
  
  # Get subset of matrix with only complete rows
  complete_matrix <- data_matrix[complete_rows, , drop = FALSE]
  
  # Calculate variance for each complete row
  row_vars <- apply(complete_matrix, 1, var)
  
  # Get indices of top n rows by variance
  n_to_select <- min(n, length(row_vars))
  top_indices <- order(row_vars, decreasing = TRUE)[1:n_to_select]
  
  # Return the subsetted matrix with only top n rows
  return(complete_matrix[top_indices, , drop = FALSE])
}

# Get top glycosites with complete data
top_glycoproteins <- select_top_complete_glycosites(protein_gly_comp_log, n = 20)

# Create heatmap with only the top glycosites
create_glyco_heatmap(
  top_glycoproteins,
  title = "Top 20 Protein Glycan Composition",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  filename = "figures/top20_protein_glycan_comp_heatmap.png"
)

# Get top glycosites with complete data
top_glycosites <- select_top_complete_glycosites(glycosite_gly_comp_log, n = 20)

# Create heatmap with only the top glycosites
create_glyco_heatmap(
  top_glycosites,
  title = "Top 20 Glycosite Glycan Composition",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  filename = "figures/top20_glycosite_glycan_comp_heatmap.png"
)

# Debug and create sialic acid count heatmap
# First, inspect the data
print("Checking protein_gly_sia_count_log data:")
print(str(protein_gly_sia_count_log))
print("Dimensions:")
print(dim(protein_gly_sia_count_log))
print("Any NA values:")
print(sum(is.na(protein_gly_sia_count_log)))
print("Any infinite values:")
print(sum(is.infinite(as.matrix(protein_gly_sia_count_log))))

# Clean the data before creating heatmap
protein_gly_sia_count_clean <- protein_gly_sia_count_log

# Convert to matrix if it isn't already
if(!is.matrix(protein_gly_sia_count_clean)) {
  protein_gly_sia_count_clean <- as.matrix(protein_gly_sia_count_clean)
}

# Remove rows with all NAs
protein_gly_sia_count_clean <- protein_gly_sia_count_clean[
  rowSums(!is.na(protein_gly_sia_count_clean)) > 0,
]

# Replace infinite values with NA
protein_gly_sia_count_clean[is.infinite(protein_gly_sia_count_clean)] <- NA

# Calculate row means and SDs for scaling
row_means <- rowMeans(protein_gly_sia_count_clean, na.rm = TRUE)
row_sds <- apply(protein_gly_sia_count_clean, 1, sd, na.rm = TRUE)

# Scale the data
protein_gly_sia_count_scaled <- sweep(protein_gly_sia_count_clean, 1, row_means, "-")
protein_gly_sia_count_scaled <- sweep(protein_gly_sia_count_scaled, 1, row_sds, "/")

# Replace any remaining NA/Inf values with 0
protein_gly_sia_count_scaled[is.na(protein_gly_sia_count_scaled)] <- 0
protein_gly_sia_count_scaled[is.infinite(protein_gly_sia_count_scaled)] <- 0

# Now try creating the heatmap with the cleaned data
create_glyco_heatmap(
  protein_gly_sia_count_scaled,
  title = "Protein Glycan Sialic Acid Count", 
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "none",  # Data is already scaled
  filename = "figures/protein_glycan_sia_count_heatmap.png"
)

#' Analyze Log-Transformed Data and Create Volcano Plots
#' 
#' This function performs statistical analysis on log-transformed data, comparing two groups 
#' and generates volcano plots to visualize the results.
#' 
#' @param log_matrix A numeric matrix containing log-transformed data. Rows represent 
#'        features and columns represent samples. Column names should include group info.
#' @param title Character string for the plot title (default: "Volcano Plot")
#' @param fc_cutoff Log2 fold change cutoff for significance (default: 1)
#' @param p_cutoff P-value cutoff for significance (default: 0.05)
#' @param group_pattern Pattern for identifying groups in column names (default: c("HC", "M"))
#' @param group_names Names for the groups (default: c("Healthy", "MECFS"))
#'
#' @return A list containing:
#'   - results: Data frame with statistical results
#'   - plot: ggplot2 volcano plot object
#'   - sample_groups: Group assignments
#'
#' @export
analyze_log_data <- function(log_matrix, 
                           title = "Volcano Plot",
                           fc_cutoff = 1,
                           p_cutoff = 0.05,
                           group_pattern = c("HC", "M"),
                           group_names = c("Healthy", "MECFS")) {
  
  # Input validation
  if (!is.matrix(log_matrix)) {
    stop("Input must be a matrix")
  }
  
  # Create sample group mapping
  sample_groups <- data.frame(
    sample = colnames(log_matrix),
    group = NA,
    stringsAsFactors = FALSE
  )
  
  # Assign groups based on patterns
  for (i in seq_along(group_pattern)) {
    sample_groups$group[grepl(paste0("^", group_pattern[i]), sample_groups$sample)] <- group_names[i]
  }
  
  # Print sample grouping for verification
  cat("\nSample Grouping for", title, ":\n")
  print(sample_groups)
  
  # Check for unassigned samples
  if(any(is.na(sample_groups$group))) {
    warning("Some samples could not be assigned to groups: ", 
            paste(sample_groups$sample[is.na(sample_groups$group)], collapse=", "))
  }
  
  # Verify sample sizes
  group_counts <- table(sample_groups$group)
  cat("\nSamples per group:\n")
  print(group_counts)
  
  if(length(group_counts) < 2 || any(group_counts < 2)) {
    stop("Need at least 2 samples in each group for comparison")
  }
  
  # Initialize results dataframe
  results <- data.frame(
    Feature = rownames(log_matrix),
    P_Value = NA_real_,
    FDR = NA_real_,
    Log2_FC = NA_real_,
    Mean_Group1 = NA_real_,
    Mean_Group2 = NA_real_,
    SD_Group1 = NA_real_,
    SD_Group2 = NA_real_,
    N_Group1 = NA_integer_,
    N_Group2 = NA_integer_,
    stringsAsFactors = FALSE
  )
  
  # Calculate statistics for each feature
  for(i in 1:nrow(log_matrix)) {
    feature_data <- log_matrix[i,]
    
    # Split data by group
    group1_data <- feature_data[sample_groups$group == group_names[1]]
    group2_data <- feature_data[sample_groups$group == group_names[2]]
    
    # Store sample counts
    results$N_Group1[i] <- sum(!is.na(group1_data))
    results$N_Group2[i] <- sum(!is.na(group2_data))
    
    # Calculate means and SDs
    results$Mean_Group1[i] <- mean(group1_data, na.rm = TRUE)
    results$Mean_Group2[i] <- mean(group2_data, na.rm = TRUE)
    results$SD_Group1[i] <- sd(group1_data, na.rm = TRUE)
    results$SD_Group2[i] <- sd(group2_data, na.rm = TRUE)
    
    # Calculate log2 fold change (difference of means for log-transformed data)
    results$Log2_FC[i] <- results$Mean_Group2[i] - results$Mean_Group1[i]
    
    # Perform t-test if we have enough samples
    if(results$N_Group1[i] >= 3 && results$N_Group2[i] >= 3) {
      t_test <- try({
        t.test(group2_data, group1_data, var.equal = FALSE)
      })
      
      if(!inherits(t_test, "try-error")) {
        results$P_Value[i] <- t_test$p.value
      }
    }
  }
  
  # Calculate FDR
  results$FDR <- p.adjust(results$P_Value, method = "BH")
  
  # Create volcano plot
  volcano_plot <- EnhancedVolcano(
    results,
    lab = results$Feature,
    x = 'Log2_FC',
    y = 'FDR',
    title = title,
    subtitle = paste0(
      'Log-transformed data\n',
      group_names[1], ' n=', group_counts[group_names[1]], 
      ', ', group_names[2], ' n=', group_counts[group_names[2]], '\n',
      'Log2FC cutoff = ±', fc_cutoff, 
      '; p-value cutoff = ', p_cutoff
    ),
    pCutoff = p_cutoff,
    FCcutoff = fc_cutoff,
    pointSize = 2.0,
    labSize = 3.0,
    col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
    colAlpha = 0.5,
    legendPosition = 'right',
    legendLabSize = 10,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    caption = paste("Total features:", nrow(results))
  )
  
  # Return results
  return(list(
    results = results,
    plot = volcano_plot,
    sample_groups = sample_groups
  ))
}

# Example usage:
# results <- analyze_log_data(protein_gly_comp_log, 
#                           title = "Glycan Composition Analysis",
#                           fc_cutoff = 1,
#                           p_cutoff = 0.05)
# 
# # Save results
# write.csv(results$results, "glycan_composition_analysis.csv")
# 
# # Save plot
# ggsave("glycan_composition_volcano.png", results$plot, width = 10, height = 8)



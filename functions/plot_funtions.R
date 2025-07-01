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




#' Analyze fucosylation relative abundance with comprehensive statistical testing
#' @param relative_abundance_df The relative abundance dataframe
#' @param feature_col The column name containing the feature categories (e.g., "contains_Fuc", "fuc_count")
#' @param output_file The output CSV file path for statistical results
#' @param plot_file The output plot file path for the visualization
#' @return A list containing statistical results and plot
analyze_fucosylation_abundance <- function(relative_abundance_df, 
                                          feature_col, 
                                          output_file, 
                                          plot_file) {
  # Load required libraries
  library(ggplot2)
  library(viridis)
  library(dplyr)
  library(stats)
  library(stringr) # for str_ends
  
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
    feature = character(),
    healthy_n = numeric(),
    mecfs_n = numeric(),
    healthy_mean = numeric(),
    mecfs_mean = numeric(),
    healthy_median = numeric(),
    mecfs_median = numeric(),
    healthy_sd = numeric(),
    mecfs_sd = numeric(),
    normality_test_p = numeric(),
    is_normal = logical(),
    test_used = character(),
    test_statistic = numeric(),
    p_value = numeric(),
    fdr_p_value = numeric(),
    significant = logical(),
    effect_size = numeric(),
    stringsAsFactors = FALSE
  )
  
  # For each feature category, perform analysis
  for (feature_value in unique_features) {
    # Filter data for this feature
    feature_data <- relative_abundance_df %>% 
      dplyr::filter(!!rlang::sym(feature_col) == feature_value)
    
    # Split data by disease status
    healthy_data <- feature_data %>% 
      dplyr::filter(disease_status == "Healthy") %>%
      dplyr::pull(relative_percentage)
    
    mecfs_data <- feature_data %>% 
      dplyr::filter(disease_status == "MECFS") %>%
      dplyr::pull(relative_percentage)
    
    # Check if we have enough samples
    if (length(healthy_data) < 3 || length(mecfs_data) < 3) {
      warning(paste("Not enough samples for feature", feature_value, 
                   "(Healthy:", length(healthy_data), 
                   ", MECFS:", length(mecfs_data), ")"))
      next
    }
    
    # Calculate descriptive statistics
    healthy_mean <- mean(healthy_data, na.rm = TRUE)
    mecfs_mean <- mean(mecfs_data, na.rm = TRUE)
    healthy_median <- median(healthy_data, na.rm = TRUE)
    mecfs_median <- median(mecfs_data, na.rm = TRUE)
    healthy_sd <- sd(healthy_data, na.rm = TRUE)
    mecfs_sd <- sd(mecfs_data, na.rm = TRUE)
    
    # Perform Shapiro-Wilk test for normality on combined data
    combined_data <- c(healthy_data, mecfs_data)
    if (length(unique(combined_data[!is.na(combined_data)])) > 1) {
      normality_test <- shapiro.test(combined_data)
      is_normal <- normality_test$p.value > 0.05
    } else {
      normality_test <- list(p.value = NA)
      is_normal <- FALSE
      warning(paste("All values are identical or not enough unique values for feature", feature_value, "- skipping normality test."))
    }
    
    # Perform appropriate statistical test
    if (is_normal) {
      # Use t-test for normally distributed data
      test_result <- t.test(mecfs_data, healthy_data, var.equal = FALSE)
      test_used <- "Welch's t-test"
      test_statistic <- test_result$statistic
      p_value <- test_result$p.value
      effect_size <- abs(mecfs_mean - healthy_mean) / sqrt((healthy_sd^2 + mecfs_sd^2) / 2)
    } else {
      # Use Wilcoxon test for non-normally distributed data
      test_result <- wilcox.test(mecfs_data, healthy_data, exact = FALSE)
      test_used <- "Wilcoxon rank-sum test"
      test_statistic <- test_result$statistic
      p_value <- test_result$p.value
      # Calculate effect size for non-parametric test
      z_stat <- qnorm(p_value / 2)
      n_total <- length(healthy_data) + length(mecfs_data)
      effect_size <- abs(z_stat) / sqrt(n_total)
    }
    
    # Add results to dataframe
    results <- rbind(results, data.frame(
      feature = feature_value,
      healthy_n = length(healthy_data),
      mecfs_n = length(mecfs_data),
      healthy_mean = healthy_mean,
      mecfs_mean = mecfs_mean,
      healthy_median = healthy_median,
      mecfs_median = mecfs_median,
      healthy_sd = healthy_sd,
      mecfs_sd = mecfs_sd,
      normality_test_p = normality_test$p.value,
      is_normal = is_normal,
      test_used = test_used,
      test_statistic = test_statistic,
      p_value = p_value,
      fdr_p_value = NA,  # Will be calculated later
      significant = p_value < 0.05,
      effect_size = effect_size,
      stringsAsFactors = FALSE
    ))
  }
  
  # Calculate FDR (Benjamini-Hochberg correction)
  if (nrow(results) > 0) {
    results$fdr_p_value <- p.adjust(results$p_value, method = "BH")
    results$significant_fdr <- results$fdr_p_value < 0.05
  }
  
  # Arrange by p-value
  results <- results %>%
    dplyr::arrange(p_value)
  
  # Save results
  write.csv(results, 
            file = output_file, 
            row.names = FALSE)
  
  # Create visualization
  if (nrow(results) > 0) {
    # Prepare data for plotting
    plot_data <- relative_abundance_df %>%
      dplyr::filter(!!rlang::sym(feature_col) %in% results$feature[1:min(10, nrow(results))]) %>%
      dplyr::mutate(feature_label = paste0(!!rlang::sym(feature_col), 
                                          "\np = ", format(results$p_value[match(!!rlang::sym(feature_col), results$feature)], 
                                                          scientific = TRUE, digits = 2)))
    
    # Create box plot with statistical annotations
    p <- ggplot(plot_data, 
                aes(x = feature_label, 
                    y = relative_percentage, 
                    fill = disease_status)) +
      geom_boxplot(outlier.shape = 21, 
                   outlier.fill = "white",
                   outlier.alpha = 0.6,
                   position = position_dodge(width = 0.8),
                   alpha = 0.7) +
      geom_point(position = position_jitterdodge(jitter.width = 0.2),
                 alpha = 0.4,
                 size = 1) +
      scale_fill_viridis(discrete = TRUE, option = "D") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "top",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "cm")
      ) +
      labs(title = paste("Fucosylation Relative Abundance Analysis:", feature_col),
           x = "Feature (with p-value)",
           y = "Relative Percentage (%)",
           fill = "Disease Status")
    
    # Save plot
    ggsave(plot_file, 
           p, 
           width = 12, 
           height = 8, 
           dpi = 300,
           bg = "white")
  }
  
  # Print summary
  cat("\nFucosylation Analysis Summary for", feature_col, ":\n")
  cat("Total features analyzed:", nrow(results), "\n")
  cat("Features with normal distribution:", sum(results$is_normal), "\n")
  cat("Features with significant difference (p < 0.05):", sum(results$significant), "\n")
  cat("Features with significant difference (FDR < 0.05):", sum(results$significant_fdr), "\n")
  
  # Print top 5 significant features
  if (sum(results$significant) > 0) {
    cat("\nTop 5 Significant Features (p < 0.05):\n")
    top_5 <- results %>% 
      dplyr::filter(significant) %>% 
      dplyr::arrange(p_value) %>% 
      head(5)
    
    for (i in 1:nrow(top_5)) {
      cat(i, ". ", top_5$feature[i], 
          " (", top_5$test_used[i],
          ", p = ", format(top_5$p_value[i], scientific = TRUE, digits = 2),
          ", FDR = ", format(top_5$fdr_p_value[i], scientific = TRUE, digits = 2),
          ")\n", sep = "")
      cat("   Healthy: n = ", top_5$healthy_n[i], 
          ", Mean = ", round(top_5$healthy_mean[i], 2), 
          "%, SD = ", round(top_5$healthy_sd[i], 2), "%\n", sep = "")
      cat("   MECFS: n = ", top_5$mecfs_n[i], 
          ", Mean = ", round(top_5$mecfs_mean[i], 2), 
          "%, SD = ", round(top_5$mecfs_sd[i], 2), "%\n", sep = "")
    }
  }
  
  return(list(results = results, plot = p))
}


library(tidyr)
library(dplyr)

# Read the wide-format CSV
glycosite_gly_fuc_RA <- read.csv("output_data/glycosite/RA/glycosite_gly_fuc_RA.csv")

# Rename the first column to glycosite_id
colnames(glycosite_gly_fuc_RA)[1] <- "glycosite_id"

# Filter to only include fucosylated glycosites (ending in TRUE)
glycosite_gly_fuc_RA_fuc_only <- glycosite_gly_fuc_RA %>%
  filter(str_ends(glycosite_id, "TRUE"))

# Reshape to long format
glycosite_gly_fuc_RA_long <- glycosite_gly_fuc_RA_fuc_only %>%
  pivot_longer(
    cols = -glycosite_id,
    names_to = "sample",
    values_to = "relative_percentage"
  ) %>%
  mutate(
    disease_status = ifelse(grepl("Healthy", sample), "Healthy", "MECFS")
  )

# Remove NA values
glycosite_gly_fuc_RA_long_clean <- glycosite_gly_fuc_RA_long %>%
  filter(!is.na(relative_percentage))

# Perform statistical analysis comparing Healthy vs ME/CFS for each fucosylated glycosite
# Initialize results dataframe
results <- data.frame(
  glycosite_id = character(),
  healthy_n = numeric(),
  mecfs_n = numeric(),
  healthy_mean = numeric(),
  mecfs_mean = numeric(),
  healthy_sd = numeric(),
  mecfs_sd = numeric(),
  test_used = character(),
  test_statistic = numeric(),
  p_value = numeric(),
  fdr_p_value = numeric(),
  significant = logical(),
  effect_size = numeric(),
  stringsAsFactors = FALSE
)

# Get unique glycosite IDs
unique_glycosites <- unique(glycosite_gly_fuc_RA_long_clean$glycosite_id)

# Analyze each glycosite
for (glyc_id in unique_glycosites) {
  # Filter data for this glycosite
  glyc_data <- glycosite_gly_fuc_RA_long_clean %>%
    filter(glycosite_id == glyc_id)
  
  # Split by disease status
  healthy_data <- glyc_data %>%
    filter(disease_status == "Healthy") %>%
    pull(relative_percentage)
  
  mecfs_data <- glyc_data %>%
    filter(disease_status == "MECFS") %>%
    pull(relative_percentage)
  
  # Check if we have enough samples
  if (length(healthy_data) < 2 || length(mecfs_data) < 2) {
    warning(paste("Not enough samples for glycosite", glyc_id))
    next
  }
  
  # Calculate descriptive statistics
  healthy_mean <- mean(healthy_data, na.rm = TRUE)
  mecfs_mean <- mean(mecfs_data, na.rm = TRUE)
  healthy_sd <- sd(healthy_data, na.rm = TRUE)
  mecfs_sd <- sd(mecfs_data, na.rm = TRUE)
  
  # Perform statistical test
  if (length(healthy_data) >= 3 && length(mecfs_data) >= 3) {
    # Try t-test first
    t_test <- try({
      t.test(mecfs_data, healthy_data, var.equal = FALSE)
    })
    
    if(!inherits(t_test, "try-error")) {
      test_used <- "Welch's t-test"
      test_statistic <- t_test$statistic
      p_value <- t_test$p.value
      effect_size <- abs(mecfs_mean - healthy_mean) / sqrt((healthy_sd^2 + mecfs_sd^2) / 2)
    } else {
      # Fall back to Wilcoxon test
      wilcox_test <- try({
        wilcox.test(mecfs_data, healthy_data, exact = FALSE)
      })
      
      if(!inherits(wilcox_test, "try-error")) {
        test_used <- "Wilcoxon rank-sum test"
        test_statistic <- wilcox_test$statistic
        p_value <- wilcox_test$p.value
        z_stat <- qnorm(p_value / 2)
        n_total <- length(healthy_data) + length(mecfs_data)
        effect_size <- abs(z_stat) / sqrt(n_total)
      } else {
        test_used <- "No valid test"
        test_statistic <- NA
        p_value <- NA
        effect_size <- NA
      }
    }
  } else {
    test_used <- "Insufficient data"
    test_statistic <- NA
    p_value <- NA
    effect_size <- NA
  }
  
  # Add results to dataframe
  results <- rbind(results, data.frame(
    glycosite_id = glyc_id,
    healthy_n = length(healthy_data),
    mecfs_n = length(mecfs_data),
    healthy_mean = healthy_mean,
    mecfs_mean = mecfs_mean,
    healthy_sd = healthy_sd,
    mecfs_sd = mecfs_sd,
    test_used = test_used,
    test_statistic = test_statistic,
    p_value = p_value,
    fdr_p_value = NA,  # Will be calculated later
    significant = ifelse(is.na(p_value), FALSE, p_value < 0.05),
    effect_size = effect_size,
    stringsAsFactors = FALSE
  ))
}

# Calculate FDR (Benjamini-Hochberg correction)
if (nrow(results) > 0) {
  results$fdr_p_value <- p.adjust(results$p_value, method = "BH")
  results$significant_fdr <- results$fdr_p_value < 0.05
}

# Arrange by p-value
results <- results %>%
  arrange(p_value)

# Save results
write.csv(results, 
          file = "output_data/glycosite_fucose_analysis_results.csv", 
          row.names = FALSE)

# Create visualization for top significant glycosites
if (nrow(results) > 0) {
  # Get top 10 significant glycosites
  top_glycosites <- results %>%
    filter(significant) %>%
    arrange(p_value) %>%
    slice_head(n = 10) %>%
    pull(glycosite_id)
  
  if (length(top_glycosites) > 0) {
    # Create individual box plots for each glycosite
    for (i in seq_along(top_glycosites)) {
      glyc_id <- top_glycosites[i]
      
      # Get the statistics for this glycosite
      glyc_p_value <- results$p_value[results$glycosite_id == glyc_id]
      glyc_fdr <- results$fdr_p_value[results$glycosite_id == glyc_id]
      glyc_test_stat <- results$test_statistic[results$glycosite_id == glyc_id]
      glyc_effect_size <- results$effect_size[results$glycosite_id == glyc_id]
      
      # Clean the glycosite ID by removing "_TRUE" suffix
      clean_glyc_id <- str_replace(glyc_id, "_TRUE$", "")
      
      # Prepare data for plotting this specific glycosite
      plot_data <- glycosite_gly_fuc_RA_long_clean %>%
        filter(glycosite_id == glyc_id)
      
      # Create individual box plot
      p <- ggplot(plot_data, 
                  aes(x = disease_status, 
                      y = relative_percentage, 
                      fill = disease_status)) +
        geom_boxplot(outlier.shape = 21, 
                     outlier.fill = "white",
                     outlier.alpha = 0.6,
                     alpha = 0.7) +
        geom_point(position = position_jitter(width = 0.2),
                   alpha = 0.6,
                   size = 2) +
        scale_fill_viridis(discrete = TRUE, option = "D") +
        theme_minimal() +
        theme(
          axis.text.x = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "none",
          panel.grid.major = element_line(color = "gray90"),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(1, 1, 1, 1), "cm")
        ) +
        labs(title = paste("Glycosite:", clean_glyc_id),
             subtitle = paste("p =", format(glyc_p_value, digits = 4),
                            ", FDR =", format(glyc_fdr, digits = 4),
                            "\nF-test p =", format(glyc_p_value, digits = 4),
                            ", Effect Size =", format(glyc_effect_size, digits = 3)),
             x = "Disease Status",
             y = "Relative Percentage (%)")
      
      # Save individual plot
      safe_glyc_name <- str_replace_all(clean_glyc_id, "[^a-zA-Z0-9]", "_")
      plot_filename <- file.path("output_data/boxplots", paste0("glycosite_", safe_glyc_name, "_boxplot.png"))
      
      ggsave(plot_filename, 
             p, 
             width = 8, 
             height = 6, 
             dpi = 300,
             bg = "white")
      
      cat("Created box plot for:", clean_glyc_id, "\n")
    }
    
    # Also create a combined plot for overview
    combined_plot_data <- glycosite_gly_fuc_RA_long_clean %>%
      filter(glycosite_id %in% top_glycosites)
    
    combined_p <- ggplot(combined_plot_data, 
                        aes(x = glycosite_id, 
                            y = relative_percentage, 
                            fill = disease_status)) +
      geom_boxplot(outlier.shape = 21, 
                   outlier.fill = "white",
                   outlier.alpha = 0.6,
                   position = position_dodge(width = 0.8),
                   alpha = 0.7) +
      geom_point(position = position_jitterdodge(jitter.width = 0.2),
                 alpha = 0.4,
                 size = 1) +
      scale_fill_viridis(discrete = TRUE, option = "D") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "top",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "cm")
      ) +
      labs(title = "Fucosylated Glycosites: Healthy vs ME/CFS Comparison (Combined View)",
           x = "Glycosite ID",
           y = "Relative Percentage (%)",
           fill = "Disease Status")
    
    # Save combined plot
    ggsave("output_data/boxplots/glycosite_fucose_analysis_combined_plot.png", 
           combined_p, 
           width = 12, 
           height = 8, 
           dpi = 300,
           bg = "white")
  }
}

# Print summary
cat("\nFucosylation Analysis Summary:\n")
cat("Total fucosylated glycosites analyzed:", nrow(results), "\n")
cat("Glycosites with significant difference (p < 0.05):", sum(results$significant, na.rm = TRUE), "\n")
cat("Glycosites with significant difference (FDR < 0.05):", sum(results$significant_fdr, na.rm = TRUE), "\n")

# Print top 5 significant glycosites
if (sum(results$significant, na.rm = TRUE) > 0) {
  cat("\nTop 5 Significant Fucosylated Glycosites (p < 0.05):\n")
  top_5 <- results %>% 
    filter(significant) %>% 
    arrange(p_value) %>% 
    head(5)
  
  for (i in 1:nrow(top_5)) {
    cat(i, ". ", top_5$glycosite_id[i], 
        " (", top_5$test_used[i],
        ", p = ", format(top_5$p_value[i], scientific = TRUE, digits = 2),
        ", FDR = ", format(top_5$fdr_p_value[i], scientific = TRUE, digits = 2),
        ")\n", sep = "")
    cat("   Healthy: n = ", top_5$healthy_n[i], 
        ", Mean = ", round(top_5$healthy_mean[i], 2), 
        "%, SD = ", round(top_5$healthy_sd[i], 2), "%\n", sep = "")
    cat("   ME/CFS: n = ", top_5$mecfs_n[i], 
        ", Mean = ", round(top_5$mecfs_mean[i], 2), 
        "%, SD = ", round(top_5$mecfs_sd[i], 2), "%\n", sep = "")
  }
}

cat("\nResults saved to output_data directory.\n")



#' Analyze protein presence for specific glycans in Healthy vs ME/CFS groups
#' @param glycan_composition_file Path to the glycan composition CSV file
#' @param output_dir Directory to save results and plots
#' @param min_abundance Minimum relative abundance threshold (default: 0)
#' @return List containing analysis results and plots
analyze_protein_glycan_overlap <- function(glycan_composition_file, 
                                          output_dir = "output_data/protein_glycan_analysis",
                                          min_abundance = 0) {
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
  library(gridExtra)
  library(stringr)
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Read the glycan composition data
  glycan_data <- read.csv(glycan_composition_file)
  colnames(glycan_data)[1] <- "glycosite_id"
  
  # Reshape to long format
  glycan_long <- glycan_data %>%
    pivot_longer(
      cols = -glycosite_id,
      names_to = "sample",
      values_to = "relative_abundance"
    ) %>%
    mutate(
      disease_status = ifelse(grepl("Healthy", sample), "Healthy", "MECFS")
    )
  
  # Extract protein and glycan information from glycosite_id
  glycan_long <- glycan_long %>%
    mutate(
      protein_id = str_extract(glycosite_id, "^[^_]+"),
      glycan_composition = str_extract(glycosite_id, "[^_]+$")
    )
  
  # Filter by minimum abundance threshold
  glycan_long_filtered <- glycan_long %>%
    filter(relative_abundance >= min_abundance)
  
  # Get unique glycan compositions
  unique_glycans <- unique(glycan_long_filtered$glycan_composition)
  
  # Initialize results list
  results <- list()
  
  # Analyze each glycan composition
  for (glycan in unique_glycans) {
    cat("Analyzing glycan:", glycan, "\n")
    
    # Filter data for this glycan
    glycan_subset <- glycan_long_filtered %>%
      filter(glycan_composition == glycan)
    
    # Get proteins present in each group
    healthy_proteins <- glycan_subset %>%
      filter(disease_status == "Healthy") %>%
      pull(protein_id) %>%
      unique()
    
    mecfs_proteins <- glycan_subset %>%
      filter(disease_status == "MECFS") %>%
      pull(protein_id) %>%
      unique()
    
    # Calculate overlaps
    common_proteins <- intersect(healthy_proteins, mecfs_proteins)
    healthy_only <- setdiff(healthy_proteins, mecfs_proteins)
    mecfs_only <- setdiff(mecfs_proteins, healthy_proteins)
    
    # Store results
    results[[glycan]] <- list(
      healthy_proteins = healthy_proteins,
      mecfs_proteins = mecfs_proteins,
      common_proteins = common_proteins,
      healthy_only = healthy_only,
      mecfs_only = mecfs_only,
      healthy_count = length(healthy_proteins),
      mecfs_count = length(mecfs_proteins),
      common_count = length(common_proteins),
      healthy_only_count = length(healthy_only),
      mecfs_only_count = length(mecfs_only)
    )
    
    # Create Venn diagram
    if (length(healthy_proteins) > 0 || length(mecfs_proteins) > 0) {
      # Prepare data for Venn diagram
      venn_data <- list(
        Healthy = healthy_proteins,
        `ME/CFS` = mecfs_proteins
      )
      
      # Create a simple text-based summary instead of Venn diagram
      # Save detailed results table
      safe_glycan_name <- str_replace_all(glycan, "[^a-zA-Z0-9]", "_")
      
      # Create detailed results table
      results_table <- data.frame(
        Category = c("Healthy only", "ME/CFS only", "Common"),
        Count = c(length(healthy_only), length(mecfs_only), length(common_proteins)),
        Proteins = c(
          paste(healthy_only, collapse = "; "),
          paste(mecfs_only, collapse = "; "),
          paste(common_proteins, collapse = "; ")
        ),
        stringsAsFactors = FALSE
      )
      
      # Save results table
      table_filename <- file.path(output_dir, paste0("results_", safe_glycan_name, ".csv"))
      write.csv(results_table, table_filename, row.names = FALSE)
      
      # Create a simple bar plot showing the counts
      plot_data <- data.frame(
        Category = c("Healthy only", "ME/CFS only", "Common"),
        Count = c(length(healthy_only), length(mecfs_only), length(common_proteins)),
        stringsAsFactors = FALSE
      )
      
      overlap_plot <- ggplot(plot_data, aes(x = Category, y = Count, fill = Category)) +
        geom_bar(stat = "identity", alpha = 0.7) +
        scale_fill_manual(values = c("Healthy only" = "#440154FF", 
                                   "ME/CFS only" = "#21908CFF", 
                                   "Common" = "#FDE725FF")) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "none"
        ) +
        labs(
          title = paste("Protein Counts for Glycan:", glycan),
          x = "Category",
          y = "Number of Proteins"
        ) +
        geom_text(aes(label = Count), vjust = -0.5, size = 4, fontface = "bold")
      
      # Save the plot
      plot_filename <- file.path(output_dir, paste0("overlap_", safe_glycan_name, ".png"))
      ggsave(plot_filename, overlap_plot, width = 8, height = 6, dpi = 300, bg = "white")
    }
  }
  
  # Create summary table
  summary_data <- data.frame(
    Glycan_Composition = names(results),
    Healthy_Proteins = sapply(results, function(x) x$healthy_count),
    MECFS_Proteins = sapply(results, function(x) x$mecfs_count),
    Common_Proteins = sapply(results, function(x) x$common_count),
    Healthy_Only = sapply(results, function(x) x$healthy_only_count),
    MECFS_Only = sapply(results, function(x) x$mecfs_only_count),
    Total_Proteins = sapply(results, function(x) x$healthy_count + x$mecfs_count - x$common_count),
    stringsAsFactors = FALSE
  )
  
  # Save summary
  write.csv(summary_data, file.path(output_dir, "protein_glycan_summary.csv"), row.names = FALSE)
  
  # Create summary plot
  summary_plot <- ggplot(summary_data, aes(x = reorder(Glycan_Composition, Total_Proteins))) +
    geom_bar(aes(y = Healthy_Proteins, fill = "Healthy"), stat = "identity", alpha = 0.7) +
    geom_bar(aes(y = MECFS_Proteins, fill = "ME/CFS"), stat = "identity", alpha = 0.7) +
    scale_fill_manual(values = c("Healthy" = "#440154FF", "ME/CFS" = "#21908CFF")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "top",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    ) +
    labs(
      title = "Protein Count by Glycan Composition",
      x = "Glycan Composition",
      y = "Number of Proteins",
      fill = "Disease Status"
    )
  
  # Save summary plot
  ggsave(file.path(output_dir, "protein_glycan_summary.png"), 
         summary_plot, 
         width = 12, 
         height = 8, 
         dpi = 300,
         bg = "white")
  
  # Print summary
  cat("\nProtein-Glycan Analysis Summary:\n")
  cat("Total glycan compositions analyzed:", length(unique_glycans), "\n")
  cat("Results saved to:", output_dir, "\n")
  
  # Print top 5 glycans by total protein count
  top_glycans <- summary_data %>%
    arrange(desc(Total_Proteins)) %>%
    head(5)
  
  cat("\nTop 5 Glycans by Total Protein Count:\n")
  for (i in 1:nrow(top_glycans)) {
    cat(i, ". ", top_glycans$Glycan_Composition[i], 
        " (Total: ", top_glycans$Total_Proteins[i],
        ", Healthy: ", top_glycans$Healthy_Proteins[i],
        ", ME/CFS: ", top_glycans$MECFS_Proteins[i],
        ", Common: ", top_glycans$Common_Proteins[i], ")\n", sep = "")
  }
  
  return(list(
    results = results,
    summary = summary_data,
    summary_plot = summary_plot
  ))
}

# Run the protein-glycan overlap analysis
protein_glycan_analysis <- analyze_protein_glycan_overlap(
  glycan_composition_file = "output_data/glycosite/RA/glycosite_gly_comp_RA.csv",
  output_dir = "output_data/protein_glycan_analysis",
  min_abundance = 0  # You can adjust this threshold as needed
)


#' Analyze proteins with fucosylated glycans to identify common vs unique proteins
#' @param fucose_data_file Path to the fucose data CSV file
#' @param output_dir Directory to save results and plots
#' @return List containing analysis results and plots
analyze_fucose_protein_distribution <- function(fucose_data_file = "output_data/glycosite/RA/glycosite_gly_fuc_RA.csv",
                                               output_dir = "output_data/fucose_protein_analysis") {
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(viridis)
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Read the fucose data
  fucose_data <- read.csv(fucose_data_file)
  colnames(fucose_data)[1] <- "glycosite_id"
  
  # Extract protein information from glycosite_id
  fucose_data <- fucose_data %>%
    mutate(
      protein_id = str_extract(glycosite_id, "^[^_]+"),
      glycosite_position = str_extract(glycosite_id, "(?<=_)\\d+(?=_)"),
      fucose_status = str_extract(glycosite_id, "[^_]+$")
    )
  
  # Filter to only include fucosylated glycosites (TRUE)
  fucose_data_fuc <- fucose_data %>%
    filter(fucose_status == "TRUE")
  
  # Reshape to long format
  fucose_long <- fucose_data_fuc %>%
    pivot_longer(
      cols = -c(glycosite_id, protein_id, glycosite_position, fucose_status),
      names_to = "sample",
      values_to = "relative_abundance"
    ) %>%
    mutate(
      disease_status = ifelse(grepl("Healthy", sample), "Healthy", "MECFS")
    )
  
  # Get proteins present in each group (no threshold filtering)
  healthy_proteins <- fucose_long %>%
    filter(disease_status == "Healthy") %>%
    pull(protein_id) %>%
    unique()
  
  mecfs_proteins <- fucose_long %>%
    filter(disease_status == "MECFS") %>%
    pull(protein_id) %>%
    unique()
  
  # Calculate overlaps
  common_proteins <- intersect(healthy_proteins, mecfs_proteins)
  healthy_only <- setdiff(healthy_proteins, mecfs_proteins)
  mecfs_only <- setdiff(mecfs_proteins, healthy_proteins)
  
  # Create summary statistics for figures
  summary_stats <- data.frame(
    Category = c("Common to Both", "Healthy Only", "ME/CFS Only", "Total Unique"),
    Count = c(length(common_proteins), length(healthy_only), length(mecfs_only), 
              length(healthy_only) + length(mecfs_only)),
    Percentage = c(
      round(length(common_proteins) / (length(healthy_proteins) + length(mecfs_proteins) - length(common_proteins)) * 100, 1),
      round(length(healthy_only) / (length(healthy_proteins) + length(mecfs_proteins) - length(common_proteins)) * 100, 1),
      round(length(mecfs_only) / (length(healthy_proteins) + length(mecfs_proteins) - length(common_proteins)) * 100, 1),
      round((length(healthy_only) + length(mecfs_only)) / (length(healthy_proteins) + length(mecfs_proteins) - length(common_proteins)) * 100, 1)
    ),
    stringsAsFactors = FALSE
  )
  
  # Create visualization
  # Bar plot of protein counts
  count_plot <- ggplot(summary_stats, aes(x = Category, y = Count, fill = Category)) +
    geom_bar(stat = "identity", alpha = 0.7) +
    scale_fill_viridis(discrete = TRUE, option = "D") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "none"
    ) +
    labs(
      title = "Distribution of Proteins with Fucosylated Glycans",
      subtitle = "All fucosylated glycosites included",
      x = "Category",
      y = "Number of Proteins"
    ) +
    geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")), 
              vjust = -0.5, size = 4, fontface = "bold")
  
  # Save the plot
  ggsave(file.path(output_dir, "fucose_protein_distribution.png"), 
         count_plot, 
         width = 10, 
         height = 8, 
         dpi = 300,
         bg = "white")
  
  # Create Venn diagram-like visualization
  venn_data <- data.frame(
    x = c(1, 2, 1.5),
    y = c(1, 1, 0.5),
    label = c("Healthy\nOnly", "ME/CFS\nOnly", "Common"),
    count = c(length(healthy_only), length(mecfs_only), length(common_proteins)),
    stringsAsFactors = FALSE
  )
  
  venn_plot <- ggplot(venn_data, aes(x = x, y = y, label = paste0(label, "\n(", count, ")"))) +
    geom_point(size = 15, alpha = 0.3, color = "blue") +
    geom_text(size = 4, fontface = "bold") +
    theme_void() +
    labs(title = "Protein Distribution with Fucosylated Glycans") +
    xlim(0.5, 2.5) +
    ylim(0, 1.5)
  
  # Save Venn plot
  ggsave(file.path(output_dir, "fucose_protein_venn.png"), 
         venn_plot, 
         width = 8, 
         height = 6, 
         dpi = 300,
         bg = "white")
  
  # Create detailed analysis of glycosites per protein
  glycosite_analysis <- fucose_long %>%
    group_by(protein_id, disease_status) %>%
    summarise(
      glycosite_count = n(),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = disease_status,
      values_from = glycosite_count,
      names_prefix = ""
    )
  
  # Add category information
  glycosite_analysis <- glycosite_analysis %>%
    mutate(
      category = case_when(
        protein_id %in% common_proteins ~ "Common",
        protein_id %in% healthy_only ~ "Healthy Only",
        protein_id %in% mecfs_only ~ "ME/CFS Only",
        TRUE ~ "Unknown"
      )
    )
  
  # Save glycosite analysis
  write.csv(glycosite_analysis, file.path(output_dir, "glycosite_analysis_by_protein.csv"), row.names = FALSE)
  
  # Print summary
  cat("\nFucose Protein Distribution Analysis:\n")
  cat("=====================================\n")
  cat("Total proteins with fucosylated glycans (Healthy):", length(healthy_proteins), "\n")
  cat("Total proteins with fucosylated glycans (ME/CFS):", length(mecfs_proteins), "\n")
  cat("Proteins common to both groups:", length(common_proteins), "\n")
  cat("Proteins unique to Healthy:", length(healthy_only), "\n")
  cat("Proteins unique to ME/CFS:", length(mecfs_only), "\n")
  cat("Total unique proteins:", length(healthy_only) + length(mecfs_only), "\n")
  
  # Create ranked protein lists based on glycosite count only
  protein_ranking <- fucose_long %>%
    group_by(protein_id, disease_status) %>%
    summarise(
      glycosite_count = n(),
      .groups = "drop"
    ) %>%
    mutate(
      category = case_when(
        protein_id %in% common_proteins ~ "Common",
        protein_id %in% healthy_only ~ "Healthy Only",
        protein_id %in% mecfs_only ~ "ME/CFS Only",
        TRUE ~ "Unknown"
      )
    )
  
  # Rank proteins by glycosite count (descending)
  top_common_ranked <- protein_ranking %>%
    filter(category == "Common") %>%
    arrange(desc(glycosite_count)) %>%
    head(10)
  
  top_healthy_ranked <- protein_ranking %>%
    filter(category == "Healthy Only") %>%
    arrange(desc(glycosite_count)) %>%
    head(10)
  
  top_mecfs_ranked <- protein_ranking %>%
    filter(category == "ME/CFS Only") %>%
    arrange(desc(glycosite_count)) %>%
    head(10)
  
  cat("\nTop 10 Common Proteins (ranked by number of fucosylated glycosites):\n")
  if (nrow(top_common_ranked) > 0) {
    for (i in 1:nrow(top_common_ranked)) {
      cat(i, ".", top_common_ranked$protein_id[i], 
          " (Glycosites:", top_common_ranked$glycosite_count[i], ")\n")
    }
  }
  
  cat("\nTop 10 Healthy-Only Proteins (ranked by number of fucosylated glycosites):\n")
  if (nrow(top_healthy_ranked) > 0) {
    for (i in 1:nrow(top_healthy_ranked)) {
      cat(i, ".", top_healthy_ranked$protein_id[i], 
          " (Glycosites:", top_healthy_ranked$glycosite_count[i], ")\n")
    }
  }
  
  cat("\nTop 10 ME/CFS-Only Proteins (ranked by number of fucosylated glycosites):\n")
  if (nrow(top_mecfs_ranked) > 0) {
    for (i in 1:nrow(top_mecfs_ranked)) {
      cat(i, ".", top_mecfs_ranked$protein_id[i], 
          " (Glycosites:", top_mecfs_ranked$glycosite_count[i], ")\n")
    }
  }
  
  # Save ranked protein lists
  write.csv(top_common_ranked, file.path(output_dir, "top_10_common_proteins_ranked.csv"), row.names = FALSE)
  write.csv(top_healthy_ranked, file.path(output_dir, "top_10_healthy_only_proteins_ranked.csv"), row.names = FALSE)
  write.csv(top_mecfs_ranked, file.path(output_dir, "top_10_mecfs_only_proteins_ranked.csv"), row.names = FALSE)
  
  cat("\nResults saved to:", output_dir, "\n")
  
  return(list(
    healthy_proteins = healthy_proteins,
    mecfs_proteins = mecfs_proteins,
    common_proteins = common_proteins,
    healthy_only = healthy_only,
    mecfs_only = mecfs_only
  ))
}

# Run the fucose protein distribution analysis
fucose_protein_analysis <- analyze_fucose_protein_distribution(
  fucose_data_file = "output_data/glycosite/RA/glycosite_gly_fuc_RA.csv",
  output_dir = "output_data/fucose_protein_analysis"
)


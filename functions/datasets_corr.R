#this script is used to correlate the datasets of the different methods
library("tidyr")
library("dplyr")
library("ggplot2")
library("readr")
library("correlation")
library(calibrate)
library(Correlplot)
library(ToolsForCoDa)
#load the data
setwd("~/Desktop/PD_GlyQuant/")
comp_map <- read.csv("input_data/Corr/struc-comp.csv")
glycomics <- read.csv("input_data/Corr/10_glycomics.csv")
pesudoglycomics_intensity <- read.csv("input_data/Corr/glycan_composition_peptidegroups_intensity_summary.csv")
pesudoglycomics_psm_count <- read.csv("input_data/Corr/relative_abundance_by_composition_PSM-count_RA.csv")


#calculate the means and standard deviations across all sample in the datasets

# Recalculate glycomics data by summing abundances for structures with same composition
glycomics_mapped <- glycomics %>%
  inner_join(comp_map, by = c("Structures" = "glycan.structures")) %>%
  # First sum abundances for each composition within each sample
  group_by(glycan.comps) %>%
  summarise(
    across(S1:S10, ~sum(.x, na.rm = TRUE)), # Sum values for each sample first
    .groups = 'drop'
  ) %>%
  # Then calculate mean and sd across samples
  mutate(
    mean_abundance = rowMeans(across(S1:S10), na.rm = TRUE),
    sd_abundance = apply(across(S1:S10), 1, sd, na.rm = TRUE),
    n_samples = rowSums(!is.na(across(S1:S10)))
  ) %>%
  # Rename the composition column and select only what we need
  rename(glycan_composition = glycan.comps) %>%
  dplyr::select(glycan_composition, mean_abundance, sd_abundance, n_samples)

# Debug: Check what columns are in glycomics_mapped
cat("\nColumns in glycomics_mapped:\n")
print(colnames(glycomics_mapped))
cat("\nFirst few rows of glycomics_mapped:\n")
print(head(glycomics_mapped))

# Print summary of compositions that were summed
composition_counts <- glycomics %>%
  inner_join(comp_map, by = c("Structures" = "glycan.structures")) %>%
  group_by(glycan.comps) %>%
  summarise(
    n_structures = n_distinct(Structures),
    structures = paste(Structures, collapse = ", "),
    .groups = 'drop'
  ) %>%
  filter(n_structures > 1)

cat("\nGlycan compositions with multiple structures:\n")
for(i in 1:nrow(composition_counts)) {
  cat(sprintf("%s: %d structures (%s)\n", 
              composition_counts$glycan.comps[i],
              composition_counts$n_structures[i],
              composition_counts$structures[i]))
}

# Join pseudoglycoproteomics intensity data with glycomics data
joined_data <- pesudoglycomics_intensity %>%
  mutate(glycan_composition = trimws(glycan_composition)) %>%  # Remove trailing spaces
  inner_join(glycomics_mapped, 
             by = "glycan_composition",
             suffix = c("_pseudo", "_glycomics")) %>%
  # Rename columns for clarity
  rename(
    mean_abundance_pseudo = mean_relative_abundance,
    sd_abundance_pseudo = se_relative_abundance
  )

# Print summary of the join
cat("\nJoined dataset summary:\n")
cat("Total compositions matched:", nrow(joined_data), "\n")
cat("Compositions in pseudoglycoproteomics:", nrow(pesudoglycomics_intensity), "\n")
cat("Compositions in glycomics:", nrow(glycomics_mapped), "\n\n")

# Print first few rows to verify join
print("First few rows of joined data:")
print(head(joined_data))

# Find compositions that didn't match between datasets
unmatched_pseudo <- pesudoglycomics_intensity %>%
  mutate(glycan_composition = trimws(glycan_composition)) %>%
  anti_join(glycomics_mapped, 
            by = "glycan_composition")

unmatched_glycomics <- glycomics_mapped %>%
  anti_join(pesudoglycomics_intensity %>% mutate(glycan_composition = trimws(glycan_composition)), 
            by = "glycan_composition")

# Print summary of unmatched compositions
cat("\nUnmatched compositions summary:\n")
cat("Compositions only in pseudoglycoproteomics:", nrow(unmatched_pseudo), "\n")
cat("Compositions only in glycomics:", nrow(unmatched_glycomics), "\n")

cat("\nCompositions only in pseudoglycoproteomics:\n")
if(nrow(unmatched_pseudo) > 0) {
  print(unmatched_pseudo %>% dplyr::select(glycan_composition, mean_relative_abundance) %>% arrange(desc(mean_relative_abundance)))
} else {
  print("None")
}

cat("\nCompositions only in glycomics:\n")
if(nrow(unmatched_glycomics) > 0) {
  print(unmatched_glycomics %>% dplyr::select(glycan_composition, mean_abundance) %>% arrange(desc(mean_abundance)))
} else {
  print("None")
}

#reshape data before processing 
# Keep only the essential columns for correlation analysis
joined_data_clean <- joined_data %>%
  dplyr::select(glycan_composition, mean_abundance_pseudo, mean_abundance) %>%
  rename(mean_abundance_glycomics = mean_abundance)

# Print summary of cleaned data
cat("\nCleaned dataset summary:\n")
cat("Total compositions:", nrow(joined_data_clean), "\n")
cat("Columns retained:", paste(colnames(joined_data_clean), collapse=", "), "\n")

#Transform the datasets from compositional data to centre log ratio transformed data

# Calculate CLR transformation for pseudoglycoproteomics data
# First ensure no zeros (replace with small value if needed)
min_nonzero <- min(joined_data_clean$mean_abundance_pseudo[joined_data_clean$mean_abundance_pseudo > 0], na.rm = TRUE)
joined_data_clean <- joined_data_clean %>%
  mutate(mean_abundance_pseudo_nozero = ifelse(mean_abundance_pseudo == 0, min_nonzero/2, mean_abundance_pseudo))

# Calculate geometric mean
geom_mean_pseudo <- exp(mean(log(joined_data_clean$mean_abundance_pseudo_nozero)))

# Apply CLR transformation
joined_data_clean <- joined_data_clean %>%
  mutate(clr_pseudo = log(mean_abundance_pseudo_nozero/geom_mean_pseudo))

# Print summary of CLR transformation
cat("\nCLR Transformation Summary:\n")
cat("Geometric mean:", geom_mean_pseudo, "\n")
cat("CLR range:", paste(range(joined_data_clean$clr_pseudo, na.rm=TRUE), collapse=" to "), "\n")

# Calculate CLR transformation for glycomics data
# First ensure no zeros (replace with small value if needed)
min_nonzero_glycomics <- min(joined_data_clean$mean_abundance_glycomics[joined_data_clean$mean_abundance_glycomics > 0], na.rm = TRUE)
joined_data_clean <- joined_data_clean %>%
  mutate(mean_abundance_glycomics_nozero = ifelse(mean_abundance_glycomics == 0, min_nonzero_glycomics/2, mean_abundance_glycomics))

# Calculate geometric mean for glycomics data
geom_mean_glycomics <- exp(mean(log(joined_data_clean$mean_abundance_glycomics_nozero)))

# Apply CLR transformation
joined_data_clean <- joined_data_clean %>%
  mutate(clr_glycomics = log(mean_abundance_glycomics_nozero/geom_mean_glycomics))

# Print summary of glycomics CLR transformation
cat("\nGlycomics CLR Transformation Summary:\n")
cat("Geometric mean:", geom_mean_glycomics, "\n")
cat("CLR range:", paste(range(joined_data_clean$clr_glycomics, na.rm=TRUE), collapse=" to "), "\n")


#caluate the correlation between the matching varaible means between the datasets

# Calculate correlations between CLR-transformed datasets
correlations <- joined_data_clean %>%
  summarise(
    pearson_cor = cor(clr_pseudo, clr_glycomics, method = "pearson"),
    pearson_p = cor.test(clr_pseudo, clr_glycomics, method = "pearson")$p.value,
    spearman_cor = cor(clr_pseudo, clr_glycomics, method = "spearman"), 
    spearman_p = cor.test(clr_pseudo, clr_glycomics, method = "spearman")$p.value
  )

# Print correlation results
cat("\nCorrelation Results:\n")
cat("Pearson correlation:", round(correlations$pearson_cor, 3), 
    "(p =", format(correlations$pearson_p, scientific = TRUE, digits = 3), ")\n")
cat("Spearman correlation:", round(correlations$spearman_cor, 3),
    "(p =", format(correlations$spearman_p, scientific = TRUE, digits = 3), ")\n")

# Add correlations to joined data
joined_data_clean <- joined_data_clean %>%
  mutate(
    pearson_cor = correlations$pearson_cor,
    pearson_p = correlations$pearson_p,
    spearman_cor = correlations$spearman_cor,
    spearman_p = correlations$spearman_p
  )

#calulate the difference bween the datasets
joined_data_clean <- joined_data_clean %>%
  mutate(
    difference = clr_pseudo - clr_glycomics
  )

#print the summary of the difference
cat("\nDifference Summary:\n")
cat("Mean difference:", mean(joined_data_clean$difference, na.rm=TRUE), "\n")
cat("SD difference:", sd(joined_data_clean$difference, na.rm=TRUE), "\n")

# Create violin plot of differences
difference_plot <- ggplot(joined_data_clean, aes(x = "", y = difference)) +
  geom_violin(fill = "lightblue", alpha = 0.5) +
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.4) +
  theme_minimal() +
  labs(
    title = "Distribution of Differences Between CLR-Transformed Datasets",
    x = "",
    y = "Difference (CLR Pseudoglycomics - CLR Glycomics)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11)
  )

# Save the difference plot
ggsave(
  "figures/clr_differences_violin.png",
  difference_plot,
  width = 8,
  height = 6,
  dpi = 300
)

library("ggpp")
library("ggpmisc")

# Create scatter plot comparing CLR-transformed values
# Calculate mean and SD of differences to identify outliers
joined_data_clean <- joined_data_clean %>%
  mutate(
    diff = clr_pseudo - clr_glycomics,
    z_score = abs(scale(diff))
  )

# Label points with z-score > 2 as outliers
outliers <- joined_data_clean %>%
  filter(z_score > 2) %>%
  pull(glycan_composition)

clr_comparison_plot <- ggplot(joined_data_clean, aes(x = clr_glycomics, y = clr_pseudo)) +
  # Regular points in gray, outliers in red
  geom_point(data = filter(joined_data_clean, !(glycan_composition %in% outliers)), 
            size = 3, alpha = 0.6, color = "black") +
  geom_point(data = filter(joined_data_clean, glycan_composition %in% outliers),
             size = 3, alpha = 0.6, color = "red") +
  geom_smooth(method = "lm", formula = y ~ x + 0, se = TRUE, color = "blue", alpha = 0.2) +
  stat_poly_eq(aes(label = after_stat(eq.label)), 
               formula = y ~ x + 0, 
               parse = TRUE,
               label.x = 0.05, 
               label.y = 0.95,
               size = 4,
               color = "blue") +
  geom_text_repel(
    data = filter(joined_data_clean, glycan_composition %in% outliers),
    aes(label = glycan_composition),
    size = 3,
    max.overlaps = 20,
    box.padding = 0.5,
    color = "red"
  ) +
  theme_minimal() +
  labs(
    title = "Comparison of CLR-Transformed Glycan Compositions",
    x = "CLR-Transformed Glycomics",
    y = "CLR-Transformed Pseudoglycomics",
    caption = paste0(
      "Pearson r = ", round(correlations$pearson_cor, 3),
      "\nOutliers labeled (z-score > 2)"
    )
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )

# Save the comparison plot
ggsave(
  "figures/clr_comparison_scatter.png",
  clr_comparison_plot,
  width = 8,
  height = 6,
  dpi = 300
)



# Calculate cosine similarity between CLR-transformed datasets
cosine_sim <- sum(joined_data_clean$clr_pseudo * joined_data_clean$clr_glycomics) /
  (sqrt(sum(joined_data_clean$clr_pseudo^2)) * sqrt(sum(joined_data_clean$clr_glycomics^2)))

# Print cosine similarity
cat("\nCosine Similarity:\n")
cat("Similarity score:", round(cosine_sim, 3), "\n")


# Calculate mean and standard error for each glycan composition in pseudoglycomics
pesudoglycomics_psm_count_summary <- pesudoglycomics_psm_count %>%
  group_by(glycan_composition) %>%
  summarise(
    mean_psm_count_RA = mean(relative_percentage, na.rm = TRUE),
    sd_psm_count_RA = sd(relative_percentage, na.rm = TRUE),
    .groups = 'drop'
  )

# Debug: Check what glycan compositions are in each dataset
cat("\nDebug: Glycan compositions in joined_data_clean:\n")
print(head(joined_data_clean$glycan_composition, 5))
cat("\nDebug: Glycan compositions in pesudoglycomics_psm_count_summary:\n")
print(head(pesudoglycomics_psm_count_summary$glycan_composition, 5))

# Check for exact matches
common_compositions <- intersect(joined_data_clean$glycan_composition, pesudoglycomics_psm_count_summary$glycan_composition)
cat("\nCommon compositions:", length(common_compositions), "\n")
if(length(common_compositions) > 0) {
  cat("First few common compositions:", paste(head(common_compositions, 3), collapse=", "), "\n")
}

# Join the PSM count summary with the CLR-transformed data
# First normalize the case for matching
joined_data_clean_normalized <- joined_data_clean %>%
  mutate(glycan_composition_norm = tolower(glycan_composition))

pesudoglycomics_psm_count_summary_normalized <- pesudoglycomics_psm_count_summary %>%
  mutate(glycan_composition_norm = tolower(glycan_composition))

# First do the join and check what columns we have
joined_temp <- joined_data_clean_normalized %>%
  left_join(pesudoglycomics_psm_count_summary_normalized, by = "glycan_composition_norm")

# Debug: Check what columns are available
cat("\nColumns after join:\n")
print(colnames(joined_temp))

joined_data_complete <- joined_temp %>%
  # Remove any rows with NA values after joining
  filter(!is.na(mean_psm_count_RA)) %>%
  # Remove the normalized column and keep original glycan_composition
  dplyr::select(-glycan_composition_norm, -glycan_composition.y) %>%
  rename(glycan_composition = glycan_composition.x)

# Print summary of joined data
cat("\nJoined Data Summary:\n")
cat("Total glycan compositions:", nrow(joined_data_complete), "\n")
cat("Compositions with complete data:", sum(complete.cases(joined_data_complete)), "\n")

# CLR transform PSM count data
# First ensure no zeros (replace with small value if needed)
min_nonzero_psm <- min(joined_data_complete$mean_psm_count_RA[joined_data_complete$mean_psm_count_RA > 0], na.rm = TRUE)
joined_data_complete <- joined_data_complete %>%
  mutate(mean_psm_count_RA_nozero = ifelse(mean_psm_count_RA == 0, min_nonzero_psm/2, mean_psm_count_RA))

# Calculate geometric mean for PSM count data
geom_mean_psm <- exp(mean(log(joined_data_complete$mean_psm_count_RA_nozero)))

# Apply CLR transformation
joined_data_complete <- joined_data_complete %>%
  mutate(clr_psm = log(mean_psm_count_RA_nozero/geom_mean_psm))

# Print summary of PSM CLR transformation
cat("\nPSM Count CLR Transformation Summary:\n")
cat("Geometric mean:", geom_mean_psm, "\n")
cat("CLR range:", paste(range(joined_data_complete$clr_psm, na.rm=TRUE), collapse=" to "), "\n")

# Calculate correlation between CLR-transformed data
clr_correlation <- cor.test(joined_data_complete$clr_pseudo, 
                          joined_data_complete$clr_psm,
                          method = "pearson")

# Print correlation results
cat("\nCorrelation between CLR-transformed pseudoglycoproteomics (intensity vs PSM count):\n")
cat("Pearson correlation:", round(clr_correlation$estimate, 3), "\n")
cat("p-value:", format(clr_correlation$p.value, scientific = TRUE, digits = 3), "\n")

# Identify outliers for the intensity vs PSM count comparison
# Calculate residuals from the regression model
intensity_psm_model <- lm(clr_psm ~ clr_pseudo + 0, data = joined_data_complete)
joined_data_complete <- joined_data_complete %>%
  mutate(
    intensity_psm_residuals = residuals(intensity_psm_model),
    intensity_psm_z_score = abs(scale(intensity_psm_residuals))
  )

# Identify outliers for this specific comparison
intensity_psm_outliers <- joined_data_complete %>%
  filter(intensity_psm_z_score > 2) %>%
  pull(glycan_composition)

cat("\nOutliers in intensity vs PSM count comparison:", length(intensity_psm_outliers), "\n")
if(length(intensity_psm_outliers) > 0) {
  cat("Outlier compositions:", paste(intensity_psm_outliers, collapse = ", "), "\n")
}

# Create scatter plot comparing CLR-transformed values
clr_scatter <- ggplot(joined_data_complete, aes(x = clr_pseudo, y = clr_psm)) +
  # Regular points in black, outliers in red
  geom_point(data = filter(joined_data_complete, !(glycan_composition %in% intensity_psm_outliers)), 
             size = 3, alpha = 0.6, color = "black") +
  geom_point(data = filter(joined_data_complete, glycan_composition %in% intensity_psm_outliers),
             size = 3, alpha = 0.6, color = "red") +
  geom_smooth(method = "lm", formula = y ~ x + 0, se = TRUE, color = "blue", alpha = 0.2) +
  stat_poly_eq(aes(label = after_stat(eq.label)), 
               formula = y ~ x + 0, 
               parse = TRUE,
               label.x = 0.05, 
               label.y = 0.95,
               size = 4,
               color = "blue") +
  geom_text_repel(
    data = filter(joined_data_complete, glycan_composition %in% intensity_psm_outliers),
    aes(label = glycan_composition),
    size = 3,
    max.overlaps = 20,
    box.padding = 0.5,
    color = "red"
  ) +
  theme_minimal() +
  labs(
    title = "Comparison of CLR-Transformed Abundance Measurements",
    x = "CLR-Transformed Intensity-based Abundance",
    y = "CLR-Transformed PSM Count-based Abundance",
    caption = paste0(
      "Pearson r = ", round(clr_correlation$estimate, 3),
      "\nOutliers labeled (z-score > 2)"
    )
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )

ggsave("figures/clr_scatter_psm_count.png", clr_scatter, width = 8, height = 6, dpi = 300)



#plot the correlations between the datasets

#save the data and plots

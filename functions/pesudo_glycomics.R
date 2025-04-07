##pseudo glycomics 

pseudo_glycomics <- glycoPSMs %>%
  group_by(disease_status) %>%
  # First calculate total intensity per disease status
  mutate(total_intensity = sum(intensity, na.rm = TRUE)) %>%
  # Then group by both disease status and glycan composition
  group_by(glycan_composition, disease_status) %>%
  summarise(
    glycan_intensity = sum(intensity, na.rm = TRUE),
    total_group_intensity = first(total_intensity),
    relative_percentage = (glycan_intensity / total_group_intensity) * 100,
    .groups = 'drop'
  ) %>%
  arrange(disease_status, desc(relative_percentage))

# Print summary
cat("\nRelative Glycan Composition Percentages by Disease Status:\n")
print(pseudo_glycomics)

# Save results
write.csv(pseudo_glycomics, "output_data/glycan_composition_percentages.csv", row.names = FALSE)

# Create summary statistics
glycan_summary <- pseudo_glycomics %>%
  group_by(disease_status) %>%
  summarise(
    total_glycans = n(),
    mean_percentage = mean(relative_percentage),
    median_percentage = median(relative_percentage),
    max_percentage = max(relative_percentage),
    min_percentage = min(relative_percentage)
  )

# Print summary statistics
cat("\nSummary Statistics:\n")
print(glycan_summary)

# Create detailed summary statistics for each glycan composition by disease status
glycan_detailed_summary <- glycoPSMs %>%
  group_by(glycan_composition, disease_status) %>%
  summarise(
    n_observations = n(),
    mean_intensity = mean(intensity, na.rm = TRUE),
    median_intensity = median(intensity, na.rm = TRUE),
    sd_intensity = sd(intensity, na.rm = TRUE),
    cv = (sd_intensity / mean_intensity) * 100,  # Coefficient of variation
    min_intensity = min(intensity, na.rm = TRUE),
    max_intensity = max(intensity, na.rm = TRUE),
    q25 = quantile(intensity, 0.25, na.rm = TRUE),
    q75 = quantile(intensity, 0.75, na.rm = TRUE),
    relative_percentage = (sum(intensity, na.rm = TRUE) / 
                         sum(glycoPSMs$intensity[glycoPSMs$disease_status == cur_group()$disease_status], 
                             na.rm = TRUE)) * 100,
    missing_values = sum(is.na(intensity)),
    .groups = 'drop'
  ) %>%
  arrange(glycan_composition, disease_status)

# Add difference between disease states
glycan_comparison <- glycan_detailed_summary %>%
  dplyr::select(glycan_composition, disease_status, relative_percentage) %>%
  tidyr::pivot_wider(names_from = disease_status, 
              values_from = relative_percentage) %>%
  dplyr::mutate(percentage_difference = MECFS - Healthy) %>%
  tidyr::pivot_longer(cols = c(Healthy, MECFS, percentage_difference),
               names_to = "metric",
               values_to = "value")

# Combine with detailed summary
glycan_detailed_summary <- glycan_detailed_summary %>%
  dplyr::left_join(
    glycan_comparison %>% 
      dplyr::filter(metric == "percentage_difference") %>%
      dplyr::select(glycan_composition, value),
    by = "glycan_composition"
  ) %>%
  dplyr::rename(percentage_diff_from_healthy = value)

# Round numeric columns to 3 decimal places
glycan_detailed_summary <- glycan_detailed_summary %>%
  dplyr::mutate(across(where(is.numeric), ~round(., 3)))

# Save detailed summary to CSV
write.csv(glycan_detailed_summary, 
          "output_data/glycan_composition_detailed_summary.csv", 
          row.names = FALSE)

# Create a simplified summary for quick review
glycan_simple_summary <- glycan_detailed_summary %>%
  dplyr::select(glycan_composition, 
         disease_status, 
         n_observations, 
         mean_intensity, 
         relative_percentage, 
         cv,
         percentage_diff_from_healthy) %>%
  arrange(desc(abs(percentage_diff_from_healthy)))

# Save simplified summary
write.csv(glycan_simple_summary, 
          "output_data/glycan_composition_simple_summary.csv", 
          row.names = FALSE)

# Print summary of the most different glycans
cat("\nTop 10 glycans with largest difference between disease states:\n")
print(glycan_simple_summary %>% 
        dplyr::filter(disease_status == "MECFS") %>% 
        dplyr::select(glycan_composition, percentage_diff_from_healthy) %>% 
        dplyr::arrange(desc(abs(percentage_diff_from_healthy))) %>% 
        head(10))

# Print overall summary
cat("\nSummary of analysis:\n")
cat("Total unique glycan compositions:", length(unique(glycoPSMs$glycan_composition)), "\n")
cat("Number of disease states:", length(unique(glycoPSMs$disease_status)), "\n")
cat("Total observations:", nrow(glycoPSMs), "\n")

# Create bar plot
glycan_plot <- pseudo_glycomics %>%
  ggplot(aes(x = reorder(glycan_composition, relative_percentage), 
             y = relative_percentage, 
             fill = disease_status)) +
  geom_bar(stat = "identity", 
           position = "dodge", 
           width = 0.8) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, 
                              hjust = 1, 
                              size = 8),
    plot.title = element_text(hjust = 0.5, 
                             size = 14),
    legend.position = "top",
    legend.title = element_text(size = 10),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Glycan Composition Relative Percentages by Disease Status",
    x = "Glycan Composition",
    y = "Relative Percentage (%)",
    fill = "Disease Status"
  ) +
  scale_fill_brewer(palette = "Set2")

# Save the plot
ggsave("figures/glycan_composition_comparison.pdf", 
       glycan_plot, 
       width = 12, 
       height = 8, 
       dpi = 300)

# Save the complete plot as both PDF and PNG
ggsave("figures/glycan_composition_comparison.png", 
       glycan_plot, 
       width = 12, 
       height = 8, 
       dpi = 300)

# Create a focused plot showing only top differential glycans
top_glycans_plot <- pseudo_glycomics %>%
  # Get top 15 glycans with biggest differences between groups
  group_by(glycan_composition) %>%
  summarise(diff = abs(diff(relative_percentage))) %>%
  arrange(desc(diff)) %>%
  head(15) %>%
  pull(glycan_composition) %>%
  # Filter original data for these glycans
  {filter(pseudo_glycomics, glycan_composition %in% .)} %>%
  ggplot(aes(x = reorder(glycan_composition, relative_percentage), 
             y = relative_percentage, 
             fill = disease_status)) +
  geom_bar(stat = "identity", 
           position = "dodge", 
           width = 0.8) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, 
                              hjust = 1, 
                              size = 10),
    plot.title = element_text(hjust = 0.5, 
                             size = 14),
    legend.position = "top",
    legend.title = element_text(size = 10)
  ) +
  labs(
    title = "Top 15 Differential Glycan Compositions",
    x = "Glycan Composition",
    y = "Relative Percentage (%)",
    fill = "Disease Status"
  ) +
  scale_fill_brewer(palette = "Set2")

# Save the focused plot
ggsave("figures/top_differential_glycans.pdf", 
       top_glycans_plot, 
       width = 10, 
       height = 6, 
       dpi = 300)

# Save the focused plot as both PDF and PNG
ggsave("figures/top_differential_glycans.png", 
       top_glycans_plot, 
       width = 10, 
       height = 6, 
       dpi = 300)

########################################################
# Calculate fucosylation percentages by sample
fucosylation_by_sample <- glycoPSMs %>%
  # Group by sample and disease status to get total intensity per sample
  group_by(sample, disease_status) %>%
  mutate(total_sample_intensity = sum(intensity, na.rm = TRUE)) %>%
  # Calculate Fuc percentage per sample using existing contains_Fuc column
  summarise(
    fuc_intensity = sum(intensity[contains_Fuc == TRUE], na.rm = TRUE),
    total_intensity = first(total_sample_intensity),
    fuc_percentage = (fuc_intensity / total_intensity) * 100,
    .groups = 'drop'
  )

# Perform t-test
t_test_result <- t.test(fuc_percentage ~ disease_status, data = fucosylation_by_sample)
t_test_label <- sprintf("p = %.3g", t_test_result$p.value)

# Create boxplot with individual points
fuc_boxplot <- ggplot(fucosylation_by_sample, 
                     aes(x = disease_status, 
                         y = fuc_percentage, 
                         fill = disease_status)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  scale_fill_viridis_d() +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  labs(
    title = "Fucosylation Percentage by Disease Status",
    x = "Disease Status",
    y = "Fucosylation Percentage (%)"
  ) +
  # Add t-test p-value
  annotate("text", 
           x = 1.5, 
           y = max(fucosylation_by_sample$fuc_percentage) * 1.05,
           label = t_test_label)

# Save the plot
ggsave("figures/fucosylation_comparison.pdf", 
       fuc_boxplot, 
       width = 8, 
       height = 6, 
       dpi = 300)

ggsave("figures/fucosylation_comparison.png", 
       fuc_boxplot, 
       width = 8, 
       height = 6, 
       dpi = 300)

# Print summary statistics
cat("\nFucosylation Summary Statistics:\n")
fucosylation_summary <- fucosylation_by_sample %>%
  group_by(disease_status) %>%
  summarise(
    mean_fuc = mean(fuc_percentage),
    sd_fuc = sd(fuc_percentage),
    median_fuc = median(fuc_percentage),
    n_samples = n()
  )
print(fucosylation_summary)

cat("\nT-test Results:\n")
print(t_test_result)


########################################################
# Calculate relative abundance of glycan classes by sample
glycan_class_by_sample <- glycoPSMs %>%
  # Group by sample, disease status, and glycan_class to get intensities
  group_by(sample, disease_status, glycan_class) %>%
  summarise(
    group_intensity = sum(intensity, na.rm = TRUE),
    .groups = 'keep'
  ) %>%
  # Calculate total intensity per sample for percentage
  group_by(sample, disease_status) %>%
  mutate(
    total_sample_intensity = sum(group_intensity),
    relative_percentage = (group_intensity / total_sample_intensity) * 100
  ) %>%
  ungroup()

# Check if the data frame was created
print("First few rows of glycan_class_by_sample:")
print(head(glycan_class_by_sample))

# Continue with the rest of your code only if the above works
if(exists("glycan_class_by_sample")) {
  # Perform t-test for each glycan class category
  t_test_results <- glycan_class_by_sample %>%
    group_by(glycan_class) %>%
    summarise(
      p_value = t.test(relative_percentage ~ disease_status)$p.value,
      .groups = 'drop'
    )
  
  # Create boxplot with individual points
  glycan_class_boxplot <- ggplot(glycan_class_by_sample, 
                       aes(x = glycan_class, 
                           y = relative_percentage, 
                           fill = disease_status)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
    scale_fill_viridis_d() +
    theme_minimal() +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    labs(
      title = "Relative Abundance of Glycan Classes by Disease Status",
      x = "Glycan Class",
      y = "Relative Abundance (%)",
      fill = "Disease Status"
    ) +
    # Add t-test p-values
    geom_text(data = t_test_results,
              aes(x = glycan_class,
                  y = max(glycan_class_by_sample$relative_percentage) + 5,
                  label = sprintf("p = %.3g", p_value)),
              inherit.aes = FALSE)
  
  # Save the plot
  ggsave("figures/glycan_class_distribution_comparison.pdf", 
         glycan_class_boxplot, 
         width = 10, 
         height = 7, 
         dpi = 300)
  
  ggsave("figures/glycan_class_distribution_comparison.png", 
         glycan_class_boxplot, 
         width = 10, 
         height = 7, 
         dpi = 300)
  
  # Print summary statistics
  cat("\nGlycan Class Distribution Summary Statistics:\n")
  glycan_class_summary <- glycan_class_by_sample %>%
    group_by(disease_status, glycan_class) %>%
    summarise(
      mean_percentage = mean(relative_percentage),
      sd_percentage = sd(relative_percentage),
      median_percentage = median(relative_percentage),
      n_samples = n(),
      .groups = 'drop'
    )
  print(glycan_class_summary)
  
  cat("\nT-test Results for each Glycan Class:\n")
  print(t_test_results)
}

# Create separate boxplots for each glycan class
# Complex glycans
complex_boxplot <- glycan_class_by_sample %>%
  filter(glycan_class == "Complex") %>%
  {
    # Calculate t-test p-value
    t_test <- t.test(relative_percentage ~ disease_status, data = .)
    p_val <- sprintf("p = %.3g", t_test$p.value)
    
    ggplot(., aes(x = disease_status, 
               y = relative_percentage, 
               fill = disease_status)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
      scale_fill_viridis_d() +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
      ) +
      labs(
        title = "Complex Glycans",
        x = "Disease Status",
        y = "Relative Abundance (%)"
      ) +
      # Add t-test p-value
      annotate("text", 
               x = 1.5, 
               y = max(.$relative_percentage) * 1.05,
               label = p_val)
  }

# Hybrid glycans
hybrid_boxplot <- glycan_class_by_sample %>%
  filter(glycan_class == "Hybrid") %>%
  {
    # Calculate t-test p-value
    t_test <- t.test(relative_percentage ~ disease_status, data = .)
    p_val <- sprintf("p = %.3g", t_test$p.value)
    
    ggplot(., aes(x = disease_status, 
               y = relative_percentage, 
               fill = disease_status)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
      scale_fill_viridis_d() +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
      ) +
      labs(
        title = "Hybrid Glycans",
        x = "Disease Status",
        y = "Relative Abundance (%)"
      ) +
      # Add t-test p-value
      annotate("text", 
               x = 1.5, 
               y = max(.$relative_percentage) * 1.05,
               label = p_val)
  }

# Oligomannose glycans
oligomannose_boxplot <- glycan_class_by_sample %>%
  filter(glycan_class == "Oligomannose") %>%
  {
    # Calculate t-test p-value
    t_test <- t.test(relative_percentage ~ disease_status, data = .)
    p_val <- sprintf("p = %.3g", t_test$p.value)
    
    ggplot(., aes(x = disease_status, 
               y = relative_percentage, 
               fill = disease_status)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
      scale_fill_viridis_d() +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
      ) +
      labs(
        title = "Oligomannose Glycans",
        x = "Disease Status",
        y = "Relative Abundance (%)"
      ) +
      # Add t-test p-value
      annotate("text", 
               x = 1.5, 
               y = max(.$relative_percentage) * 1.05,
               label = p_val)
  }

# Rest of the code remains the same
library(patchwork)
combined_plot <- complex_boxplot + hybrid_boxplot + oligomannose_boxplot +
  plot_layout(ncol = 3)

# Save individual plots
ggsave("figures/complex_glycans_boxplot.pdf", complex_boxplot, width = 6, height = 5, dpi = 300)
ggsave("figures/complex_glycans_boxplot.png", complex_boxplot, width = 6, height = 5, dpi = 300)

ggsave("figures/hybrid_glycans_boxplot.pdf", hybrid_boxplot, width = 6, height = 5, dpi = 300)
ggsave("figures/hybrid_glycans_boxplot.png", hybrid_boxplot, width = 6, height = 5, dpi = 300)

ggsave("figures/oligomannose_glycans_boxplot.pdf", oligomannose_boxplot, width = 6, height = 5, dpi = 300)
ggsave("figures/oligomannose_glycans_boxplot.png", oligomannose_boxplot, width = 6, height = 5, dpi = 300)

# Save combined plot
ggsave("figures/all_glycan_classes_boxplots.pdf", combined_plot, width = 18, height = 5, dpi = 300)
ggsave("figures/all_glycan_classes_boxplots.png", combined_plot, width = 18, height = 5, dpi = 300)

# Print summary statistics for each class
for(class in unique(glycan_class_by_sample$glycan_class)) {
  cat(sprintf("\nSummary Statistics for %s Glycans:\n", class))
  print(glycan_class_by_sample %>%
    filter(glycan_class == class) %>%
    group_by(disease_status) %>%
    summarise(
      mean_percentage = mean(relative_percentage),
      sd_percentage = sd(relative_percentage),
      median_percentage = median(relative_percentage),
      n_samples = n(),
      .groups = 'drop'
    ))
  
  # Print t-test results
  t_test <- t.test(relative_percentage ~ disease_status, 
                   data = filter(glycan_class_by_sample, glycan_class == class))
  cat(sprintf("\nT-test Results for %s Glycans:\n", class))
  print(t_test)
}

########################################################
# Calculate NeuAc percentages by sample
neuac_by_sample <- glycoPSMs %>%
  # Group by sample and disease status to get total intensity per sample
  group_by(sample, disease_status) %>%
  mutate(total_sample_intensity = sum(intensity, na.rm = TRUE)) %>%
  # Calculate NeuAc percentage per sample using existing contains_NeuAc column
  summarise(
    neuac_intensity = sum(intensity[contains_NeuAc == TRUE], na.rm = TRUE),
    total_intensity = first(total_sample_intensity),
    neuac_percentage = (neuac_intensity / total_intensity) * 100,
    .groups = 'drop'
  )

# Perform t-test
t_test_result_neuac <- t.test(neuac_percentage ~ disease_status, data = neuac_by_sample)
t_test_label_neuac <- sprintf("p = %.3g", t_test_result_neuac$p.value)

# Create boxplot with individual points
neuac_boxplot <- ggplot(neuac_by_sample, 
                       aes(x = disease_status, 
                           y = neuac_percentage, 
                           fill = disease_status)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  scale_fill_viridis_d() +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  labs(
    title = "Sialylation (NeuAc) Percentage by Disease Status",
    x = "Disease Status",
    y = "Sialylation Percentage (%)"
  ) +
  # Add t-test p-value
  annotate("text", 
           x = 1.5, 
           y = max(neuac_by_sample$neuac_percentage) * 1.05,
           label = t_test_label_neuac)

# Save the plot
ggsave("figures/sialylation_comparison.pdf", 
       neuac_boxplot, 
       width = 8, 
       height = 6, 
       dpi = 300)

ggsave("figures/sialylation_comparison.png", 
       neuac_boxplot, 
       width = 8, 
       height = 6, 
       dpi = 300)

# Print summary statistics
cat("\nSialylation Summary Statistics:\n")
neuac_summary <- neuac_by_sample %>%
  group_by(disease_status) %>%
  summarise(
    mean_neuac = mean(neuac_percentage),
    sd_neuac = sd(neuac_percentage),
    median_neuac = median(neuac_percentage),
    n_samples = n()
  )
print(neuac_summary)

cat("\nT-test Results:\n")
print(t_test_result_neuac)

########################################################
# Calculate relative abundance of sialic acid counts by sample
sia_by_sample <- glycoPSMs %>%
  # Group by sample, disease status, and sia_count to get intensities
  group_by(sample, disease_status, sia_count) %>%
  summarise(
    group_intensity = sum(intensity, na.rm = TRUE),
    .groups = 'keep'
  ) %>%
  # Calculate total intensity per sample for percentage
  group_by(sample, disease_status) %>%
  mutate(
    total_sample_intensity = sum(group_intensity),
    relative_percentage = (group_intensity / total_sample_intensity) * 100
  ) %>%
  ungroup()

# Perform t-test for each sia_count category
t_test_results <- sia_by_sample %>%
  group_by(sia_count) %>%
  summarise(
    p_value = t.test(relative_percentage ~ disease_status)$p.value,
    .groups = 'drop'
  )

# Create boxplot with individual points
sia_boxplot <- ggplot(sia_by_sample, 
                     aes(x = factor(sia_count), 
                         y = relative_percentage, 
                         fill = disease_status)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  scale_fill_viridis_d() +
  theme_minimal() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  labs(
    title = "Relative Abundance of Sialic Acid Counts by Disease Status",
    x = "Number of Sialic Acids (NeuAc)",
    y = "Relative Abundance (%)",
    fill = "Disease Status"
  ) +
  # Add t-test p-values
  geom_text(data = t_test_results,
            aes(x = factor(sia_count),
                y = max(sia_by_sample$relative_percentage) + 5,
                label = sprintf("p = %.3g", p_value)),
            inherit.aes = FALSE)

# Save the plot
ggsave("figures/sialic_acid_distribution_comparison.pdf", 
       sia_boxplot, 
       width = 10, 
       height = 7, 
       dpi = 300)

ggsave("figures/sialic_acid_distribution_comparison.png", 
       sia_boxplot, 
       width = 10, 
       height = 7, 
       dpi = 300)

# Print summary statistics
cat("\nSialic Acid Distribution Summary Statistics:\n")
sia_summary <- sia_by_sample %>%
  group_by(disease_status, sia_count) %>%
  summarise(
    mean_percentage = mean(relative_percentage),
    sd_percentage = sd(relative_percentage),
    median_percentage = median(relative_percentage),
    n_samples = n(),
    .groups = 'drop'
  )
print(sia_summary)

cat("\nT-test Results for each Sialic Acid Count:\n")
print(t_test_results)

# Calculate relative abundance of glycan classes by sample
glycan_class_by_sample <- glycoPSMs %>%
  # Group by sample, disease status, and glycan_class to get intensities
  group_by(sample, disease_status, glycan_class) %>%
  summarise(
    group_intensity = sum(intensity, na.rm = TRUE),
    .groups = 'keep'
  ) %>%
  # Calculate total intensity per sample for percentage
  group_by(sample, disease_status) %>%
  mutate(
    total_sample_intensity = sum(group_intensity),
    relative_percentage = (group_intensity / total_sample_intensity) * 100
  ) %>%
  ungroup()


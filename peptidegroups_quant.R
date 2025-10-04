library(tidyverse)
library(here)
library(janitor)
library(stringr)
library(readr)
library(dplyr)
library(ggplot2)
library(broom)
library(viridis)
library(patchwork)


# Set working directory to the project root
setwd("/Users/thomasreilly/Desktop/PD_GlyQuant")

StudyInformation <- read_tsv("input_data/10S_MECFS_GPEPS_250125_StudyInformation.txt")
PeptideGroups <- read_tsv("input_data/10S_MECFS_GPEPS_250125_PeptideGroups.txt")

PeptideGroups <- PeptideGroups %>%
  clean_names(case = "snake") %>%
    mutate(
        pep_glycosite = str_extract(modifications, "(?<=\\[N)\\d+(?=\\])"),
        protein_names = str_extract(master_protein_descriptions, "^[^O]+(?= OS=)"),
        gene_name = str_extract(master_protein_descriptions, "(?<=GN=)[^\\s]+(?= PE=)"),
        contains_Fuc = str_detect(glycan_composition, "Fuc"),
        contains_NeuAc = str_detect(glycan_composition, "NeuAc"),
        pep_glycosite = as.numeric(pep_glycosite),
        glycan_composition = str_remove(glycan_composition, "@ n \\| rare1$"),
        protein_glycosite = position_in_protein + pep_glycosite - 1
               )

glyco_peptide_groups <- PeptideGroups %>% filter(!is.na(glycan_composition))
glyco_peptide_groups <- glyco_peptide_groups %>% filter(pep_2d_by_search_engine_a2_pmi_byonic < 0.001)

StudyInformation <- StudyInformation %>%
    clean_names(case = "snake") 

glycan_class_map <- read_csv(file = "input_data/glycan_class_map.csv")

match_glycan_class <- function(input_df, reference_df) {
               # Ensure required columns exist
               if (!"glycan_composition" %in% colnames(input_df)) {
                 stop("Column 'glycan_composition' not found in input dataframe")
               }
               if (!all(c("glycans", "glycan_class") %in% colnames(reference_df))) {
                 stop("Columns 'glycans' and 'glycan_class' not found in reference dataframe")
               }
               
               # Create a copy of the input dataframe to avoid modifying the original
               result_df <- input_df
               
               # Clean glycan compositions by removing the "@ N | rare1 [NXX]" part
               cleaned_glycans <- str_replace(input_df$glycan_composition, " @ N \\| rare1( \\[N\\d+\\])?", "")
               
               # Convert to character and trim whitespace
               cleaned_glycans <- trimws(tolower(as.character(cleaned_glycans)))
               reference_df$glycans <- trimws(tolower(as.character(reference_df$glycans)))
               
               # Initialize an empty vector to store the matched glycan classes
               matched_classes <- vector("character", length = nrow(input_df))
               
              # Iterate over each glycan_composition in the input dataframe
               for (i in 1:nrow(input_df)) {
                 glycan_comp <- cleaned_glycans[i]
                 
                 # Find the row in the reference dataframe that matches the glycan composition
                 matched_row <- reference_df[reference_df$glycans == glycan_comp, ]
                 
                 # If a match is found, store the glycan_class, otherwise store NA
                 if (nrow(matched_row) > 0) {
                   matched_classes[i] <- matched_row$glycan_class[1]
                 } else {
                   matched_classes[i] <- NA
                 }
               }
               
               # Add the matched glycan_class to the input dataframe
               result_df$glycan_class <- matched_classes
               
               return(result_df)
             }

glyco_peptide_groups <- match_glycan_class(glyco_peptide_groups, glycan_class_map)

glyco_peptide_groups <- glyco_peptide_groups %>%
               mutate(gsite_ID = paste(protein_accessions, protein_glycosite, sep = "_"))
 

# Count unique glycopeptides
glycopeptides <- find_unique_values(glyco_peptide_groups, "peptide_groups_peptide_group_id")
glycopeptide_count <- length(glycopeptides)
cat("Number of unique glycopeptides:", glycopeptide_count, "\n")

# Find common glycopeptides across all samples
common_glycopeptides <- glyco_peptide_groups_long %>%
  group_by(peptide_groups_peptide_group_id) %>%
  summarise(
    num_samples = n_distinct(sample),
    .groups = 'drop'
  ) %>%
  filter(num_samples == total_samples) %>%
  pull(peptide_groups_peptide_group_id)

# Count common glycopeptides
common_glycopeptide_count <- length(common_glycopeptides)
cat("Number of glycopeptides common across all samples:", common_glycopeptide_count, "\n")
 

# Pivot longer by abundance columns
glyco_peptide_groups_long <- glyco_peptide_groups %>%
  select(-any_of(starts_with("abundances_")), -any_of(starts_with("abundance_ratio")), -any_of(starts_with("found_in_"))) %>%
  pivot_longer(
    cols = starts_with("abundance_f"),
    names_to = "sample",
    values_to = "abundance"
  ) %>%
  filter(!is.na(abundance))

 find_unique_values <- function(df, column_name) {
               # Check if the column exists
               if (!(column_name %in% names(df))) {
                 stop("Column not found in the data frame")
               }
               
               # Extract unique values
               unique_values <- unique(df[[column_name]])
               
               # Convert to list and return
               return(as.list(unique_values))
             }

gsites <- find_unique_values(glyco_peptide_groups_long, "gsite_ID")
# Count gsites
gsite_count <- length(gsites)
cat("Number of unique glycosites:", gsite_count, "\n")

# Count glycoprots
glycoprots <- find_unique_values(glyco_peptide_groups_long, "protein_accessions")
glycoprot_count <- length(glycoprots)
cat("Number of unique glycoproteins:", glycoprot_count, "\n")

# Find features present across all samples
all_samples <- unique(glyco_peptide_groups_long$sample)
total_samples <- length(all_samples)

# Find common protein_accessions
common_proteins <- glyco_peptide_groups_long %>%
  group_by(protein_accessions) %>%
  summarise(
    num_samples = n_distinct(sample),
    .groups = 'drop'
  ) %>%
  filter(num_samples == total_samples) %>%
  pull(protein_accessions)

# Find common glycosites
common_gsites <- glyco_peptide_groups_long %>%
  group_by(gsite_ID) %>%
  summarise(
    num_samples = n_distinct(sample),
    .groups = 'drop'
  ) %>%
  filter(num_samples == total_samples) %>%
  pull(gsite_ID)

# Print results
cat("Number of proteins common across all samples:", length(common_proteins), "\n")
cat("Number of glycosites common across all samples:", length(common_gsites), "\n")



#mutate the sample column to remove evverything before and including the 4th _ 
glyco_peptide_groups_long <- glyco_peptide_groups_long %>%
  mutate(sample = str_remove(sample, "^([^_]*_){4}")) %>%
  mutate(sample = str_remove(sample, "_[^_]*$"))

# Count unique gene names
gene_names <- find_unique_values(glyco_peptide_groups_long, "gene_name")
gene_count <- length(gene_names)
cat("Number of unique genes:", gene_count, "\n")

# Find gene names present across all samples
common_genes <- glyco_peptide_groups_long %>%
  group_by(gene_name) %>%
  summarise(
    num_samples = n_distinct(sample),
    .groups = 'drop'
  ) %>%
  filter(num_samples == total_samples) %>%
  pull(gene_name)

cat("Number of genes common across all samples:", length(common_genes), "\n")


##################################################################################################
# Overview  volcano plots
#glycoproteins not RA normailiesd 
# Sum abundance by gene name and sample
gene_abundance <- glyco_peptide_groups_long %>%
  group_by(gene_name, sample, group) %>%
  summarise(
    total_abundance = sum(abundance, na.rm = TRUE),
    .groups = 'drop'
  ) 

# Calculate statistics for volcano plot
volcano_stats <- gene_abundance %>%
  group_by(gene_name) %>%
  summarise(
    # Calculate mean abundance for each group
    mean_healthy = mean(total_abundance[group == "Healthy"], na.rm = TRUE),
    mean_mecfs = mean(total_abundance[group == "MECFS"], na.rm = TRUE),
    # Calculate log2 fold change
    log2FC = log2(mean_mecfs / mean_healthy),
    # Perform t-test with error handling
    pvalue = tryCatch({
      t.test(
        total_abundance[group == "MECFS"],
        total_abundance[group == "Healthy"]
      )$p.value
    }, error = function(e) NA_real_),
    .groups = 'drop'
  ) %>%
  # Calculate adjusted p-values (handle NA values)
  mutate(
    adj_pvalue = ifelse(is.na(pvalue), NA_real_, p.adjust(pvalue, method = "BH")),
    sig = !is.na(adj_pvalue) & adj_pvalue < 0.05,
    # Add labels for significant genes
    label = ifelse(sig & abs(log2FC) > 1, gene_name, "")
  )

# Load EnhancedVolcano package
library(EnhancedVolcano)
# Create enhanced volcano plot
volcano_plot <- EnhancedVolcano(
  volcano_stats,
  lab = volcano_stats$gene_name,
  x = 'log2FC',
  y = 'adj_pvalue',
  title = 'Differential glycoprotein aubundance in ME/CFS vs Healthy',
  subtitle = 'ME/CFS vs Healthy Controls',
  caption = 'FC cutoff: 2; p-value cutoff: 0.05',
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 3.0,
  labSize = 3.0,
  labCol = 'black',
  labFace = 'bold',
  boxedLabels = TRUE,
  parseLabels = FALSE,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  colAlpha = 0.75,
  colConnectors = 'grey50',
  gridlines.major = TRUE,
  gridlines.minor = FALSE,
  border = 'partial',
  borderWidth = 0.5,
  borderColour = 'black',
  xlim = c(min(volcano_stats$log2FC, na.rm = TRUE) - 0.5, 
           max(volcano_stats$log2FC, na.rm = TRUE) + 0.5),
  ylim = c(0, max(-log10(volcano_stats$adj_pvalue), na.rm = TRUE) + 1),
  xlab = bquote(~Log[2]~ 'fold change'),
  ylab = bquote(~-Log[10]~ 'adjusted p-value'),
  axisLabSize = 16,
  titleLabSize = 18,
  subtitleLabSize = 14,
  captionLabSize = 12,
  legendPosition = 'right',
  legendLabSize = 12,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  typeConnectors = 'closed',
  endsConnectors = 'first',
  lengthConnectors = unit(0.01, 'npc'),
  hline = c(10e-8, 10e-4, 10e-2),
  hlineCol = c('red', 'red', 'red'),
  hlineType = c('twodash', 'twodash', 'twodash'),
  hlineWidth = c(0.5, 0.5, 0.5),
  vline = c(-1, 1),
  vlineCol = c('red', 'red'),
  vlineType = c('twodash', 'twodash'),
  vlineWidth = c(0.5, 0.5),
  selectLab = volcano_stats$gene_name[volcano_stats$sig | abs(volcano_stats$log2FC) > 1],
  raster = FALSE,
  max.overlaps = 40 # Increased max.overlaps to show more labels
) 

# Save the plot
ggsave(
  "figures/peptidegroups_intensity/gene_volcano_plot.png",
  volcano_plot,
  width = 10,
  height = 8,
  dpi = 300
)

##################################################################################################
# Protein-Glycan Composition Volcano Plot
# Sum abundance by protein and glycan composition for each sample
protein_glycan_abundance <- glyco_peptide_groups_long %>%
  group_by(gene_name, glycan_composition, sample, group) %>%
  summarise(
    total_abundance = sum(abundance, na.rm = TRUE),
    .groups = 'drop'
  )

# Calculate statistics for protein-glycan volcano plot
protein_glycan_volcano_stats <- protein_glycan_abundance %>%
  group_by(gene_name, glycan_composition) %>%
  summarise(
    mean_mecfs = mean(total_abundance[group == "MECFS"], na.rm = TRUE),
    mean_healthy = mean(total_abundance[group == "Healthy"], na.rm = TRUE),
    pvalue = tryCatch({
      t.test(
        total_abundance[group == "MECFS"],
        total_abundance[group == "Healthy"]
      )$p.value
    }, error = function(e) NA_real_),
    .groups = 'drop'
  ) %>%
  mutate(
    log2FC = log2(mean_mecfs / mean_healthy),
    adj_pvalue = ifelse(is.na(pvalue), NA_real_, p.adjust(pvalue, method = "BH")),
    sig = !is.na(adj_pvalue) & adj_pvalue < 0.05,
    # Create combined labels for protein-glycan pairs
    protein_glycan_label = paste0(gene_name, ":", glycan_composition),
    label = ifelse(sig & abs(log2FC) > 1, protein_glycan_label, "")
  )

# Create enhanced volcano plot for protein-glycan combinations
protein_glycan_volcano_plot <- EnhancedVolcano(
  protein_glycan_volcano_stats,
  lab = protein_glycan_volcano_stats$protein_glycan_label,
  x = 'log2FC',
  y = 'adj_pvalue',
  title = 'Differential Protein-Glycan Abundance in ME/CFS vs Healthy',
  subtitle = 'ME/CFS vs Healthy Controls - Protein:Glycan combinations',
  caption = 'FC cutoff: 2; p-value cutoff: 0.05',
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.0,
  labSize = 2.5,
  labCol = 'black',
  labFace = 'bold',
  boxedLabels = TRUE,
  parseLabels = FALSE,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  colAlpha = 0.75,
  colConnectors = 'grey50',
  gridlines.major = TRUE,
  gridlines.minor = FALSE,
  border = 'partial',
  borderWidth = 0.5,
  borderColour = 'black',
  xlim = c(min(protein_glycan_volcano_stats$log2FC, na.rm = TRUE) - 0.5, 
           max(protein_glycan_volcano_stats$log2FC, na.rm = TRUE) + 0.5),
  ylim = c(0, max(-log10(protein_glycan_volcano_stats$adj_pvalue), na.rm = TRUE) + 1),
  xlab = bquote(~Log[2]~ 'fold change'),
  ylab = bquote(~-Log[10]~ 'adjusted p-value'),
  axisLabSize = 16,
  titleLabSize = 18,
  subtitleLabSize = 14,
  captionLabSize = 12,
  legendPosition = 'right',
  legendLabSize = 12,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  typeConnectors = 'closed',
  endsConnectors = 'first',
  lengthConnectors = unit(0.01, 'npc'),
  hline = c(10e-8, 10e-4, 10e-2),
  hlineCol = c('red', 'red', 'red'),
  hlineType = c('twodash', 'twodash', 'twodash'),
  hlineWidth = c(0.5, 0.5, 0.5),
  vline = c(-1, 1),
  vlineCol = c('red', 'red'),
  vlineType = c('twodash', 'twodash'),
  vlineWidth = c(0.5, 0.5),
  selectLab = protein_glycan_volcano_stats$protein_glycan_label[startsWith(protein_glycan_volcano_stats$protein_glycan_label, "FN1")],
  raster = FALSE,
  max.overlaps = 50 # Increased max.overlaps for protein-glycan combinations
)

# Save the protein-glycan volcano plot
ggsave(
  "figures/peptidegroups_intensity/protein_glycan_volcano_plot.png",
  protein_glycan_volcano_plot,
  width = 12,
  height = 10,
  dpi = 300
)

# Print summary statistics
cat("Protein-Glycan Volcano Plot Summary:\n")
cat("Total protein-glycan combinations:", nrow(protein_glycan_volcano_stats), "\n")
cat("Significant combinations (p < 0.05):", sum(protein_glycan_volcano_stats$sig, na.rm = TRUE), "\n")
cat("High fold change combinations (|log2FC| > 1):", sum(abs(protein_glycan_volcano_stats$log2FC) > 1, na.rm = TRUE), "\n")
cat("Significant AND high fold change:", sum(protein_glycan_volcano_stats$sig & abs(protein_glycan_volcano_stats$log2FC) > 1, na.rm = TRUE), "\n")
# Print top 10 highest absolute fold change data points
cat("\nTop 10 highest absolute fold change protein-glycan combinations:\n")
protein_glycan_volcano_stats %>%
  arrange(desc(abs(log2FC))) %>%
  head(10) %>%
  select(protein_glycan_label, log2FC, adj_pvalue) %>%
  print(n = 10)


##################################################################################################
# Protein-Glycosite-Glycan Composition Volcano Plot
# Sum abundance by protein, glycosite, and glycan composition for each sample
protein_glycosite_glycan_abundance <- glyco_peptide_groups_long %>%
  group_by(gene_name, protein_glycosite, glycan_composition, sample, group) %>%
  summarise(
    total_abundance = sum(abundance, na.rm = TRUE),
    .groups = 'drop'
  )

# Calculate statistics for protein-glycosite-glycan volcano plot
protein_glycosite_glycan_volcano_stats <- protein_glycosite_glycan_abundance %>%
  group_by(gene_name, protein_glycosite, glycan_composition) %>%
  summarise(
    mean_mecfs = mean(total_abundance[group == "MECFS"], na.rm = TRUE),
    mean_healthy = mean(total_abundance[group == "Healthy"], na.rm = TRUE),
    pvalue = tryCatch({
      t.test(
        total_abundance[group == "MECFS"],
        total_abundance[group == "Healthy"]
      )$p.value
    }, error = function(e) NA_real_),
    .groups = 'drop'
  ) %>%
  mutate(
    log2FC = log2(mean_mecfs / mean_healthy),
    adj_pvalue = ifelse(is.na(pvalue), NA_real_, p.adjust(pvalue, method = "BH")),
    sig = !is.na(adj_pvalue) & adj_pvalue < 0.05,
    # Create combined labels for protein-glycosite-glycan triplets
    protein_glycosite_glycan_label = paste0(gene_name, ":", protein_glycosite, ":", glycan_composition),
    label = ifelse(sig & abs(log2FC) > 1, protein_glycosite_glycan_label, "")
  )

# Debug: Check the top 10 data points and their labels
cat("Debugging top 10 data points:\n")
top_10_indices <- order(abs(protein_glycosite_glycan_volcano_stats$log2FC), decreasing = TRUE)[1:min(10, nrow(protein_glycosite_glycan_volcano_stats))]
top_10_data <- protein_glycosite_glycan_volcano_stats[top_10_indices, ]
cat("Top 10 by absolute fold change:\n")
print(top_10_data[, c("gene_name", "protein_glycosite", "glycan_composition", "log2FC", "adj_pvalue", "protein_glycosite_glycan_label")])
cat("Number of valid labels:", sum(!is.na(top_10_data$protein_glycosite_glycan_label) & top_10_data$protein_glycosite_glycan_label != ""), "\n")

# Create enhanced volcano plot for protein-glycosite-glycan combinations
protein_glycosite_glycan_volcano_plot <- EnhancedVolcano(
  protein_glycosite_glycan_volcano_stats,
  lab = protein_glycosite_glycan_volcano_stats$protein_glycosite_glycan_label,
  x = 'log2FC',
  y = 'adj_pvalue',
  title = 'Differential Protein-Glycosite-Glycan Abundance in ME/CFS vs Healthy',
  subtitle = 'ME/CFS vs Healthy Controls - Protein:Glycosite:Glycan combinations',
  caption = 'FC cutoff: 2; p-value cutoff: 0.05',
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 1.5,
  labSize = 4.0,
  labCol = 'black',
  labFace = 'bold',
  boxedLabels = TRUE,
  parseLabels = FALSE,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  colAlpha = 0.75,
  colConnectors = 'grey50',
  gridlines.major = TRUE,
  gridlines.minor = FALSE,
  border = 'partial',
  borderWidth = 0.5,
  borderColour = 'black',
  xlim = c(min(protein_glycosite_glycan_volcano_stats$log2FC, na.rm = TRUE) - 0.5, 
           max(protein_glycosite_glycan_volcano_stats$log2FC, na.rm = TRUE) + 0.5),
  ylim = c(0, max(-log10(protein_glycosite_glycan_volcano_stats$adj_pvalue), na.rm = TRUE) + 0.7),
  xlab = bquote(~Log[2]~ 'fold change'),
  ylab = bquote(~-Log[10]~ 'adjusted p-value'),
  axisLabSize = 16,
  titleLabSize = 18,
  subtitleLabSize = 14,
  captionLabSize = 12,
  legendPosition = 'right',
  legendLabSize = 12,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.3,
  typeConnectors = 'closed',
  endsConnectors = 'first',
  lengthConnectors = unit(0.005, 'npc'),
  hline = c(10e-8, 10e-4, 10e-2),
  hlineCol = c('red', 'red', 'red'),
  hlineType = c('twodash', 'twodash', 'twodash'),
  hlineWidth = c(0.5, 0.5, 0.5),
  vline = c(-1, 1),
  vlineCol = c('red', 'red'),
  vlineType = c('twodash', 'twodash'),
  vlineWidth = c(0.5, 0.5),
selectLab = protein_glycosite_glycan_volcano_stats$protein_glycosite_glycan_label[startsWith(protein_glycosite_glycan_volcano_stats$protein_glycosite_glycan_label, "FN1")],  raster = FALSE,
  max.overlaps = 100
)

# Save the protein-glycosite-glycan volcano plot
ggsave(
  "figures/peptidegroups_intensity/protein_glycosite_glycan_volcano_plot.png",
  protein_glycosite_glycan_volcano_plot,
  width = 14,
  height = 12,
  dpi = 600
)

# Print summary statistics for protein-glycosite-glycan combinations
cat("\nProtein-Glycosite-Glycan Volcano Plot Summary:\n")
cat("Total protein-glycosite-glycan combinations:", nrow(protein_glycosite_glycan_volcano_stats), "\n")
cat("Significant combinations (p < 0.05):", sum(protein_glycosite_glycan_volcano_stats$sig, na.rm = TRUE), "\n")
cat("High fold change combinations (|log2FC| > 1):", sum(abs(protein_glycosite_glycan_volcano_stats$log2FC) > 1, na.rm = TRUE), "\n")
cat("Significant AND high fold change:", sum(protein_glycosite_glycan_volcano_stats$sig & abs(protein_glycosite_glycan_volcano_stats$log2FC) > 1, na.rm = TRUE), "\n")

# Show top significant protein-glycosite-glycan combinations
top_significant <- protein_glycosite_glycan_volcano_stats %>%
  filter(sig & abs(log2FC) > 1) %>%
  arrange(desc(abs(log2FC))) %>%
  head(20)

if(nrow(top_significant) > 0) {
  cat("\nTop 20 most significant protein-glycosite-glycan combinations:\n")
  print(top_significant[, c("gene_name", "protein_glycosite", "glycan_composition", "log2FC", "adj_pvalue")])
}



##################################################################################################
# Calculate relative abundance by glycan composition and sample
relative_abundance <- glyco_peptide_groups_long %>%
  group_by(sample, glycan_composition) %>%
  summarise(
    total_abundance = sum(abundance, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  group_by(sample) %>%
  mutate(
    sample_total = sum(total_abundance, na.rm = TRUE),
    relative_abundance = (total_abundance / sample_total) * 100
  ) %>%
  ungroup() %>%
  # Add group classification based on sample names
  mutate(
    group = case_when(
      str_starts(sample, "hc") ~ "Healthy",
      str_starts(sample, "m") ~ "MECFS",
      TRUE ~ "Unknown"
    )
  )

# Calculate mean and standard error for each glycan composition
glycan_summary <- relative_abundance %>%
  group_by(glycan_composition) %>%
  summarise(
    mean_relative_abundance = mean(relative_abundance, na.rm = TRUE),
    se_relative_abundance = sd(relative_abundance, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )
#save the glycan summary
write_csv(glycan_summary, "input_data/Corr/glycan_composition_peptidegroups_intensity_summary.csv")

# Create bar plot with error bars
library(ggplot2)
library(patchwork)

glycan_plot <- ggplot(glycan_summary, aes(x = reorder(glycan_composition, -mean_relative_abundance), 
                                         y = mean_relative_abundance)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_relative_abundance - se_relative_abundance, 
                    ymax = mean_relative_abundance + se_relative_abundance), 
                width = 0.2, position = position_dodge(0.9)) +
  labs(title = "Relative Abundance of Glycan Compositions",
       x = "Glycan Composition",
       y = "Relative Abundance (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(hjust = 0.5))

# Save the plot
ggsave("figures/glycan_composition_relative_abundance.png", glycan_plot, 
       width = 12, height = 8, dpi = 300)

cat("Glycan composition plot saved to: figures/glycan_composition_relative_abundance.png\n")

# Calculate the sum of top 5 glycans for each sample
top_5_by_sample <- relative_abundance %>%
  group_by(sample) %>%
  arrange(desc(relative_abundance), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  summarise(
    top_5_sum = sum(relative_abundance),
    top_5_glycans = paste(glycan_composition, collapse = ", ")
  )


# Print summary
cat("\nTop 5 Glycan Compositions Summary:\n")


# Calculate the overall top 5 glycan compositions across all samples
overall_top_5 <- relative_abundance %>%
  group_by(glycan_composition) %>%
  summarise(
    mean_abundance = mean(relative_abundance, na.rm = TRUE),
    sd_abundance = sd(relative_abundance, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 5)

# Print the overall top 5 glycan compositions
cat("\nOverall Top 5 Most Abundant Glycan Compositions:\n")
for(i in 1:nrow(overall_top_5)) {
  cat(sprintf("%d. %s (Mean: %.2f%% ± %.2f%% SD)\n",
              i,
              overall_top_5$glycan_composition[i],
              overall_top_5$mean_abundance[i],
              overall_top_5$sd_abundance[i]))
}

# Perform t-test and F-test for each glycan composition between disease groups
glycan_ttests <- relative_abundance %>%
  group_by(glycan_composition) %>%
  summarise(
    # F-test for variance
    f_stat = var.test(relative_abundance[group == "MECFS"], 
                      relative_abundance[group == "Healthy"])$statistic,
    f_pvalue = var.test(relative_abundance[group == "MECFS"],
                        relative_abundance[group == "Healthy"])$p.value,
    # t-test 
    t_stat = t.test(relative_abundance ~ group)$statistic,
    p_value = t.test(relative_abundance ~ group)$p.value,
    mean_mecfs = mean(relative_abundance[group == "MECFS"], na.rm = TRUE),
    mean_healthy = mean(relative_abundance[group == "Healthy"], na.rm = TRUE), 
    sd_mecfs = sd(relative_abundance[group == "MECFS"], na.rm = TRUE),
    sd_healthy = sd(relative_abundance[group == "Healthy"], na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  # Add BH correction
  mutate(
    p_value_adj = p.adjust(p_value, method = "BH"),
    significant = p_value_adj < 0.05,
    fold_change = mean_mecfs / mean_healthy,
    variance_equal = f_pvalue >= 0.05
  ) %>%
  arrange(p_value_adj)

# Save t-test results
write.csv(glycan_ttests, "output_data/peptidegroups_intensity/glycan_composition_ttests.csv", row.names = TRUE)

# Print significant results
cat("\nSignificant differences in glycan compositions (p < 0.05):\n")
significant_glycans <- glycan_ttests %>% filter(significant)

if(nrow(significant_glycans) > 0) {
  for(i in 1:nrow(significant_glycans)) {
    cat(sprintf("\n%d. %s\n", i, significant_glycans$glycan_composition[i]))
    cat(sprintf("   MECFS: %.2f%% ± %.2f%% SD\n", 
                significant_glycans$mean_mecfs[i],
                significant_glycans$sd_mecfs[i]))
    cat(sprintf("   Healthy: %.2f%% ± %.2f%% SD\n", 
                significant_glycans$mean_healthy[i],
                significant_glycans$sd_healthy[i]))
    cat(sprintf("   p-value: %.3e\n", significant_glycans$p_value[i]))
    cat(sprintf("   Fold change: %.2f\n", significant_glycans$fold_change[i]))
  }
} else {
  cat("No significant differences found between MECFS and Healthy groups.\n")
}

# Find glycan compositions with fold change > 2 or < 0.5 (1/2)
high_fold_changes <- glycan_ttests %>%
  filter(fold_change > 2 | fold_change < (1/2)) %>%
  arrange(desc(fold_change))

# Print summary
cat("\nGlycan compositions with >2-fold change between groups:\n")
cat("Total count:", nrow(high_fold_changes), "\n")

if(nrow(high_fold_changes) > 0) {
  for(i in 1:nrow(high_fold_changes)) {
    cat(sprintf("\n%d. %s\n", i, high_fold_changes$glycan_composition[i]))
    cat(sprintf("   Fold change: %.2f\n", high_fold_changes$fold_change[i]))
    cat(sprintf("   MECFS: %.2f%% ± %.2f%% SD\n", 
                high_fold_changes$mean_mecfs[i],
                high_fold_changes$sd_mecfs[i]))
    cat(sprintf("   Healthy: %.2f%% ± %.2f%% SD\n", 
                high_fold_changes$mean_healthy[i],
                high_fold_changes$sd_healthy[i]))
    cat(sprintf("   Adjusted p-value: %.3e\n", high_fold_changes$p_value_adj[i]))
  }
} else {
    cat("No glycan compositions had >2-fold change between groups.\n")
}


# Create a bar plot of the top 5 glycan compositions
top_5_plot <- ggplot(overall_top_5, 
                     aes(x = reorder(glycan_composition, mean_abundance), 
                         y = mean_abundance)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_abundance - sd_abundance,
                    ymax = mean_abundance + sd_abundance),
                width = 0.2) +
  labs(title = "Top 5 Most Abundant Glycan Compositions",
       x = "Glycan Composition",
       y = "Mean Relative Abundance (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

# Save the plot
ggsave("figures/top_5_glycan_compositions.png", top_5_plot,
       width = 8, height = 6, dpi = 300)

cat("\nTop 5 glycan compositions plot saved to: figures/top_5_glycan_compositions.png\n")


# Calculate relative abundance by glycan class and sample
relative_abundance_class <- glyco_peptide_groups_long %>%
  group_by(sample, glycan_class) %>%
  summarise(
    total_abundance = sum(abundance, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  group_by(sample) %>%
  mutate(
    sample_total = sum(total_abundance, na.rm = TRUE),
    relative_abundance = (total_abundance / sample_total) * 100
  ) %>%
  ungroup()

# Create boxplot for glycan classes
glycan_class_plot <- ggplot(relative_abundance_class, aes(x = glycan_class, y = relative_abundance)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7, outlier.shape = 1) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  labs(title = "Relative Abundance by Glycan Class",
       x = "Glycan Class",
       y = "Relative Abundance (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

# Save the boxplot
ggsave("figures/glycan_class_relative_abundance_boxplot.png", glycan_class_plot, 
       width = 10, height = 8, dpi = 300)

cat("Glycan class boxplot saved to: figures/glycan_class_relative_abundance_boxplot.png\n")


# Create histogram of peptide scores
score_hist <- ggplot(glyco_peptide_groups_long, aes(x = log2(pep_2d_by_search_engine_a2_pmi_byonic))) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, color = "black") +
  labs(title = "Distribution of log2 Peptide Pep2d Scores",
       x = "log2 Peptide Pep2d Score", 
       y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Save the histogram
ggsave("figures/peptidegroups_intensity/peptide_log2_pep2d_histogram.png", score_hist,
       width = 8, height = 6, dpi = 300)

cat("Peptide scores histogram saved to: figures/peptide_log2_pep2d_histogram.png\n")

# Create boxplot of pep2d scores by sample
pep2d_sample_plot <- ggplot(glyco_peptide_groups_long, 
                           aes(x = sample, y = log2(pep_2d_by_search_engine_a2_pmi_byonic))) +
  geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.shape = 1) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  labs(title = "Peptide log2 Pep2d Score Distribution by Sample",
       x = "Sample", 
       y = "log2 Pep2d Score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

# Save the pep2d score boxplot
ggsave("figures/peptidegroups_intensity/log2_pep2d_score_by_sample_boxplot.png",
       pep2d_sample_plot,
       width = 12,
       height = 8,
       dpi = 300)

cat("Pep2d scores by sample boxplot saved to: figures/peptidegroups_intensity/log2_pep2d_score_by_sample_boxplot.png\n")


# Create boxplot of abundances by sample
sample_abundance_plot <- ggplot(glyco_peptide_groups_long, aes(x = sample, y = log2(abundance))) +
  geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.shape = 1) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  labs(title = "Abundance log2 Distribution by Sample",
       x = "Sample",
       y = "Abundance log2") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

# Save the sample abundance boxplot
ggsave("figures/peptidegroups_intensity/sample_abundance_boxplot.png", 
       sample_abundance_plot,
       width = 12, 
       height = 8, 
       dpi = 300)

cat("Sample abundance boxplot saved to: figures/peptidegroups_intensity/sample_abundance_boxplot.png\n")

# Calculate summary statistics for glycopeptide groups per sample
sample_stats <- glyco_peptide_groups_long %>%
  group_by(sample) %>%
  summarise(
    glycopeptide_count = n_distinct(peptide_groups_peptide_group_id),
    mean_abundance = mean(abundance, na.rm = TRUE),
    sd_abundance = sd(abundance, na.rm = TRUE)
  )

# Calculate overall mean and SD of glycopeptide counts
overall_stats <- sample_stats %>%
  summarise(
    mean_glycopeptides = mean(glycopeptide_count),
    sd_glycopeptides = sd(glycopeptide_count)
  )

# Print results
cat("\nGlycopeptide groups per sample summary:\n")
print(sample_stats)

cat("\nOverall statistics:\n")
cat(sprintf("Mean glycopeptide groups per sample: %.1f\n", overall_stats$mean_glycopeptides))
cat(sprintf("Standard deviation: %.1f\n", overall_stats$sd_glycopeptides))
# Calculate the number of unique peptide groups
unique_peptide_groups <- n_distinct(glyco_peptide_groups_long$peptide_groups_peptide_group_id)

cat("\nTotal number of unique peptide groups:", unique_peptide_groups, "\n")

# Save results to file
write.csv(sample_stats, "output_data/glycopeptide_groups_per_sample.csv", row.names = FALSE)
cat("\nDetailed per-sample statistics saved to: output_data/glycopeptide_groups_per_sample.csv\n")


# Calculate descriptive statistics for log2-transformed glycopeptide abundances
log2_glyco_stats <- data.frame(
  Log2_Abundance = c(
    Mean = mean(log2(glyco_peptide_groups_long$abundance), na.rm = TRUE),
    Median = median(log2(glyco_peptide_groups_long$abundance), na.rm = TRUE), 
    SD = sd(log2(glyco_peptide_groups_long$abundance), na.rm = TRUE),
    Min = min(log2(glyco_peptide_groups_long$abundance), na.rm = TRUE),
    Max = max(log2(glyco_peptide_groups_long$abundance), na.rm = TRUE),
    Q1 = quantile(log2(glyco_peptide_groups_long$abundance), 0.25, na.rm = TRUE),
    Q3 = quantile(log2(glyco_peptide_groups_long$abundance), 0.75, na.rm = TRUE),
    N = sum(!is.na(log2(glyco_peptide_groups_long$abundance))),
    N_missing = sum(is.na(log2(glyco_peptide_groups_long$abundance)))
  )
)

# Print log2-transformed statistics
cat("\nLog2-transformed abundance statistics:\n")
print(log2_glyco_stats)

# Save log2-transformed statistics to file
write.csv(log2_glyco_stats, "output_data/peptidegroups_intensity/log2_abundance_statistics.csv")
cat("\nLog2-transformed statistics saved to: output_data/peptidegroups_intensity/log2_abundance_statistics.csv\n")


# Create histogram of log2 transformed abundance across all samples
abundance_histogram <- ggplot(glyco_peptide_groups_long, aes(x = log2(abundance))) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of Log2 Abundance Values",
       x = "Log2 Abundance",
       y = "Count") +
  theme_minimal()

# Save the histogram
ggsave("figures/peptidegroups_intensity/abundance_distribution_histogram.png",
       abundance_histogram,
       width = 10,
       height = 6,
       dpi = 300)

# Perform normality tests for abundance data
cat("\n=== NORMALITY TESTING FOR ABUNDANCE DATA ===\n")

# Get log2 transformed abundance data
log_abundance <- log2(glyco_peptide_groups_long$abundance)
n_obs <- length(log_abundance)

cat("Total observations:", n_obs, "\n")

# Initialize results dataframe
normality_results <- data.frame()

# 1. Shapiro-Wilk test (if sample size allows)
if(n_obs <= 5000) {
  shapiro_test <- shapiro.test(log_abundance)
  cat("\nShapiro-Wilk normality test results:\n")
  cat("W statistic:", shapiro_test$statistic, "\n")
  cat("p-value:", shapiro_test$p.value, "\n")
  cat("Interpretation:", ifelse(shapiro_test$p.value < 0.05,
                              "Data significantly deviates from normal distribution",
                              "Data appears to be normally distributed"), "\n")
  
  normality_results <- rbind(normality_results, data.frame(
    test = "Shapiro-Wilk",
    statistic = shapiro_test$statistic,
    p_value = shapiro_test$p.value,
    sample_size = n_obs,
    stringsAsFactors = FALSE
  ))
} else {
  cat("\nShapiro-Wilk test skipped: sample size (", n_obs, ") exceeds limit of 5000\n")
  
  # Perform Shapiro-Wilk on a random sample
  set.seed(123)  # For reproducibility
  sample_data <- sample(log_abundance, 5000)
  shapiro_sample <- shapiro.test(sample_data)
  
  cat("Shapiro-Wilk test on random sample of 5000 observations:\n")
  cat("W statistic:", shapiro_sample$statistic, "\n")
  cat("p-value:", shapiro_sample$p.value, "\n")
  cat("Interpretation:", ifelse(shapiro_sample$p.value < 0.05,
                              "Sample significantly deviates from normal distribution",
                              "Sample appears to be normally distributed"), "\n")
  
  normality_results <- rbind(normality_results, data.frame(
    test = "Shapiro-Wilk (sample)",
    statistic = shapiro_sample$statistic,
    p_value = shapiro_sample$p.value,
    sample_size = 5000,
    stringsAsFactors = FALSE
  ))
}

# 2. Kolmogorov-Smirnov test (no sample size limit)
ks_test <- ks.test(log_abundance, "pnorm", mean = mean(log_abundance, na.rm = TRUE), 
                   sd = sd(log_abundance, na.rm = TRUE))
cat("\nKolmogorov-Smirnov normality test results:\n")
cat("D statistic:", ks_test$statistic, "\n")
cat("p-value:", ks_test$p.value, "\n")
cat("Interpretation:", ifelse(ks_test$p.value < 0.05,
                            "Data significantly deviates from normal distribution",
                            "Data appears to be normally distributed"), "\n")

normality_results <- rbind(normality_results, data.frame(
  test = "Kolmogorov-Smirnov",
  statistic = ks_test$statistic,
  p_value = ks_test$p.value,
  sample_size = n_obs,
  stringsAsFactors = FALSE
))

# 3. Anderson-Darling test (if nortest package is available)
if(require(nortest, quietly = TRUE)) {
  ad_test <- ad.test(log_abundance)
  cat("\nAnderson-Darling normality test results:\n")
  cat("A statistic:", ad_test$statistic, "\n")
  cat("p-value:", ad_test$p.value, "\n")
  cat("Interpretation:", ifelse(ad_test$p.value < 0.05,
                              "Data significantly deviates from normal distribution",
                              "Data appears to be normally distributed"), "\n")
  
  normality_results <- rbind(normality_results, data.frame(
    test = "Anderson-Darling",
    statistic = ad_test$statistic,
    p_value = ad_test$p.value,
    sample_size = n_obs,
    stringsAsFactors = FALSE
  ))
} else {
  cat("\nAnderson-Darling test skipped: 'nortest' package not available\n")
  cat("To install: install.packages('nortest')\n")
}

# 4. Summary statistics for normality assessment
cat("\nSummary statistics for log2(abundance):\n")
summary_stats <- summary(log_abundance)
print(summary_stats)

# Calculate skewness and kurtosis
if(require(moments, quietly = TRUE)) {
  skewness_val <- skewness(log_abundance)
  kurtosis_val <- kurtosis(log_abundance)
  
  cat("\nSkewness:", round(skewness_val, 4), "\n")
  cat("Kurtosis:", round(kurtosis_val, 4), "\n")
  cat("Interpretation:\n")
  cat("  Skewness: 0 = symmetric, >0 = right-skewed, <0 = left-skewed\n")
  cat("  Kurtosis: 3 = normal, >3 = heavy-tailed, <3 = light-tailed\n")
  
  # Add skewness and kurtosis to results
  normality_results <- rbind(normality_results, data.frame(
    test = c("Skewness", "Kurtosis"),
    statistic = c(skewness_val, kurtosis_val),
    p_value = c(NA, NA),
    sample_size = c(n_obs, n_obs),
    stringsAsFactors = FALSE
  ))
} else {
  cat("\nSkewness and kurtosis skipped: 'moments' package not available\n")
  cat("To install: install.packages('moments')\n")
}

# Save test results to file
write.csv(normality_results, "output_data/peptidegroups_intensity/abundance_normality_test.csv", row.names = FALSE)
cat("\nNormality test results saved to: output_data/peptidegroups_intensity/abundance_normality_test.csv\n")

# Create a histogram with normal curve overlay for visual assessment
normality_plot <- ggplot(data.frame(log_abundance = log_abundance), 
                        aes(x = log_abundance)) +
  geom_histogram(aes(y = ..density..), bins = 50, alpha = 0.7, 
                 fill = "#9DD4CC", color = "black") +
  stat_function(fun = dnorm, 
                args = list(mean = mean(log_abundance, na.rm = TRUE), 
                           sd = sd(log_abundance, na.rm = TRUE)),
                color = "#E49CB1", size = 1.2) +
  labs(title = "Distribution of log2(Abundance) with Normal Curve Overlay",
       x = "log2(Abundance)",
       y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Save the plot
ggsave("figures/peptidegroups_intensity/abundance_normality_plot.png", 
       normality_plot, width = 10, height = 6, dpi = 300)
cat("Normality assessment plot saved to: figures/peptidegroups_intensity/abundance_normality_plot.png\n")



# Calculate descriptive statistics for glycopeptide group abundances
glyco_stats <- data.frame(
  Raw_Abundance = c(
    Mean = mean(glyco_peptide_groups_long$abundance, na.rm = TRUE),
    Median = median(glyco_peptide_groups_long$abundance, na.rm = TRUE),
    SD = sd(glyco_peptide_groups_long$abundance, na.rm = TRUE),
    Min = min(glyco_peptide_groups_long$abundance, na.rm = TRUE),
    Max = max(glyco_peptide_groups_long$abundance, na.rm = TRUE),
    Q1 = quantile(glyco_peptide_groups_long$abundance, 0.25, na.rm = TRUE),
    Q3 = quantile(glyco_peptide_groups_long$abundance, 0.75, na.rm = TRUE),
    N = sum(!is.na(glyco_peptide_groups_long$abundance)),
    N_missing = sum(is.na(glyco_peptide_groups_long$abundance))
  ),
  Log2_Abundance = c(
    Mean = mean(log2(glyco_peptide_groups_long$abundance), na.rm = TRUE),
    Median = median(log2(glyco_peptide_groups_long$abundance), na.rm = TRUE),
    SD = sd(log2(glyco_peptide_groups_long$abundance), na.rm = TRUE),
    Min = min(log2(glyco_peptide_groups_long$abundance), na.rm = TRUE),
    Max = max(log2(glyco_peptide_groups_long$abundance), na.rm = TRUE),
    Q1 = quantile(log2(glyco_peptide_groups_long$abundance), 0.25, na.rm = TRUE),
    Q3 = quantile(log2(glyco_peptide_groups_long$abundance), 0.75, na.rm = TRUE),
    N = sum(!is.na(log2(glyco_peptide_groups_long$abundance))),
    N_missing = sum(is.na(log2(glyco_peptide_groups_long$abundance)))
  ),
  row.names = c("Mean", "Median", "SD", "Min", "Max", "Q1", "Q3", "N", "N_missing")
)

# Print the statistics
print("Descriptive Statistics for Glycopeptide Group Abundances:")
print(round(glyco_stats, 3))

# Calculate statistics by sample using tapply for each metric
sample_means <- tapply(glyco_peptide_groups_long$abundance, glyco_peptide_groups_long$sample, mean, na.rm = TRUE)
sample_medians <- tapply(glyco_peptide_groups_long$abundance, glyco_peptide_groups_long$sample, median, na.rm = TRUE)
sample_sds <- tapply(glyco_peptide_groups_long$abundance, glyco_peptide_groups_long$sample, sd, na.rm = TRUE)
sample_n <- tapply(glyco_peptide_groups_long$abundance, glyco_peptide_groups_long$sample, length)
sample_missing <- tapply(glyco_peptide_groups_long$abundance, glyco_peptide_groups_long$sample, function(x) sum(is.na(x)))

# Combine into a data frame
sample_stats_df <- data.frame(
  Sample = names(sample_means),
  Mean = as.numeric(sample_means),
  Median = as.numeric(sample_medians),
  SD = as.numeric(sample_sds),
  N = as.numeric(sample_n),
  Missing = as.numeric(sample_missing)
)

# Print sample statistics
print("Sample-level statistics:")
# Round only numeric columns, excluding the Sample column
sample_stats_df_rounded <- sample_stats_df
sample_stats_df_rounded[, c("Mean", "Median", "SD", "N", "Missing")] <- 
  round(sample_stats_df[, c("Mean", "Median", "SD", "N", "Missing")], 3)
print(sample_stats_df_rounded)

# Save statistics to file
write.csv(glyco_stats, "output_data/peptidegroups_intensity/glycopeptide_groups_overall_stats.csv", row.names = TRUE)
write.csv(sample_stats_df, "output_data/peptidegroups_intensity/glycopeptide_groups_sample_stats.csv", row.names = FALSE)

cat("\nDescriptive statistics saved to:\n")
cat("- output_data/glycopeptide_groups_overall_stats.csv\n")
cat("- output_data/glycopeptide_groups_sample_stats.csv\n")


###########################################################################################################
##pseudo glycomics 

# Add group classification based on sample names if not already present
if (!"group" %in% colnames(glyco_peptide_groups_long)) {
  glyco_peptide_groups_long <- glyco_peptide_groups_long %>%
    mutate(group = case_when(
      str_starts(sample, "hc") ~ "Healthy",
      str_starts(sample, "m") ~ "MECFS",
      TRUE ~ "Unknown"
    ))
}

pseudo_glycomics <- glyco_peptide_groups_long %>%
  group_by(group) %>%
  # First calculate total abundance per group
  mutate(total_abundance = sum(abundance, na.rm = TRUE)) %>%
  # Then group by both group and glycan composition
  group_by(glycan_composition, group) %>%
  summarise(
    glycan_abundance = sum(abundance, na.rm = TRUE),
    total_group_abundance = first(total_abundance),
    relative_percentage = (glycan_abundance / total_group_abundance) * 100,
    .groups = 'drop'
  ) %>%
  arrange(group, desc(relative_percentage))

# Print summary
cat("\nRelative Glycan Composition Percentages by Group:\n")
print(pseudo_glycomics)

# Save results
write.csv(pseudo_glycomics, "output_data/peptidegroups_intensity/glycan_composition_percentages.csv", row.names = FALSE)

# Create summary statistics
glycan_summary <- pseudo_glycomics %>%
  group_by(group) %>%
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

# Create detailed summary statistics for each glycan composition by group
glycan_detailed_summary <- glyco_peptide_groups_long %>%
  group_by(glycan_composition, group) %>%
  summarise(
    n_observations = n(),
    mean_abundance = mean(abundance, na.rm = TRUE),
    median_abundance = median(abundance, na.rm = TRUE),
    sd_abundance = sd(abundance, na.rm = TRUE),
    cv = (sd_abundance / mean_abundance) * 100,  # Coefficient of variation
    min_abundance = min(abundance, na.rm = TRUE),
    max_abundance = max(abundance, na.rm = TRUE),
    q25 = quantile(abundance, 0.25, na.rm = TRUE),
    q75 = quantile(abundance, 0.75, na.rm = TRUE),
    relative_percentage = (sum(abundance, na.rm = TRUE) / 
                         sum(glyco_peptide_groups_long$abundance[glyco_peptide_groups_long$group == cur_group()$group], 
                             na.rm = TRUE)) * 100,
    missing_values = sum(is.na(abundance)),
    .groups = 'drop'
  ) %>%
  arrange(glycan_composition, group)

# Add difference between groups
glycan_comparison <- glycan_detailed_summary %>%
  dplyr::select(glycan_composition, group, relative_percentage) %>%
  tidyr::pivot_wider(names_from = group, 
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
          "output_data/peptidegroups_intensity/glycan_composition_detailed_summary.csv", 
          row.names = FALSE)

# Create a simplified summary for quick review
glycan_simple_summary <- glycan_detailed_summary %>%
  dplyr::select(glycan_composition, 
         group, 
         n_observations, 
         mean_abundance, 
         relative_percentage, 
         cv,
         percentage_diff_from_healthy) %>%
  arrange(desc(abs(percentage_diff_from_healthy)))

# Save simplified summary
write.csv(glycan_simple_summary, 
          "output_data/peptidegroups_intensity/glycan_composition_simple_summary.csv", 
          row.names = FALSE)

# Print summary of the most different glycans
cat("\nTop 10 glycans with largest difference between groups:\n")
print(glycan_simple_summary %>% 
        dplyr::filter(group == "MECFS") %>% 
        dplyr::select(glycan_composition, percentage_diff_from_healthy) %>% 
        dplyr::arrange(desc(abs(percentage_diff_from_healthy))) %>% 
        head(10))

# Print overall summary
cat("\nSummary of analysis:\n")
cat("Total unique glycan compositions:", length(unique(glyco_peptide_groups_long$glycan_composition)), "\n")
cat("Number of groups:", length(unique(glyco_peptide_groups_long$group)), "\n")
cat("Total observations:", nrow(glyco_peptide_groups_long), "\n")

# Calculate summary statistics for error bars from original data
# We need to calculate relative percentages for each sample first, then get the mean and SE
glycan_summary <- glyco_peptide_groups_long %>%
  # Calculate total abundance per sample
  group_by(sample) %>%
  mutate(sample_total = sum(abundance, na.rm = TRUE)) %>%
  # Calculate relative percentage for each glycan composition per sample
  group_by(sample, glycan_composition) %>%
  summarise(
    glycan_abundance = sum(abundance, na.rm = TRUE),
    sample_total = first(sample_total),
    relative_percentage = (glycan_abundance / sample_total) * 100,
    .groups = 'drop'
  ) %>%
  # Add group classification
  mutate(
    group = case_when(
      str_starts(sample, "hc") ~ "Healthy",
      str_starts(sample, "m") ~ "MECFS",
      TRUE ~ "Unknown"
    )
  ) %>%
  # Now calculate summary statistics with proper standard error
  group_by(glycan_composition, group) %>%
  summarise(
    mean_percentage = mean(relative_percentage, na.rm = TRUE),
    sd_percentage = sd(relative_percentage, na.rm = TRUE),
    std.error = sd_percentage / sqrt(n()),
    n_samples = n(),
    .groups = 'drop'
  )

# Print summary to verify error bars are calculated
cat("\nGlycan Summary with Standard Errors:\n")
print(head(glycan_summary, 10))

# Create bar plot with error bars
glycan_plot <- glycan_summary %>%
  ggplot(aes(x = reorder(glycan_composition, mean_percentage), 
             y = mean_percentage, 
             fill = group)) +
  geom_bar(stat = "identity", 
           position = "dodge", 
           width = 0.8) +
  geom_errorbar(aes(ymin = mean_percentage - std.error, 
                    ymax = mean_percentage + std.error),
                position = position_dodge(width = 0.8),
                width = 0.3,
                color = "black",
                size = 0.5) +
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
    title = "Glycan Composition Relative Percentages by Group (with Standard Error)",
    x = "Glycan Composition",
    y = "Relative Percentage (%)",
    fill = "Group"
  ) +
  scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1"))

# Save the plot
ggsave("figures/peptidegroups_intensity/glycan_composition_comparison.pdf", 
       glycan_plot, 
       width = 12, 
       height = 8, 
       dpi = 300)

# Save the complete plot as both PDF and PNG
ggsave("figures/peptidegroups_intensity/glycan_composition_comparison.png", 
       glycan_plot, 
       width = 12, 
       height = 8, 
       dpi = 300)

# Create a focused plot showing only top differential glycans with error bars
top_glycans_data <- pseudo_glycomics %>%
  # Get top 15 glycans with biggest differences between groups
  group_by(glycan_composition) %>%
  summarise(diff = abs(diff(relative_percentage))) %>%
  arrange(desc(diff)) %>%
  head(15) %>%
  pull(glycan_composition)

# Calculate summary statistics for top differential glycans from original data
top_glycans_summary <- glyco_peptide_groups_long %>%
  # Filter to only the top differential glycans
  filter(glycan_composition %in% top_glycans_data) %>%
  # Calculate total abundance per sample
  group_by(sample) %>%
  mutate(sample_total = sum(abundance, na.rm = TRUE)) %>%
  # Calculate relative percentage for each glycan composition per sample
  group_by(sample, glycan_composition) %>%
  summarise(
    glycan_abundance = sum(abundance, na.rm = TRUE),
    sample_total = first(sample_total),
    relative_percentage = (glycan_abundance / sample_total) * 100,
    .groups = 'drop'
  ) %>%
  # Add group classification
  mutate(
    group = case_when(
      str_starts(sample, "hc") ~ "Healthy",
      str_starts(sample, "m") ~ "MECFS",
      TRUE ~ "Unknown"
    )
  ) %>%
  # Now calculate summary statistics with proper standard error
  group_by(glycan_composition, group) %>%
  summarise(
    mean_percentage = mean(relative_percentage, na.rm = TRUE),
    sd_percentage = sd(relative_percentage, na.rm = TRUE),
    std.error = sd_percentage / sqrt(n()),
    n_samples = n(),
    .groups = 'drop'
  )

# Print summary to verify error bars are calculated for top glycans
cat("\nTop Glycans Summary with Standard Errors:\n")
print(head(top_glycans_summary, 10))

top_glycans_plot <- top_glycans_summary %>%
  ggplot(aes(x = reorder(glycan_composition, mean_percentage), 
             y = mean_percentage, 
             fill = group)) +
  geom_bar(stat = "identity", 
           position = "dodge", 
           width = 0.8) +
  geom_errorbar(aes(ymin = mean_percentage - std.error, 
                    ymax = mean_percentage + std.error),
                position = position_dodge(width = 0.8),
                width = 0.3,
                color = "black",
                size = 0.5) +
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
    title = "Top 15 Differential Glycan Compositions (with Standard Error)",
    x = "Glycan Composition",
    y = "Relative Percentage (%)",
    fill = "Group"
  ) +
  scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1"))

# Save the focused plot
ggsave("figures/peptidegroups_intensity/top_differential_glycans.pdf", 
       top_glycans_plot, 
       width = 10, 
       height = 6, 
       dpi = 300)

# Save the focused plot as both PDF and PNG
ggsave("figures/peptidegroups_intensity/top_differential_glycans.png", 
       top_glycans_plot, 
       width = 10, 
       height = 6, 
       dpi = 300)


########################################################################################

# Calculate fucosylation percentages by sample
fucosylation_by_sample <- glyco_peptide_groups_long %>%
  # Group by sample and group to get total abundance per sample
  group_by(sample, group) %>%
  mutate(total_sample_abundance = sum(abundance, na.rm = TRUE)) %>%
  # Calculate Fuc percentage per sample using existing contains_Fuc column
  summarise(
    fuc_abundance = sum(abundance[contains_Fuc == TRUE], na.rm = TRUE),
    total_abundance = first(total_sample_abundance),
    fuc_percentage = (fuc_abundance / total_abundance) * 100,
    .groups = 'drop'
  )

# Perform t-test
t_test_result <- t.test(fuc_percentage ~ group, data = fucosylation_by_sample)

# Calculate F-test for variance (homogeneity of variance test)
healthy_var <- var(fucosylation_by_sample$fuc_percentage[fucosylation_by_sample$group == "Healthy"])
mecfs_var <- var(fucosylation_by_sample$fuc_percentage[fucosylation_by_sample$group == "MECFS"])

# F-test for variance (larger variance in numerator)
if (healthy_var >= mecfs_var) {
  f_ratio <- healthy_var / mecfs_var
  df1 <- sum(fucosylation_by_sample$group == "Healthy") - 1
  df2 <- sum(fucosylation_by_sample$group == "MECFS") - 1
} else {
  f_ratio <- mecfs_var / healthy_var
  df1 <- sum(fucosylation_by_sample$group == "MECFS") - 1
  df2 <- sum(fucosylation_by_sample$group == "Healthy") - 1
}

f_p_value <- 2 * (1 - pf(f_ratio, df1, df2))  # Two-tailed test

# Calculate Cohen's d effect size
library(effectsize)
cohens_d_result <- cohens_d(fuc_percentage ~ group, data = fucosylation_by_sample)
cohens_d_value <- cohens_d_result$Cohens_d

# Create labels for the plot
t_test_label <- sprintf("t-test: p = %.3g", t_test_result$p.value)
f_test_label <- sprintf("F-var = %.3f, p = %.3g", f_ratio, f_p_value)
cohens_d_label <- sprintf("Cohen's d = %.3f", cohens_d_value)

# Create boxplot with individual points
fuc_boxplot <- ggplot(fucosylation_by_sample, 
                     aes(x = group, 
                         y = fuc_percentage, 
                         fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  labs(
    title = "Fucosylation Percentage by Group",
    x = "Group",
    y = "Fucosylation Percentage (%)"
  ) +
  # Add statistical annotations
  annotate("text", 
           x = 1.5, 
           y = max(fucosylation_by_sample$fuc_percentage) * 1.15,
           label = t_test_label,
           size = 3.5) +
  annotate("text", 
           x = 1.5, 
           y = max(fucosylation_by_sample$fuc_percentage) * 1.05,
           label = f_test_label,
           size = 3.5) +
  annotate("text", 
           x = 1.5, 
           y = max(fucosylation_by_sample$fuc_percentage) * 0.95,
           label = cohens_d_label,
           size = 3.5)

# Save the plot
ggsave("figures/peptidegroups_intensity/fucosylation_comparison.pdf", 
       fuc_boxplot, 
       width = 8, 
       height = 6, 
       dpi = 300)

ggsave("figures/peptidegroups_intensity/fucosylation_comparison.png", 
       fuc_boxplot, 
       width = 8, 
       height = 6, 
       dpi = 300)

# Calculate standard error and create error plots
########################################################################################

# Calculate summary statistics including standard error
fucosylation_summary <- fucosylation_by_sample %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean_fuc_percentage = mean(fuc_percentage, na.rm = TRUE),
    sd_fuc_percentage = sd(fuc_percentage, na.rm = TRUE),
    se_fuc_percentage = sd_fuc_percentage / sqrt(n),
    ci_lower = mean_fuc_percentage - (1.96 * se_fuc_percentage),
    ci_upper = mean_fuc_percentage + (1.96 * se_fuc_percentage),
    .groups = 'drop'
  )

# Print the summary statistics
cat("Fucosylation Percentage Summary Statistics:\n")
print(fucosylation_summary)

# Create bar plot with error bars (mean ± SE)
fuc_barplot_se <- ggplot(fucosylation_summary, 
                         aes(x = group, 
                             y = mean_fuc_percentage, 
                             fill = group)) +
  geom_col(alpha = 0.7, width = 0.6) +
  geom_errorbar(aes(ymin = mean_fuc_percentage - se_fuc_percentage,
                    ymax = mean_fuc_percentage + se_fuc_percentage),
                width = 0.2, size = 1, color = "black") +
  scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  labs(
    title = "Fucosylation Percentage by Group (Mean ± SE)",
    x = "Group",
    y = "Fucosylation Percentage (%)"
  ) +
  # Add sample size annotations
  geom_text(aes(label = paste("n =", n)), 
            vjust = -0.5, size = 3.5) +
  # Add mean value annotations
  geom_text(aes(label = sprintf("Mean = %.2f%%", mean_fuc_percentage)), 
            vjust = -1.5, size = 3.5) +
  # Add SE annotations
  geom_text(aes(label = sprintf("SE = %.2f%%", se_fuc_percentage)), 
            vjust = -2.5, size = 3.5)

# Create bar plot with 95% confidence intervals
fuc_barplot_ci <- ggplot(fucosylation_summary, 
                         aes(x = group, 
                             y = mean_fuc_percentage, 
                             fill = group)) +
  geom_col(alpha = 0.7, width = 0.6) +
  geom_errorbar(aes(ymin = ci_lower,
                    ymax = ci_upper),
                width = 0.2, size = 1, color = "black") +
  scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  labs(
    title = "Fucosylation Percentage by Group (Mean ± 95% CI)",
    x = "Group",
    y = "Fucosylation Percentage (%)"
  ) +
  # Add sample size annotations
  geom_text(aes(label = paste("n =", n)), 
            vjust = -0.5, size = 3.5) +
  # Add mean value annotations
  geom_text(aes(label = sprintf("Mean = %.2f%%", mean_fuc_percentage)), 
            vjust = -1.5, size = 3.5) +
  # Add CI annotations
  geom_text(aes(label = sprintf("95%% CI = [%.2f, %.2f]", ci_lower, ci_upper)), 
            vjust = -2.5, size = 3)

# Update the existing boxplot to include standard error information
fuc_boxplot_with_se <- ggplot(fucosylation_by_sample, 
                             aes(x = group, 
                                 y = fuc_percentage, 
                                 fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  labs(
    title = "Fucosylation Percentage by Group (Boxplot with SE Info)",
    x = "Group",
    y = "Fucosylation Percentage (%)"
  ) +
  # Add statistical annotations including SE
  annotate("text", 
           x = 1.5, 
           y = max(fucosylation_by_sample$fuc_percentage) * 1.15,
           label = t_test_label,
           size = 3.5) +
  annotate("text", 
           x = 1.5, 
           y = max(fucosylation_by_sample$fuc_percentage) * 1.05,
           label = f_test_label,
           size = 3.5) +
  annotate("text", 
           x = 1.5, 
           y = max(fucosylation_by_sample$fuc_percentage) * 0.95,
           label = cohens_d_label,
           size = 3.5) +
  # Add SE information
  annotate("text", 
           x = 1, 
           y = max(fucosylation_by_sample$fuc_percentage) * 0.85,
           label = sprintf("SE = %.2f%%", fucosylation_summary$se_fuc_percentage[fucosylation_summary$group == "Healthy"]),
           size = 3.5) +
  annotate("text", 
           x = 2, 
           y = max(fucosylation_by_sample$fuc_percentage) * 0.85,
           label = sprintf("SE = %.2f%%", fucosylation_summary$se_fuc_percentage[fucosylation_summary$group == "MECFS"]),
           size = 3.5)

# Save the new error plots
ggsave("figures/peptidegroups_intensity/fucosylation_mean_se.pdf", 
       fuc_barplot_se, 
       width = 8, 
       height = 6, 
       dpi = 300)

ggsave("figures/peptidegroups_intensity/fucosylation_mean_se.png", 
       fuc_barplot_se, 
       width = 8, 
       height = 6, 
       dpi = 300)

ggsave("figures/peptidegroups_intensity/fucosylation_mean_ci.pdf", 
       fuc_barplot_ci, 
       width = 8, 
       height = 6, 
       dpi = 300)

ggsave("figures/peptidegroups_intensity/fucosylation_mean_ci.png", 
       fuc_barplot_ci, 
       width = 8, 
       height = 6, 
       dpi = 300)

ggsave("figures/peptidegroups_intensity/fucosylation_boxplot_with_se.pdf", 
       fuc_boxplot_with_se, 
       width = 8, 
       height = 6, 
       dpi = 300)

ggsave("figures/peptidegroups_intensity/fucosylation_boxplot_with_se.png", 
       fuc_boxplot_with_se, 
       width = 8, 
       height = 6, 
       dpi = 300)

# Create a combined summary table for export
fucosylation_detailed_summary <- fucosylation_by_sample %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean = mean(fuc_percentage, na.rm = TRUE),
    median = median(fuc_percentage, na.rm = TRUE),
    sd = sd(fuc_percentage, na.rm = TRUE),
    se = sd / sqrt(n),
    ci_lower = mean - (1.96 * se),
    ci_upper = mean + (1.96 * se),
    min = min(fuc_percentage, na.rm = TRUE),
    max = max(fuc_percentage, na.rm = TRUE),
    .groups = 'drop'
  )

# Save the detailed summary
write.csv(fucosylation_detailed_summary, 
          "output_data/peptidegroups_intensity/fucosylation_detailed_summary_with_se.csv", 
          row.names = FALSE)

cat("\nStandard Error Analysis Complete!\n")
cat("Created plots with error bars and saved detailed summary statistics.\n")

# Print summary statistics
cat("\nFucosylation Summary Statistics:\n")
fucosylation_summary <- fucosylation_by_sample %>%
  group_by(group) %>%
  summarise(
    mean_fuc = mean(fuc_percentage),
    sd_fuc = sd(fuc_percentage),
    median_fuc = median(fuc_percentage),
    n_samples = n()
  )
print(fucosylation_summary)

cat("\nT-test Results:\n")
print(t_test_result)

cat("\nF-test for Variance (Homogeneity of Variance):\n")
cat(sprintf("Healthy group variance: %.3f\n", healthy_var))
cat(sprintf("MECFS group variance: %.3f\n", mecfs_var))
cat(sprintf("F-ratio: %.3f\n", f_ratio))
cat(sprintf("Degrees of freedom: %d, %d\n", df1, df2))
cat(sprintf("p-value: %.3g\n", f_p_value))

cat("\nCohen's d Effect Size:\n")
print(cohens_d_result)

# Interpret Cohen's d
cohens_d_interpretation <- case_when(
  abs(cohens_d_value) < 0.2 ~ "negligible",
  abs(cohens_d_value) < 0.5 ~ "small",
  abs(cohens_d_value) < 0.8 ~ "medium",
  TRUE ~ "large"
)

cat(sprintf("\nEffect Size Interpretation: Cohen's d = %.3f (%s effect)\n", 
            cohens_d_value, cohens_d_interpretation))

########################################################################################

# Calculate NeuAc percentages by sample
neuac_by_sample <- glyco_peptide_groups_long %>%
  # Group by sample and group to get total abundance per sample
  group_by(sample, group) %>%
  mutate(total_sample_abundance = sum(abundance, na.rm = TRUE)) %>%
  # Calculate NeuAc percentage per sample using existing contains_NeuAc column
  summarise(
    neuac_abundance = sum(abundance[contains_NeuAc == TRUE], na.rm = TRUE),
    total_abundance = first(total_sample_abundance),
    neuac_percentage = (neuac_abundance / total_abundance) * 100,
    .groups = 'drop'
  )

# Perform t-test
t_test_result_neuac <- t.test(neuac_percentage ~ group, data = neuac_by_sample)

# Calculate F-test for variance (homogeneity of variance test)
healthy_var_neuac <- var(neuac_by_sample$neuac_percentage[neuac_by_sample$group == "Healthy"])
mecfs_var_neuac <- var(neuac_by_sample$neuac_percentage[neuac_by_sample$group == "MECFS"])

# F-test for variance (larger variance in numerator)
if (healthy_var_neuac >= mecfs_var_neuac) {
  f_ratio_neuac <- healthy_var_neuac / mecfs_var_neuac
  df1_neuac <- sum(neuac_by_sample$group == "Healthy") - 1
  df2_neuac <- sum(neuac_by_sample$group == "MECFS") - 1
} else {
  f_ratio_neuac <- mecfs_var_neuac / healthy_var_neuac
  df1_neuac <- sum(neuac_by_sample$group == "MECFS") - 1
  df2_neuac <- sum(neuac_by_sample$group == "Healthy") - 1
}

f_p_value_neuac <- 2 * (1 - pf(f_ratio_neuac, df1_neuac, df2_neuac))  # Two-tailed test

# Calculate Cohen's d effect size
cohens_d_result_neuac <- cohens_d(neuac_percentage ~ group, data = neuac_by_sample)
cohens_d_value_neuac <- cohens_d_result_neuac$Cohens_d

# Create labels for the plot
t_test_label_neuac <- sprintf("t-test: p = %.3g", t_test_result_neuac$p.value)
f_test_label_neuac <- sprintf("F-var = %.3f, p = %.3g", f_ratio_neuac, f_p_value_neuac)
cohens_d_label_neuac <- sprintf("Cohen's d = %.3f", cohens_d_value_neuac)

# Create boxplot with individual points
neuac_boxplot <- ggplot(neuac_by_sample, 
                       aes(x = group, 
                           y = neuac_percentage, 
                           fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  labs(
    title = "Sialylation (NeuAc) Percentage by Group",
    x = "Group",
    y = "Sialylation Percentage (%)"
  ) +
  # Add statistical annotations
  annotate("text", 
           x = 1.5, 
           y = max(neuac_by_sample$neuac_percentage) * 1.15,
           label = t_test_label_neuac,
           size = 3.5) +
  annotate("text", 
           x = 1.5, 
           y = max(neuac_by_sample$neuac_percentage) * 1.05,
           label = f_test_label_neuac,
           size = 3.5) +
  annotate("text", 
           x = 1.5, 
           y = max(neuac_by_sample$neuac_percentage) * 0.95,
           label = cohens_d_label_neuac,
           size = 3.5)

# Save the plot
ggsave("figures/peptidegroups_intensity/sialylation_comparison.pdf", 
       neuac_boxplot, 
       width = 8, 
       height = 6, 
       dpi = 300)

ggsave("figures/peptidegroups_intensity/sialylation_comparison.png", 
       neuac_boxplot, 
       width = 8, 
       height = 6, 
       dpi = 300)

# Calculate standard error and create error plots for sialylation
########################################################################################

# Calculate summary statistics including standard error for sialylation
neuac_summary <- neuac_by_sample %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean_neuac_percentage = mean(neuac_percentage, na.rm = TRUE),
    sd_neuac_percentage = sd(neuac_percentage, na.rm = TRUE),
    se_neuac_percentage = sd_neuac_percentage / sqrt(n),
    ci_lower = mean_neuac_percentage - (1.96 * se_neuac_percentage),
    ci_upper = mean_neuac_percentage + (1.96 * se_neuac_percentage),
    .groups = 'drop'
  )

# Print the summary statistics
cat("Sialylation Percentage Summary Statistics:\n")
print(neuac_summary)

# Create bar plot with error bars (mean ± SE) for sialylation
neuac_barplot_se <- ggplot(neuac_summary, 
                           aes(x = group, 
                               y = mean_neuac_percentage, 
                               fill = group)) +
  geom_col(alpha = 0.7, width = 0.6) +
  geom_errorbar(aes(ymin = mean_neuac_percentage - se_neuac_percentage,
                    ymax = mean_neuac_percentage + se_neuac_percentage),
                width = 0.2, size = 1, color = "black") +
  scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  labs(
    title = "Sialylation Percentage by Group (Mean ± SE)",
    x = "Group",
    y = "Sialylation Percentage (%)"
  ) +
  # Add sample size annotations
  geom_text(aes(label = paste("n =", n)), 
            vjust = -0.5, size = 3.5) +
  # Add mean value annotations
  geom_text(aes(label = sprintf("Mean = %.2f%%", mean_neuac_percentage)), 
            vjust = -1.5, size = 3.5) +
  # Add SE annotations
  geom_text(aes(label = sprintf("SE = %.2f%%", se_neuac_percentage)), 
            vjust = -2.5, size = 3.5)

# Create bar plot with 95% confidence intervals for sialylation
neuac_barplot_ci <- ggplot(neuac_summary, 
                           aes(x = group, 
                               y = mean_neuac_percentage, 
                               fill = group)) +
  geom_col(alpha = 0.7, width = 0.6) +
  geom_errorbar(aes(ymin = ci_lower,
                    ymax = ci_upper),
                width = 0.2, size = 1, color = "black") +
  scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  labs(
    title = "Sialylation Percentage by Group (Mean ± 95% CI)",
    x = "Group",
    y = "Sialylation Percentage (%)"
  ) +
  # Add sample size annotations
  geom_text(aes(label = paste("n =", n)), 
            vjust = -0.5, size = 3.5) +
  # Add mean value annotations
  geom_text(aes(label = sprintf("Mean = %.2f%%", mean_neuac_percentage)), 
            vjust = -1.5, size = 3.5) +
  # Add CI annotations
  geom_text(aes(label = sprintf("95%% CI = [%.2f, %.2f]", ci_lower, ci_upper)), 
            vjust = -2.5, size = 3)

# Save the sialylation error plots
ggsave("figures/peptidegroups_intensity/sialylation_mean_se.pdf", 
       neuac_barplot_se, 
       width = 8, 
       height = 6, 
       dpi = 300)

ggsave("figures/peptidegroups_intensity/sialylation_mean_se.png", 
       neuac_barplot_se, 
       width = 8, 
       height = 6, 
       dpi = 300)

ggsave("figures/peptidegroups_intensity/sialylation_mean_ci.pdf", 
       neuac_barplot_ci, 
       width = 8, 
       height = 6, 
       dpi = 300)

ggsave("figures/peptidegroups_intensity/sialylation_mean_ci.png", 
       neuac_barplot_ci, 
       width = 8, 
       height = 6, 
       dpi = 300)

# Create a combined summary table for sialylation export
neuac_detailed_summary <- neuac_by_sample %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean = mean(neuac_percentage, na.rm = TRUE),
    median = median(neuac_percentage, na.rm = TRUE),
    sd = sd(neuac_percentage, na.rm = TRUE),
    se = sd / sqrt(n),
    ci_lower = mean - (1.96 * se),
    ci_upper = mean + (1.96 * se),
    min = min(neuac_percentage, na.rm = TRUE),
    max = max(neuac_percentage, na.rm = TRUE),
    .groups = 'drop'
  )

# Save the detailed sialylation summary
write.csv(neuac_detailed_summary, 
          "output_data/peptidegroups_intensity/sialylation_detailed_summary_with_se.csv", 
          row.names = FALSE)

cat("\nSialylation Standard Error Analysis Complete!\n")
cat("Created sialylation plots with error bars and saved detailed summary statistics.\n")

# Print summary statistics
cat("\nSialylation Summary Statistics:\n")
neuac_summary <- neuac_by_sample %>%
  group_by(group) %>%
  summarise(
    mean_neuac = mean(neuac_percentage),
    sd_neuac = sd(neuac_percentage),
    median_neuac = median(neuac_percentage),
    n_samples = n()
  )
print(neuac_summary)

cat("\nT-test Results:\n")
print(t_test_result_neuac)

cat("\nF-test for Variance (Homogeneity of Variance):\n")
cat(sprintf("Healthy group variance: %.3f\n", healthy_var_neuac))
cat(sprintf("MECFS group variance: %.3f\n", mecfs_var_neuac))
cat(sprintf("F-ratio: %.3f\n", f_ratio_neuac))
cat(sprintf("Degrees of freedom: %d, %d\n", df1_neuac, df2_neuac))
cat(sprintf("p-value: %.3g\n", f_p_value_neuac))

cat("\nCohen's d Effect Size:\n")
print(cohens_d_result_neuac)

# Interpret Cohen's d
cohens_d_interpretation_neuac <- case_when(
  abs(cohens_d_value_neuac) < 0.2 ~ "negligible",
  abs(cohens_d_value_neuac) < 0.5 ~ "small",
  abs(cohens_d_value_neuac) < 0.8 ~ "medium",
  TRUE ~ "large"
)

cat(sprintf("\nEffect Size Interpretation: Cohen's d = %.3f (%s effect)\n", 
            cohens_d_value_neuac, cohens_d_interpretation_neuac))

########################################################################################

# Calculate relative abundance of glycan classes by sample
glycan_class_by_sample <- glyco_peptide_groups_long %>%
  # Group by sample, group, and glycan_class to get abundances
  group_by(sample, group, glycan_class) %>%
  summarise(
    group_abundance = sum(abundance, na.rm = TRUE),
    .groups = 'keep'
  ) %>%
  # Calculate total abundance per sample for percentage
  group_by(sample, group) %>%
  mutate(
    total_sample_abundance = sum(group_abundance),
    relative_percentage = (group_abundance / total_sample_abundance) * 100
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
      p_value = t.test(relative_percentage ~ group)$p.value,
      .groups = 'drop'
    )
  
  # Create boxplot with individual points
  glycan_class_boxplot <- ggplot(glycan_class_by_sample, 
                       aes(x = glycan_class, 
                           y = relative_percentage, 
                           fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
    scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
    theme_minimal() +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    labs(
      title = "Relative Abundance of Glycan Classes by Group",
      x = "Glycan Class",
      y = "Relative Abundance (%)",
      fill = "Group"
    ) +
    # Add t-test p-values
    geom_text(data = t_test_results,
              aes(x = glycan_class,
                  y = max(glycan_class_by_sample$relative_percentage) + 5,
                  label = sprintf("p = %.3g", p_value)),
              inherit.aes = FALSE)
  
  # Save the plot
  ggsave("figures/peptidegroups_intensity/glycan_class_distribution_comparison.pdf", 
         glycan_class_boxplot, 
         width = 10, 
         height = 7, 
         dpi = 300)
  
  ggsave("figures/peptidegroups_intensity/glycan_class_distribution_comparison.png", 
         glycan_class_boxplot, 
         width = 10, 
         height = 7, 
         dpi = 300)
  
  # Print summary statistics
  cat("\nGlycan Class Distribution Summary Statistics:\n")
  glycan_class_summary <- glycan_class_by_sample %>%
    group_by(group, glycan_class) %>%
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

# Create separate boxplots for each glycan class with comprehensive statistics
# Complex glycans
complex_boxplot <- glycan_class_by_sample %>%
  filter(glycan_class == "Complex") %>%
  {
    # Calculate comprehensive statistics
    t_test <- t.test(relative_percentage ~ group, data = .)
    
    # F-test for variance
    healthy_var <- var(.$relative_percentage[.$group == "Healthy"])
    mecfs_var <- var(.$relative_percentage[.$group == "MECFS"])
    f_ratio <- max(healthy_var, mecfs_var) / min(healthy_var, mecfs_var)
    df1 <- sum(.$group == ifelse(healthy_var >= mecfs_var, "Healthy", "MECFS")) - 1
    df2 <- sum(.$group == ifelse(healthy_var >= mecfs_var, "MECFS", "Healthy")) - 1
    f_p_value <- 2 * (1 - pf(f_ratio, df1, df2))
    
    # Cohen's d
    cohens_d_result <- cohens_d(relative_percentage ~ group, data = .)
    cohens_d_value <- cohens_d_result$Cohens_d
    
    # Create labels
    t_test_label <- sprintf("t-test: p = %.3g", t_test$p.value)
    f_test_label <- sprintf("F-var = %.3f, p = %.3g", f_ratio, f_p_value)
    cohens_d_label <- sprintf("Cohen's d = %.3f", cohens_d_value)
    
    ggplot(., aes(x = group, 
               y = relative_percentage, 
               fill = group)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
      scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
      ) +
      labs(
        title = "Complex Glycans",
        x = "Group",
        y = "Relative Abundance (%)"
      ) +
      # Add statistical annotations
      annotate("text", 
               x = 1.5, 
               y = max(.$relative_percentage) * 1.15,
               label = t_test_label,
               size = 3) +
      annotate("text", 
               x = 1.5, 
               y = max(.$relative_percentage) * 1.05,
               label = f_test_label,
               size = 3) +
      annotate("text", 
               x = 1.5, 
               y = max(.$relative_percentage) * 0.95,
               label = cohens_d_label,
               size = 3)
  }

# Hybrid glycans
hybrid_boxplot <- glycan_class_by_sample %>%
  filter(glycan_class == "Hybrid") %>%
  {
    # Calculate comprehensive statistics
    t_test <- t.test(relative_percentage ~ group, data = .)
    
    # F-test for variance
    healthy_var <- var(.$relative_percentage[.$group == "Healthy"])
    mecfs_var <- var(.$relative_percentage[.$group == "MECFS"])
    f_ratio <- max(healthy_var, mecfs_var) / min(healthy_var, mecfs_var)
    df1 <- sum(.$group == ifelse(healthy_var >= mecfs_var, "Healthy", "MECFS")) - 1
    df2 <- sum(.$group == ifelse(healthy_var >= mecfs_var, "MECFS", "Healthy")) - 1
    f_p_value <- 2 * (1 - pf(f_ratio, df1, df2))
    
    # Cohen's d
    cohens_d_result <- cohens_d(relative_percentage ~ group, data = .)
    cohens_d_value <- cohens_d_result$Cohens_d
    
    # Create labels
    t_test_label <- sprintf("t-test: p = %.3g", t_test$p.value)
    f_test_label <- sprintf("F-var = %.3f, p = %.3g", f_ratio, f_p_value)
    cohens_d_label <- sprintf("Cohen's d = %.3f", cohens_d_value)
    
    ggplot(., aes(x = group, 
               y = relative_percentage, 
               fill = group)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
      scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
      ) +
      labs(
        title = "Hybrid Glycans",
        x = "Group",
        y = "Relative Abundance (%)"
      ) +
      # Add statistical annotations
      annotate("text", 
               x = 1.5, 
               y = max(.$relative_percentage) * 1.15,
               label = t_test_label,
               size = 3) +
      annotate("text", 
               x = 1.5, 
               y = max(.$relative_percentage) * 1.05,
               label = f_test_label,
               size = 3) +
      annotate("text", 
               x = 1.5, 
               y = max(.$relative_percentage) * 0.95,
               label = cohens_d_label,
               size = 3)
  }

# Oligomannose glycans
oligomannose_boxplot <- glycan_class_by_sample %>%
  filter(glycan_class == "Oligomannose") %>%
  {
    # Calculate comprehensive statistics
    t_test <- t.test(relative_percentage ~ group, data = .)
    
    # F-test for variance
    healthy_var <- var(.$relative_percentage[.$group == "Healthy"])
    mecfs_var <- var(.$relative_percentage[.$group == "MECFS"])
    f_ratio <- max(healthy_var, mecfs_var) / min(healthy_var, mecfs_var)
    df1 <- sum(.$group == ifelse(healthy_var >= mecfs_var, "Healthy", "MECFS")) - 1
    df2 <- sum(.$group == ifelse(healthy_var >= mecfs_var, "MECFS", "Healthy")) - 1
    f_p_value <- 2 * (1 - pf(f_ratio, df1, df2))
    
    # Cohen's d
    cohens_d_result <- cohens_d(relative_percentage ~ group, data = .)
    cohens_d_value <- cohens_d_result$Cohens_d
    
    # Create labels
    t_test_label <- sprintf("t-test: p = %.3g", t_test$p.value)
    f_test_label <- sprintf("F-var = %.3f, p = %.3g", f_ratio, f_p_value)
    cohens_d_label <- sprintf("Cohen's d = %.3f", cohens_d_value)
    
    ggplot(., aes(x = group, 
               y = relative_percentage, 
               fill = group)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
      scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
      ) +
      labs(
        title = "Oligomannose Glycans",
        x = "Group",
        y = "Relative Abundance (%)"
      ) +
      # Add statistical annotations
      annotate("text", 
               x = 1.5, 
               y = max(.$relative_percentage) * 1.15,
               label = t_test_label,
               size = 3) +
      annotate("text", 
               x = 1.5, 
               y = max(.$relative_percentage) * 1.05,
               label = f_test_label,
               size = 3) +
      annotate("text", 
               x = 1.5, 
               y = max(.$relative_percentage) * 0.95,
               label = cohens_d_label,
               size = 3)
  }

# Rest of the code remains the same
library(patchwork)
combined_plot <- complex_boxplot + hybrid_boxplot + oligomannose_boxplot +
  plot_layout(ncol = 3)

# Save individual plots
ggsave("figures/peptidegroups_intensity/complex_glycans_boxplot.pdf", complex_boxplot, width = 6, height = 5, dpi = 300)
ggsave("figures/peptidegroups_intensity/complex_glycans_boxplot.png", complex_boxplot, width = 6, height = 5, dpi = 300)

ggsave("figures/peptidegroups_intensity/hybrid_glycans_boxplot.pdf", hybrid_boxplot, width = 6, height = 5, dpi = 300)
ggsave("figures/peptidegroups_intensity/hybrid_glycans_boxplot.png", hybrid_boxplot, width = 6, height = 5, dpi = 300)

ggsave("figures/peptidegroups_intensity/oligomannose_glycans_boxplot.pdf", oligomannose_boxplot, width = 6, height = 5, dpi = 300)
ggsave("figures/peptidegroups_intensity/oligomannose_glycans_boxplot.png", oligomannose_boxplot, width = 6, height = 5, dpi = 300)

# Save combined plot
ggsave("figures/peptidegroups_intensity/all_glycan_classes_boxplots.pdf", combined_plot, width = 18, height = 5, dpi = 300)
ggsave("figures/peptidegroups_intensity/all_glycan_classes_boxplots.png", combined_plot, width = 18, height = 5, dpi = 300)

# Print comprehensive summary statistics for each class
for(class in unique(glycan_class_by_sample$glycan_class)) {
  cat(sprintf("\n=== Summary Statistics for %s Glycans ===\n", class))
  
  # Basic summary statistics
  class_summary <- glycan_class_by_sample %>%
    filter(glycan_class == class) %>%
    group_by(group) %>%
    summarise(
      mean_percentage = mean(relative_percentage),
      sd_percentage = sd(relative_percentage),
      median_percentage = median(relative_percentage),
      n_samples = n(),
      .groups = 'drop'
    )
  print(class_summary)
  
  # T-test results
  class_data <- filter(glycan_class_by_sample, glycan_class == class)
  t_test <- t.test(relative_percentage ~ group, data = class_data)
  cat(sprintf("\nT-test Results for %s Glycans:\n", class))
  print(t_test)
  
  # F-test for variance
  healthy_var <- var(class_data$relative_percentage[class_data$group == "Healthy"])
  mecfs_var <- var(class_data$relative_percentage[class_data$group == "MECFS"])
  f_ratio <- max(healthy_var, mecfs_var) / min(healthy_var, mecfs_var)
  df1 <- sum(class_data$group == ifelse(healthy_var >= mecfs_var, "Healthy", "MECFS")) - 1
  df2 <- sum(class_data$group == ifelse(healthy_var >= mecfs_var, "MECFS", "Healthy")) - 1
  f_p_value <- 2 * (1 - pf(f_ratio, df1, df2))
  
  cat(sprintf("\nF-test for Variance (Homogeneity of Variance):\n"))
  cat(sprintf("Healthy group variance: %.3f\n", healthy_var))
  cat(sprintf("MECFS group variance: %.3f\n", mecfs_var))
  cat(sprintf("F-ratio: %.3f\n", f_ratio))
  cat(sprintf("Degrees of freedom: %d, %d\n", df1, df2))
  cat(sprintf("p-value: %.3g\n", f_p_value))
  
  # Cohen's d effect size
  cohens_d_result <- cohens_d(relative_percentage ~ group, data = class_data)
  cohens_d_value <- cohens_d_result$Cohens_d
  
  cat(sprintf("\nCohen's d Effect Size:\n"))
  print(cohens_d_result)
  
  # Interpret Cohen's d
  cohens_d_interpretation <- case_when(
    abs(cohens_d_value) < 0.2 ~ "negligible",
    abs(cohens_d_value) < 0.5 ~ "small",
    abs(cohens_d_value) < 0.8 ~ "medium",
    TRUE ~ "large"
  )
  
  cat(sprintf("Effect Size Interpretation: Cohen's d = %.3f (%s effect)\n", 
              cohens_d_value, cohens_d_interpretation))
  cat("\n" + strrep("=", 50) + "\n")
}

# Create multi-panel figure using patchwork
library(patchwork)
# Create each row separately first
row1 <- top_glycans_plot + plot_spacer()
row2 <- fuc_boxplot + neuac_boxplot
row3 <- complex_boxplot + hybrid_boxplot + oligomannose_boxplot

# Then combine them
combined_plot <- row1 / row2 / row3 +
  plot_annotation(
    title = "Glycosylation Patterns in ME/CFS vs Healthy Controls",
    theme = theme(
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
      plot.margin = margin(20, 20, 20, 20)
    )
  ) +
  plot_layout(
    ncol = 2,
    nrow = 3,
    heights = c(1.2, 1, 1),
    widths = c(1, 1),
    guides = "collect"
  )



ggsave(
  "figures/peptidegroups_intensity/glycosylation_patterns_combined.png", 
  combined_plot,
  width = 15,
  height = 18,
  dpi = 300
)



###########################################################################################################
## Protein-level fucosylation and sialylation analysis

# Function to analyze protein-level glycosylation patterns
# Uses Isometric Log Ratio (ILR) transformation for compositional data analysis
# ILR transformation: log(proportion_fucosylated / proportion_non_fucosylated)
# This transformation makes the data suitable for standard statistical tests
analyze_protein_glycosylation <- function(data, output_dir = "output_data/peptidegroups_intensity/protein_level", 
                                         figures_dir = "figures/protein_level") {
  
  # Load required libraries
  library(dplyr)
  library(ggplot2)
  
  # Create output directories if they don't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Split data by protein accession
  protein_list <- split(data, data$protein_accessions)
  
  # Initialize results storage
  fuc_results <- list()
  neuac_results <- list()
  summary_stats <- data.frame()
  
  # Track proteins for summary
  proteins_analyzed <- 0
  proteins_insufficient_data <- 0
  proteins_insufficient_groups <- 0
  proteins_with_fuc_data <- 0
  proteins_with_neuac_data <- 0
  fuc_tests_performed <- 0
  neuac_tests_performed <- 0
  
  # Track failed tests
  failed_fuc_tests <- data.frame(
    protein_accession = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
  
  failed_neuac_tests <- data.frame(
    protein_accession = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
  
  cat("Analyzing", length(protein_list), "proteins for fucosylation and sialylation patterns...\n")
  
  # Process each protein
  for(i in seq_along(protein_list)) {
    protein_name <- names(protein_list)[i]
    protein_data <- protein_list[[i]]
    
    cat(sprintf("\nProcessing protein %d/%d: %s\n", i, length(protein_list), protein_name))
    
    # Skip proteins with insufficient data
    if(n_distinct(protein_data$sample) < 3) {
      cat("  Skipping - insufficient total samples\n")
      proteins_insufficient_data <- proteins_insufficient_data + 1
      next
    }
    
    # Check if we have at least 3 samples in each group
    sample_counts <- protein_data %>%
      group_by(group) %>%
      summarise(n_samples = n_distinct(sample), .groups = 'drop')
    
    if(any(sample_counts$n_samples < 3)) {
      cat("  Skipping - insufficient samples in each group (need ≥3 per group)\n")
      cat("    Sample counts:", paste(sample_counts$group, "=", sample_counts$n_samples, collapse = ", "), "\n")
      proteins_insufficient_groups <- proteins_insufficient_groups + 1
      next
    }
    
    proteins_analyzed <- proteins_analyzed + 1
    
    # Calculate fucosylation percentages by sample
    fuc_by_sample <- protein_data %>%
      group_by(sample, group) %>%
      summarise(
        total_abundance = sum(abundance, na.rm = TRUE),
        fuc_abundance = sum(abundance[contains_Fuc == TRUE], na.rm = TRUE),
        fuc_percentage = (fuc_abundance / total_abundance) * 100,
        .groups = 'drop'
      )
    
    # Apply ILR transformation to fucosylation data
    # ILR requires compositional data (proportions that sum to 1)
    fuc_by_sample <- fuc_by_sample %>%
      mutate(
        # Convert percentages to proportions
        fuc_proportion = fuc_percentage / 100,
        non_fuc_proportion = 1 - fuc_proportion,
        # Apply ILR transformation: log(fuc_proportion / non_fuc_proportion)
        fuc_ilr = log(fuc_proportion / non_fuc_proportion)
      )
    
    # Calculate sialylation percentages by sample
    neuac_by_sample <- protein_data %>%
      group_by(sample, group) %>%
      summarise(
        total_abundance = sum(abundance, na.rm = TRUE),
        neuac_abundance = sum(abundance[contains_NeuAc == TRUE], na.rm = TRUE),
        neuac_percentage = (neuac_abundance / total_abundance) * 100,
        .groups = 'drop'
      )
    
    # Apply ILR transformation to sialylation data
    neuac_by_sample <- neuac_by_sample %>%
      mutate(
        # Convert percentages to proportions
        neuac_proportion = neuac_percentage / 100,
        non_neuac_proportion = 1 - neuac_proportion,
        # Apply ILR transformation: log(neuac_proportion / non_neuac_proportion)
        # Handle edge cases for ILR transformation
        neuac_proportion_adj = ifelse(neuac_proportion == 0, 0.001, 
                                     ifelse(neuac_proportion == 1, 0.999, neuac_proportion)),
        non_neuac_proportion_adj = 1 - neuac_proportion_adj,
        neuac_ilr = log(neuac_proportion_adj / non_neuac_proportion_adj)
      ) %>%
      # Remove infinite or NaN values
      filter(is.finite(neuac_ilr))
    

    
    # Check if protein has sufficient fucosylation data (≥3 samples per group with fucosylation)
    fuc_sample_counts <- fuc_by_sample %>%
      filter(fuc_percentage > 0) %>%  # Only samples with some fucosylation
      group_by(group) %>%
      summarise(n_samples = n(), .groups = 'drop')
    
    has_sufficient_fuc <- nrow(fuc_sample_counts) >= 2 && all(fuc_sample_counts$n_samples >= 3)
    
    # Check if protein has sufficient sialylation data (≥3 samples per group with sialylation)
    neuac_sample_counts <- neuac_by_sample %>%
      filter(neuac_percentage > 0) %>%  # Only samples with some sialylation
      group_by(group) %>%
      summarise(n_samples = n(), .groups = 'drop')
    
    has_sufficient_neuac <- nrow(neuac_sample_counts) >= 2 && all(neuac_sample_counts$n_samples >= 3)
    
    # Store results only if sufficient data exists
    if(has_sufficient_fuc) {
      fuc_results[[protein_name]] <- fuc_by_sample
      proteins_with_fuc_data <- proteins_with_fuc_data + 1
    } else {
      cat("    No sufficient fucosylation data (need ≥3 samples per group with fucosylation)\n")
    }
    if(has_sufficient_neuac) {
      neuac_results[[protein_name]] <- neuac_by_sample
      proteins_with_neuac_data <- proteins_with_neuac_data + 1
    } else {
      cat("    No sufficient sialylation data (need ≥3 samples per group with sialylation)\n")
    }
    
    # Perform statistical tests for fucosylation
    if(has_sufficient_fuc) {
      
      # Check for constant data within groups
      healthy_fuc_data <- fuc_by_sample$fuc_percentage[fuc_by_sample$group == "Healthy"]
      mecfs_fuc_data <- fuc_by_sample$fuc_percentage[fuc_by_sample$group == "MECFS"]
      
      # Check if data is constant within groups
      if(length(unique(healthy_fuc_data)) <= 1 && length(unique(mecfs_fuc_data)) <= 1) {
        failed_fuc_tests <- rbind(failed_fuc_tests, data.frame(
          protein_accession = protein_name,
          reason = "Constant data within both groups",
          stringsAsFactors = FALSE
        ))
        cat("    Fucosylation test failed: Constant data within both groups\n")
      } else if(length(unique(healthy_fuc_data)) <= 1) {
        failed_fuc_tests <- rbind(failed_fuc_tests, data.frame(
          protein_accession = protein_name,
          reason = "Constant data in Healthy group",
          stringsAsFactors = FALSE
        ))
        cat("    Fucosylation test failed: Constant data in Healthy group\n")
      } else if(length(unique(mecfs_fuc_data)) <= 1) {
        failed_fuc_tests <- rbind(failed_fuc_tests, data.frame(
          protein_accession = protein_name,
          reason = "Constant data in MECFS group",
          stringsAsFactors = FALSE
        ))
        cat("    Fucosylation test failed: Constant data in MECFS group\n")
      } else {
        # F-test for variance
        fuc_f_test <- try({
          var.test(fuc_percentage ~ group, data = fuc_by_sample)
        }, silent = TRUE)
        
        # T-test
        fuc_t_test <- try({
          t.test(fuc_percentage ~ group, data = fuc_by_sample)
        }, silent = TRUE)
        
        if(!inherits(fuc_t_test, "try-error")) {
          fuc_stats <- data.frame(
            protein_accession = protein_name,
            feature = "Fucosylation",
            mean_healthy = mean(healthy_fuc_data, na.rm = TRUE),
            mean_mecfs = mean(mecfs_fuc_data, na.rm = TRUE),
            sd_healthy = sd(healthy_fuc_data, na.rm = TRUE),
            sd_mecfs = sd(mecfs_fuc_data, na.rm = TRUE),
            t_stat = fuc_t_test$statistic,
            p_value = fuc_t_test$p.value,
            f_stat = ifelse(!inherits(fuc_f_test, "try-error"), fuc_f_test$statistic, NA),
            f_p_value = ifelse(!inherits(fuc_f_test, "try-error"), fuc_f_test$p.value, NA),
            n_healthy = sum(fuc_by_sample$group == "Healthy"),
            n_mecfs = sum(fuc_by_sample$group == "MECFS"),
            stringsAsFactors = FALSE
          )
          summary_stats <- rbind(summary_stats, fuc_stats)
          fuc_tests_performed <- fuc_tests_performed + 1
        } else {
          failed_fuc_tests <- rbind(failed_fuc_tests, data.frame(
            protein_accession = protein_name,
            reason = paste("T-test error:", fuc_t_test[1]),
            stringsAsFactors = FALSE
          ))
          cat("    Fucosylation t-test failed:", fuc_t_test[1], "\n")
        }
      }
    }
    
    # Perform statistical tests for sialylation
    if(has_sufficient_neuac) {
      
      # Check for constant data within groups
      healthy_neuac_data <- neuac_by_sample$neuac_percentage[neuac_by_sample$group == "Healthy"]
      mecfs_neuac_data <- neuac_by_sample$neuac_percentage[neuac_by_sample$group == "MECFS"]
      
      # Check if data is constant within groups
      if(length(unique(healthy_neuac_data)) <= 1 && length(unique(mecfs_neuac_data)) <= 1) {
        failed_neuac_tests <- rbind(failed_neuac_tests, data.frame(
          protein_accession = protein_name,
          reason = "Constant data within both groups",
          stringsAsFactors = FALSE
        ))
        cat("    Sialylation test failed: Constant data within both groups\n")
      } else if(length(unique(healthy_neuac_data)) <= 1) {
        failed_neuac_tests <- rbind(failed_neuac_tests, data.frame(
          protein_accession = protein_name,
          reason = "Constant data in Healthy group",
          stringsAsFactors = FALSE
        ))
        cat("    Sialylation test failed: Constant data in Healthy group\n")
      } else if(length(unique(mecfs_neuac_data)) <= 1) {
        failed_neuac_tests <- rbind(failed_neuac_tests, data.frame(
          protein_accession = protein_name,
          reason = "Constant data in MECFS group",
          stringsAsFactors = FALSE
        ))
        cat("    Sialylation test failed: Constant data in MECFS group\n")
      } else {
        # F-test for variance
        neuac_f_test <- try({
          var.test(neuac_percentage ~ group, data = neuac_by_sample)
        }, silent = TRUE)
        
        # T-test
        neuac_t_test <- try({
          t.test(neuac_percentage ~ group, data = neuac_by_sample)
        }, silent = TRUE)
        
        if(!inherits(neuac_t_test, "try-error")) {
          neuac_stats <- data.frame(
            protein_accession = protein_name,
            feature = "Sialylation",
            mean_healthy = mean(healthy_neuac_data, na.rm = TRUE),
            mean_mecfs = mean(mecfs_neuac_data, na.rm = TRUE),
            sd_healthy = sd(healthy_neuac_data, na.rm = TRUE),
            sd_mecfs = sd(mecfs_neuac_data, na.rm = TRUE),
            t_stat = neuac_t_test$statistic,
            p_value = neuac_t_test$p.value,
                       f_stat = ifelse(!inherits(neuac_f_test, "try-error"), neuac_f_test$statistic, NA),
             f_p_value = ifelse(!inherits(neuac_f_test, "try-error"), neuac_f_test$p.value, NA),
            n_healthy = sum(neuac_by_sample$group == "Healthy"),
            n_mecfs = sum(neuac_by_sample$group == "MECFS"),
            stringsAsFactors = FALSE
          )
          summary_stats <- rbind(summary_stats, neuac_stats)
          neuac_tests_performed <- neuac_tests_performed + 1
        } else {
          failed_neuac_tests <- rbind(failed_neuac_tests, data.frame(
            protein_accession = protein_name,
            reason = paste("T-test error:", neuac_t_test[1]),
            stringsAsFactors = FALSE
          ))
          cat("    Sialylation t-test failed:", neuac_t_test[1], "\n")
        }
      }
    }
  }
  
  # Apply BH correction to all p-values
  if(nrow(summary_stats) > 0) {
    summary_stats$p_value_adj <- p.adjust(summary_stats$p_value, method = "BH")
    summary_stats$significant <- summary_stats$p_value_adj < 0.05
    summary_stats$fold_change <- summary_stats$mean_mecfs / summary_stats$mean_healthy
  }
  
  # Create separate boxplots for fucosylation and sialylation
  cat("\nCreating separate boxplots for fucosylation and sialylation...\n")
  
  for(protein_name in names(protein_list)) {
    safe_protein_name <- gsub("[^a-zA-Z0-9]", "_", protein_name)
    
    # Create fucosylation boxplot if sufficient data exists
    if(protein_name %in% names(fuc_results)) {
      fuc_plot <- ggplot(fuc_results[[protein_name]], 
                        aes(x = group, y = fuc_percentage, fill = group)) +
        geom_boxplot(alpha = 0.7, outlier.shape = 1) +
        geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
        scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
        labs(title = paste("Fucosylation -", protein_name),
             x = "Group", 
             y = "Fucosylation (%)") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 14),
              legend.position = "top",
              axis.text = element_text(size = 12))
      
      # Save fucosylation plot
      ggsave(paste0(figures_dir, "/", safe_protein_name, "_fucosylation.png"), 
             fuc_plot, width = 8, height = 6, dpi = 300)
    }
    
    # Create sialylation boxplot if sufficient data exists
    if(protein_name %in% names(neuac_results)) {
      neuac_plot <- ggplot(neuac_results[[protein_name]], 
                          aes(x = group, y = neuac_percentage, fill = group)) +
        geom_boxplot(alpha = 0.7, outlier.shape = 1) +
        geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
        scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
        labs(title = paste("Sialylation -", protein_name),
             x = "Group", 
             y = "Sialylation (%)") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 14),
              legend.position = "top",
              axis.text = element_text(size = 12))
      
      # Save sialylation plot
      ggsave(paste0(figures_dir, "/", safe_protein_name, "_sialylation.png"), 
             neuac_plot, width = 8, height = 6, dpi = 300)
    }
  }
  
  # Save results
  saveRDS(fuc_results, paste0(output_dir, "/fucosylation_by_protein.rds"))
  saveRDS(neuac_results, paste0(output_dir, "/sialylation_by_protein.rds"))
  write.csv(summary_stats, paste0(output_dir, "/protein_glycosylation_statistics.csv"), row.names = FALSE)
  
  # Create separate summary files for fucosylation and sialylation
  fuc_summary <- summary_stats %>% filter(feature == "Fucosylation")
  neuac_summary <- summary_stats %>% filter(feature == "Sialylation")
  
  write.csv(fuc_summary, paste0(output_dir, "/fucosylation_statistics.csv"), row.names = FALSE)
  write.csv(neuac_summary, paste0(output_dir, "/sialylation_statistics.csv"), row.names = FALSE)
  
  # Create detailed summary dataframes for each feature
  fuc_detailed <- bind_rows(lapply(names(fuc_results), function(protein) {
    data <- fuc_results[[protein]]
    data$protein_accession <- protein
    return(data)
  }))
  
  neuac_detailed <- bind_rows(lapply(names(neuac_results), function(protein) {
    data <- neuac_results[[protein]]
    data$protein_accession <- protein
    return(data)
  }))
  
  write.csv(fuc_detailed, paste0(output_dir, "/fucosylation_detailed_data.csv"), row.names = FALSE)
  write.csv(neuac_detailed, paste0(output_dir, "/sialylation_detailed_data.csv"), row.names = FALSE)
  
  # Save failed test information
  write.csv(failed_fuc_tests, paste0(output_dir, "/failed_fucosylation_tests.csv"), row.names = FALSE)
  write.csv(failed_neuac_tests, paste0(output_dir, "/failed_sialylation_tests.csv"), row.names = FALSE)
  
  # Print summary
  cat("\n=== PROTEIN-LEVEL GLYCOSYLATION ANALYSIS SUMMARY ===\n")
  cat("Total proteins in dataset:", length(protein_list), "\n")
  cat("Proteins with sufficient data (≥3 samples per group):", proteins_analyzed, "\n")
  cat("Proteins excluded - insufficient total samples:", proteins_insufficient_data, "\n")
  cat("Proteins excluded - insufficient samples per group:", proteins_insufficient_groups, "\n")
  cat("Proteins with sufficient fucosylation data (≥3 samples per group with fucosylation):", proteins_with_fuc_data, "\n")
  cat("Proteins with sufficient sialylation data (≥3 samples per group with sialylation):", proteins_with_neuac_data, "\n")
  cat("Fucosylation statistical tests performed:", fuc_tests_performed, "\n")
  cat("Sialylation statistical tests performed:", neuac_tests_performed, "\n")
  cat("Total statistical comparisons with results:", nrow(summary_stats), "\n")
  
  if(nrow(summary_stats) > 0) {
    cat("Significant differences (FDR < 0.05):", sum(summary_stats$significant, na.rm = TRUE), "\n")
    
    # Fucosylation summary
    fuc_stats <- summary_stats %>% filter(feature == "Fucosylation")
    cat(sprintf("\n=== FUCOSYLATION ANALYSIS ===\n"))
    cat("Proteins with sufficient fucosylation data:", proteins_with_fuc_data, "\n")
    cat("Proteins with successful fucosylation statistical tests:", nrow(fuc_stats), "\n")
    cat("Significant fucosylation differences:", sum(fuc_stats$significant, na.rm = TRUE), "\n")
    
    if(nrow(fuc_stats) > 0) {
      cat("Mean fucosylation by group:\n")
      fuc_means <- fuc_stats %>%
        summarise(
          mean_healthy = mean(mean_healthy, na.rm = TRUE),
          mean_mecfs = mean(mean_mecfs, na.rm = TRUE),
          .groups = 'drop'
        )
      cat(sprintf("  Healthy: %.2f%%\n", fuc_means$mean_healthy))
      cat(sprintf("  MECFS: %.2f%%\n", fuc_means$mean_mecfs))
    }
    
    # Significant fucosylation differences
    sig_fuc <- fuc_stats %>% filter(significant)
    if(nrow(sig_fuc) > 0) {
      cat("\nSignificant fucosylation differences:\n")
      for(i in 1:nrow(sig_fuc)) {
        cat(sprintf("  %s: p = %.3e, FC = %.2f\n", 
                   sig_fuc$protein_accession[i], 
                   sig_fuc$p_value_adj[i], 
                   sig_fuc$fold_change[i]))
      }
    }
    
    # Sialylation summary
    neuac_stats <- summary_stats %>% filter(feature == "Sialylation")
    cat(sprintf("\n=== SIALYLATION ANALYSIS ===\n"))
    cat("Proteins with sufficient sialylation data:", proteins_with_neuac_data, "\n")
    cat("Proteins with successful sialylation statistical tests:", nrow(neuac_stats), "\n")
    cat("Significant sialylation differences:", sum(neuac_stats$significant, na.rm = TRUE), "\n")
    
    if(nrow(neuac_stats) > 0) {
      cat("Mean sialylation by group:\n")
      neuac_means <- neuac_stats %>%
        summarise(
          mean_healthy = mean(mean_healthy, na.rm = TRUE),
          mean_mecfs = mean(mean_mecfs, na.rm = TRUE),
          .groups = 'drop'
        )
      cat(sprintf("  Healthy: %.2f%%\n", neuac_means$mean_healthy))
      cat(sprintf("  MECFS: %.2f%%\n", neuac_means$mean_mecfs))
    }
    
    # Significant sialylation differences
    sig_neuac <- neuac_stats %>% filter(significant)
    if(nrow(sig_neuac) > 0) {
      cat("\nSignificant sialylation differences:\n")
      for(i in 1:nrow(sig_neuac)) {
        cat(sprintf("  %s: p = %.3e, FC = %.2f\n", 
                   sig_neuac$protein_accession[i], 
                   sig_neuac$p_value_adj[i], 
                   sig_neuac$fold_change[i]))
      }
    }
  }
  
  cat("\nResults saved to:", output_dir, "\n")
  cat("  - fucosylation_statistics.csv: Fucosylation statistical results\n")
  cat("  - sialylation_statistics.csv: Sialylation statistical results\n")
  cat("  - fucosylation_detailed_data.csv: Detailed fucosylation data by protein and sample\n")
  cat("  - sialylation_detailed_data.csv: Detailed sialylation data by protein and sample\n")
  cat("  - protein_glycosylation_statistics.csv: Combined statistical results\n")
  cat("  - fucosylation_by_protein.rds: Fucosylation data as R object\n")
  cat("  - sialylation_by_protein.rds: Sialylation data as R object\n")
  cat("  - failed_fucosylation_tests.csv: Proteins with failed fucosylation tests\n")
  cat("  - failed_sialylation_tests.csv: Proteins with failed sialylation tests\n")
  cat("\nFigures saved to:", figures_dir, "\n")
  cat("  - *_fucosylation.png: Individual fucosylation boxplots\n")
  cat("  - *_sialylation.png: Individual sialylation boxplots\n")
  
  # Report failed tests
  if(nrow(failed_fuc_tests) > 0) {
    cat("\n=== FAILED FUCOSYLATION TESTS ===\n")
    cat("Total failed fucosylation tests:", nrow(failed_fuc_tests), "\n")
    
    # Count reasons
    fuc_reason_counts <- table(failed_fuc_tests$reason)
    cat("Reasons for failure:\n")
    for(i in 1:length(fuc_reason_counts)) {
      cat(sprintf("  %s: %d proteins\n", names(fuc_reason_counts)[i], fuc_reason_counts[i]))
    }
    
    cat("\nProteins with failed fucosylation tests:\n")
    for(i in 1:nrow(failed_fuc_tests)) {
      cat(sprintf("  %s: %s\n", failed_fuc_tests$protein_accession[i], failed_fuc_tests$reason[i]))
    }
  }
  
  if(nrow(failed_neuac_tests) > 0) {
    cat("\n=== FAILED SIALYLATION TESTS ===\n")
    cat("Total failed sialylation tests:", nrow(failed_neuac_tests), "\n")
    
    # Count reasons
    neuac_reason_counts <- table(failed_neuac_tests$reason)
    cat("Reasons for failure:\n")
    for(i in 1:length(neuac_reason_counts)) {
      cat(sprintf("  %s: %d proteins\n", names(neuac_reason_counts)[i], neuac_reason_counts[i]))
    }
    
    cat("\nProteins with failed sialylation tests:\n")
    for(i in 1:nrow(failed_neuac_tests)) {
      cat(sprintf("  %s: %s\n", failed_neuac_tests$protein_accession[i], failed_neuac_tests$reason[i]))
    }
  }
  
  return(list(
    fucosylation_data = fuc_results,
    sialylation_data = neuac_results,
    statistics = summary_stats,
    failed_fucosylation_tests = failed_fuc_tests,
    failed_sialylation_tests = failed_neuac_tests
  ))
}

# Add group classification to glyco_peptide_groups_long if not already present
if (!"group" %in% colnames(glyco_peptide_groups_long)) {
  glyco_peptide_groups_long <- glyco_peptide_groups_long %>%
    mutate(group = case_when(
      str_starts(sample, "hc") ~ "Healthy",
      str_starts(sample, "m") ~ "MECFS",
      TRUE ~ "Unknown"
    ))
}

# Run the protein-level analysis
protein_glycosylation_results <- analyze_protein_glycosylation(glyco_peptide_groups_long)

# Additional ILR-based analysis function
analyze_protein_glycosylation_ilr <- function(data, original_results = NULL, output_dir = "output_data/peptidegroups_intensity/protein_level_ilr", 
                                             figures_dir = "figures/protein_level_ilr") {
  
  # Load required libraries
  library(dplyr)
  library(ggplot2)
  
  # Create output directories if they don't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Split data by protein accession
  protein_list <- split(data, data$protein_accessions)
  
  # Initialize results storage
  ilr_summary_stats <- data.frame()
  
  cat("Performing ILR-based statistical analysis...\n")
  
  # Analyze each protein with ILR transformation
  for(protein_name in names(protein_list)) {
    protein_data <- protein_list[[protein_name]]
    
    # Skip if insufficient data
    if(n_distinct(protein_data$sample) < 3) next
    
    # Check group sample sizes
    sample_counts <- protein_data %>%
      group_by(group) %>%
      summarise(n_samples = n_distinct(sample), .groups = 'drop')
    
    if(any(sample_counts$n_samples < 3)) next
    
    # Calculate fucosylation and sialylation with ILR
    fuc_by_sample <- protein_data %>%
      group_by(sample, group) %>%
      summarise(
        total_abundance = sum(abundance, na.rm = TRUE),
        fuc_abundance = sum(abundance[contains_Fuc == TRUE], na.rm = TRUE),
        fuc_percentage = (fuc_abundance / total_abundance) * 100,
        .groups = 'drop'
      ) %>%
      mutate(
        fuc_proportion = fuc_percentage / 100,
        non_fuc_proportion = 1 - fuc_proportion,
        # Handle edge cases for ILR transformation
        fuc_proportion_adj = ifelse(fuc_proportion == 0, 0.001, 
                                   ifelse(fuc_proportion == 1, 0.999, fuc_proportion)),
        non_fuc_proportion_adj = 1 - fuc_proportion_adj,
        fuc_ilr = log(fuc_proportion_adj / non_fuc_proportion_adj)
      ) %>%
      # Remove infinite or NaN values
      filter(is.finite(fuc_ilr))
    
    neuac_by_sample <- protein_data %>%
      group_by(sample, group) %>%
      summarise(
        total_abundance = sum(abundance, na.rm = TRUE),
        neuac_abundance = sum(abundance[contains_NeuAc == TRUE], na.rm = TRUE),
        neuac_percentage = (neuac_abundance / total_abundance) * 100,
        .groups = 'drop'
      ) %>%
      mutate(
        neuac_proportion = neuac_percentage / 100,
        non_neuac_proportion = 1 - neuac_proportion,
        # Handle edge cases for ILR transformation
        neuac_proportion_adj = ifelse(neuac_proportion == 0, 0.001, 
                                     ifelse(neuac_proportion == 1, 0.999, neuac_proportion)),
        non_neuac_proportion_adj = 1 - neuac_proportion_adj,
        neuac_ilr = log(neuac_proportion_adj / non_neuac_proportion_adj)
      ) %>%
      # Remove infinite or NaN values
      filter(is.finite(neuac_ilr))
    
    # Check if sufficient data for ILR analysis
    fuc_sample_counts <- fuc_by_sample %>%
      group_by(group) %>%
      summarise(n_samples = n(), .groups = 'drop')
    
    neuac_sample_counts <- neuac_by_sample %>%
      group_by(group) %>%
      summarise(n_samples = n(), .groups = 'drop')
    
    has_sufficient_fuc_ilr <- nrow(fuc_sample_counts) >= 2 && all(fuc_sample_counts$n_samples >= 3)
    has_sufficient_neuac_ilr <- nrow(neuac_sample_counts) >= 2 && all(neuac_sample_counts$n_samples >= 3)
    
    # ILR-based fucosylation analysis
    if(has_sufficient_fuc_ilr) {
      healthy_fuc_ilr_data <- fuc_by_sample$fuc_ilr[fuc_by_sample$group == "Healthy"]
      mecfs_fuc_ilr_data <- fuc_by_sample$fuc_ilr[fuc_by_sample$group == "MECFS"]
      
      if(length(unique(healthy_fuc_ilr_data)) > 1 || length(unique(mecfs_fuc_ilr_data)) > 1) {
        fuc_ilr_t_test <- try({
          t.test(fuc_ilr ~ group, data = fuc_by_sample)
        }, silent = TRUE)
        
        if(!inherits(fuc_ilr_t_test, "try-error")) {
          fuc_ilr_stats <- data.frame(
            protein_accession = protein_name,
            feature = "Fucosylation_ILR",
            mean_healthy = mean(healthy_fuc_ilr_data, na.rm = TRUE),
            mean_mecfs = mean(mecfs_fuc_ilr_data, na.rm = TRUE),
            sd_healthy = sd(healthy_fuc_ilr_data, na.rm = TRUE),
            sd_mecfs = sd(mecfs_fuc_ilr_data, na.rm = TRUE),
            t_stat = fuc_ilr_t_test$statistic,
            p_value = fuc_ilr_t_test$p.value,
            n_healthy = sum(fuc_by_sample$group == "Healthy"),
            n_mecfs = sum(fuc_by_sample$group == "MECFS"),
            stringsAsFactors = FALSE
          )
          ilr_summary_stats <- rbind(ilr_summary_stats, fuc_ilr_stats)
        }
      }
    }
    
    # ILR-based sialylation analysis
    if(has_sufficient_neuac_ilr) {
      healthy_neuac_ilr_data <- neuac_by_sample$neuac_ilr[neuac_by_sample$group == "Healthy"]
      mecfs_neuac_ilr_data <- neuac_by_sample$neuac_ilr[neuac_by_sample$group == "MECFS"]
      
      if(length(unique(healthy_neuac_ilr_data)) > 1 || length(unique(mecfs_neuac_ilr_data)) > 1) {
        neuac_ilr_t_test <- try({
          t.test(neuac_ilr ~ group, data = neuac_by_sample)
        }, silent = TRUE)
        
        if(!inherits(neuac_ilr_t_test, "try-error")) {
          neuac_ilr_stats <- data.frame(
            protein_accession = protein_name,
            feature = "Sialylation_ILR",
            mean_healthy = mean(healthy_neuac_ilr_data, na.rm = TRUE),
            mean_mecfs = mean(mecfs_neuac_ilr_data, na.rm = TRUE),
            sd_healthy = sd(healthy_neuac_ilr_data, na.rm = TRUE),
            sd_mecfs = sd(mecfs_neuac_ilr_data, na.rm = TRUE),
            t_stat = neuac_ilr_t_test$statistic,
            p_value = neuac_ilr_t_test$p.value,
            n_healthy = sum(neuac_by_sample$group == "Healthy"),
            n_mecfs = sum(neuac_by_sample$group == "MECFS"),
            stringsAsFactors = FALSE
          )
          ilr_summary_stats <- rbind(ilr_summary_stats, neuac_ilr_stats)
        }
      }
    }
  }
  
  # Apply BH correction to ILR p-values
  if(nrow(ilr_summary_stats) > 0) {
    ilr_summary_stats$p_value_adj <- p.adjust(ilr_summary_stats$p_value, method = "BH")
    ilr_summary_stats$significant <- ilr_summary_stats$p_value_adj < 0.05
    ilr_summary_stats$fold_change <- ilr_summary_stats$mean_mecfs / ilr_summary_stats$mean_healthy
  }
  
  # Save ILR results
  write.csv(ilr_summary_stats, file.path(output_dir, "protein_glycosylation_ilr_results.csv"), row.names = FALSE)
  
  # Create ILR-based plots
  cat("\nCreating ILR-based plots...\n")
  
  # Store ILR data for plotting
  ilr_fuc_results <- list()
  ilr_neuac_results <- list()
  
  # Re-calculate ILR data for plotting
  for(protein_name in names(protein_list)) {
    protein_data <- protein_list[[protein_name]]
    
    # Skip if insufficient data
    if(n_distinct(protein_data$sample) < 3) next
    
    # Check group sample sizes
    sample_counts <- protein_data %>%
      group_by(group) %>%
      summarise(n_samples = n_distinct(sample), .groups = 'drop')
    
    if(any(sample_counts$n_samples < 3)) next
    
    # Calculate fucosylation with ILR
    fuc_by_sample <- protein_data %>%
      group_by(sample, group) %>%
      summarise(
        total_abundance = sum(abundance, na.rm = TRUE),
        fuc_abundance = sum(abundance[contains_Fuc == TRUE], na.rm = TRUE),
        fuc_percentage = (fuc_abundance / total_abundance) * 100,
        .groups = 'drop'
      ) %>%
      mutate(
        fuc_proportion = fuc_percentage / 100,
        non_fuc_proportion = 1 - fuc_proportion,
        fuc_ilr = log(fuc_proportion / non_fuc_proportion)
      )
    
    # Calculate sialylation with ILR
    neuac_by_sample <- protein_data %>%
      group_by(sample, group) %>%
      summarise(
        total_abundance = sum(abundance, na.rm = TRUE),
        neuac_abundance = sum(abundance[contains_NeuAc == TRUE], na.rm = TRUE),
        neuac_percentage = (neuac_abundance / total_abundance) * 100,
        .groups = 'drop'
      ) %>%
      mutate(
        neuac_proportion = neuac_percentage / 100,
        non_neuac_proportion = 1 - neuac_proportion,
        # Handle edge cases for ILR transformation
        neuac_proportion_adj = ifelse(neuac_proportion == 0, 0.001, 
                                     ifelse(neuac_proportion == 1, 0.999, neuac_proportion)),
        non_neuac_proportion_adj = 1 - neuac_proportion_adj,
        neuac_ilr = log(neuac_proportion_adj / non_neuac_proportion_adj)
      ) %>%
      # Remove infinite or NaN values
      filter(is.finite(neuac_ilr))
    
    # Check if sufficient data for plotting
    fuc_sample_counts <- fuc_by_sample %>%
      group_by(group) %>%
      summarise(n_samples = n(), .groups = 'drop')
    
    neuac_sample_counts <- neuac_by_sample %>%
      group_by(group) %>%
      summarise(n_samples = n(), .groups = 'drop')
    
    has_sufficient_fuc_ilr <- nrow(fuc_sample_counts) >= 2 && all(fuc_sample_counts$n_samples >= 3)
    has_sufficient_neuac_ilr <- nrow(neuac_sample_counts) >= 2 && all(neuac_sample_counts$n_samples >= 3)
    
    # Store data for plotting
    if(has_sufficient_fuc_ilr) {
      ilr_fuc_results[[protein_name]] <- fuc_by_sample
    }
    if(has_sufficient_neuac_ilr) {
      ilr_neuac_results[[protein_name]] <- neuac_by_sample
    }
  }
  
  # Create ILR-based boxplots
  for(protein_name in names(protein_list)) {
    safe_protein_name <- gsub("[^a-zA-Z0-9]", "_", protein_name)
    
    # Create fucosylation ILR boxplot if sufficient data exists
    if(protein_name %in% names(ilr_fuc_results)) {
      fuc_ilr_plot <- ggplot(ilr_fuc_results[[protein_name]], 
                            aes(x = group, y = fuc_ilr, fill = group)) +
        geom_boxplot(alpha = 0.7, outlier.shape = 1) +
        geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
        scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
        labs(title = paste("Fucosylation ILR -", protein_name),
             x = "Group", 
             y = "Fucosylation ILR") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 14),
              legend.position = "top",
              axis.text = element_text(size = 12))
      
      # Save fucosylation ILR plot
      ggsave(paste0(figures_dir, "/", safe_protein_name, "_fucosylation_ilr.png"), 
             fuc_ilr_plot, width = 8, height = 6, dpi = 300)
    }
    
    # Create sialylation ILR boxplot if sufficient data exists
    if(protein_name %in% names(ilr_neuac_results)) {
      neuac_ilr_plot <- ggplot(ilr_neuac_results[[protein_name]], 
                              aes(x = group, y = neuac_ilr, fill = group)) +
        geom_boxplot(alpha = 0.7, outlier.shape = 1) +
        geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
        scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
        labs(title = paste("Sialylation ILR -", protein_name),
             x = "Group", 
             y = "Sialylation ILR") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 14),
              legend.position = "top",
              axis.text = element_text(size = 12))
      
      # Save sialylation ILR plot
      ggsave(paste0(figures_dir, "/", safe_protein_name, "_sialylation_ilr.png"), 
             neuac_ilr_plot, width = 8, height = 6, dpi = 300)
    }
  }
  
  # Create comparison plots (percentage vs ILR side by side)
  cat("\nCreating comparison plots (Percentage vs ILR)...\n")
  
  # Get original results if available
  if(!is.null(original_results)) {
    fuc_results <- original_results$fucosylation_data
    neuac_results <- original_results$sialylation_data
  } else {
    # If no original results provided, skip comparison plots
    cat("No original results provided, skipping comparison plots...\n")
    fuc_results <- list()
    neuac_results <- list()
  }
  
  for(protein_name in names(protein_list)) {
    safe_protein_name <- gsub("[^a-zA-Z0-9]", "_", protein_name)
    
    # Check if we have both percentage and ILR data
    has_percentage_fuc <- protein_name %in% names(fuc_results)
    has_ilr_fuc <- protein_name %in% names(ilr_fuc_results)
    has_percentage_neuac <- protein_name %in% names(neuac_results)
    has_ilr_neuac <- protein_name %in% names(ilr_neuac_results)
    
    # Create fucosylation comparison plot
    if(has_percentage_fuc && has_ilr_fuc) {
      # Combine data for comparison
      fuc_comparison_data <- fuc_results[[protein_name]] %>%
        select(sample, group, fuc_percentage) %>%
        mutate(data_type = "Percentage") %>%
        rename(value = fuc_percentage) %>%
        rbind(
          ilr_fuc_results[[protein_name]] %>%
            select(sample, group, fuc_ilr) %>%
            mutate(data_type = "ILR") %>%
            rename(value = fuc_ilr)
        )
      
      fuc_comparison_plot <- ggplot(fuc_comparison_data, 
                                   aes(x = group, y = value, fill = group)) +
        geom_boxplot(alpha = 0.7, outlier.shape = 1) +
        geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
        scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
        facet_wrap(~data_type, scales = "free_y") +
        labs(title = paste("Fucosylation Comparison -", protein_name),
             x = "Group", 
             y = "Value") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 14),
              legend.position = "top",
              axis.text = element_text(size = 12),
              strip.text = element_text(size = 12, face = "bold"))
      
      # Save fucosylation comparison plot
      ggsave(paste0(figures_dir, "/", safe_protein_name, "_fucosylation_comparison.png"), 
             fuc_comparison_plot, width = 12, height = 6, dpi = 300)
    }
    
    # Create sialylation comparison plot
    if(has_percentage_neuac && has_ilr_neuac) {
      # Combine data for comparison
      neuac_comparison_data <- neuac_results[[protein_name]] %>%
        select(sample, group, neuac_percentage) %>%
        mutate(data_type = "Percentage") %>%
        rename(value = neuac_percentage) %>%
        rbind(
          ilr_neuac_results[[protein_name]] %>%
            select(sample, group, neuac_ilr) %>%
            mutate(data_type = "ILR") %>%
            rename(value = neuac_ilr)
        )
      
      neuac_comparison_plot <- ggplot(neuac_comparison_data, 
                                     aes(x = group, y = value, fill = group)) +
        geom_boxplot(alpha = 0.7, outlier.shape = 1) +
        geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
        scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
        facet_wrap(~data_type, scales = "free_y") +
        labs(title = paste("Sialylation Comparison -", protein_name),
             x = "Group", 
             y = "Value") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 14),
              legend.position = "top",
              axis.text = element_text(size = 12),
              strip.text = element_text(size = 12, face = "bold"))
      
      # Save sialylation comparison plot
      ggsave(paste0(figures_dir, "/", safe_protein_name, "_sialylation_comparison.png"), 
             neuac_comparison_plot, width = 12, height = 6, dpi = 300)
    }
  }
  
  cat("ILR analysis completed. Results saved to:", file.path(output_dir, "protein_glycosylation_ilr_results.csv"), "\n")
  cat("ILR plots saved to:", figures_dir, "\n")
  
  return(list(
    ilr_statistics = ilr_summary_stats,
    fucosylation_data_ilr = ilr_fuc_results,
    sialylation_data_ilr = ilr_neuac_results
  ))
}

# Run the ILR-based analysis
protein_glycosylation_ilr_results <- analyze_protein_glycosylation_ilr(glyco_peptide_groups_long, protein_glycosylation_results)

###########################################################################################################
## Glycosite-level fucosylation and sialylation analysis

# Function to analyze glycosite-level glycosylation patterns
analyze_glycosite_glycosylation <- function(data, output_dir = "output_data/peptidegroups_intensity/glycosite_level", 
                                           figures_dir = "figures/glycosite_level") {
  
  # Create output directories if they don't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Split data by glycosite ID
  glycosite_list <- split(data, data$gsite_ID)
  
  # Initialize results storage
  fuc_results <- list()
  neuac_results <- list()
  summary_stats <- data.frame()
  
  # Track glycosites for summary
  glycosites_analyzed <- 0
  glycosites_insufficient_data <- 0
  glycosites_insufficient_groups <- 0
  glycosites_with_fuc_data <- 0
  glycosites_with_neuac_data <- 0
  fuc_tests_performed <- 0
  neuac_tests_performed <- 0
  
  # Track failed tests
  failed_fuc_tests <- data.frame(
    glycosite_ID = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
  
  failed_neuac_tests <- data.frame(
    glycosite_ID = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
  
  cat("Analyzing", length(glycosite_list), "glycosites for fucosylation and sialylation patterns...\n")
  
  # Process each glycosite
  for(i in seq_along(glycosite_list)) {
    glycosite_name <- names(glycosite_list)[i]
    glycosite_data <- glycosite_list[[i]]
    
    cat(sprintf("\nProcessing glycosite %d/%d: %s\n", i, length(glycosite_list), glycosite_name))
    
    # Skip glycosites with insufficient data
    if(n_distinct(glycosite_data$sample) < 3) {
      cat("  Skipping - insufficient total samples\n")
      glycosites_insufficient_data <- glycosites_insufficient_data + 1
      next
    }
    
    # Check if we have at least 3 samples in each group
    sample_counts <- glycosite_data %>%
      group_by(group) %>%
      summarise(n_samples = n_distinct(sample), .groups = 'drop')
    
    if(any(sample_counts$n_samples < 3)) {
      cat("  Skipping - insufficient samples in each group (need ≥3 per group)\n")
      cat("    Sample counts:", paste(sample_counts$group, "=", sample_counts$n_samples, collapse = ", "), "\n")
      glycosites_insufficient_groups <- glycosites_insufficient_groups + 1
      next
    }
    
    glycosites_analyzed <- glycosites_analyzed + 1
    
    # Calculate fucosylation percentages by sample
    fuc_by_sample <- glycosite_data %>%
      group_by(sample, group) %>%
      summarise(
        total_abundance = sum(abundance, na.rm = TRUE),
        fuc_abundance = sum(abundance[contains_Fuc == TRUE], na.rm = TRUE),
        fuc_percentage = (fuc_abundance / total_abundance) * 100,
        .groups = 'drop'
      )
    
    # Calculate sialylation percentages by sample
    neuac_by_sample <- glycosite_data %>%
      group_by(sample, group) %>%
      summarise(
        total_abundance = sum(abundance, na.rm = TRUE),
        neuac_abundance = sum(abundance[contains_NeuAc == TRUE], na.rm = TRUE),
        neuac_percentage = (neuac_abundance / total_abundance) * 100,
        .groups = 'drop'
      )
    
    # Check if glycosite has sufficient fucosylation data (≥3 samples per group with fucosylation)
    fuc_sample_counts <- fuc_by_sample %>%
      filter(fuc_percentage > 0) %>%  # Only samples with some fucosylation
      group_by(group) %>%
      summarise(n_samples = n(), .groups = 'drop')
    
    has_sufficient_fuc <- nrow(fuc_sample_counts) >= 2 && all(fuc_sample_counts$n_samples >= 3)
    
    # Check if glycosite has sufficient sialylation data (≥3 samples per group with sialylation)
    neuac_sample_counts <- neuac_by_sample %>%
      filter(neuac_percentage > 0) %>%  # Only samples with some sialylation
      group_by(group) %>%
      summarise(n_samples = n(), .groups = 'drop')
    
    has_sufficient_neuac <- nrow(neuac_sample_counts) >= 2 && all(neuac_sample_counts$n_samples >= 3)
    
    # Store results only if sufficient data exists
    if(has_sufficient_fuc) {
      fuc_results[[glycosite_name]] <- fuc_by_sample
      glycosites_with_fuc_data <- glycosites_with_fuc_data + 1
    } else {
      cat("    No sufficient fucosylation data (need ≥3 samples per group with fucosylation)\n")
    }
    if(has_sufficient_neuac) {
      neuac_results[[glycosite_name]] <- neuac_by_sample
      glycosites_with_neuac_data <- glycosites_with_neuac_data + 1
    } else {
      cat("    No sufficient sialylation data (need ≥3 samples per group with sialylation)\n")
    }
    
    # Perform statistical tests for fucosylation
    if(has_sufficient_fuc) {
      
      # Check for constant data within groups
      healthy_fuc_data <- fuc_by_sample$fuc_percentage[fuc_by_sample$group == "Healthy"]
      mecfs_fuc_data <- fuc_by_sample$fuc_percentage[fuc_by_sample$group == "MECFS"]
      
      # Check if data is constant within groups
      if(length(unique(healthy_fuc_data)) <= 1 && length(unique(mecfs_fuc_data)) <= 1) {
        failed_fuc_tests <- rbind(failed_fuc_tests, data.frame(
          glycosite_ID = glycosite_name,
          reason = "Constant data within both groups",
          stringsAsFactors = FALSE
        ))
        cat("    Fucosylation test failed: Constant data within both groups\n")
      } else if(length(unique(healthy_fuc_data)) <= 1) {
        failed_fuc_tests <- rbind(failed_fuc_tests, data.frame(
          glycosite_ID = glycosite_name,
          reason = "Constant data in Healthy group",
          stringsAsFactors = FALSE
        ))
        cat("    Fucosylation test failed: Constant data in Healthy group\n")
      } else if(length(unique(mecfs_fuc_data)) <= 1) {
        failed_fuc_tests <- rbind(failed_fuc_tests, data.frame(
          glycosite_ID = glycosite_name,
          reason = "Constant data in MECFS group",
          stringsAsFactors = FALSE
        ))
        cat("    Fucosylation test failed: Constant data in MECFS group\n")
      } else {
        # F-test for variance
        fuc_f_test <- try({
          var.test(fuc_percentage ~ group, data = fuc_by_sample)
        }, silent = TRUE)
        
        # T-test
        fuc_t_test <- try({
          t.test(fuc_percentage ~ group, data = fuc_by_sample)
        }, silent = TRUE)
        
        if(!inherits(fuc_t_test, "try-error")) {
          fuc_stats <- data.frame(
            glycosite_ID = glycosite_name,
            feature = "Fucosylation",
            mean_healthy = mean(healthy_fuc_data, na.rm = TRUE),
            mean_mecfs = mean(mecfs_fuc_data, na.rm = TRUE),
            sd_healthy = sd(healthy_fuc_data, na.rm = TRUE),
            sd_mecfs = sd(mecfs_fuc_data, na.rm = TRUE),
            t_stat = fuc_t_test$statistic,
            p_value = fuc_t_test$p.value,
            f_stat = ifelse(!inherits(fuc_f_test, "try-error"), fuc_f_test$statistic, NA),
            f_p_value = ifelse(!inherits(fuc_f_test, "try-error"), fuc_f_test$p.value, NA),
            n_healthy = sum(fuc_by_sample$group == "Healthy"),
            n_mecfs = sum(fuc_by_sample$group == "MECFS"),
            stringsAsFactors = FALSE
          )
          summary_stats <- rbind(summary_stats, fuc_stats)
          fuc_tests_performed <- fuc_tests_performed + 1
        } else {
          failed_fuc_tests <- rbind(failed_fuc_tests, data.frame(
            glycosite_ID = glycosite_name,
            reason = paste("T-test error:", fuc_t_test[1]),
            stringsAsFactors = FALSE
          ))
          cat("    Fucosylation t-test failed:", fuc_t_test[1], "\n")
        }
      }
    }
    
    # Perform statistical tests for sialylation
    if(has_sufficient_neuac) {
      
      # Check for constant data within groups
      healthy_neuac_data <- neuac_by_sample$neuac_percentage[neuac_by_sample$group == "Healthy"]
      mecfs_neuac_data <- neuac_by_sample$neuac_percentage[neuac_by_sample$group == "MECFS"]
      
      # Check if data is constant within groups
      if(length(unique(healthy_neuac_data)) <= 1 && length(unique(mecfs_neuac_data)) <= 1) {
        failed_neuac_tests <- rbind(failed_neuac_tests, data.frame(
          glycosite_ID = glycosite_name,
          reason = "Constant data within both groups",
          stringsAsFactors = FALSE
        ))
        cat("    Sialylation test failed: Constant data within both groups\n")
      } else if(length(unique(healthy_neuac_data)) <= 1) {
        failed_neuac_tests <- rbind(failed_neuac_tests, data.frame(
          glycosite_ID = glycosite_name,
          reason = "Constant data in Healthy group",
          stringsAsFactors = FALSE
        ))
        cat("    Sialylation test failed: Constant data in Healthy group\n")
      } else if(length(unique(mecfs_neuac_data)) <= 1) {
        failed_neuac_tests <- rbind(failed_neuac_tests, data.frame(
          glycosite_ID = glycosite_name,
          reason = "Constant data in MECFS group",
          stringsAsFactors = FALSE
        ))
        cat("    Sialylation test failed: Constant data in MECFS group\n")
      } else {
        # F-test for variance
        neuac_f_test <- try({
          var.test(neuac_percentage ~ group, data = neuac_by_sample)
        }, silent = TRUE)
        
        # T-test
        neuac_t_test <- try({
          t.test(neuac_percentage ~ group, data = neuac_by_sample)
        }, silent = TRUE)
        
        if(!inherits(neuac_t_test, "try-error")) {
          neuac_stats <- data.frame(
            glycosite_ID = glycosite_name,
            feature = "Sialylation",
            mean_healthy = mean(healthy_neuac_data, na.rm = TRUE),
            mean_mecfs = mean(mecfs_neuac_data, na.rm = TRUE),
            sd_healthy = sd(healthy_neuac_data, na.rm = TRUE),
            sd_mecfs = sd(mecfs_neuac_data, na.rm = TRUE),
            t_stat = neuac_t_test$statistic,
            p_value = neuac_t_test$p.value,
            f_stat = ifelse(!inherits(neuac_f_test, "try-error"), neuac_f_test$statistic, NA),
            f_p_value = ifelse(!inherits(neuac_f_test, "try-error"), neuac_f_test$p.value, NA),
            n_healthy = sum(neuac_by_sample$group == "Healthy"),
            n_mecfs = sum(neuac_by_sample$group == "MECFS"),
            stringsAsFactors = FALSE
          )
          summary_stats <- rbind(summary_stats, neuac_stats)
          neuac_tests_performed <- neuac_tests_performed + 1
        } else {
          failed_neuac_tests <- rbind(failed_neuac_tests, data.frame(
            glycosite_ID = glycosite_name,
            reason = paste("T-test error:", neuac_t_test[1]),
            stringsAsFactors = FALSE
          ))
          cat("    Sialylation t-test failed:", neuac_t_test[1], "\n")
        }
      }
    }
  }
  
  # Apply BH correction to all p-values
  if(nrow(summary_stats) > 0) {
    summary_stats$p_value_adj <- p.adjust(summary_stats$p_value, method = "BH")
    summary_stats$significant <- summary_stats$p_value_adj < 0.05
    summary_stats$fold_change <- summary_stats$mean_mecfs / summary_stats$mean_healthy
  }
  
  # Create separate boxplots for fucosylation and sialylation
  cat("\nCreating separate boxplots for fucosylation and sialylation...\n")
  
  for(glycosite_name in names(glycosite_list)) {
    safe_glycosite_name <- gsub("[^a-zA-Z0-9]", "_", glycosite_name)
    
    # Create fucosylation boxplot if sufficient data exists
    if(glycosite_name %in% names(fuc_results)) {
      fuc_plot <- ggplot(fuc_results[[glycosite_name]], 
                        aes(x = group, y = fuc_percentage, fill = group)) +
        geom_boxplot(alpha = 0.7, outlier.shape = 1) +
        geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
        scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
        labs(title = paste("Fucosylation -", glycosite_name),
             x = "Group", 
             y = "Fucosylation (%)") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 10),
              legend.position = "top")
      
      # Save fucosylation plot
      ggsave(paste0(figures_dir, "/", safe_glycosite_name, "_fucosylation.png"), 
             fuc_plot, width = 8, height = 6, dpi = 300)
    }
    
    # Create sialylation boxplot if sufficient data exists
    if(glycosite_name %in% names(neuac_results)) {
      neuac_plot <- ggplot(neuac_results[[glycosite_name]], 
                          aes(x = group, y = neuac_percentage, fill = group)) +
        geom_boxplot(alpha = 0.7, outlier.shape = 1) +
        geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
        scale_fill_manual(values = c("Healthy" = "#9DD4CC", "MECFS" = "#E49CB1")) +
        labs(title = paste("Sialylation -", glycosite_name),
             x = "Group", 
             y = "Sialylation (%)") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 10),
              legend.position = "top")
      
      # Save sialylation plot
      ggsave(paste0(figures_dir, "/", safe_glycosite_name, "_sialylation.png"), 
             neuac_plot, width = 8, height = 6, dpi = 300)
    }
  }
  
  # Save results
  saveRDS(fuc_results, paste0(output_dir, "/fucosylation_by_glycosite.rds"))
  saveRDS(neuac_results, paste0(output_dir, "/sialylation_by_glycosite.rds"))
  write.csv(summary_stats, paste0(output_dir, "/glycosite_glycosylation_statistics.csv"), row.names = FALSE)
  
  # Create separate summary files for fucosylation and sialylation
  fuc_summary <- summary_stats %>% filter(feature == "Fucosylation")
  neuac_summary <- summary_stats %>% filter(feature == "Sialylation")
  
  write.csv(fuc_summary, paste0(output_dir, "/fucosylation_statistics.csv"), row.names = FALSE)
  write.csv(neuac_summary, paste0(output_dir, "/sialylation_statistics.csv"), row.names = FALSE)
  
  # Create detailed summary dataframes for each feature
  fuc_detailed <- bind_rows(lapply(names(fuc_results), function(glycosite) {
    data <- fuc_results[[glycosite]]
    data$glycosite_ID <- glycosite
    return(data)
  }))
  
  neuac_detailed <- bind_rows(lapply(names(neuac_results), function(glycosite) {
    data <- neuac_results[[glycosite]]
    data$glycosite_ID <- glycosite
    return(data)
  }))
  
  write.csv(fuc_detailed, paste0(output_dir, "/fucosylation_detailed_data.csv"), row.names = FALSE)
  write.csv(neuac_detailed, paste0(output_dir, "/sialylation_detailed_data.csv"), row.names = FALSE)
  
  # Save failed test information
  write.csv(failed_fuc_tests, paste0(output_dir, "/failed_fucosylation_tests.csv"), row.names = FALSE)
  write.csv(failed_neuac_tests, paste0(output_dir, "/failed_sialylation_tests.csv"), row.names = FALSE)
  
  # Print summary
  cat("\n=== GLYCOSITE-LEVEL GLYCOSYLATION ANALYSIS SUMMARY ===\n")
  cat("Total glycosites in dataset:", length(glycosite_list), "\n")
  cat("Glycosites with sufficient data (≥3 samples per group):", glycosites_analyzed, "\n")
  cat("Glycosites excluded - insufficient total samples:", glycosites_insufficient_data, "\n")
  cat("Glycosites excluded - insufficient samples per group:", glycosites_insufficient_groups, "\n")
  cat("Glycosites with sufficient fucosylation data (≥3 samples per group with fucosylation):", glycosites_with_fuc_data, "\n")
  cat("Glycosites with sufficient sialylation data (≥3 samples per group with sialylation):", glycosites_with_neuac_data, "\n")
  cat("Fucosylation statistical tests performed:", fuc_tests_performed, "\n")
  cat("Sialylation statistical tests performed:", neuac_tests_performed, "\n")
  cat("Total statistical comparisons with results:", nrow(summary_stats), "\n")
  
  if(nrow(summary_stats) > 0) {
    cat("Significant differences (FDR < 0.05):", sum(summary_stats$significant, na.rm = TRUE), "\n")
    
    # Fucosylation summary
    fuc_stats <- summary_stats %>% filter(feature == "Fucosylation")
    cat(sprintf("\n=== FUCOSYLATION ANALYSIS ===\n"))
    cat("Glycosites with sufficient fucosylation data:", glycosites_with_fuc_data, "\n")
    cat("Glycosites with successful fucosylation statistical tests:", nrow(fuc_stats), "\n")
    cat("Significant fucosylation differences:", sum(fuc_stats$significant, na.rm = TRUE), "\n")
    
    if(nrow(fuc_stats) > 0) {
      cat("Mean fucosylation by group:\n")
      fuc_means <- fuc_stats %>%
        summarise(
          mean_healthy = mean(mean_healthy, na.rm = TRUE),
          mean_mecfs = mean(mean_mecfs, na.rm = TRUE),
          .groups = 'drop'
        )
      cat(sprintf("  Healthy: %.2f%%\n", fuc_means$mean_healthy))
      cat(sprintf("  MECFS: %.2f%%\n", fuc_means$mean_mecfs))
    }
    
    # Significant fucosylation differences
    sig_fuc <- fuc_stats %>% filter(significant)
    if(nrow(sig_fuc) > 0) {
      cat("\nSignificant fucosylation differences:\n")
      for(i in 1:nrow(sig_fuc)) {
        cat(sprintf("  %s: p = %.3e, FC = %.2f\n", 
                   sig_fuc$glycosite_ID[i], 
                   sig_fuc$p_value_adj[i], 
                   sig_fuc$fold_change[i]))
      }
    }
    
    # Sialylation summary
    neuac_stats <- summary_stats %>% filter(feature == "Sialylation")
    cat(sprintf("\n=== SIALYLATION ANALYSIS ===\n"))
    cat("Glycosites with sufficient sialylation data:", glycosites_with_neuac_data, "\n")
    cat("Glycosites with successful sialylation statistical tests:", nrow(neuac_stats), "\n")
    cat("Significant sialylation differences:", sum(neuac_stats$significant, na.rm = TRUE), "\n")
    
    if(nrow(neuac_stats) > 0) {
      cat("Mean sialylation by group:\n")
      neuac_means <- neuac_stats %>%
        summarise(
          mean_healthy = mean(mean_healthy, na.rm = TRUE),
          mean_mecfs = mean(mean_mecfs, na.rm = TRUE),
          .groups = 'drop'
        )
      cat(sprintf("  Healthy: %.2f%%\n", neuac_means$mean_healthy))
      cat(sprintf("  MECFS: %.2f%%\n", neuac_means$mean_mecfs))
    }
    
    # Significant sialylation differences
    sig_neuac <- neuac_stats %>% filter(significant)
    if(nrow(sig_neuac) > 0) {
      cat("\nSignificant sialylation differences:\n")
      for(i in 1:nrow(sig_neuac)) {
        cat(sprintf("  %s: p = %.3e, FC = %.2f\n", 
                   sig_neuac$glycosite_ID[i], 
                   sig_neuac$p_value_adj[i], 
                   sig_neuac$fold_change[i]))
      }
    }
  }
  
  cat("\nResults saved to:", output_dir, "\n")
  cat("  - fucosylation_statistics.csv: Fucosylation statistical results\n")
  cat("  - sialylation_statistics.csv: Sialylation statistical results\n")
  cat("  - fucosylation_detailed_data.csv: Detailed fucosylation data by glycosite and sample\n")
  cat("  - sialylation_detailed_data.csv: Detailed sialylation data by glycosite and sample\n")
  cat("  - glycosite_glycosylation_statistics.csv: Combined statistical results\n")
  cat("  - fucosylation_by_glycosite.rds: Fucosylation data as R object\n")
  cat("  - sialylation_by_glycosite.rds: Sialylation data as R object\n")
  cat("  - failed_fucosylation_tests.csv: Glycosites with failed fucosylation tests\n")
  cat("  - failed_sialylation_tests.csv: Glycosites with failed sialylation tests\n")
  cat("\nFigures saved to:", figures_dir, "\n")
  cat("  - *_fucosylation.png: Individual fucosylation boxplots\n")
  cat("  - *_sialylation.png: Individual sialylation boxplots\n")
  
  # Report failed tests
  if(nrow(failed_fuc_tests) > 0) {
    cat("\n=== FAILED FUCOSYLATION TESTS ===\n")
    cat("Total failed fucosylation tests:", nrow(failed_fuc_tests), "\n")
    
    # Count reasons
    fuc_reason_counts <- table(failed_fuc_tests$reason)
    cat("Reasons for failure:\n")
    for(i in 1:length(fuc_reason_counts)) {
      cat(sprintf("  %s: %d glycosites\n", names(fuc_reason_counts)[i], fuc_reason_counts[i]))
    }
    
    cat("\nGlycosites with failed fucosylation tests:\n")
    for(i in 1:nrow(failed_fuc_tests)) {
      cat(sprintf("  %s: %s\n", failed_fuc_tests$glycosite_ID[i], failed_fuc_tests$reason[i]))
    }
  }
  
  if(nrow(failed_neuac_tests) > 0) {
    cat("\n=== FAILED SIALYLATION TESTS ===\n")
    cat("Total failed sialylation tests:", nrow(failed_neuac_tests), "\n")
    
    # Count reasons
    neuac_reason_counts <- table(failed_neuac_tests$reason)
    cat("Reasons for failure:\n")
    for(i in 1:length(neuac_reason_counts)) {
      cat(sprintf("  %s: %d glycosites\n", names(neuac_reason_counts)[i], neuac_reason_counts[i]))
    }
    
    cat("\nGlycosites with failed sialylation tests:\n")
    for(i in 1:nrow(failed_neuac_tests)) {
      cat(sprintf("  %s: %s\n", failed_neuac_tests$glycosite_ID[i], failed_neuac_tests$reason[i]))
    }
  }
  
  return(list(
    fucosylation_data = fuc_results,
    sialylation_data = neuac_results,
    statistics = summary_stats,
    failed_fucosylation_tests = failed_fuc_tests,
    failed_sialylation_tests = failed_neuac_tests
  ))
}

# Run the glycosite-level analysis
glycosite_glycosylation_results <- analyze_glycosite_glycosylation(glyco_peptide_groups_long)

###########################################################################################################
library(dplyr)
library(tidyr)
library(stringr)
library(EnhancedVolcano)

# --- Summarise glycopeptide data to gene level ---
gene_level_data <- glyco_peptide_groups_long %>%
  group_by(gene_name, sample) %>%
  summarise(
    total_abundance = sum(abundance, na.rm = TRUE),
    .groups = 'drop'
  )

# --- Map samples to groups ---
gene_stats <- gene_level_data %>%
  mutate(
    group = case_when(
      str_detect(sample, "^hc") ~ "Healthy",
      str_detect(sample, "^m")  ~ "MECFS",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(gene_name, group) %>%
  summarise(
    mean_abundance = mean(total_abundance, na.rm = TRUE),
    sd_abundance   = sd(total_abundance, na.rm = TRUE),
    n              = n(),
    .groups = 'drop'
  ) %>%
  pivot_wider(
    names_from = group,
    values_from = c(mean_abundance, sd_abundance, n),
    names_glue = "{.value}_{group}"
  )

# Debug: Check what columns were actually created
print("Column names after pivot_wider:")
print(colnames(gene_stats))
print("First few rows:")
print(head(gene_stats))

# --- Differential statistics ---
gene_stats <- gene_stats %>%
  mutate(
    log2FC = log2(mean_abundance_MECFS / mean_abundance_Healthy),
    se_diff = sqrt((sd_abundance_Healthy^2 / n_Healthy) +
                   (sd_abundance_MECFS^2 / n_MECFS)),
    df = n_Healthy + n_MECFS - 2,
    t_stat = (mean_abundance_MECFS - mean_abundance_Healthy) / se_diff,
    p_value = 2 * pt(-abs(t_stat), df),
    adj_p_value = p.adjust(p_value, method = "BH")
  )

# --- Volcano plot with EnhancedVolcano ---
EnhancedVolcano(
  gene_stats,
  lab = gene_stats$gene_name,
  x = "log2FC",
  y = "adj_p_value",
  xlab = "log2(Fold Change) MECFS / Healthy",
  ylab = "-log10(Adjusted P-value)",
  pCutoff = 0.08,
  FCcutoff = 0.5,      # |log2FC| > 0.5
  pointSize = 2.5,
  labSize = 3.5,
  colAlpha = 0.6,
  title = "Gene-Level Glycopeptide Changes",
  subtitle = "MECFS vs Healthy",
  caption = paste("Total genes:", nrow(gene_stats)),
  legendLabels = c("NS", "Log2FC", "Adj.P", "Adj.P & Log2FC"),
  legendPosition = "top",
  drawConnectors = TRUE,
  widthConnectors = 0.5
)

# --- Save plot ---
ggsave("figures/gene_level_volcano_plot_enhanced.png", width = 10, height = 8, dpi = 300)

# --- Summary ---
cat("\nGene-level differential abundance analysis:\n")
cat("Total genes analyzed:", nrow(gene_stats), "\n")
cat("Significantly different genes:", sum(gene_stats$adj_p_value < 0.05, na.rm = TRUE), "\n")
cat("Upregulated in MECFS:", sum(gene_stats$adj_p_value < 0.05 & gene_stats$log2FC > 0, na.rm = TRUE), "\n")
cat("Downregulated in MECFS:", sum(gene_stats$adj_p_value < 0.05 & gene_stats$log2FC < 0, na.rm = TRUE), "\n")

# Debug: Check what sample names we actually have
print("Sample names in gene_level_data:")
print(unique(gene_level_data$sample))

# Debug: Check the group assignment
gene_stats_debug <- gene_level_data %>%
  mutate(
    group = case_when(
      str_detect(sample, "^hc") ~ "Healthy",
      str_detect(sample, "^m")  ~ "MECFS",
      TRUE ~ NA_character_
    )
  )

print("Group assignments:")
print(gene_stats_debug %>% select(sample, group) %>% distinct())

print("Unique groups after assignment:")
print(unique(gene_stats_debug$group))

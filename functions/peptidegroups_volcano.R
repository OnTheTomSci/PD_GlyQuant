# Peptide Groups Volcano Analysis Functions
# This module contains functions for creating volcano plots at different levels of analysis

library(tidyverse)
library(EnhancedVolcano)
library(ggplot2)

#' Calculate volcano plot statistics for differential abundance analysis
#' 
#' @param data Long format data with abundance, sample, and group columns
#' @param group_cols Columns to group by for aggregation (e.g., gene_name, or c(gene_name, glycan_composition))
#' @param group1_name Name of first group (e.g., "Healthy")
#' @param group2_name Name of second group (e.g., "MECFS")
#' @return Data frame with volcano plot statistics
calculate_volcano_stats <- function(data, group_cols, group1_name = "Healthy", group2_name = "MECFS") {
  
  # Ensure group column exists
  if (!"group" %in% colnames(data)) {
    stop("Column 'group' not found in data")
  }
  
  # Calculate abundance by specified groups
  abundance_data <- data %>%
    group_by(across(c(all_of(group_cols), sample, group))) %>%
    summarise(
      total_abundance = sum(abundance, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Calculate statistics for volcano plot
  volcano_stats <- abundance_data %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      # Calculate mean abundance for each group
      mean_group1 = mean(total_abundance[group == group1_name], na.rm = TRUE),
      mean_group2 = mean(total_abundance[group == group2_name], na.rm = TRUE),
      # Perform t-test with error handling
      pvalue = tryCatch({
        t.test(
          total_abundance[group == group2_name],
          total_abundance[group == group1_name]
        )$p.value
      }, error = function(e) NA_real_),
      .groups = 'drop'
    ) %>%
    # Calculate log2 fold change
    mutate(
      log2FC = log2(mean_group2 / mean_group1),
      # Calculate adjusted p-values (handle NA values)
      adj_pvalue = ifelse(is.na(pvalue), NA_real_, p.adjust(pvalue, method = "BH")),
      sig = !is.na(adj_pvalue) & adj_pvalue < 0.05
    )
  
  return(volcano_stats)
}

#' Create gene-level volcano plot
#' 
#' @param data Long format data with abundance, sample, group columns
#' @param output_path Path to save the volcano plot
#' @param group1_name Name of first group (default: "Healthy")
#' @param group2_name Name of second group (default: "MECFS")
#' @param p_cutoff P-value cutoff for significance (default: 0.05)
#' @param fc_cutoff Fold change cutoff (default: 1)
#' @return List containing volcano plot and statistics
create_gene_volcano <- function(data, output_path = "figures/peptidegroups_intensity/gene_volcano_plot.png",
                               group1_name = "Healthy", group2_name = "MECFS",
                               p_cutoff = 0.05, fc_cutoff = 1) {
  
  # Calculate volcano statistics
  volcano_stats <- calculate_volcano_stats(data, "gene_name", group1_name, group2_name)
  
  # Add labels for significant genes
  volcano_stats <- volcano_stats %>%
    mutate(label = ifelse(sig & abs(log2FC) > fc_cutoff, gene_name, ""))
  
  # Create enhanced volcano plot
  volcano_plot <- EnhancedVolcano(
    volcano_stats,
    lab = volcano_stats$gene_name,
    x = 'log2FC',
    y = 'adj_pvalue',
    title = 'Differential glycoprotein abundance in ME/CFS vs Healthy',
    subtitle = 'ME/CFS vs Healthy Controls',
    caption = paste0('FC cutoff: ', 2^fc_cutoff, '; p-value cutoff: ', p_cutoff),
    pCutoff = p_cutoff,
    FCcutoff = fc_cutoff,
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
    vline = c(-fc_cutoff, fc_cutoff),
    vlineCol = c('red', 'red'),
    vlineType = c('twodash', 'twodash'),
    vlineWidth = c(0.5, 0.5),
    selectLab = volcano_stats$gene_name[volcano_stats$sig | abs(volcano_stats$log2FC) > fc_cutoff],
    raster = FALSE,
    max.overlaps = 40
  )
  
  # Create output directory if it doesn't exist
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  
  # Save the plot
  ggsave(
    output_path,
    volcano_plot,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  cat("Gene-level volcano plot saved to:", output_path, "\n")
  
  # Print summary statistics
  cat("\nGene-level volcano plot summary:\n")
  cat("Total genes analyzed:", nrow(volcano_stats), "\n")
  cat("Significant genes (p <", p_cutoff, "):", sum(volcano_stats$sig, na.rm = TRUE), "\n")
  cat("High fold change genes (|log2FC| >", fc_cutoff, "):", sum(abs(volcano_stats$log2FC) > fc_cutoff, na.rm = TRUE), "\n")
  cat("Significant AND high fold change:", sum(volcano_stats$sig & abs(volcano_stats$log2FC) > fc_cutoff, na.rm = TRUE), "\n")
  
  return(list(
    plot = volcano_plot,
    statistics = volcano_stats
  ))
}

#' Create protein-glycan combination volcano plot
#' 
#' @param data Long format data with abundance, sample, group columns
#' @param output_path Path to save the volcano plot
#' @param group1_name Name of first group (default: "Healthy")
#' @param group2_name Name of second group (default: "MECFS")
#' @param p_cutoff P-value cutoff for significance (default: 0.05)
#' @param fc_cutoff Fold change cutoff (default: 1)
#' @return List containing volcano plot and statistics
create_protein_glycan_volcano <- function(data, output_path = "figures/peptidegroups_intensity/protein_glycan_volcano_plot.png",
                                         group1_name = "Healthy", group2_name = "MECFS",
                                         p_cutoff = 0.05, fc_cutoff = 1) {
  
  # Calculate volcano statistics for protein-glycan combinations
  volcano_stats <- calculate_volcano_stats(data, c("gene_name", "glycan_composition"), group1_name, group2_name)
  
  # Create combined labels and add significance labels
  volcano_stats <- volcano_stats %>%
    mutate(
      protein_glycan_label = paste0(gene_name, ":", glycan_composition),
      label = ifelse(sig & abs(log2FC) > fc_cutoff, protein_glycan_label, "")
    )
  
  # Create enhanced volcano plot for protein-glycan combinations
  volcano_plot <- EnhancedVolcano(
    volcano_stats,
    lab = volcano_stats$protein_glycan_label,
    x = 'log2FC',
    y = 'adj_pvalue',
    title = 'Differential Protein-Glycan Abundance in ME/CFS vs Healthy',
    subtitle = 'ME/CFS vs Healthy Controls - Protein:Glycan combinations',
    caption = paste0('FC cutoff: ', 2^fc_cutoff, '; p-value cutoff: ', p_cutoff),
    pCutoff = p_cutoff,
    FCcutoff = fc_cutoff,
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
    vline = c(-fc_cutoff, fc_cutoff),
    vlineCol = c('red', 'red'),
    vlineType = c('twodash', 'twodash'),
    vlineWidth = c(0.5, 0.5),
    selectLab = volcano_stats$protein_glycan_label[startsWith(volcano_stats$protein_glycan_label, "FN1")],
    raster = FALSE,
    max.overlaps = 50
  )
  
  # Create output directory if it doesn't exist
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  
  # Save the protein-glycan volcano plot
  ggsave(
    output_path,
    volcano_plot,
    width = 12,
    height = 10,
    dpi = 300
  )
  
  cat("Protein-glycan volcano plot saved to:", output_path, "\n")
  
  # Print summary statistics
  cat("\nProtein-glycan volcano plot summary:\n")
  cat("Total protein-glycan combinations:", nrow(volcano_stats), "\n")
  cat("Significant combinations (p <", p_cutoff, "):", sum(volcano_stats$sig, na.rm = TRUE), "\n")
  cat("High fold change combinations (|log2FC| >", fc_cutoff, "):", sum(abs(volcano_stats$log2FC) > fc_cutoff, na.rm = TRUE), "\n")
  cat("Significant AND high fold change:", sum(volcano_stats$sig & abs(volcano_stats$log2FC) > fc_cutoff, na.rm = TRUE), "\n")
  
  # Print top 10 highest absolute fold change data points
  cat("\nTop 10 highest absolute fold change protein-glycan combinations:\n")
  volcano_stats %>%
    arrange(desc(abs(log2FC))) %>%
    head(10) %>%
    select(protein_glycan_label, log2FC, adj_pvalue) %>%
    print(n = 10)
  
  return(list(
    plot = volcano_plot,
    statistics = volcano_stats
  ))
}

#' Create protein-glycosite-glycan combination volcano plot
#' 
#' @param data Long format data with abundance, sample, group columns
#' @param output_path Path to save the volcano plot
#' @param group1_name Name of first group (default: "Healthy")
#' @param group2_name Name of second group (default: "MECFS")
#' @param p_cutoff P-value cutoff for significance (default: 0.05)
#' @param fc_cutoff Fold change cutoff (default: 1)
#' @return List containing volcano plot and statistics
create_protein_glycosite_glycan_volcano <- function(data, output_path = "figures/peptidegroups_intensity/protein_glycosite_glycan_volcano_plot.png",
                                                   group1_name = "Healthy", group2_name = "MECFS",
                                                   p_cutoff = 0.05, fc_cutoff = 1) {
  
  # Calculate volcano statistics for protein-glycosite-glycan combinations
  volcano_stats <- calculate_volcano_stats(data, c("gene_name", "protein_glycosite", "glycan_composition"), group1_name, group2_name)
  
  # Create combined labels and add significance labels
  volcano_stats <- volcano_stats %>%
    mutate(
      protein_glycosite_glycan_label = paste0(gene_name, ":", protein_glycosite, ":", glycan_composition),
      label = ifelse(sig & abs(log2FC) > fc_cutoff, protein_glycosite_glycan_label, "")
    )
  
  # Create enhanced volcano plot for protein-glycosite-glycan combinations
  volcano_plot <- EnhancedVolcano(
    volcano_stats,
    lab = volcano_stats$protein_glycosite_glycan_label,
    x = 'log2FC',
    y = 'adj_pvalue',
    title = 'Differential Protein-Glycosite-Glycan Abundance in ME/CFS vs Healthy',
    subtitle = 'ME/CFS vs Healthy Controls - Protein:Glycosite:Glycan combinations',
    caption = paste0('FC cutoff: ', 2^fc_cutoff, '; p-value cutoff: ', p_cutoff),
    pCutoff = p_cutoff,
    FCcutoff = fc_cutoff,
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
    xlim = c(min(volcano_stats$log2FC, na.rm = TRUE) - 0.5, 
             max(volcano_stats$log2FC, na.rm = TRUE) + 0.5),
    ylim = c(0, max(-log10(volcano_stats$adj_pvalue), na.rm = TRUE) + 0.7),
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
    vline = c(-fc_cutoff, fc_cutoff),
    vlineCol = c('red', 'red'),
    vlineType = c('twodash', 'twodash'),
    vlineWidth = c(0.5, 0.5),
    selectLab = volcano_stats$protein_glycosite_glycan_label[startsWith(volcano_stats$protein_glycosite_glycan_label, "FN1")],
    raster = FALSE,
    max.overlaps = 100
  )
  
  # Create output directory if it doesn't exist
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  
  # Save the protein-glycosite-glycan volcano plot
  ggsave(
    output_path,
    volcano_plot,
    width = 14,
    height = 12,
    dpi = 600
  )
  
  cat("Protein-glycosite-glycan volcano plot saved to:", output_path, "\n")
  
  # Print summary statistics
  cat("\nProtein-glycosite-glycan volcano plot summary:\n")
  cat("Total protein-glycosite-glycan combinations:", nrow(volcano_stats), "\n")
  cat("Significant combinations (p <", p_cutoff, "):", sum(volcano_stats$sig, na.rm = TRUE), "\n")
  cat("High fold change combinations (|log2FC| >", fc_cutoff, "):", sum(abs(volcano_stats$log2FC) > fc_cutoff, na.rm = TRUE), "\n")
  cat("Significant AND high fold change:", sum(volcano_stats$sig & abs(volcano_stats$log2FC) > fc_cutoff, na.rm = TRUE), "\n")
  
  # Show top significant protein-glycosite-glycan combinations
  top_significant <- volcano_stats %>%
    filter(sig & abs(log2FC) > fc_cutoff) %>%
    arrange(desc(abs(log2FC))) %>%
    head(20)
  
  if(nrow(top_significant) > 0) {
    cat("\nTop 20 most significant protein-glycosite-glycan combinations:\n")
    print(top_significant[, c("gene_name", "protein_glycosite", "glycan_composition", "log2FC", "adj_pvalue")])
  }
  
  return(list(
    plot = volcano_plot,
    statistics = volcano_stats
  ))
}
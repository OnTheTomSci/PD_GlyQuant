

#' Create gene-level volcano plot (EnhancedVolcano)
#'
#' @param data Long format data with abundance, sample, group columns
#' @param output_path Path to save the volcano plot
#' @param group1_name Name of first group (default: "Healthy")
#' @param group2_name Name of second group (default: "MECFS")
#' @param p_cutoff P-value cutoff for significance (default: 0.05, adjusted p OK)
#' @param fc_cutoff Linear fold-change cutoff (default: 1; e.g., 2 = 2-fold)
#' @param highlight_genes Vector of gene names to highlight/label (e.g., DIA_sig_genes)
#' @param use_adjusted_pvalue If TRUE, use adjusted p-values; if FALSE, use raw p-values (default: TRUE)
#' @return list(plot, statistics)
create_gene_volcano_sel_lab <- function(
    data,
    output_path = "figures/peptidegroups_intensity/gene_volcano_plot.png",
    group1_name = "Healthy",
    group2_name = "MECFS",
    p_cutoff = 0.05,
    fc_cutoff = 2,
    highlight_genes = NULL,
    use_adjusted_pvalue = TRUE
) {
  # Packages
  require(EnhancedVolcano)
  require(ggplot2)
  require(dplyr)
  require(grid)
  require(ggrepel)
  
  # Compute stats (expects columns: gene_name, log2FC, pvalue or adj_pvalue, sig)
  volcano_stats <- calculate_volcano_stats(data, "gene_name", group1_name, group2_name)
  
  # Ensure we have the necessary p-value columns (calculate_volcano_stats returns 'pvalue' and 'adj_pvalue')
  if (!"pvalue" %in% names(volcano_stats)) {
    stop("volcano_stats must contain 'pvalue' column.")
  }
  if (!"adj_pvalue" %in% names(volcano_stats)) {
    # If adj_pvalue doesn't exist, create it from pvalue
    volcano_stats <- mutate(volcano_stats, adj_pvalue = p.adjust(pvalue, method = "BH"))
  }
  
  # Choose which p-value to use for plotting and significance
  if (use_adjusted_pvalue) {
    volcano_stats <- volcano_stats %>%
      mutate(
        pval_for_plot = adj_pvalue,
        pval_type = "adjusted"
      )
  } else {
    # Use raw p-values
    volcano_stats <- volcano_stats %>%
      mutate(
        pval_for_plot = pvalue,  # Note: it's 'pvalue' not 'p_value'
        pval_type = "raw"
      )
  }
  
  # Recalculate significance based on selected p-value type
  volcano_stats <- volcano_stats %>%
    mutate(sig = !is.na(pval_for_plot) & pval_for_plot < p_cutoff & abs(log2FC) > log2(fc_cutoff))
  
  # Labels for significant + large-effect genes
  volcano_stats <- volcano_stats %>%
    mutate(label = ifelse(sig & abs(log2FC) > log2(fc_cutoff), gene_name, ""))
  
  # Determine which genes to label in the plot
  if (!is.null(highlight_genes) && length(highlight_genes) > 0) {
    # Label genes from highlight_genes list AND genes that meet FC cutoff
    genes_from_list <- highlight_genes[highlight_genes %in% volcano_stats$gene_name]
    genes_high_fc <- volcano_stats$gene_name[abs(volcano_stats$log2FC) > log2(fc_cutoff)]
    genes_to_label <- unique(c(genes_from_list, genes_high_fc))
    
    message("Labeling ", length(genes_to_label), " genes total:")
    message("  - ", length(genes_from_list), " from custom list (", length(highlight_genes), " provided)")
    message("  - ", length(genes_high_fc), " meeting FC cutoff (|log2FC| > ", round(log2(fc_cutoff), 2), ")")
  } else {
    # Default: label genes that meet FC cutoff (regardless of p-value)
    genes_to_label <- volcano_stats$gene_name[abs(volcano_stats$log2FC) > log2(fc_cutoff)]
    message("Labeling ", length(genes_to_label), " genes meeting FC cutoff (|log2FC| > ", round(log2(fc_cutoff), 2), ")")
  }
  
  # Axis ranges
  xpad <- 0.5
  xmax <- max(abs(volcano_stats$log2FC), na.rm = TRUE)
  xlim_use <- c(-xmax - xpad, xmax + xpad)
  ylim_use <- c(0, max(-log10(pmax(volcano_stats$pval_for_plot, .Machine$double.xmin)), na.rm = TRUE) + 1)
  
  # Lines (remember EnhancedVolcano uses -log10(p) units on y)
  hlines_p <- c(1e-7, 1e-4, 1e-2)
  
  # Determine labels based on p-value type (can't use ifelse with bquote)
  if (use_adjusted_pvalue) {
    y_label <- bquote(~-Log[10]~"adjusted p-value")
  } else {
    y_label <- bquote(~-Log[10]~"p-value")
  }
  
  caption_text <- sprintf("Linear FC cutoff: %s; %s p cutoff: %s", 
                         fc_cutoff, 
                         ifelse(use_adjusted_pvalue, "adjusted", "raw"),
                         p_cutoff)
  
  volcano_plot <- EnhancedVolcano(
    volcano_stats,
    lab         = volcano_stats$gene_name,
    x           = "log2FC",
    y           = "pval_for_plot",
    title       = "Differential glycoprotein abundance: ME/CFS vs Healthy",
    caption     = caption_text,
    
    pCutoff     = p_cutoff,                   # p-value threshold
    FCcutoff    = log2(fc_cutoff),           # log2(FC) threshold
    
    pointSize   = 3.0,
    labSize     = 3.0,
    labCol      = "black",
    labFace     = "bold",
    boxedLabels = TRUE,
    parseLabels = FALSE,
    
    col         = c("grey30", "forestgreen", "royalblue", "red2"),
    colAlpha    = 0.75,
    
    gridlines.major = TRUE,
    gridlines.minor = FALSE,
    border      = "partial",
    borderWidth = 0.5,
    borderColour= "black",
    
    xlim        = xlim_use,
    ylim        = ylim_use,
    xlab        = bquote(~Log[2]~"fold change"),
    ylab        = y_label,
    axisLabSize = 16,
    titleLabSize= 18,
    subtitleLabSize = 14,
    captionLabSize  = 12,
    legendPosition  = "top",
    legendLabSize   = 12,
    legendIconSize  = 4.0,
    
    drawConnectors   = TRUE,
    widthConnectors  = 0.5,
    typeConnectors   = "closed",
    endsConnectors   = "first",
    lengthConnectors = grid::unit(0.01, "npc"),
    
    hline       = -log10(hlines_p),
    hlineCol    = rep("red", length(hlines_p)),
    hlineType   = rep("twodash", length(hlines_p)),
    hlineWidth  = rep(0.5, length(hlines_p)),
    
    vline       = c(-log2(fc_cutoff), log2(fc_cutoff)),
    vlineCol    = c("red", "red"),
    vlineType   = c("twodash", "twodash"),
    vlineWidth  = c(0.5, 0.5),
    
    selectLab   = genes_to_label,
    raster      = FALSE,
    max.overlaps= 40
  )
  
  # Ensure output dir exists and save
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = output_path, plot = volcano_plot, width = 10, height = 8, dpi = 300, device = "png")
  
  message("Gene-level volcano plot saved to: ", output_path)
  message("Total genes: ", nrow(volcano_stats))
  message("Significant (p < ", p_cutoff, "): ", sum(volcano_stats$sig, na.rm = TRUE))
  message("|log2FC| > ", round(log2(fc_cutoff), 3), ": ", sum(abs(volcano_stats$log2FC) > log2(fc_cutoff), na.rm = TRUE))
  message("Sig & high |log2FC|: ",
          sum(volcano_stats$sig & abs(volcano_stats$log2FC) > log2(fc_cutoff), na.rm = TRUE))
  message("Total genes labeled: ", length(genes_to_label))
  
  invisible(list(plot = volcano_plot, statistics = volcano_stats, labeled_genes = genes_to_label))
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





# Example: Define DIA_sig_genes (replace with your actual gene list)
DIA_sig_genes <- c("FN1",
                   "LYVE1",
                   "SHBG",
                   "CFHR5",
                   "IGHV3-74",
                   "PIGR",
                   "THBS1",
                   "VASN",
                   "EFEMP1",
                   "IGLL1",
                   "DNAJC2",
                   "PROC",
                   "B2M",
                   "CDH5",
                   "IGKV6-21",
                   "PRG4",
                   "PIEZO1",
                   "PCYOX1"
  # Add your significant genes from DIA analysis here
  # Example: "FN1", "APOA1", "SERPINA1", "HP", "C3", etc.
)

# Create volcano plots at different levels

# Example 1: With adjusted p-values (default)
gene_volcano_results_sel_lab <- create_gene_volcano_sel_lab(
  data = glyco_peptide_groups_long,
  output_path = "figures/peptidegroups_intensity/gene_volcano_plot_selected_genes.png",
  group1_name = "Healthy",
  group2_name = "MECFS",
  p_cutoff = 0.05,
  fc_cutoff = 2,
  highlight_genes = DIA_sig_genes,
  use_adjusted_pvalue = TRUE  # Use adjusted p-values
)

# Example 2: With raw (non-adjusted) p-values
gene_volcano_results_sel_lab_raw <- create_gene_volcano_sel_lab(
  data = glyco_peptide_groups_long,
  output_path = "figures/peptidegroups_intensity/gene_volcano_plot_selected_genes_raw_pval.png",
  group1_name = "Healthy",
  group2_name = "MECFS",
  p_cutoff = 0.05,
  fc_cutoff = 2,
  use_adjusted_pvalue = FALSE  # Use raw p-values
)
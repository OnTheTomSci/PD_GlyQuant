

             #' @param gpeps_dataframe the input glycopsm dataframe with all the glyco features anotated to the data frame 
             #' @param top_lev_group the top level grouping is t how you subset and grup data to peform glycan feature anlysse at eg. protein or glycosites
             #' @param value_col Column name for the measurement values to analyze
             #' @param sample_col Column name to use for pivoting to wide format
             #' @param group_values Vector of expected group values (default: c("Healthy", "MECFS"))
             #' @param min_samples Minimum number of samples required per group (default: 3)
             #' @param file_prefix Prefix for output files (default: "Analysis_")
            
             #' @return a matrix of relative aubundances for each top level groupingings and for each sample 
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
    file_prefix = "Analysis_"
             ) {
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
               final_matrix <- as.matrix(result_df[, -c(1, 2)])  # Remove the first two columns (top_level_group and glycofeature)
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

             # Add sia_count column if it doesn't exist
             if (!"sia_count" %in% colnames(glyco_peptide_groups_long)) {
               glyco_peptide_groups_long$sia_count <- sapply(glyco_peptide_groups_long$glycan_composition, function(x) {
                 match <- regmatches(x, regexpr("NeuAc\\((\\d+)\\)", x, perl = TRUE))
                 if (length(match) > 0) {
                   paste("NeuAc", sub("NeuAc\\((\\d+)\\)", "\\1", match))
                 } else {
                   "NeuAc 0"
                 }
               })
             }
             
             protein_gly_class <- glyco_matrix(
               gpeps_dataframe = glyco_peptide_groups_long,
               top_lev_group = "protein_accessions",
               glycofeature_group = "glycan_class",
               value_col = "abundance",
               sample_col = "sample",
               group_col = "group",
               file_prefix = "protein_gly_class"
             )            
             write.csv(protein_gly_class, file = "output_data/protein_gly_class_RA.csv")

             protein_gly_sia <- glyco_matrix(
               gpeps_dataframe = glyco_peptide_groups_long,
               top_lev_group = "protein_accessions",
               glycofeature_group = "contains_NeuAc",
               value_col = "abundance",
               sample_col = "sample",
               group_col = "group",
               file_prefix = "protein_gly_sia"
             )            
             write.csv(protein_gly_sia, file = "output_data/protein_gly_sia_RA.csv")

             protein_gly_fuc <- glyco_matrix(
               gpeps_dataframe = glyco_peptide_groups_long,
               top_lev_group = "protein_accessions",
               glycofeature_group = "contains_Fuc",
               value_col = "abundance",
               sample_col = "sample",
               group_col = "group",
               file_prefix = "protein_gly_fuc"
             )       
             write.csv(protein_gly_fuc, file = "output_data/protein_gly_fuc_RA.csv")

             

             protein_gly_comp <- glyco_matrix(
               gpeps_dataframe = glyco_peptide_groups_long,
               top_lev_group = "protein_accessions",
               glycofeature_group = "glycan_composition",
               value_col = "abundance",
               sample_col = "sample",
               group_col = "group",
               file_prefix = "protein_gly_comp"
             )                         
             write.csv(protein_gly_comp, file = "output_data/protein_gly_comp_RA.csv")

             glycosite_gly_class <- glyco_matrix(
               gpeps_dataframe = glyco_peptide_groups_long,
               top_lev_group = "gsite_ID",
               glycofeature_group = "glycan_class",
               value_col = "abundance",
               sample_col = "sample",
               group_col = "group",
               file_prefix = "glycosite_gly_class"
             )            
             write.csv(glycosite_gly_class, file = "output_data/glycosite_gly_class_RA.csv")

             glycosite_gly_sia <- glyco_matrix(
               gpeps_dataframe = glyco_peptide_groups_long,
               top_lev_group = "gsite_ID",
               glycofeature_group = "contains_NeuAc",
               value_col = "abundance",
               sample_col = "sample",
               group_col = "group",
               file_prefix = "glycosite_gly_sia"
             )            
             write.csv(glycosite_gly_sia, file = "output_data/glycosite_gly_sia_RA.csv")

             glycosite_gly_fuc <- glyco_matrix(
               gpeps_dataframe = glyco_peptide_groups_long,
               top_lev_group = "gsite_ID",
               glycofeature_group = "contains_Fuc",
               value_col = "abundance",
               sample_col = "sample",
               group_col = "group",
               file_prefix = "glycosite_gly_fuc"
             )       
             write.csv(glycosite_gly_fuc, file = "output_data/glycosite_gly_fuc_RA.csv")

            

             glycosite_gly_comp <- glyco_matrix(
               gpeps_dataframe = glyco_peptide_groups_long,
               top_lev_group = "gsite_ID",
               glycofeature_group = "glycan_composition",
               value_col = "abundance",
               sample_col = "sample",
               group_col = "group",
               file_prefix = "glycosite_gly_comp"
             ) 
             write.csv(glycosite_gly_comp, file = "output_data/glycosite_gly_comp_RA.csv")        

            # Load required libraries for volcano plots
            library(ggplot2)
            library(patchwork)
            library(ggrepel)  # For non-overlapping labels
             
            # Function to create volcano plot from relative abundance matrix
            #' @param ra_matrix Relative abundance matrix
            #' @param title Plot title
            #' @param output_file Output file path
            #' @param use_adjusted_pvalue If TRUE, use adjusted p-values (default: TRUE)
            #' @param highlight_proteins Vector of protein accessions to highlight/label
            #' @param fc_cutoff Fold change cutoff for significance (default: 2)
            #' @param p_cutoff P-value cutoff for significance (default: 0.05)
            create_volcano_from_matrix <- function(ra_matrix, title, output_file, 
                                                  use_adjusted_pvalue = TRUE,
                                                  highlight_proteins = NULL,
                                                  fc_cutoff = 2,
                                                  p_cutoff = 0.05) {
               # Convert matrix to long format for analysis
               ra_long <- as.data.frame(ra_matrix) %>%
                 tibble::rownames_to_column("feature") %>%
                 tidyr::pivot_longer(
                   cols = -feature,
                   names_to = "sample",
                   values_to = "relative_abundance"
                 ) %>%
                 # Parse sample names to extract group information
                 mutate(
                   sample_clean = str_remove(sample, "_.*$"),
                   group = case_when(
                     str_detect(sample, "_Healthy$") ~ "Healthy",
                     str_detect(sample, "_MECFS$") ~ "MECFS",
                     TRUE ~ "Unknown"
                   ),
                   # Convert relative abundance to numeric
                   relative_abundance = as.numeric(relative_abundance)
                 ) %>%
                 filter(!is.na(relative_abundance), relative_abundance > 0)
               
             # Calculate volcano statistics
              volcano_stats <- ra_long %>%
                group_by(feature) %>%
                summarise(
                  mean_healthy = mean(relative_abundance[group == "Healthy"], na.rm = TRUE),
                  mean_mecfs = mean(relative_abundance[group == "MECFS"], na.rm = TRUE),
                  pvalue = tryCatch({
                    t.test(
                      relative_abundance[group == "MECFS"],
                      relative_abundance[group == "Healthy"]
                    )$p.value
                  }, error = function(e) NA_real_),
                  .groups = 'drop'
                ) %>%
                mutate(
                  log2FC = log2(mean_mecfs / mean_healthy),
                  adj_pvalue = ifelse(is.na(pvalue), NA_real_, p.adjust(pvalue, method = "BH"))
                )
              
              # Choose which p-value to use for plotting
              if (use_adjusted_pvalue) {
                volcano_stats <- volcano_stats %>%
                  mutate(
                    pval_for_plot = adj_pvalue,
                    pval_type = "adjusted"
                  )
              } else {
                volcano_stats <- volcano_stats %>%
                  mutate(
                    pval_for_plot = pvalue,
                    pval_type = "raw"
                  )
              }
              
              # Continue with volcano stats calculation
              volcano_stats <- volcano_stats %>%
                mutate(
                  neg_log10_pval = -log10(pval_for_plot),
                  direction = case_when(
                    is.na(log2FC) ~ "NS",
                    log2FC > log2(fc_cutoff) ~ "Up in ME/CFS",
                    log2FC < -log2(fc_cutoff) ~ "Down in ME/CFS",
                    TRUE ~ "NS"
                  ),
                  sig = !is.na(pval_for_plot) & pval_for_plot < p_cutoff & direction != "NS"
                )
              
              # Extract protein accessions from feature names (before underscore)
              volcano_stats <- volcano_stats %>%
                mutate(protein_accession = str_extract(feature, "^[^_]+"))
              
              # Clean feature names by removing "@ N | rare1" or similar patterns
              volcano_stats <- volcano_stats %>%
                mutate(
                  feature_clean = str_remove(feature, "\\s*@\\s*[^|]*\\|\\s*rare\\d+\\s*$"),
                  feature_clean = str_trim(feature_clean)
                )
              
              # Determine which features to label
              if (!is.null(highlight_proteins) && length(highlight_proteins) > 0) {
                # Label features matching highlight_proteins OR meeting FC cutoff
                volcano_stats <- volcano_stats %>%
                  mutate(
                    is_highlighted = protein_accession %in% highlight_proteins,
                    meets_fc_cutoff = abs(log2FC) > log2(fc_cutoff),
                    label = ifelse(is_highlighted | meets_fc_cutoff, feature_clean, "")
                  )
                n_highlighted <- sum(volcano_stats$is_highlighted, na.rm = TRUE)
                n_fc_cutoff <- sum(volcano_stats$meets_fc_cutoff, na.rm = TRUE)
                message("Labeling ", sum(volcano_stats$label != "", na.rm = TRUE), " features total:")
                message("  - ", n_highlighted, " from custom protein list (", length(highlight_proteins), " proteins provided)")
                message("  - ", n_fc_cutoff, " meeting FC cutoff (|log2FC| > ", round(log2(fc_cutoff), 2), ")")
              } else {
                # Label features that meet FC cutoff
                volcano_stats <- volcano_stats %>%
                  mutate(
                    is_highlighted = FALSE,
                    meets_fc_cutoff = abs(log2FC) > log2(fc_cutoff),
                    label = ifelse(meets_fc_cutoff, feature_clean, "")
                  )
                message("Labeling ", sum(volcano_stats$label != "", na.rm = TRUE), 
                       " features meeting FC cutoff (|log2FC| > ", round(log2(fc_cutoff), 2), ")")
              }
              
              volcano_stats <- volcano_stats %>%
                filter(!is.na(log2FC), is.finite(log2FC), !is.na(neg_log10_pval), is.finite(neg_log10_pval))
               
              # Create volcano plot with ggplot2
              # Determine y-axis label based on p-value type
              y_label <- ifelse(use_adjusted_pvalue,
                               expression(-log[10] ~ "adjusted p-value"),
                               expression(-log[10] ~ "p-value"))
              
              # Determine caption based on settings
              caption_text <- sprintf("FC cutoff: %.1f; %s p-value cutoff: %.3f", 
                                     fc_cutoff, 
                                     ifelse(use_adjusted_pvalue, "adjusted", "raw"),
                                     p_cutoff)
              
             volcano_plot <- ggplot(volcano_stats, aes(x = log2FC, y = neg_log10_pval)) +
                # Add background grid
                theme_bw() +
                theme(
                  panel.grid.major = element_line(color = "grey90", size = 0.3),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(color = "black", size = 0.5),
                  legend.position = "none",
                  plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 7, hjust = 0.5),
                  plot.caption = element_text(size = 6, hjust = 0.5),
                  axis.title = element_text(size = 8),
                  axis.text = element_text(size = 7)
                ) +
                # Add scatter points with color coding
               geom_point(aes(color = direction, alpha = sig), size = 1.5) +
               scale_color_manual(values = c(
                 "Down in ME/CFS" = "#1f78b4",  # blue
                 "NS" = "#bdbdbd",              # grey
                 "Up in ME/CFS" = "#e31a1c"      # red
               )) +
               scale_alpha_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.4), guide = "none") +
                # Add significance thresholds
                geom_vline(xintercept = c(-log2(fc_cutoff), log2(fc_cutoff)), linetype = "twodash", color = "blue", size = 0.3) +
                geom_hline(yintercept = -log10(p_cutoff), linetype = "twodash", color = "red", size = 0.3) +
                # Add non-overlapping labels using ggrepel
                geom_text_repel(
                  aes(label = label),
                  size = 1.8,
                  color = "black",
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.color = "grey50",
                  segment.size = 0.2,
                  max.overlaps = 50,
                  min.segment.length = 0,
                  force = 2,
                  force_pull = 0.5
                ) +
                # Set axis labels and title
                labs(
                  title = title,
                  subtitle = "ME/CFS vs Healthy Controls",
                  caption = caption_text,
                  x = expression(log[2] ~ "fold change"),
                  y = y_label
                ) +
                # Set axis limits
                xlim(min(volcano_stats$log2FC, na.rm = TRUE) - 0.5, 
                     max(volcano_stats$log2FC, na.rm = TRUE) + 0.5) +
                ylim(0, max(volcano_stats$neg_log10_pval, na.rm = TRUE) + 1)
               
               # Save plot with larger dimensions for better visibility
               ggsave(output_file, volcano_plot, width = 120, height = 120, units = "mm", dpi = 300)
               
               # Return both plot and statistics
               return(list(plot = volcano_plot, stats = volcano_stats))
             }
             
             # Create volcano plots for each relative abundance matrix
             
             # 1. Protein Glycan Class Volcano Plot
             cat("Creating protein glycan class volcano plot...\n")
            protein_gly_class_volcano <- create_volcano_from_matrix(
              protein_gly_class,
              "Protein glycan class",
               "figures/protein_glycan_class_volcano.png",
              use_adjusted_pvalue = FALSE,

             )
             write.csv(protein_gly_class_volcano$stats, 
                      file = "output_data/protein_glycan_class_volcano_stats.csv", 
                      row.names = FALSE)
             
            # 2. Protein Glycan Sialic Acid Volcano Plot
            cat("Creating protein glycan sialic acid volcano plot...\n")
            # Filter out rows ending in "FALSE" (non-sialylated)
            protein_gly_sia_filtered <- protein_gly_sia[!grepl("FALSE$", rownames(protein_gly_sia)), ]
            cat("  Filtered out ", nrow(protein_gly_sia) - nrow(protein_gly_sia_filtered), 
                " non-sialylated features (ending in FALSE)\n")
           protein_gly_sia_volcano <- create_volcano_from_matrix(
             protein_gly_sia_filtered,
             "Protein sialylation",
              "figures/protein_glycan_sia_volcano.png",
              use_adjusted_pvalue = FALSE,
            )
             write.csv(protein_gly_sia_volcano$stats, 
                      file = "output_data/protein_glycan_sia_volcano_stats.csv", 
                      row.names = FALSE)
             
            # 3. Protein Glycan Fucose Volcano Plot
            cat("Creating protein glycan fucose volcano plot...\n")
            # Filter out rows ending in "FALSE" (non-fucosylated)
            protein_gly_fuc_filtered <- protein_gly_fuc[!grepl("FALSE$", rownames(protein_gly_fuc)), ]
            cat("  Filtered out ", nrow(protein_gly_fuc) - nrow(protein_gly_fuc_filtered), 
                " non-fucosylated features (ending in FALSE)\n")
           protein_gly_fuc_volcano <- create_volcano_from_matrix(
             protein_gly_fuc_filtered,
             "Protein fucosylation",
              "figures/protein_glycan_fuc_volcano.png",
              use_adjusted_pvalue = FALSE,
            )
             write.csv(protein_gly_fuc_volcano$stats, 
                      file = "output_data/protein_glycan_fuc_volcano_stats.csv", 
                      row.names = FALSE)
             
             # 4. Protein Glycan Composition Volcano Plot
             cat("Creating protein glycan composition volcano plot...\n")
            protein_gly_comp_volcano <- create_volcano_from_matrix(
              protein_gly_comp,
              "Protein glycan composition",
               "figures/protein_glycan_composition_volcano.png",
               use_adjusted_pvalue = TRUE,
             )
             write.csv(protein_gly_comp_volcano$stats, 
                      file = "output_data/protein_glycan_composition_volcano_stats.csv", 
                      row.names = TRUE)
             
             # 5. Glycosite Glycan Class Volcano Plot
             cat("Creating glycosite glycan class volcano plot...\n")
            glycosite_gly_class_volcano <- create_volcano_from_matrix(
              glycosite_gly_class,
              "Glycosite glycan class",
               "figures/glycosite_glycan_class_volcano.png",
               use_adjusted_pvalue = FALSE,
             )
             write.csv(glycosite_gly_class_volcano$stats, 
                      file = "output_data/glycosite_glycan_class_volcano_stats.csv", 
                      row.names = FALSE)
             
            # 6. Glycosite Glycan Sialic Acid Volcano Plot
            cat("Creating glycosite glycan sialic acid volcano plot...\n")
            # Filter out rows ending in "FALSE" (non-sialylated)
            glycosite_gly_sia_filtered <- glycosite_gly_sia[!grepl("FALSE$", rownames(glycosite_gly_sia)), ]
            cat("  Filtered out ", nrow(glycosite_gly_sia) - nrow(glycosite_gly_sia_filtered), 
                " non-sialylated features (ending in FALSE)\n")
           glycosite_gly_sia_volcano <- create_volcano_from_matrix(
             glycosite_gly_sia_filtered,
             "Glycosite sialylation",
              "figures/glycosite_glycan_sia_volcano.png",
              use_adjusted_pvalue = FALSE,
            )
             write.csv(glycosite_gly_sia_volcano$stats, 
                      file = "output_data/glycosite_glycan_sia_volcano_stats.csv", 
                      row.names = FALSE)
             
            # 7. Glycosite Glycan Fucose Volcano Plot
            cat("Creating glycosite glycan fucose volcano plot...\n")
            # Filter out rows ending in "FALSE" (non-fucosylated)
            glycosite_gly_fuc_filtered <- glycosite_gly_fuc[!grepl("FALSE$", rownames(glycosite_gly_fuc)), ]
            cat("  Filtered out ", nrow(glycosite_gly_fuc) - nrow(glycosite_gly_fuc_filtered), 
                " non-fucosylated features (ending in FALSE)\n")
           glycosite_gly_fuc_volcano <- create_volcano_from_matrix(
             glycosite_gly_fuc_filtered,
             "Glycosite fucosylation",
              "figures/glycosite_glycan_fuc_volcano.png",
              use_adjusted_pvalue = FALSE,
            )
             write.csv(glycosite_gly_fuc_volcano$stats, 
                      file = "output_data/glycosite_glycan_fuc_volcano_stats.csv", 
                      row.names = FALSE)
             
             # 8. Glycosite Glycan Composition Volcano Plot
             cat("Creating glycosite glycan composition volcano plot...\n")
            glycosite_gly_comp_volcano <- create_volcano_from_matrix(
              glycosite_gly_comp,
              "Glycosite glycan composition",
              "figures/glycosite_glycan_composition_volcano.png",
               use_adjusted_pvalue = TRUE,
             )
             write.csv(glycosite_gly_comp_volcano$stats, 
                      file = "output_data/glycosite_glycan_composition_volcano_stats.csv", 
                      row.names = FALSE)
             
             # Create summary table of significant results
             cat("Creating summary table of significant results...\n")
             
             # Combine all significant results
             all_significant_results <- rbind(
               protein_gly_class_volcano$stats %>% 
                 filter(sig) %>% 
                 mutate(analysis_type = "Protein_Glycan_Class", 
                        feature_type = "Glycan_Class"),
               protein_gly_sia_volcano$stats %>% 
                 filter(sig) %>% 
                 mutate(analysis_type = "Protein_Glycan_Sialic_Acid", 
                        feature_type = "Sialic_Acid"),
               protein_gly_fuc_volcano$stats %>% 
                 filter(sig) %>% 
                 mutate(analysis_type = "Protein_Glycan_Fucose", 
                        feature_type = "Fucose"),
               protein_gly_comp_volcano$stats %>% 
                 filter(sig) %>% 
                 mutate(analysis_type = "Protein_Glycan_Composition", 
                        feature_type = "Glycan_Composition"),
               glycosite_gly_class_volcano$stats %>% 
                 filter(sig) %>% 
                 mutate(analysis_type = "Glycosite_Glycan_Class", 
                        feature_type = "Glycan_Class"),
               glycosite_gly_sia_volcano$stats %>% 
                 filter(sig) %>% 
                 mutate(analysis_type = "Glycosite_Glycan_Sialic_Acid", 
                        feature_type = "Sialic_Acid"),
               glycosite_gly_fuc_volcano$stats %>% 
                 filter(sig) %>% 
                 mutate(analysis_type = "Glycosite_Glycan_Fucose", 
                        feature_type = "Fucose"),
               glycosite_gly_comp_volcano$stats %>% 
                 filter(sig) %>% 
                 mutate(analysis_type = "Glycosite_Glycan_Composition", 
                        feature_type = "Glycan_Composition")
            ) %>%
            select(analysis_type, feature_type, feature, mean_healthy, mean_mecfs, 
                   log2FC, pvalue, pval_for_plot, sig) %>%
            arrange(analysis_type, desc(abs(log2FC)))
             
             # Save summary results
             write.csv(all_significant_results, 
                      file = "output_data/all_significant_volcano_results.csv", 
                      row.names = FALSE)
             
             # Print summary statistics
             cat("\n=== VOLCANO PLOT ANALYSIS SUMMARY ===\n")
             cat("Total significant features found:", nrow(all_significant_results), "\n")
             cat("\nSignificant features by analysis type:\n")
             print(table(all_significant_results$analysis_type))
             
             cat("\nSignificant features by feature type:\n")
             print(table(all_significant_results$feature_type))
             
            cat("\nTop 10 most significant features (by p-value):\n")
            top_significant <- all_significant_results %>%
              arrange(pval_for_plot) %>%
              head(10) %>%
              select(feature, analysis_type, log2FC, pvalue, pval_for_plot)
            print(top_significant)
             
             # Create combined figure panels using patchwork
             cat("Creating combined figure panels...\n")
             
             # Create protein-level combined panel (2x2 grid)
             protein_combined <- (protein_gly_class_volcano$plot + 
                                 protein_gly_sia_volcano$plot) / 
                                (protein_gly_fuc_volcano$plot + 
                                 protein_gly_comp_volcano$plot) +
               plot_annotation(
                 title = "Protein-Level Glycan Analysis: ME/CFS vs Healthy Controls",
                 subtitle = "Volcano plots showing differential glycan abundance patterns",
                 theme = theme(plot.title = element_text(size = 10, face = "bold"),
                             plot.subtitle = element_text(size = 7))
               )
             
             # Create glycosite-level combined panel (2x2 grid)
             glycosite_combined <- (glycosite_gly_class_volcano$plot + 
                                   glycosite_gly_sia_volcano$plot) / 
                                  (glycosite_gly_fuc_volcano$plot + 
                                   glycosite_gly_comp_volcano$plot) +
               plot_annotation(
                 title = "Glycosite-Level Glycan Analysis: ME/CFS vs Healthy Controls",
                 subtitle = "Volcano plots showing differential glycan abundance patterns at specific glycosites",
                 theme = theme(plot.title = element_text(size = 10, face = "bold"),
                             plot.subtitle = element_text(size = 7))
               )
             
             # Save combined panels with larger dimensions for better spacing
             # A4 dimensions: 210mm x 297mm
             # For 2x2 grid of 120mm plots, we need approximately 250mm width and 250mm height
             # This fits well within A4 page with margins
             
             cat("Saving protein-level combined panel...\n")
             ggsave("figures/protein_level_volcano_combined.png", 
                   protein_combined, 
                   width = 250, height = 250, units = "mm", dpi = 300)
             
             cat("Saving glycosite-level combined panel...\n")
             ggsave("figures/glycosite_level_volcano_combined.png", 
                   glycosite_combined, 
                   width = 250, height = 250, units = "mm", dpi = 300)
             
             # Create a single combined panel with all 8 plots (4x2 grid)
             # This might be too large for A4, so we'll create it but note the size
             all_combined <- (protein_gly_class_volcano$plot + 
                             protein_gly_sia_volcano$plot + 
                             protein_gly_fuc_volcano$plot + 
                             protein_gly_comp_volcano$plot) /
                            (glycosite_gly_class_volcano$plot + 
                             glycosite_gly_sia_volcano$plot + 
                             glycosite_gly_fuc_volcano$plot + 
                             glycosite_gly_comp_volcano$plot) +
               plot_annotation(
                 title = "Complete Glycan Analysis: Protein and Glycosite Levels",
                 subtitle = "Volcano plots showing differential glycan abundance patterns at protein and glycosite levels",
                 theme = theme(plot.title = element_text(size = 12, face = "bold"),
                             plot.subtitle = element_text(size = 8))
               )
             
             # Save all-combined panel (larger size for 4x2 grid)
             cat("Saving all-combined panel (4x2 grid)...\n")
             ggsave("figures/all_volcano_combined.png", 
                   all_combined, 
                   width = 500, height = 250, units = "mm", dpi = 300)
             
            cat("\n=== ANALYSIS COMPLETE ===\n")
            cat("Individual volcano plots saved to: figures/\n")
            cat("Combined panels saved:\n")
            cat("  - figures/protein_level_volcano_combined.png (2x2, 250x250mm)\n")
            cat("  - figures/glycosite_level_volcano_combined.png (2x2, 250x250mm)\n")
            cat("  - figures/all_volcano_combined.png (4x2, 500x250mm)\n")
            cat("Statistical results saved to: output_data/\n")
            cat("Summary results saved to: output_data/all_significant_volcano_results.csv\n")

            # ============================
            # FULL SET: Custom volcano plots with non-adjusted p-values and DIA protein highlighting
            # ============================
            
            # Define proteins of interest (e.g., DIA significant genes)
            DIA_sig_proteins <- c(
              "FN1", "LYVE1", "SHBG", "CFHR5", "IGHV3-74", "PIGR", 
              "THBS1", "VASN", "EFEMP1", "IGLL1", "DNAJC2", "PROC", 
              "B2M", "CDH5", "IGKV6-21", "PRG4", "PIEZO1", "PCYOX1"
            )
            
            # Create all 8 volcano plots with:
            # - Raw (non-adjusted) p-values
            # - Highlighting DIA significant proteins
            # - FC cutoff: 2
            # - P-value cutoff: 0.05
            cat("\n\n=== CREATING FULL SET OF CUSTOM VOLCANO PLOTS (Raw p-values + DIA proteins) ===\n")
            
            # Create output directories if they don't exist
            dir.create("figures/DIA_custom", recursive = TRUE, showWarnings = FALSE)
            dir.create("output_data/DIA_custom", recursive = TRUE, showWarnings = FALSE)
            cat("Created output directories for DIA custom analysis\n")
            
            # Check the format of protein accessions in the data
            cat("\nChecking protein accession format in data...\n")
            sample_features <- head(rownames(protein_gly_class), 10)
            cat("Sample features from data:\n")
            print(sample_features)
            sample_accessions <- str_extract(sample_features, "^[^_]+")
            cat("\nExtracted protein accessions:\n")
            print(sample_accessions)
            cat("\nDIA proteins to highlight:\n")
            print(DIA_sig_proteins)
            cat("\n")
            
            # 1. Protein Glycan Class - Custom
            cat("Creating protein glycan class volcano plot (DIA custom)...\n")
            protein_gly_class_volcano_DIA <- create_volcano_from_matrix(
              protein_gly_class,
              "Protein glycan class (DIA proteins)",
              "figures/DIA_custom/protein_glycan_class_volcano_DIA.png",
              use_adjusted_pvalue = FALSE,
              highlight_proteins = DIA_sig_proteins,
              fc_cutoff = 2,
              p_cutoff = 0.05
            )
            write.csv(protein_gly_class_volcano_DIA$stats, 
                     file = "output_data/DIA_custom/protein_glycan_class_volcano_DIA_stats.csv", 
                     row.names = FALSE)
            
            # 2. Protein Glycan Sialic Acid - Custom
            cat("Creating protein glycan sialic acid volcano plot (DIA custom)...\n")
            protein_gly_sia_volcano_DIA <- create_volcano_from_matrix(
              protein_gly_sia_filtered,  # Use filtered matrix
              "Protein sialylation (DIA proteins)",
              "figures/DIA_custom/protein_glycan_sia_volcano_DIA.png",
              use_adjusted_pvalue = FALSE,
              highlight_proteins = DIA_sig_proteins,
              fc_cutoff = 2,
              p_cutoff = 0.05
            )
            write.csv(protein_gly_sia_volcano_DIA$stats, 
                     file = "output_data/DIA_custom/protein_glycan_sia_volcano_DIA_stats.csv", 
                     row.names = FALSE)
            
            # 3. Protein Glycan Fucose - Custom
            cat("Creating protein glycan fucose volcano plot (DIA custom)...\n")
            protein_gly_fuc_volcano_DIA <- create_volcano_from_matrix(
              protein_gly_fuc_filtered,  # Use filtered matrix
              "Protein fucosylation (DIA proteins)",
              "figures/DIA_custom/protein_glycan_fuc_volcano_DIA.png",
              use_adjusted_pvalue = FALSE,
              highlight_proteins = DIA_sig_proteins,
              fc_cutoff = 2,
              p_cutoff = 0.05
            )
            write.csv(protein_gly_fuc_volcano_DIA$stats, 
                     file = "output_data/DIA_custom/protein_glycan_fuc_volcano_DIA_stats.csv", 
                     row.names = FALSE)
            
            # 4. Protein Glycan Composition - Custom
            cat("Creating protein glycan composition volcano plot (DIA custom)...\n")
            protein_gly_comp_volcano_DIA <- create_volcano_from_matrix(
              protein_gly_comp,
              "Protein glycan composition (DIA proteins)",
              "figures/DIA_custom/protein_glycan_comp_volcano_DIA.png",
              use_adjusted_pvalue = FALSE,
              highlight_proteins = DIA_sig_proteins,
              fc_cutoff = 2,
              p_cutoff = 0.05
            )
            write.csv(protein_gly_comp_volcano_DIA$stats, 
                     file = "output_data/DIA_custom/protein_glycan_comp_volcano_DIA_stats.csv", 
                     row.names = FALSE)
            
            # 5. Glycosite Glycan Class - Custom
            cat("Creating glycosite glycan class volcano plot (DIA custom)...\n")
            glycosite_gly_class_volcano_DIA <- create_volcano_from_matrix(
              glycosite_gly_class,
              "Glycosite glycan class (DIA proteins)",
              "figures/DIA_custom/glycosite_glycan_class_volcano_DIA.png",
              use_adjusted_pvalue = FALSE,
              highlight_proteins = DIA_sig_proteins,
              fc_cutoff = 2,
              p_cutoff = 0.05
            )
            write.csv(glycosite_gly_class_volcano_DIA$stats, 
                     file = "output_data/DIA_custom/glycosite_glycan_class_volcano_DIA_stats.csv", 
                     row.names = FALSE)
            
            # 6. Glycosite Glycan Sialic Acid - Custom
            cat("Creating glycosite glycan sialic acid volcano plot (DIA custom)...\n")
            glycosite_gly_sia_volcano_DIA <- create_volcano_from_matrix(
              glycosite_gly_sia_filtered,  # Use filtered matrix
              "Glycosite sialylation (DIA proteins)",
              "figures/DIA_custom/glycosite_glycan_sia_volcano_DIA.png",
              use_adjusted_pvalue = FALSE,
              highlight_proteins = DIA_sig_proteins,
              fc_cutoff = 2,
              p_cutoff = 0.05
            )
            write.csv(glycosite_gly_sia_volcano_DIA$stats, 
                     file = "output_data/DIA_custom/glycosite_glycan_sia_volcano_DIA_stats.csv", 
                     row.names = FALSE)
            
            # 7. Glycosite Glycan Fucose - Custom
            cat("Creating glycosite glycan fucose volcano plot (DIA custom)...\n")
            glycosite_gly_fuc_volcano_DIA <- create_volcano_from_matrix(
              glycosite_gly_fuc_filtered,  # Use filtered matrix
              "Glycosite fucosylation (DIA proteins)",
              "figures/DIA_custom/glycosite_glycan_fuc_volcano_DIA.png",
              use_adjusted_pvalue = FALSE,
              highlight_proteins = DIA_sig_proteins,
              fc_cutoff = 2,
              p_cutoff = 0.05
            )
            write.csv(glycosite_gly_fuc_volcano_DIA$stats, 
                     file = "output_data/DIA_custom/glycosite_glycan_fuc_volcano_DIA_stats.csv", 
                     row.names = FALSE)
            
            # 8. Glycosite Glycan Composition - Custom
            cat("Creating glycosite glycan composition volcano plot (DIA custom)...\n")
            glycosite_gly_comp_volcano_DIA <- create_volcano_from_matrix(
              glycosite_gly_comp,
              "Glycosite glycan composition (DIA proteins)",
              "figures/DIA_custom/glycosite_glycan_comp_volcano_DIA.png",
              use_adjusted_pvalue = FALSE,
              highlight_proteins = DIA_sig_proteins,
              fc_cutoff = 2,
              p_cutoff = 0.05
            )
            write.csv(glycosite_gly_comp_volcano_DIA$stats, 
                     file = "output_data/DIA_custom/glycosite_glycan_comp_volcano_DIA_stats.csv", 
                     row.names = FALSE)
            
            # Create combined DIA custom panels
            cat("Creating combined DIA custom figure panels...\n")
            
            # Protein-level DIA combined panel (2x2 grid)
            protein_combined_DIA <- (protein_gly_class_volcano_DIA$plot + 
                                    protein_gly_sia_volcano_DIA$plot) / 
                                   (protein_gly_fuc_volcano_DIA$plot + 
                                    protein_gly_comp_volcano_DIA$plot) +
              plot_annotation(
                title = "Protein-Level Glycan Analysis: DIA Significant Proteins Highlighted (Raw p-values)",
                subtitle = "ME/CFS vs Healthy Controls - Highlighting DIA proteomics hits",
                theme = theme(plot.title = element_text(size = 10, face = "bold"),
                            plot.subtitle = element_text(size = 7))
              )
            
            # Glycosite-level DIA combined panel (2x2 grid)
            glycosite_combined_DIA <- (glycosite_gly_class_volcano_DIA$plot + 
                                      glycosite_gly_sia_volcano_DIA$plot) / 
                                     (glycosite_gly_fuc_volcano_DIA$plot + 
                                      glycosite_gly_comp_volcano_DIA$plot) +
              plot_annotation(
                title = "Glycosite-Level Glycan Analysis: DIA Significant Proteins Highlighted (Raw p-values)",
                subtitle = "ME/CFS vs Healthy Controls - Highlighting DIA proteomics hits at specific glycosites",
                theme = theme(plot.title = element_text(size = 10, face = "bold"),
                            plot.subtitle = element_text(size = 7))
              )
            
            # Save DIA combined panels
            cat("Saving DIA protein-level combined panel...\n")
            ggsave("figures/DIA_custom/protein_level_volcano_combined_DIA.png", 
                  protein_combined_DIA, 
                  width = 250, height = 250, units = "mm", dpi = 300)
            
            cat("Saving DIA glycosite-level combined panel...\n")
            ggsave("figures/DIA_custom/glycosite_level_volcano_combined_DIA.png", 
                  glycosite_combined_DIA, 
                  width = 250, height = 250, units = "mm", dpi = 300)
            
            # All DIA plots combined (4x2 grid)
            all_combined_DIA <- (protein_gly_class_volcano_DIA$plot + 
                                protein_gly_sia_volcano_DIA$plot + 
                                protein_gly_fuc_volcano_DIA$plot + 
                                protein_gly_comp_volcano_DIA$plot) /
                               (glycosite_gly_class_volcano_DIA$plot + 
                                glycosite_gly_sia_volcano_DIA$plot + 
                                glycosite_gly_fuc_volcano_DIA$plot + 
                                glycosite_gly_comp_volcano_DIA$plot) +
              plot_annotation(
                title = "Complete Glycan Analysis: DIA Significant Proteins Highlighted (Raw p-values)",
                subtitle = "Volcano plots with raw p-values highlighting DIA proteomics significant hits",
                theme = theme(plot.title = element_text(size = 12, face = "bold"),
                            plot.subtitle = element_text(size = 8))
              )
            
            cat("Saving DIA all-combined panel (4x2 grid)...\n")
            ggsave("figures/DIA_custom/all_volcano_combined_DIA.png", 
                  all_combined_DIA, 
                  width = 500, height = 250, units = "mm", dpi = 300)
            
            # Create summary of DIA protein features
            cat("Creating summary of DIA protein significant results...\n")
            all_DIA_results <- rbind(
              protein_gly_class_volcano_DIA$stats %>% 
                filter(sig) %>% 
                mutate(analysis_type = "Protein_Glycan_Class", 
                       feature_type = "Glycan_Class"),
              protein_gly_sia_volcano_DIA$stats %>% 
                filter(sig) %>% 
                mutate(analysis_type = "Protein_Glycan_Sialic_Acid", 
                       feature_type = "Sialic_Acid"),
              protein_gly_fuc_volcano_DIA$stats %>% 
                filter(sig) %>% 
                mutate(analysis_type = "Protein_Glycan_Fucose", 
                       feature_type = "Fucose"),
              protein_gly_comp_volcano_DIA$stats %>% 
                filter(sig) %>% 
                mutate(analysis_type = "Protein_Glycan_Composition", 
                       feature_type = "Glycan_Composition"),
              glycosite_gly_class_volcano_DIA$stats %>% 
                filter(sig) %>% 
                mutate(analysis_type = "Glycosite_Glycan_Class", 
                       feature_type = "Glycan_Class"),
              glycosite_gly_sia_volcano_DIA$stats %>% 
                filter(sig) %>% 
                mutate(analysis_type = "Glycosite_Glycan_Sialic_Acid", 
                       feature_type = "Sialic_Acid"),
              glycosite_gly_fuc_volcano_DIA$stats %>% 
                filter(sig) %>% 
                mutate(analysis_type = "Glycosite_Glycan_Fucose", 
                       feature_type = "Fucose"),
              glycosite_gly_comp_volcano_DIA$stats %>% 
                filter(sig) %>% 
                mutate(analysis_type = "Glycosite_Glycan_Composition", 
                       feature_type = "Glycan_Composition")
            ) %>%
            select(analysis_type, feature_type, feature, protein_accession, is_highlighted, meets_fc_cutoff,
                   mean_healthy, mean_mecfs, log2FC, pvalue, pval_for_plot, sig) %>%
            arrange(analysis_type, desc(abs(log2FC)))
            
            # Save DIA summary results
            write.csv(all_DIA_results, 
                     file = "output_data/DIA_custom/all_significant_volcano_DIA_results.csv", 
                     row.names = FALSE)
            
            # Print DIA summary statistics
            cat("\n=== DIA CUSTOM VOLCANO PLOT ANALYSIS SUMMARY ===\n")
            cat("Total significant features found (raw p < 0.05, FC > 2):", nrow(all_DIA_results), "\n")
            cat("Features from DIA proteins:", sum(all_DIA_results$is_highlighted, na.rm = TRUE), "\n")
            cat("Features meeting FC cutoff (|log2FC| > 1):", sum(all_DIA_results$meets_fc_cutoff, na.rm = TRUE), "\n")
            cat("\nSignificant features by analysis type:\n")
            print(table(all_DIA_results$analysis_type))
            
            cat("\nTop 10 DIA protein features (by raw p-value):\n")
            top_DIA <- all_DIA_results %>%
              filter(is_highlighted) %>%
              arrange(pvalue) %>%
              head(10) %>%
              select(feature, analysis_type, log2FC, pvalue)
            print(top_DIA)
            
            cat("\n=== DIA CUSTOM ANALYSIS COMPLETE ===\n")
            cat("Individual DIA volcano plots saved to: figures/DIA_custom/\n")
            cat("Combined DIA panels saved:\n")
            cat("  - figures/DIA_custom/protein_level_volcano_combined_DIA.png (2x2, 250x250mm)\n")
            cat("  - figures/DIA_custom/glycosite_level_volcano_combined_DIA.png (2x2, 250x250mm)\n")
            cat("  - figures/DIA_custom/all_volcano_combined_DIA.png (4x2, 500x250mm)\n")
            cat("Statistical results saved to: output_data/DIA_custom/\n")
            cat("DIA summary results saved to: output_data/DIA_custom/all_significant_volcano_DIA_results.csv\n")
 
             
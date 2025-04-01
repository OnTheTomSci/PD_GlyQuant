#' ---
#' title: Proteome Discover exploratory data Analysis of 10 sample ME/CFS Plasma glycoproteomics
#'   cohort
#' output: html_notebook
#'
#' ---
#' #
#' The mass spec files of raw data are searched and analysised using proteome discoverer 
#' with byonic as a peptide search node. This analysis workflow uses a LFQ style 
#' quantification method using peak intensities even when a sample does not have an MS2 PSM 
#' corresponding to the peak. Instead the peak is still quantified based on a PSM found 
#' in another sample that matches the peak's retention times.
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Get started by loading the packages
library('tidyr')
library("knitr")
library("wrMisc")
library("wrProteo")
library("wrGraph")
library("readr")
library("ggplot2")
library("dplyr")
library("stringr")
library("purrr")
library("broom")
library("ggsignif")
library("viridis")
library("hrbrthemes")
library("janitor")
library(openxlsx2)
library(here)
library(scales)
library("EnhancedVolcano")
library(ggrepel)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
here()

# Read as tab-separated file
Proteins <- read_tsv("input_data/10S_MECFS_GPEPS_250125_Proteins.txt")
ProteinGroups <- read_tsv("input_data/10S_MECFS_GPEPS_250125_ProteinGroups.txt")
PeptideGroups <- read_tsv("input_data/10S_MECFS_GPEPS_250125_PeptideGroups.txt")
PSMs <- read_tsv("input_data/10S_MECFS_GPEPS_250125_PSMs.txt")
ConsensusFeatures <- read_tsv("input_data/10S_MECFS_GPEPS_250125_ConsensusFeatures.txt")
InputFiles <- read_tsv("input_data/10S_MECFS_GPEPS_250125_InputFiles.txt")
PathwayProteinGroups <- read_tsv("input_data/10S_MECFS_GPEPS_250125_PathwayProteinGroups.txt")
StudyInformation <- read_tsv("input_data/10S_MECFS_GPEPS_250125_StudyInformation.txt")

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
glycoPSMs <- PSMs %>% filter(!is.na(`Glycan Composition`)) # I think i may have forgot pep2D score filtering but it shoul have been done in PD using byonic as a node

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
glycoPSMs <- glycoPSMs %>%
  mutate(`Glycan Composition` = str_remove(`Glycan Composition`, "@ N \\| rare1$"))

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Glycans <- glycoPSMs %>%
  group_by(`Glycan Composition`) %>%
  summarise(count = n(), .groups = "drop")

png(
  "figures/glycan_barplot.png",
  width = 1600,
  height = 1200,
  res = 150
)  # High-res image
par(mar = c(5, 15, 4, 2))
barplot(
  Glycans$count,
  names.arg = Glycans$`Glycan Composition`,
  horiz = TRUE,
  las = 1,
  col = "steelblue",
  main = "Glycan Composition Counts",
  xlab = "Count",
  cex.names = 0.9
)
dev.off()  # Save the file


#' The following barchart of glycan PSM counts shows the plasma glycome is predominated 
#' by a few abundant glycans, thus skewing the distribution of abundances. This is what 
#' we would expect and is similarly reflected in the released glycan analysis.
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             glycoPSMs <- glycoPSMs %>%
               mutate(
                 Sample = str_extract(
                   `Spectrum File`,
                   "(?<=20250116_OE_TR_10S_MECFS_GPEP_)(.*?)(?=\\.raw)"
                 ),
                 Disease_Status = ifelse(
                   str_detect(Sample, "^HC"),
                   "Healthy",
                   ifelse(str_detect(Sample, "M"), "MECFS", NA)
                 ),
                 Protein_Names = str_extract(`Master Protein Descriptions`, "^[^O]*(?=OS=)")
               )
             
             
          
         
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             glycoPSMs <- glycoPSMs %>%
               clean_names() %>%
               mutate(
                 pep_glycosite = str_extract(modifications, "(?<=N)\\d+(?=\\()"),
                 protein_names = str_extract(master_protein_descriptions, "^[^O]*(?=OS=)"),
                 gene_name = str_extract(master_protein_descriptions, "(?<=GN=).*?(?= PE=)"),
                 contains_Fuc = str_detect(glycan_composition, "Fuc"),
                 contains_NeuAc = str_detect(glycan_composition, "NeuAc"),
                 pep_glycosite = as.numeric(pep_glycosite),
                 protein_glycosite = position_in_protein + pep_glycosite - 1
               ) %>%
               mutate(
                 sample = str_extract(spectrum_file, "(?<=GPEP_)[^\\.]+"),
                 disease_status = ifelse(
                   str_detect(sample, "HC"),
                   "Healthy",
                   ifelse(str_detect(sample, "M"), "MECFS", NA_character_)
                 )
               )
             
             
             
             #'
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             glycan_class_map <- read_csv(file = "input_data/glycan_class_map.csv")
             
             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             match_glycan_class <- function(input_df, reference_df) {
               # Ensure required columns exist
               if (!"glycan_composition" %in% colnames(input_df)) {
                 stop("Column 'glycan_composition' not found in input dataframe")
               }
               if (!all(c("glycans", "glycan_class") %in% colnames(reference_df))) {
                 stop("Columns 'glycans' and 'glycan_class' not found in reference dataframe")
               }
               
               # Convert to character and trim whitespace
               input_df$glycan_composition <- trimws(tolower(as.character(input_df$glycan_composition)))
               reference_df$glycans <- trimws(tolower(as.character(reference_df$glycans)))
               
               # Initialize an empty vector to store the matched glycan classes
               matched_classes <- vector("character", length = nrow(input_df))
               
              # Iterate over each glycan_composition in the input dataframe
               for (i in 1:nrow(input_df)) {
                 glycan_comp <- input_df$glycan_composition[i]
                 
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
               input_df$glycan_class <- matched_classes
               
               return(input_df)
             }
             
    #---------------------------------------------------------------------------------------------------------------------------------------
             count_sia <- function(input_df) {
               # Ensure required column exists
               if (!"glycan_composition" %in% colnames(input_df)) {
                 stop("Column 'glycan_composition' not found in input dataframe")
               }
               
               # Extract number of NeuAc occurrences and add as new column
               input_df$sia_count <- sapply(input_df$glycan_composition, function(x) {
                 match <- regmatches(x, regexpr("NeuAc\\((\\d+)\\)", x, perl = TRUE))
                 if (length(match) > 0) {
                   paste("NeuAc", sub("NeuAc\\((\\d+)\\)", "\\1", match))
                 } else {
                   "NeuAc 0"
                 }
               })
               
               return(input_df)
             }
             
             count_fuc <- function(input_df) {
               # Ensure required column exists
               if (!"glycan_composition" %in% colnames(input_df)) {
                 stop("Column 'glycan_composition' not found in input dataframe")
               }
               
               # Extract number of Fuc occurrences and add as new column
               input_df$fuc_count <- sapply(input_df$glycan_composition, function(x) {
                 match <- regmatches(x, regexpr("Fuc\\((\\d+)\\)", x, perl = TRUE))
                 if (length(match) > 0) {
                   paste("Fucose", sub("Fuc\\((\\d+)\\)", "\\1", match))
                 } else {
                   "Fucose 0"
                 }
               })
               
               return(input_df)
             }
             
             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             glycoPSMs <- glycoPSMs %>%
               count_sia() %>%
               count_fuc() %>%
               match_glycan_class(., glycan_class_map)  # Explicitly passing the second argument
             
             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             # Ensure gsite_ID is created
             glycoPSMs <- glycoPSMs %>%
               mutate(gsite_ID = paste(protein_accessions, protein_glycosite, sep = "_"))
             
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
             
               gsites <- find_unique_values(glycoPSMs, "gsite_ID")
               glycoprots <- find_unique_values(glycoPSMs, "protein_accessions")
               
       #-----------------------------------------------------------------------------------------------        
              
               
               #' Perform Statistical Analysis on Multiple Datasets
               #'
               #' @param data_list A list of data frames to analyze
               #' @param group_col Column name for grouping (e.g., disease status)
               #' @param value_col Column name for the measurement values to analyze
               #' @param sample_col Column name to use for pivoting to wide format
               #' @param group_values Vector of expected group values 
               #'        (default: c("Healthy", "MECFS"))
               #' @param min_samples Minimum number of samples required per group (default: 3)
               #' @param file_prefix Prefix for output files (default: "Analysis_")
               #' @param perform_tests List of tests to perform (default: all available)
               #'
               #' @return A list containing all analysis results
               #'
               #' @export
               perform_analysis_on_list <- function(data_list,
                                                    group_col,
                                                    value_col,
                                                    sample_col,
                                                    group_values = c("Healthy", "MECFS"),
                                                    min_samples = 3,
                                                    file_prefix = "Analysis_",
                                                    perform_tests = c("descriptive", "ttest", "anova", "boxplot")) {
                 
                 # Validate inputs
                 if (!is.list(data_list)) {
                   stop("data_list must be a list of data frames")
                 }
                 
                 if (length(data_list) == 0) {
                   warning("Empty data list provided. No analysis performed.")
                   return(NULL)
                 }
                 
                 # Create output filename and workbook
                 output_file <- paste0(file_prefix, "analysis_results.xlsx")
                 wb <- wb_workbook()
                 
                 # Initialize results list to return
                 all_results <- list()
                 
                 # Process each dataset
                 for (i in seq_along(data_list)) {
                   tryCatch({
                     data <- data_list[[i]]
                     dataset_name <- paste0("Dataset_", i)
                     
                     # Check if required columns exist
                     if (!all(c(group_col, value_col) %in% colnames(data))) {
                       warning(paste("Skipping", dataset_name, ": Missing required columns"))
                       next
                     }
                     
                     if (!is.null(sample_col) && !(sample_col %in% colnames(data))) {
                       warning(paste("Skipping wide format for", dataset_name, ": Missing sample column"))
                       sample_col <- NULL
                     }
                     
                     # Count samples in each group
                     disease_counts <- data %>%
                       group_by(!!sym(group_col)) %>%
                       summarise(count = n(), .groups = "drop")
                     
                     # Store counts in results
                     all_results[[paste0(dataset_name, "_Group_Counts")]] <- disease_counts
                     
                     # Check if we have all required groups with minimum sample sizes
                     has_required_groups <- all(group_values %in% unique(data[[group_col]]))
                     has_min_samples <- all(sapply(group_values, function(g) {
                       sum(data[[group_col]] == g, na.rm = TRUE) >= min_samples
                     }))
                     
                     if (!has_required_groups || !has_min_samples) {
                       message(
                         "Skipping analysis for ", dataset_name, 
                         ": Not enough samples in each required group (", 
                         paste(group_values, collapse = ", "), 
                         "). Need at least ", min_samples, " per group."
                       )
                       next
                     }
                     
                     # Initialize results for this dataset
                     dataset_results <- list()
                     
                     # Descriptive statistics
                     if ("descriptive" %in% perform_tests) {
                       descriptive_stats <- data %>%
                         group_by(!!sym(group_col)) %>%
                         summarise(
                           mean = mean(!!sym(value_col), na.rm = TRUE),
                           sd = sd(!!sym(value_col), na.rm = TRUE),
                           median = median(!!sym(value_col), na.rm = TRUE),
                           min = min(!!sym(value_col), na.rm = TRUE),
                           max = max(!!sym(value_col), na.rm = TRUE),
                           count = n(),
                           .groups = "drop"
                         )
                       
                       dataset_results[["Descriptive_Stats"]] <- descriptive_stats
                     }
                     
                     # T-test
                     if ("ttest" %in% perform_tests) {
                       # Safe t-test that handles potential errors
                       t_test_result <- tryCatch({
                         # Filter only the groups we're interested in
                         test_data <- data %>% 
                           filter(!!sym(group_col) %in% group_values)
                         
                         result <- t.test(
                           formula = as.formula(paste(value_col, "~", group_col)),
                           data = test_data,
                           var.equal = FALSE
                         )
                         
                         data.frame(
                           Statistic = result$statistic,
                           P_Value = result$p.value,
                           DF = result$parameter,
                           Confidence_Lower = result$conf.int[1],
                           Confidence_Upper = result$conf.int[2],
                           Mean_Diff = diff(result$estimate)
                         )
                       }, 
                       error = function(e) {
                         warning(paste("T-test failed for", dataset_name, ":", e$message))
                         data.frame(
                           Error = e$message
                         )
                       })
                       
                       dataset_results[["T_test_Result"]] <- t_test_result
                     }
                     
                     # ANOVA
                     if ("anova" %in% perform_tests) {
                       anova_result <- tryCatch({
                         model <- aov(as.formula(paste(value_col, "~", group_col)), data = data)
                         as.data.frame(summary(model)[[1]])
                       },
                       error = function(e) {
                         warning(paste("ANOVA failed for", dataset_name, ":", e$message))
                         data.frame(
                           Error = e$message
                         )
                       })
                       
                       dataset_results[["ANOVA_Result"]] <- anova_result
                     }
                     
                     # Wide format data 
                     if (!is.null(sample_col)) {
                       tryCatch({
                         wide_data <- data %>%
                           pivot_wider(
                             id_cols = setdiff(names(data), c(sample_col, value_col)),
                             names_from = !!sym(sample_col),
                             values_from = !!sym(value_col),
                             values_fill = NA
                           )
                         
                         dataset_results[["Wide_Data"]] <- wide_data
                       },
                       error = function(e) {
                         warning(paste("Pivot to wide format failed for", dataset_name, ":", e$message))
                       })
                     }
                     
                     # Boxplot for significant differences
                     if ("boxplot" %in% perform_tests && 
                         "T_test_Result" %in% names(dataset_results) && 
                         !("Error" %in% names(dataset_results[["T_test_Result"]])) &&
                         dataset_results[["T_test_Result"]]$P_Value < 0.05) {
                       
                       plot_filename <- paste0(file_prefix, "boxplot_", i, ".png")
                       
                       tryCatch({
                         boxplot <- ggplot(data, aes(
                           x = !!sym(group_col),
                           y = !!sym(value_col)
                         )) +
                           geom_boxplot(aes(fill = !!sym(group_col))) +
                           geom_jitter(width = 0.2, alpha = 0.5) +  # Add individual points
                           theme_minimal() +
                           labs(
                             title = paste("Boxplot for", dataset_name),
                             subtitle = paste("p-value =", round(dataset_results[["T_test_Result"]]$P_Value, 4)),
                             x = "Disease Status",
                             y = "Measurement Value"
                           ) +
                           scale_fill_manual(values = setNames(
                             c("blue", "red", rep("gray", length(unique(data[[group_col]])) - 2)), 
                             unique(data[[group_col]])
                           ))
                         
                         ggsave(
                           plot_filename,
                           plot = boxplot,
                           width = 6,
                           height = 4
                         )
                         
                         dataset_results[["Boxplot_Path"]] <- plot_filename
                       },
                       error = function(e) {
                         warning(paste("Boxplot creation failed for", dataset_name, ":", e$message))
                       })
                     }
                     
                     # Save results to workbook
                     for (name in names(dataset_results)) {
                       sheet_name <- paste0(dataset_name, "_", name)
                       tryCatch({
                         wb <- wb_add_worksheet(wb, sheet_name)
                         wb <- wb_add_data(wb, sheet = sheet_name, x = dataset_results[[name]])
                       },
                       error = function(e) {
                         warning(paste("Failed to add worksheet", sheet_name, ":", e$message))
                       })
                     }
                     
                     # Add this dataset's results to the overall results
                     all_results[[dataset_name]] <- dataset_results
                     
                   }, error = function(e) {
                     warning(paste("Error processing", paste0("Dataset_", i), ":", e$message))
                   })
                 }
                 
                 # Save the workbook
                 tryCatch({
                   wb_save(wb, output_file)
                   message("Analysis complete. Results saved to: ", output_file)
                 }, 
                 error = function(e) {
                   warning(paste("Failed to save workbook:", e$message))
                 })
                 
                 return(all_results)
               }
               
               
               
               
               
               
               
    #------------------------------------------------------------------------------------------------------           
               # Group by gsite_ID and split into a named list
             glycosite_list <- glycoPSMs %>%
               group_by(gsite_ID) %>%
               group_split() %>%
               setNames(unique(glycoPSMs$gsite_ID))  # Extract unique gsite_IDs in the same order as group_split()
             
             
             
              gly_GN_list <- glycoPSMs %>%
               dplyr::select(gene_name) %>%  # Ensure we are using dplyr's select
               distinct() %>%
               pull(gene_name) %>%
               list()
                # Convert the vector to a list
             
             
             
             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             # Function to split glycoPSMs into a list of data frames by gene_name
             split_by_gene <- function(data, gene_list) {
               gene_data_list <- lapply(gene_list, function(gene) {
                 data %>% filter(gene_name == gene)
               })
               
               # Naming the list elements with gene names for easy identification
               names(gene_data_list) <- gene_list
               
               return(gene_data_list)
             }
             
             
             
             glycoPSM_list <- split_by_gene(glycoPSMs, gene_list)
             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             Gpeps_compostion <- glycoPSMs %>%
               group_by(sample, glycan_composition) %>%
               summarise(sum_intensity = sum(intensity)) %>%
               pivot_wider(names_from = sample, values_from = sum_intensity)
             
             write_csv(Gpeps_class, "Gpeps_class.csv")
             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             Gpeps_Fuc <- glycoPSMs %>%
               group_by(sample, contains_fuc) %>%
               summarise(sum_intensity = sum(intensity)) %>%
               pivot_wider(names_from = sample, values_from = sum_intensity)
             
             write_csv(Gpeps_Fuc, "Gpeps_Fuc.csv")
             
             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             Gpeps_Sia <- glycoPSMs %>%
               group_by(sample, contains_neu_ac) %>%
               summarise(sum_intensity = sum(intensity)) %>%
               pivot_wider(names_from = sample, values_from = sum_intensity)
             
             write_csv(Gpeps_Sia, "Gpeps_Sia.csv")
             
             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             
             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             Sia_list <- lapply(glycoPSM_list, Sia_table)
             fuc_list <- lapply(glycoPSM_list, fuc_table)
             glycomp_list <- lapply(glycoPSM_list, glycomp_table)
             glyc_class_list <- lapply(glycoPSM_list, class_table)
             fuc_count_list <- lapply(glycoPSM_list, count_fuc_table)
             sia_count_list <- lapply(glycoPSM_list, count_sia_table)
             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             fix_list_names <- function(data_list) {
               names(data_list) <- ifelse(
                 is.na(names(data_list)) | names(data_list) == "",
                 paste0("Sheet_", seq_along(data_list)),
                 names(data_list)
               )
               return(data_list)
             }
             
             # Fix names for all lists
             Sia_list <- fix_list_names(Sia_list)
             glycomp_list <- fix_list_names(glycomp_list)
             glyc_class_list <- fix_list_names(glyc_class_list)
             fuc_count_list <- fix_list_names(fuc_count_list)
             sia_count_list <- fix_list_names(sia_count_list)
             fuc_list <- fix_list_names(fuc_list)
             
             
             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             
             # Function to save a list of data frames into an Excel workbook
             save_list_to_excel <- function(data_list, file_name) {
               # Ensure all elements have valid names
               names(data_list) <- ifelse(
                 is.na(names(data_list)) | names(data_list) == "",
                 paste0("Sheet_", seq_along(data_list)),
                 names(data_list)
               )
               
               # Create a new workbook
               wb <- wb_workbook()
               
               # Add each data frame to a separate sheet in the workbook
               for (name in names(data_list)) {
                 wb$add_worksheet(name)  # Add worksheet with corrected name
                 wb$add_data(sheet = name, x = data_list[[name]])  # Write data
               }
               
               # Save the workbook
               wb$save(file_name)
             }
             
             # Apply the function to multiple lists
             save_list_to_excel(Sia_list, "glyprot_sia.xlsx")
             save_list_to_excel(glycomp_list, "glyprot_glycomp.xlsx")
             save_list_to_excel(glyc_class_list, "glyprot_class.xlsx")
             save_list_to_excel(fuc_count_list, "glyprot_fuc_count.xlsx") # bug with the output of this
             save_list_to_excel(sia_count_list, "glyprot_sia_count.xlsx") # bug with the output of this
             save_list_to_excel(fuc_list, "glyprot_fuc.xlsx")
             
             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             
             unique_glycosites <- glycoPSMs %>%
               group_by(gene_name) %>%
               summarise(unique_glycosites = n_distinct(protein_glycosite))
             
             #'
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             unique_glycoproteins <- glycoPSMs %>%
               group_by(protein_accessions) %>%
               summarise(unique_glycoprot = n_distinct(protein_accessions))
             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             # Ensure gsite_ID is created
             glycoPSMs <- glycoPSMs %>%
               mutate(gsite_ID = paste(protein_accessions, protein_glycosite, sep = "_"))
             
             # Group by gsite_ID and split into a named list
             glycosite_list <- glycoPSMs %>%
               group_by(gsite_ID) %>%
               group_split() %>%
               setNames(unique(glycoPSMs$gsite_ID))  # Extract unique gsite_IDs in the same order as group_split()
             
             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             Sia_glycosite_list <- lapply(glycosite_list, Sia_table)
             fuc_glycosite_list <- lapply(glycosite_list, fuc_table)
             glycomp_glycosite_list <- lapply(glycosite_list, glycomp_table)
             glyc_class_glycosite_list <- lapply(glycosite_list, class_table)
             fuc_count_glycosite_list <- lapply(glycosite_list, count_fuc_table)
             sia_count_glycosite_list <- lapply(glycosite_list, count_sia_table)
             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             # Apply the function to multiple lists
             save_list_to_excel(Sia_glycosite_list, "glycosite_sia.xlsx")
             save_list_to_excel(glycomp_glycosite_list, "glycosite_glycomp.xlsx")
             save_list_to_excel(glyc_class_glycosite_list, "glycosite_class.xlsx")
             save_list_to_excel(fuc_count_glycosite_list, "glycosite_fuc_count.xlsx") # bug with the output of this
             save_list_to_excel(sia_count_glycosite_list, "glycosite_sia_count.xlsx") # bug with the output of this
             save_list_to_excel(fuc_glycosite_list, "glycosite_fuc.xlsx")
             
             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             count_entries_per_df <- function(df_list) {
               if (!is.list(df_list)) {
                 stop("Input must be a list of data frames")
               }
               
               # Ensure all elements have names
               if (is.null(names(df_list))) {
                 names(df_list) <- paste0("DF", seq_along(df_list))
               } else {
                 missing_names <- names(df_list) == "" | is.na(names(df_list))
                 names(df_list)[missing_names] <- paste0("DF", which(missing_names))
               }
               
               # Create a data frame with names and row counts
               result <- data.frame(
                 Name = names(df_list),
                 Row_Count = sapply(df_list, nrow),
                 stringsAsFactors = FALSE
               )
               
               return(result)
             }
             
             # Example usage:
             # df1 <- data.frame(a = 1:5, b = 6:10)
             # df2 <- data.frame(x = 1:3, y = 4:6)
             # df3 <- data.frame(z = 7:9)
             # df_list <- list(df1 = df1, df2 = df2, df3)
             # count_entries_per_df(df_list)
             
             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             gprot_psm_count <- count_entries_per_df(glycoPSM_list)
             gsite_psm_count <- count_entries_per_df(glycosite_list)
             
             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             
            
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             
             # Ensure gsite_list is a vector of distinct gsite_IDs
             gsite_list <- glycoPSMs %>%
               dplyr::select(gsite_ID) %>%  # Ensure we are using dplyr's select
               distinct() %>%
               pull(gsite_ID)  # This returns a vector, not a list
             
             # Function to split glycoPSMs into a list of data frames by gsite_ID
             split_by_gsite <- function(data, gsite_list) {
               gsite_data_list <- lapply(gsite_list, function(gsite) {
                 data %>% filter(gsite_ID == gsite)
               })
               
               # Naming the list elements with gsite_IDs for easy identification
               names(gsite_data_list) <- as.character(gsite_list)  # Ensure that names are character strings
               
               return(gsite_data_list)
             }
             
             # Example usage
             gsitePSM_list <- split_by_gsite(glycoPSMs, gsite_list)
             
             
             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             glycomp_table <- function(input_df) {
               input_df %>%
                 group_by(sample, glycan_composition) %>%
                 summarise(sum_intensity = sum(intensity), .groups = "drop") %>%
                 # Calculate relative abundance after summarization
                 group_by(sample) %>%
                 mutate(
                   total_intensity = sum(sum_intensity),
                   relative_abundance = (sum_intensity / total_intensity) * 100
                 ) %>%
                 ungroup() %>%
                 # Join disease_status to each row based on sample
                 left_join(input_df %>% select(sample, disease_status) %>% distinct(),
                           by = "sample") %>%
                 select(sample,
                        disease_status,
                        glycan_composition,
                        relative_abundance)
             }
             fuc_table <- function(input_df) {
               input_df %>%
                 group_by(sample, contains_fuc) %>%
                 summarise(sum_intensity = sum(intensity), .groups = "drop") %>%
                 # Calculate relative abundance after summarization
                 group_by(sample) %>%
                 mutate(
                   total_intensity = sum(sum_intensity),
                   relative_abundance = (sum_intensity / total_intensity) * 100
                 ) %>%
                 ungroup() %>%
                 # Join disease_status to each row based on sample
                 left_join(input_df %>% select(sample, disease_status) %>% distinct(),
                           by = "sample") %>%
                 select(sample, disease_status, contains_fuc, relative_abundance)
             }
             Sia_table <- function(input_df) {
               input_df %>%
                 group_by(sample, contains_neu_ac) %>%
                 summarise(sum_intensity = sum(intensity), .groups = "drop") %>%
                 # Calculate relative abundance after summarization
                 group_by(sample) %>%
                 mutate(
                   total_intensity = sum(sum_intensity),
                   relative_abundance = (sum_intensity / total_intensity) * 100
                 ) %>%
                 ungroup() %>%
                 # Join disease_status to each row based on sample
                 left_join(input_df %>% select(sample, disease_status) %>% distinct(),
                           by = "sample") %>%
                 select(sample, disease_status, contains_neu_ac, relative_abundance)
             }
             class_table <- function(input_df) {
               input_df %>%
                 group_by(sample, glycan_class) %>%
                 summarise(sum_intensity = sum(intensity), .groups = "drop") %>%
                 # Calculate relative abundance after summarization
                 group_by(sample) %>%
                 mutate(
                   total_intensity = sum(sum_intensity),
                   relative_abundance = (sum_intensity / total_intensity) * 100
                 ) %>%
                 ungroup() %>%
                 # Join disease_status to each row based on sample
                 left_join(input_df %>% select(sample, disease_status) %>% distinct(),
                           by = "sample") %>%
                 select(sample, disease_status, glycan_class, relative_abundance)
             }
             count_sia_table <- function(input_df) {
               input_df %>%
                 group_by(sample, sia_count) %>%
                 summarise(sum_intensity = sum(intensity), .groups = "drop") %>%
                 # Calculate relative abundance after summarization
                 group_by(sample) %>%
                 mutate(
                   total_intensity = sum(sum_intensity),
                   relative_abundance = (sum_intensity / total_intensity) * 100
                 ) %>%
                 ungroup() %>%
                 # Join disease_status to each row based on sample
                 left_join(input_df %>% select(sample, disease_status) %>% distinct(),
                           by = "sample") %>%
                 select(sample, disease_status, sia_count, relative_abundance)
             }
             count_fuc_table <- function(input_df) {
               input_df %>%
                 group_by(sample, fuc_count) %>%
                 summarise(sum_intensity = sum(intensity), .groups = "drop") %>%
                 # Calculate relative abundance after summarization
                 group_by(sample) %>%
                 mutate(
                   total_intensity = sum(sum_intensity),
                   relative_abundance = (sum_intensity / total_intensity) * 100
                 ) %>%
                 ungroup() %>%
                 # Join disease_status to each row based on sample
                 left_join(input_df %>% select(sample, disease_status) %>% distinct(),
                           by = "sample") %>%
                 select(sample, disease_status, fuc_count, relative_abundance)
             }
             
             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             Sia_glycoprot_list <- lapply(glycoPSM_list, Sia_table)
             fuc_glycoprot_list <- lapply(glycoPSM_list, fuc_table)
             glycomp_glycoprot_list <- lapply(glycoPSM_list, glycomp_table)
             glyc_class_glycoprot_list <- lapply(glycoPSM_list, class_table)
             fuc_count_glycoprot_list <- lapply(glycoPSM_list, count_fuc_table)
             sia_count_glycoprot_list <- lapply(glycoPSM_list, count_sia_table)
             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             Sia_glycosite_list <- lapply(glycosite_list, Sia_table)
             fuc_glycosite_list <- lapply(glycosite_list, fuc_table)
             glycomp_glycosite_list <- lapply(glycosite_list, glycomp_table)
             glyc_class_glycosite_list <- lapply(glycosite_list, class_table)
             fuc_count_glycosite_list <- lapply(glycosite_list, count_fuc_table)
             sia_count_glycosite_list <- lapply(glycosite_list, count_sia_table)
             
             #'
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             
            
             
             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             perform_analysis_on_list(Sia_glycoprot_list,
                                      "disease_status",
                                      "relative_abundance",
                                      "Sia_gprot_")
             perform_analysis_on_list(glycomp_glycoprot_list,
                                      "disease_status",
                                      "relative_abundance",
                                      "Glycomp_gprot_")
             perform_analysis_on_list(
               glyc_class_glycoprot_list,
               "disease_status",
               "relative_abundance",
               "glyClass_gprot_"
             )
             perform_analysis_on_list(
               fuc_count_glycoprot_list,
               "disease_status",
               "relative_abundance",
               "Fuc_count_gprot_"
             )
             perform_analysis_on_list(
               sia_count_glycoprot_list,
               "disease_status",
               "relative_abundance",
               "Sia_count_gprot_"
             )
             perform_analysis_on_list(fuc_glycoprot_list,
                                      "disease_status",
                                      "relative_abundance",
                                      "Fuc_gprot_")
             
             perform_analysis_on_list(Sia_glycosite_list,
                                      "disease_status",
                                      "relative_abundance",
                                      "Sia_gsite_")
             perform_analysis_on_list(fuc_glycosite_list,
                                      "disease_status",
                                      "relative_abundance",
                                      "Fuc_gsite_")
             perform_analysis_on_list(glycomp_glycosite_list,
                                      "disease_status",
                                      "relative_abundance",
                                      "Glycomp_gsite_")
             perform_analysis_on_list(
               glyc_class_glycosite_list,
               "disease_status",
               "relative_abundance",
               "glyClass_gsite_"
             )
             perform_analysis_on_list(
               fuc_count_glycosite_list,
               "disease_status",
               "relative_abundance",
               "Fuc_count_gsite_"
             )
             perform_analysis_on_list(
               sia_count_glycosite_list,
               "disease_status",
               "relative_abundance",
               "Sia_count_gsite_"
             )
 
             
             
#-----------------------------------------------------------------------
             

relative_abundance <- function(df, protein_list, sample_list, protein_col, sample_col, intensity_col) {
             # Load necessary library
               library(dplyr)
               library(tidyr)
               
               # Initialize an empty list to store results
               results_list <- list()
               
               # Loop over each protein in the protein list
               for (protein in protein_list) {
                 # Subset data for the current protein
                 protein_data <- df %>% filter(.data[[protein_col]] == protein)
                 
                 # Initialize an empty data frame for storing sample results
                 sample_results <- data.frame(Sample = sample_list, Relative_Abundance = NA)
                 
                 # Loop over each sample in the sample list
                 for (sample in sample_list) {
                   # Subset data for the current sample
                   sample_data <- protein_data %>% filter(.data[[sample_col]] == sample)
                   
                   # If no data for the sample, keep NA
                   if (nrow(sample_data) == 0) {
                     next
                   }
                   
                   # Sum intensities for each unique grouping
                   group_sums <- sample_data %>%
                     group_by(across(-all_of(intensity_col))) %>%  # Group by all except intensity
                     summarise(Summed_Intensity = sum(.data[[intensity_col]], na.rm = TRUE), .groups = "drop")
                   
                   # Total intensity for this sample and protein
                   total_intensity <- sum(group_sums$Summed_Intensity, na.rm = TRUE)
                   
                   # Compute relative abundance
                   if (total_intensity > 0) {
                     relative_abundance_value <- (sum(group_sums$Summed_Intensity) / total_intensity) * 100
                     sample_results$Relative_Abundance[sample_results$Sample == sample] <- relative_abundance_value
                   }
                 }
                 
                 # Store results for this protein
                 results_list[[protein]] <- sample_results$Relative_Abundance
               }
               
               # Convert the list into a wide-format data frame
               results_df <- as.data.frame(do.call(rbind, results_list))
               colnames(results_df) <- sample_list  # Set column names as samples
               results_df$Protein <- rownames(results_df)  # Add protein names as a column
               rownames(results_df) <- NULL  # Reset row names
               
               # Reorder columns to have Protein as the first column
               results_df <- results_df %>% relocate(Protein, .before = everything())
               
               return(results_df)
             }
             
             # Example usage
             df <- data.frame(
               Protein = c("P1", "P1", "P1", "P2", "P2", "P3"),
               Sample = c("S1", "S2", "S3", "S1", "S2", "S1"),
               Intensity = c(10, 20, 30, 40, 50, 60)
             )
             
             protein_list <- c("P1", "P2", "P3")
             sample_list <- c("S1", "S2", "S3")
             
             relative_abundance(df, protein_list, sample_list, "Protein", "Sample", "Intensity")
             
             
           
             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             DDA_DAPs <- ProtAbun_t_test_results %>%
               filter(Log2_FC > 0.5 & `-Log10_p.value` < 0.05)
             
             write.csv(DIA_DAPs, file = "DDA_DAPs.csv")
             
             
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             
             # Create a basic volcano plot
             ggplot(data = ProtAbun_t_test_results, aes(x = Log2_FC, y = `-Log10_p.value`)) +
               geom_vline(
                 xintercept = c(-1.2, 1.2),
                 col = "gray",
                 linetype = 'dashed'
               ) +
               geom_hline(
                 yintercept = -log10(0.3),
                 col = "gray",
                 linetype = 'dashed'
               ) +
               geom_point() +
               theme()
             
             
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             
             ProtAbun_t_test_results <- ProtAbun_t_test_results %>%
               filter(!is.na(Log2_FC) & !is.na(p.adjusted))
             
             ProtAbun_t_test_results$p.adjusted[ProtAbun_t_test_results$p.adjusted == 0] <- 1e-300
             
             Gprot_volcano <- EnhancedVolcano(
               ProtAbun_t_test_results,
               lab = ProtAbun_t_test_results$Protein_Names,
               x = 'Log2_FC',
               y = 'p.adjusted',
               xlim = c(-5, 5),
               # Adjust the range based on your data distribution
               ylim = c(
                 0,
                 max(ProtAbun_t_test_results$`-Log10_p.value`, na.rm = TRUE) + 0.08
               ),
               # Set y-axis limits
               pCutoff = 0.3,
               # Optional: Set a cutoff for significance
               FCcutoff = 1,
               # Optional: Set a cutoff for fold change
               labSize = 3.0,
               # Optional: Adjust label size
               labCol = 'black',
               colAlpha = 1,
               legendPosition = 'none',
               legendLabSize = 12,
               legendIconSize = 4.0,
               drawConnectors = TRUE,
               widthConnectors = 0.75
             )
             
             ggsave("Gprot_volcano.png", width = 10, height = 8)
             knitr::include_graphics("Gprot_volcano.png")
             
             #'
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             
             
             #' Have done differenial protein aubance analysis of the Charlie's DIA protemics data, as this is more correct dataset to perform analysis on then compared to my own glycoproteomics data set. this is because the glycoproteomics datasets has been biased by glycopeptide enrichment, therefore this may lead to inaccurate quants due to possible greater llels of missing values . the glycoproteomics data set also suffer from being collected by a DDA exmerimental method as this is best for glycoproteoms as the chimeric spectra would lead to signifcant issues for glycan asigments to peptides. hthe DDA method further creates issues of missing values at random due to ticinal featues inherant in the method.
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             wb <- wb_load(
               "/Users/thomasreilly/Desktop/MRes Data/10S_PD_quant/Charlie_DIA/MRCFS Proteomics Results.xlsx"
             )
             DIA_Prot_DiffAb <- wb %>% wb_to_df(sheet = "Comp1")
             library(janitor)
             DIA_Prot_DiffAb <- DIA_Prot_DiffAb %>% clean_names()
             
             DIA_Prot_DiffAb <- DIA_Prot_DiffAb %>%
               select(-starts_with("na"))
             
             DIA_Prot_DiffAb <- DIA_Prot_DiffAb %>%
               mutate(Log2_FC = log2(fc_proteins))
             
             
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             DIA_prot_volcano <- EnhancedVolcano(
               DIA_Prot_DiffAb,
               lab = DIA_Prot_DiffAb$protein,
               title = "DIA proteomics differentially aubundant proteins",
               x = 'Log2_FC',
               y = 'pval_proteins',
               ylim = c(0, 6),
               # Adjust this based on your data
               xlim = c(-1.2, 1.2),
               # Adjust the range based on your data distribution
               pCutoff = 0.05,
               # Optional: Set a cutoff for significance
               FCcutoff = 0.5,
               # Optional: Set a cutoff for fold change
               labSize = 3.0,
               # Optional: Adjust label size
               labCol = 'black',
               colAlpha = 1,
               legendPosition = 'right',
               legendLabSize = 12,
               legendIconSize = 4.0,
               drawConnectors = TRUE,
               widthConnectors = 0.75
             )
             
             ggsave("DIA_prot_volcano.png", width = 10, height = 8)
             knitr::include_graphics("DIA_prot_volcano.png")
             
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             knitr::kable(DIA_Prot_DiffAb,
                          caption = "DIA Protein Abundances",
                          digits = 2,
                          na = 'NA')
             
             
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             DIA_DAPs <- DIA_Prot_DiffAb %>%
               filter(Log2_FC > 0.5 & pval_proteins < 0.05)
             
             write.csv(DIA_DAPs, file = "DIA_DAPs.csv")
             
             knitr::kable(DIA_DAPs,
                          caption = "DIA plasma proteomics - Differentaly Aubundant Proteins",
                          digits = 2,
                          na = 'NA')
             
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             DIA_DAPs %>%
               summarise(n_distinct(protein))
             
             #'
             #-------------------------------------------------------------------------------------------------------
             
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
             
  #' Log Transform Matrix Data with NA handling
#' 
#' @param data_matrix A numeric matrix containing the data to be transformed
#' @param pseudo_count A small number to add before log transformation to handle zeros (default: 1)
#' @param base The logarithm base to use (default: 2)
#' @param replace_inf Whether to replace infinite values with NA (default: TRUE)
#' @param na_handling Strategy for handling NA values: 
#'        "keep" (default) - keep NAs as is
#'        "remove" - remove rows with any NAs
#'        "impute_min" - replace NAs with minimum non-NA value in dataset
#'        "impute_zero" - replace NAs with zero before adding pseudo count
#' @param min_non_na Minimum number of non-NA values required per row to keep row (default: 1)
#' 
#' @return A matrix of the same dimensions as input with log-transformed values
#' 
log_transform_matrix <- function(data_matrix, 
                               pseudo_count = 1, 
                               base = 2,
                               replace_inf = TRUE,
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
  transformed_matrix <- data_matrix + pseudo_count
  
  # Perform log transformation
  transformed_matrix <- switch(as.character(base),
    "2" = log2(transformed_matrix),
    "10" = log10(transformed_matrix),
    "2.718282" = log(transformed_matrix),  # exp(1)
    log(transformed_matrix, base)
  )
  
  # Replace infinite values with NA if requested
  if (replace_inf) {
    transformed_matrix[is.infinite(transformed_matrix)] <- NA
  }
  
  return(transformed_matrix)
}

# Log transform and write protein matrices to CSV
protein_gly_comp_log <- log_transform_matrix(protein_gly_comp)
write.csv(protein_gly_comp_log, file = "output_data/protein_gly_comp_log.csv")

protein_gly_fuc_log <- log_transform_matrix(protein_gly_fuc)  # Changed from protein_fuc to protein_gly_fuc
write.csv(protein_gly_fuc_log, file = "output_data/protein_gly_fuc_log.csv")

protein_gly_sia_log <- log_transform_matrix(protein_gly_sia)  # Changed from protein_sia to protein_gly_sia
write.csv(protein_gly_sia_log, file = "output_data/protein_gly_sia_log.csv")

protein_gly_sia_count_log <- log_transform_matrix(protein_gly_sia_count)
write.csv(protein_gly_sia_count_log, file = "output_data/protein_gly_sia_count_log.csv")

# Log transform and write glycosite matrices to CSV
glycosite_gly_comp_log <- log_transform_matrix(glycosite_gly_comp)
write.csv(glycosite_gly_comp_log, file = "output_data/glycosite_gly_comp_log.csv")

glycosite_fuc_log <- log_transform_matrix(glycosite_gly_fuc)
write.csv(glycosite_fuc_log, file = "output_data/glycosite_gly_fuc_log.csv")

glycosite_sia_log <- log_transform_matrix(glycosite_gly_sia)
write.csv(glycosite_sia_log, file = "output_data/glycosite_gly_sia_log.csv")

glycosite_sia_count_log <- log_transform_matrix(glycosite_gly_sia_count)
write.csv(glycosite_sia_count_log, file = "output_data/glycosite_gly_sia_count_log.csv")

#' Centered Log-Ratio (CLR) Transform Matrix Data
#' 
#' @param data_matrix A numeric matrix containing the data to be transformed
#' @param pseudo_count A small number to add before log transformation to handle zeros (default: 1)
#' @param base The logarithm base to use (default: 2)
#' @param na_handling Strategy for handling NA values: 
#'        "keep" (default) - keep NAs as is
#'        "remove" - remove rows with any NAs
#'        "impute_min" - replace NAs with minimum non-NA value in dataset
#'        "impute_zero" - replace NAs with zero before adding pseudo count
#' @param min_non_na Minimum number of non-NA values required per row to keep row (default: 1)
#' 
#' @return A matrix of the same dimensions as input with CLR-transformed values
#' 
clr_transform_matrix <- function(data_matrix,
                               pseudo_count = 1,
                               base = 2,
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
  
  # Function to calculate geometric mean, handling NAs
  geometric_mean <- function(x) {
    exp(mean(log(x), na.rm = TRUE))
  }
  
  # Initialize matrix for CLR transformed values
  clr_matrix <- matrix(NA, 
                      nrow = nrow(data_matrix), 
                      ncol = ncol(data_matrix),
                      dimnames = dimnames(data_matrix))
  
  # Perform CLR transformation row by row
  for (i in 1:nrow(data_matrix)) {
    row_data <- data_matrix[i, ]
    
    # Skip rows with all NAs
    if (all(is.na(row_data))) {
      next
    }
    
    # Calculate geometric mean for non-NA values
    geom_mean <- geometric_mean(row_data)
    
    # Perform CLR transformation
    if (base == 2) {
      clr_matrix[i, ] <- log2(row_data / geom_mean)
    } else if (base == 10) {
      clr_matrix[i, ] <- log10(row_data / geom_mean)
    } else if (base == exp(1)) {
      clr_matrix[i, ] <- log(row_data / geom_mean)
    } else {
      clr_matrix[i, ] <- log(row_data / geom_mean, base)
    }
  }
  
  return(clr_matrix)
}

# Using your existing matrices directly:
protein_gly_comp_clr <- clr_transform_matrix(protein_gly_comp)
write.csv(protein_gly_comp_clr, file = "output_data/protein_gly_comp_clr.csv")
glycosite_gly_comp_clr <- clr_transform_matrix(glycosite_gly_comp)
write.csv(glycosite_gly_comp_clr, file = "output_data/glycosite_gly_comp_clr.csv")

#------------------------------------------------------------------------------------------------


#' Analyze CLR-Transformed Data and Create Volcano Plots
#' 
#' This function performs statistical analysis on CLR-transformed data, comparing two groups 
#' (typically Healthy vs ME/CFS) and generates volcano plots to visualize the results.
#' 
#' @param clr_matrix A numeric matrix containing CLR-transformed data. Rows represent 
#'        features and columns represent samples. Column names should start with "HC" 
#'        for healthy controls or "M" for ME/CFS samples.
#' @param title Character string for the plot title (default: "Volcano Plot")
#' @param fc_cutoff Fold change cutoff for significance (default: 0.5)
#' @param p_cutoff P-value cutoff for significance (default: 0.05)
#'
#' @details The function performs the following steps:
#' \itemize{
#'   \item Automatically assigns samples to groups based on column name prefixes
#'   \item Calculates sample counts and verifies minimum group sizes
#'   \item Computes means for each group
#'   \item Performs t-tests when sufficient samples are available
#'   \item Calculates fold changes and adjusts p-values for multiple testing
#'   \item Generates an enhanced volcano plot
#' }
#'
#' @return A list containing three elements:
#' \itemize{
#'   \item results: Data frame with statistical results including:
#'     \itemize{
#'       \item Feature: Feature identifiers
#'       \item P_Value: Raw p-values from t-tests
#'       \item FDR: False Discovery Rate adjusted p-values
#'       \item Log2_FC: Log2 fold changes
#'       \item Mean_Healthy: Mean values for healthy group
#'       \item Mean_MECFS: Mean values for ME/CFS group
#'       \item N_Healthy: Number of non-NA samples in healthy group
#'       \item N_MECFS: Number of non-NA samples in ME/CFS group
#'     }
#'   \item plot: ggplot2 object containing the volcano plot
#'   \item sample_groups: Data frame showing sample to group assignments
#' }
#'
#' @note 
#' - Requires at least 2 samples per group for statistical analysis
#' - NA values are handled using na.rm = TRUE in calculations
#' - P-values are adjusted using Benjamini-Hochberg method
#'
#' @examples
#' \dontrun{
#' # Create example matrix with proper column names
#' test_matrix <- matrix(rnorm(100), nrow=10)
#' colnames(test_matrix) <- c("HC1", "HC2", "HC3", "M1", "M2", "M3")
#' 
#' # Run analysis
#' results <- analyze_clr_data(test_matrix, 
#'                            title = "Test Analysis",
#'                            fc_cutoff = 0.5,
#'                            p_cutoff = 0.05)
#' 
#' # View results
#' head(results$results)
#' print(results$plot)
#' }
#'
#' @importFrom stats t.test p.adjust
#' @importFrom EnhancedVolcano EnhancedVolcano
#'
#' @export
analyze_clr_data <- function(clr_matrix, 
                           title = "Volcano Plot",
                           fc_cutoff = 0.5,
                           p_cutoff = 0.05) {
  
  # Create sample group mapping based on your column names
  sample_groups <- data.frame(
    sample = colnames(clr_matrix),
    group = ifelse(grepl("^HC", colnames(clr_matrix)), "Healthy",
                  ifelse(grepl("^M", colnames(clr_matrix)), "MECFS", NA))
  )
  
  # Print sample grouping for verification
  cat("\nSample Grouping for", title, ":\n")
  print(sample_groups)
  
  # Check for any unassigned samples
  if(any(is.na(sample_groups$group))) {
    warning("Some samples could not be assigned to groups: ", 
            paste(sample_groups$sample[is.na(sample_groups$group)], collapse=", "))
  }
  
  # Verify we have samples in both groups
  group_counts <- table(sample_groups$group)
  cat("\nSamples per group:\n")
  print(group_counts)
  
  if(length(group_counts) < 2 || any(group_counts < 2)) {
    stop("Need at least 2 samples in each group for comparison")
  }
  
  # Initialize results dataframe
  results <- data.frame(
    Feature = rownames(clr_matrix),
    P_Value = NA,
    FDR = NA,
    Log2_FC = NA,
    Mean_Healthy = NA,
    Mean_MECFS = NA,
    N_Healthy = NA,
    N_MECFS = NA,
    stringsAsFactors = FALSE
  )
  
  # Calculate statistics for each feature
  for(i in 1:nrow(clr_matrix)) {
    # Extract data for current feature
    feature_data <- clr_matrix[i,]
    
    # Split data by group
    mecfs_data <- feature_data[sample_groups$group == "MECFS"]
    healthy_data <- feature_data[sample_groups$group == "Healthy"]
    
    # Store sample counts
    results$N_Healthy[i] <- sum(!is.na(healthy_data))
    results$N_MECFS[i] <- sum(!is.na(mecfs_data))
    
    # Calculate means
    results$Mean_Healthy[i] <- mean(healthy_data, na.rm = TRUE)
    results$Mean_MECFS[i] <- mean(mecfs_data, na.rm = TRUE)
    
    # Perform t-test if we have enough samples
    if(results$N_Healthy[i] >= 3 && results$N_MECFS[i] >= 3) {
      t_test <- try({
        t.test(mecfs_data, healthy_data)
      })
      
      if(!inherits(t_test, "try-error")) {
        results$P_Value[i] <- t_test$p.value
        results$Log2_FC[i] <- results$Mean_MECFS[i] - results$Mean_Healthy[i]
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
      'CLR-transformed data\n',
      'Healthy n=', group_counts["Healthy"], 
      ', ME/CFS n=', group_counts["MECFS"], '\n',
      'FC cutoff = ', fc_cutoff, 
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
    gridlines.minor = FALSE
  )
  
  # Return both results and plot
  return(list(
    results = results,
    plot = volcano_plot,
    sample_groups = sample_groups
  ))
}

# Function to verify sample names match expected pattern
verify_sample_names <- function(matrix_list) {
  for(matrix_name in names(matrix_list)) {
    cat("\nChecking sample names for", matrix_name, ":\n")
    sample_names <- colnames(matrix_list[[matrix_name]])
    
    # Check if sample names follow expected pattern
    valid_names <- grepl("^(HC|M)", sample_names)
    
    if(!all(valid_names)) {
      warning("Invalid sample names found in ", matrix_name, ":\n",
              paste(sample_names[!valid_names], collapse=", "))
    }
    
    # Print sample grouping
    groups <- data.frame(
      Sample = sample_names,
      Group = ifelse(grepl("^HC", sample_names), "Healthy",
                    ifelse(grepl("^M", sample_names), "MECFS", "Unknown"))
    )
    print(groups)
  }
}

# Verify sample names before analysis
verify_sample_names(transformed_results$clr_transformed)

# Run analysis
results <- process_all_clr_data(transformed_results$clr_transformed)

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

# Analyze protein glycan composition
protein_comp_results <- analyze_log_data(protein_gly_comp_log,
                                       title = "Protein Glycan Composition")

# Analyze glycosite composition
glycosite_comp_results <- analyze_log_data(glycosite_gly_comp_log,
                                         title = "Glycosite Composition")


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




#' Analyze Glycan Relative Abundance Matrix Data
#' 
#' This function performs statistical analysis on glycan relative abundance data,
#' comparing two groups and generating appropriate visualizations.
#' 
#' @param abundance_matrix A numeric matrix containing relative abundance data (in percentages).
#'        Rows represent glycofeatures and columns represent samples.
#' @param title Character string for the plot title (default: "Glycan Analysis")
#' @param fc_cutoff Fold change cutoff for significance (default: 1.5)
#' @param p_cutoff P-value cutoff for significance (default: 0.05)
#' @param group_pattern Pattern for identifying groups in column names (default: c("HC", "M"))
#' @param group_names Names for the groups (default: c("Healthy", "MECFS"))
#' @param min_abundance Minimum relative abundance threshold (default: 1%)
#' @param min_samples Minimum number of samples required per group (default: 3)
#'
#' @return A list containing:
#'   - results: Data frame with statistical results
#'   - volcano_plot: EnhancedVolcano plot
#'   - boxplot: ggplot2 boxplot for significant features
#'   - heatmap: pheatmap object showing abundance patterns
#'   - sample_groups: Group assignments
#'
#' @export
analyze_glycomatrix <- function(abundance_matrix,
                              title = "Glycan Analysis",
                              fc_cutoff = 1.5,
                              p_cutoff = 0.05,
                              group_pattern = c("HC", "M"),
                              group_names = c("Healthy", "MECFS"),
                              min_abundance = 1,
                              min_samples = 3) {
  
  # Required packages
  require(ggplot2)
  require(EnhancedVolcano)
  require(pheatmap)
  require(dplyr)
  
  # Input validation and conversion
  if (!is.matrix(abundance_matrix)) {
    abundance_matrix <- as.matrix(abundance_matrix)
  }
  
  # Add check for empty matrix
  if (nrow(abundance_matrix) == 0 || ncol(abundance_matrix) == 0) {
    stop("Empty matrix provided")
  }
  
  # Convert matrix to numeric, preserving row and column names
  original_rownames <- rownames(abundance_matrix)
  original_colnames <- colnames(abundance_matrix)
  
  # More robust numeric conversion with error handling
  abundance_matrix <- tryCatch({
    apply(abundance_matrix, 2, function(x) {
      as.numeric(gsub("[^0-9.-]", "", as.character(x)))
    })
  }, error = function(e) {
    stop("Failed to convert matrix values to numeric: ", e$message)
  })
  
  # Restore dimension names
  rownames(abundance_matrix) <- original_rownames
  colnames(abundance_matrix) <- original_colnames
  
  # Check for and report conversion issues
  na_count <- sum(is.na(abundance_matrix))
  if (na_count > 0) {
    warning(sprintf("%d values could not be converted to numeric and were set to NA", na_count))
  }
  
  # Modified t-test function that handles constant data
  safe_ttest <- function(group2_data, group1_data) {
    # Remove NA values
    group2_data <- group2_data[!is.na(group2_data)]
    group1_data <- group1_data[!is.na(group1_data)]
    
    # Check for constant data
    if (length(unique(c(group1_data, group2_data))) <= 1) {
      return(list(p.value = NA, warning = "Constant data"))
    }
    
    # Check for sufficient data points
    if (length(group1_data) < 2 || length(group2_data) < 2) {
      return(list(p.value = NA, warning = "Insufficient data points"))
    }
    
    # Perform t-test
    tryCatch({
      test_result <- t.test(group2_data, group1_data, var.equal = FALSE)
      return(list(p.value = test_result$p.value, warning = NULL))
    }, error = function(e) {
      return(list(p.value = NA, warning = e$message))
    })
  }
  
  # Create sample group mapping
  sample_groups <- data.frame(
    sample = colnames(abundance_matrix),
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
  
  if(length(group_counts) < 2 || any(group_counts < min_samples)) {
    stop(paste("Need at least", min_samples, "samples in each group for comparison"))
  }
  
  # Filter low abundance features
  row_means <- rowMeans(abundance_matrix, na.rm = TRUE)
  abundance_matrix <- abundance_matrix[row_means >= min_abundance, ]
  
  if(nrow(abundance_matrix) == 0) {
    stop("No features remain after abundance filtering")
  }
  
  # Initialize results dataframe
  results <- data.frame(
    Feature = rownames(abundance_matrix),
    P_Value = NA_real_,
    FDR = NA_real_,
    Fold_Change = NA_real_,
    Log2_FC = NA_real_,
    Mean_Group1 = NA_real_,
    Mean_Group2 = NA_real_,
    SD_Group1 = NA_real_,
    SD_Group2 = NA_real_,
    N_Group1 = NA_integer_,
    N_Group2 = NA_integer_,
    Significant = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Modified statistics calculation
  for(i in 1:nrow(abundance_matrix)) {
    feature_data <- abundance_matrix[i,]
    
    # Split data by group
    group1_data <- feature_data[sample_groups$group == group_names[1]]
    group2_data <- feature_data[sample_groups$group == group_names[2]]
    
    # Store sample counts (excluding NA values)
    results$N_Group1[i] <- sum(!is.na(group1_data))
    results$N_Group2[i] <- sum(!is.na(group2_data))
    
    # Calculate means and SDs with NA handling
    results$Mean_Group1[i] <- mean(group1_data, na.rm = TRUE)
    results$Mean_Group2[i] <- mean(group2_data, na.rm = TRUE)
    results$SD_Group1[i] <- sd(group1_data, na.rm = TRUE)
    results$SD_Group2[i] <- sd(group2_data, na.rm = TRUE)
    
    # Calculate fold change and log2 fold change with improved safety checks
    if (!is.na(results$Mean_Group1[i]) && 
        !is.na(results$Mean_Group2[i]) && 
        results$Mean_Group1[i] > 0 && 
        results$Mean_Group2[i] > 0) {
      results$Fold_Change[i] <- results$Mean_Group2[i] / results$Mean_Group1[i]
      results$Log2_FC[i] <- log2(results$Fold_Change[i])
    } else {
      results$Fold_Change[i] <- NA
      results$Log2_FC[i] <- NA
    }
    
    # Perform t-test if we have enough samples
    if(!is.na(results$N_Group1[i]) && 
       !is.na(results$N_Group2[i]) && 
       results$N_Group1[i] >= min_samples && 
       results$N_Group2[i] >= min_samples) {
      t_test_result <- safe_ttest(group2_data, group1_data)
      results$P_Value[i] <- t_test_result$p.value
      if (!is.null(t_test_result$warning)) {
        warning(sprintf("Feature %s: %s", results$Feature[i], t_test_result$warning))
      }
    }
  }
  
  # Calculate FDR
  results$FDR <- p.adjust(results$P_Value, method = "BH")
  
  # Mark significant features
  results$Significant <- results$FDR < p_cutoff & abs(results$Log2_FC) > log2(fc_cutoff)
  
  # Create volcano plot
  volcano_plot <- EnhancedVolcano(
    results,
    lab = results$Feature,
    x = 'Log2_FC',
    y = 'FDR',
    title = paste(title, "- Volcano Plot"),
    subtitle = paste0(
      'Relative abundance analysis\n',
      group_names[1], ' n=', group_counts[group_names[1]], 
      ', ', group_names[2], ' n=', group_counts[group_names[2]], '\n',
      'FC cutoff = ±', fc_cutoff, 
      '; p-value cutoff = ', p_cutoff
    ),
    pCutoff = p_cutoff,
    FCcutoff = log2(fc_cutoff),
    pointSize = 2.0,
    labSize = 3.0,
    col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
    colAlpha = 0.5,
    legendPosition = 'right',
    legendLabSize = 10,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 0.5
  )
  
  # Create boxplot for significant features
  sig_features <- results$Feature[results$Significant]
  if(length(sig_features) > 0) {
    plot_data <- data.frame()
    for(feature in sig_features) {
      # Add check to ensure feature exists in matrix
      if(!feature %in% rownames(abundance_matrix)) {
        warning(sprintf("Feature %s not found in abundance matrix - skipping", feature))
        next
      }
      
      tryCatch({
        feature_data <- data.frame(
          Feature = feature,
          Abundance = as.numeric(abundance_matrix[feature,]),
          Sample = colnames(abundance_matrix)
        )
        feature_data$Group <- sample_groups$group[match(feature_data$Sample, sample_groups$sample)]
        plot_data <- rbind(plot_data, feature_data)
      }, error = function(e) {
        warning(sprintf("Error processing feature %s: %s", feature, e$message))
      })
    }
    
    if(nrow(plot_data) > 0) {
      boxplot <- ggplot(plot_data, aes(x = Group, y = Abundance, fill = Group)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        facet_wrap(~Feature, scales = "free_y") +
        theme_bw() +
        labs(title = paste(title, "- Significant Features"),
             y = "Relative Abundance (%)",
             x = "Group") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    } else {
      boxplot <- NULL
      warning("No valid data available for boxplot")
    }
  } else {
    boxplot <- NULL
    warning("No significant features found for boxplot")
  }
  
  # Create heatmap with additional error handling
  tryCatch({
    # Remove rows with all NAs
    valid_rows <- rowSums(!is.na(abundance_matrix)) > 0
    heatmap_matrix <- abundance_matrix[valid_rows, , drop = FALSE]
    
    # Remove rows with infinite values
    valid_rows <- apply(heatmap_matrix, 1, function(x) !any(is.infinite(x)))
    heatmap_matrix <- heatmap_matrix[valid_rows, , drop = FALSE]
    
    # Check if we have any data left
    if(nrow(heatmap_matrix) == 0) {
      warning("No valid data remaining for heatmap after filtering NA/Inf values")
      heatmap <- NULL
    } else {
      # Create annotation
      annotation_col <- data.frame(
        Group = sample_groups$group,
        row.names = sample_groups$sample
      )
      
      # Create heatmap with try-catch
      heatmap <- try({
        pheatmap(
          heatmap_matrix,
          annotation_col = annotation_col,
          scale = "row",
          clustering_method = "ward.D2",
          show_colnames = FALSE,
          main = paste(title, "- Heatmap"),
          silent = TRUE,
          na_col = "grey"  # Color NA values as grey
        )
      })
      
      # Check if heatmap creation failed
      if(inherits(heatmap, "try-error")) {
        warning("Failed to create heatmap: ", heatmap)
        heatmap <- NULL
      }
    }
  }, error = function(e) {
    warning("Error creating heatmap: ", e$message)
    heatmap <- NULL
  })
  
  # Return results
  return(list(
    results = results,
    volcano_plot = volcano_plot,
    boxplot = boxplot,
    heatmap = heatmap,
    sample_groups = sample_groups
  ))
}

# Example usage:
# results <- analyze_glycomatrix(protein_gly_comp,
#                              title = "Protein Glycan Composition",
#                              fc_cutoff = 1.5,
#                              p_cutoff = 0.05,
#                              min_abundance = 1,
#                              min_samples = 3
# )
#
# # Save results
# write.csv(results$results, "glycan_composition_analysis.csv")
#
# # Save plots
# ggsave("glycan_volcano.png", results$volcano_plot, width = 10, height = 8)
# ggsave("glycan_boxplot.png", results$boxplot, width = 12, height = 8)
# pdf("glycan_heatmap.pdf", width = 10, height = 12)
# print(results$heatmap)
# dev.off()


# Analyze protein glycan composition
protein_glycan_results <- analyze_glycomatrix(
  protein_gly_comp,
  title = "Protein Glycan Composition",
  fc_cutoff = 1.5,
  p_cutoff = 0.05,
  min_abundance = 1,
  min_samples = 3
)

# Analyze glycosite composition
glycosite_results <- analyze_glycomatrix(
  glycosite_gly_comp,
  title = "Glycosite Composition",
  fc_cutoff = 1.5,
  p_cutoff = 0.05,
  min_abundance = 1,
  min_samples = 3
)

# Save protein glycan results and plots
write.csv(protein_glycan_results$results, 
          "output_data/protein_glycan_comp_analysis.csv")

ggsave("figures/protein_glycan_comp_volcano.png", 
       protein_glycan_results$volcano_plot, 
       width = 10, 
       height = 8)

if (!is.null(protein_glycan_results$boxplot)) {
  ggsave("figures/protein_glycan_comp_boxplot.png", 
         protein_glycan_results$boxplot, 
         width = 12, 
         height = 8)
}

# Save heatmap with error checking
if (!is.null(protein_glycan_results$heatmap)) {
  pdf("figures/protein_glycan_comp_heatmap.pdf", 
      width = 10, 
      height = 12)
  print(protein_glycan_results$heatmap)
  dev.off()
} else {
  warning("Heatmap could not be created for protein glycan composition")
}

# Save glycosite results and plots
write.csv(glycosite_results$results, 
          "output_data/glycosite_comp_analysis.csv")

ggsave("figures/glycosite_comp_volcano.png", 
       glycosite_results$volcano_plot, 
       width = 10, 
       height = 8)

if (!is.null(glycosite_results$boxplot)) {
  ggsave("figures/glycosite_comp_boxplot.png", 
         glycosite_results$boxplot, 
         width = 12, 
         height = 8)
}

# Save heatmap with error checking
if (!is.null(glycosite_results$heatmap)) {
  pdf("figures/glycosite_comp_heatmap.pdf", 
      width = 10, 
      height = 12)
  print(glycosite_results$heatmap)
  dev.off()
} else {
  warning("Heatmap could not be created for glycosite composition")
}


# Analyze protein glycan composition
protein_sia_count_results <- analyze_glycomatrix(
  protein_gly_sia_count,
  title = "Protein Glycan Sia Count",
  fc_cutoff = 1.5,
  p_cutoff = 0.05,
  min_abundance = 1,
  min_samples = 3
)
# Save protein glycan results and plots
write.csv(protein_sia_count_results$results, 
          "output_data/protein_glycan_sia_count_analysis.csv")

ggsave("figures/protein_glycan_sia_count_volcano.png", 
       protein_sia_count_results$volcano_plot, 
       width = 10, 
       height = 8)

if (!is.null(protein_glycan_results$boxplot)) {
  ggsave("figures/protein_glycan_sia_count_boxplot.png", 
         protein_sia_count_results$boxplot, 
         width = 12, 
         height = 8)
}

# Save heatmap with error checking
if (!is.null(protein_sia_count_results$heatmap)) {
  pdf("figures/protein_glycan_sia_count_heatmap.pdf", 
      width = 10, 
      height = 12)
  print(protein_sia_count_results$heatmap)
  dev.off()
} else {
  warning("Heatmap could not be created for protein glycan sia count")
}



# Analyze glycosite fucose
glycosite_fuc_results <- analyze_glycomatrix(
  glycosite_gly_fuc,
  title = "Glycosite Fucose",
  fc_cutoff = 1.5,
  p_cutoff = 0.05,
  min_abundance = 1,
  min_samples = 3
)

# Save glycosite fucose results and plots
write.csv(glycosite_fuc_results$results, 
          "output_data/glycosite_fucose_analysis.csv")

ggsave("figures/glycosite_fucose_volcano.png", 
       glycosite_fuc_results$volcano_plot, 
       width = 10, 
       height = 8)

if (!is.null(glycosite_fuc_results$boxplot)) {
  ggsave("figures/glycosite_fucose_boxplot.png", 
         glycosite_fuc_results$boxplot, 
         width = 12, 
         height = 8)
}

# Save heatmap with error checking
if (!is.null(glycosite_fuc_results$heatmap)) {
  pdf("figures/glycosite_fucose_heatmap.pdf", 
      width = 10, 
      height = 12)
  print(glycosite_fuc_results$heatmap)
  dev.off()
} else {
  warning("Heatmap could not be created for glycosite fucose")
}



#---------------------------------------------------
#' Create Heatmap from Log-Transformed Data
#' 
#' @param log_matrix A numeric matrix containing log-transformed data
#' @param title Character string for the plot title
#' @param group_pattern Pattern for identifying groups in column names (default: c("HC", "M"))
#' @param group_names Names for the groups (default: c("Healthy", "MECFS"))
#' @param scale Whether to scale the rows ("row"), columns ("column") or none ("none")
#' @param cluster_rows Boolean indicating whether to cluster rows (default: TRUE)
#' @param cluster_cols Boolean indicating whether to cluster columns (default: TRUE)
#' @param show_rownames Boolean indicating whether to show row names (default: TRUE)
#' @param show_colnames Boolean indicating whether to show column names (default: FALSE)
#' @param filename Output filename for the PDF (default: NULL, no file output)
#'
#' @return pheatmap object
#' @export
create_glyco_heatmap <- function(log_matrix,
                                title = "Glycan Analysis Heatmap",
                                group_pattern = c("HC", "M"),
                                group_names = c("Healthy", "MECFS"),
                                scale = "row",
                                cluster_rows = TRUE,
                                cluster_cols = TRUE,
                                show_rownames = TRUE,
                                show_colnames = FALSE,
                                filename = NULL) {
  
  # Required packages
  require(pheatmap)
  require(viridis)
  require(grid)
  
  # Input validation
  if (!is.matrix(log_matrix)) {
    stop("Input must be a matrix")
  }
  
  # Print initial dimensions
  cat("Initial matrix dimensions:", dim(log_matrix), "\n")
  
  # Remove rows with any NA values
  complete_rows <- complete.cases(log_matrix)
  log_matrix <- log_matrix[complete_rows, , drop = FALSE]
  
  # Print dimensions after NA removal
  cat("Matrix dimensions after removing NA rows:", dim(log_matrix), "\n")
  
  # Remove any rows with infinite values
  valid_rows <- apply(log_matrix, 1, function(x) !any(is.infinite(x)))
  log_matrix <- log_matrix[valid_rows, , drop = FALSE]
  
  # Print dimensions after removing infinite values
  cat("Final matrix dimensions:", dim(log_matrix), "\n")
  
  if (nrow(log_matrix) == 0) {
    stop("No valid data remaining after filtering NA/Inf values")
  }
  
  # Calculate dimensions based on number of rows and length of row names
  max_row_name_length <- max(nchar(rownames(log_matrix)))
  
  # Calculate height and width based on data dimensions
  height_per_row <- 15  # pixels per row
  width_per_char <- 7   # pixels per character
  
  # Calculate total dimensions in pixels
  height_px <- max(800, nrow(log_matrix) * height_per_row)
  width_px <- max(800, ncol(log_matrix) * 50 + max_row_name_length * width_per_char)
  
  # Convert to inches for PNG (assuming 72 DPI)
  height_in <- height_px / 72
  width_in <- width_px / 72
  
  # Create sample annotations
  annotation_col <- data.frame(
    Group = rep(NA, ncol(log_matrix)),
    row.names = colnames(log_matrix)
  )
  
  # Assign groups based on patterns
  for (i in seq_along(group_pattern)) {
    annotation_col$Group[grepl(paste0("^", group_pattern[i]), rownames(annotation_col))] <- group_names[i]
  }
  
  # Define colors for groups
  group_colors <- c("Healthy" = "#1f77b4", "MECFS" = "#d62728")
  ann_colors <- list(Group = group_colors)
  
  # Create color gradient using viridis
  color_gradient <- viridis(100)
  
  # Modify filename to use PNG instead of PDF if provided
  if (!is.null(filename)) {
    filename <- sub("\\.pdf$", ".png", filename)
  }
  
  # Create heatmap with error handling
  tryCatch({
    # Create base heatmap parameters
    heatmap_params <- list(
      mat = log_matrix,
      main = title,
      annotation_col = annotation_col,
      annotation_colors = ann_colors,
      scale = scale,
      cluster_rows = cluster_rows,
      cluster_cols = cluster_cols,
      show_rownames = show_rownames,
      show_colnames = show_colnames,
      color = color_gradient,
      na_col = "grey",
      fontsize_row = 10,
      fontsize_col = 10,
      cellwidth = 15,
      cellheight = 15,
      border_color = NA,
      treeheight_row = 50,
      treeheight_col = 50
    )
    
    # If filename is provided, save as PNG
    if (!is.null(filename)) {
      png(filename, 
          width = width_px, 
          height = height_px, 
          res = 72)
      heatmap <- do.call(pheatmap, heatmap_params)
      dev.off()
      
      # Print dimensions used
      cat(sprintf("Saved heatmap to %s\nDimensions: %d x %d pixels (%0.1f x %0.1f inches)\n",
                 filename, width_px, height_px, width_in, height_in))
    }
    
    # Create and return heatmap object for display in R
    heatmap <- do.call(pheatmap, heatmap_params)
    return(heatmap)
    
  }, error = function(e) {
    warning("Error creating heatmap: ", e$message)
    return(NULL)
  })
}


# Create and save heatmaps with clustering enabled
# Protein glycan composition heatmap
create_glyco_heatmap(
  protein_gly_comp_log,
  title = "Protein Glycan Composition",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  filename = "figures/protein_glycan_comp_heatmap.png"
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


#' Create Statistical Test and Boxplot for Group Comparisons
#' 
#' @param data A numeric vector of measurements
#' @param groups A vector of group labels (e.g., "HC" or "M")
#' @param feature_name Name of the feature being tested (for plot title)
#' @param output_dir Directory to save plots (default: "figures/boxplots")
#' @param min_samples Minimum number of non-NA samples per group (default: 3)
#' @param alpha Significance level for statistical tests (default: 0.05)
#'
#' @return List containing test results and plot object
create_stat_boxplot <- function(data, 
                              groups, 
                              feature_name,
                              output_dir = "figures/boxplots",
                              min_samples = 3,
                              alpha = 0.05) {
  
  require(ggplot2)
  require(rstatix)
  require(ggpubr)
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Create data frame
  df <- data.frame(
    Value = data,
    Group = groups
  )
  
  # Remove NA values
  df <- df[complete.cases(df), ]
  
  # Count samples per group
  sample_counts <- table(df$Group)
  
  # Initialize results list
  results <- list(
    feature = feature_name,
    sufficient_samples = FALSE,
    test_performed = FALSE,
    test_type = NA,
    p_value = NA,
    significant = FALSE,
    sample_counts = sample_counts,
    plot = NULL
  )
  
  # Check if we have enough samples
  if (all(sample_counts >= min_samples)) {
    results$sufficient_samples <- TRUE
    
    # Perform Shapiro-Wilk test for normality
    shapiro_test <- tapply(df$Value, df$Group, shapiro.test)
    all_normal <- all(sapply(shapiro_test, function(x) x$p.value > 0.05))
    
    # Perform statistical test
    if (all_normal) {
      # Use t-test if data is normal
      test_result <- t.test(Value ~ Group, data = df)
      results$test_type <- "t-test"
    } else {
      # Use Wilcoxon test if data is not normal
      test_result <- wilcox.test(Value ~ Group, data = df)
      results$test_type <- "Wilcoxon"
    }
    
    results$test_performed <- TRUE
    results$p_value <- test_result$p.value
    results$significant <- test_result$p.value < alpha
  }
  
  # Create boxplot
  p <- ggplot(df, aes(x = Group, y = Value, fill = Group)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    scale_fill_manual(values = c("HC" = "#1f77b4", "M" = "#d62728")) +
    theme_bw() +
    labs(
      title = feature_name,
      subtitle = if(results$test_performed) {
        sprintf("%s: p = %.3g%s", 
                results$test_type,
                results$p_value,
                if(results$significant) " *" else "")
      } else {
        "Insufficient samples for statistical testing"
      },
      caption = sprintf("Samples: %s", 
                       paste(names(sample_counts), "=", sample_counts, 
                             collapse = ", "))
    ) +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(size = 8),
      legend.position = "none"
    )
  
  results$plot <- p
  
  # Save plot if statistical test was performed or if explicitly requested
  if (results$test_performed || !is.null(output_dir)) {
    filename <- file.path(output_dir, 
                         paste0(make.names(feature_name), "_boxplot.png"))
    ggsave(filename, p, width = 6, height = 8, dpi = 300)
  }
  
  return(results)
}

# Example usage:
# For a single feature
test_single <- function(feature_data, feature_name) {
  results <- create_stat_boxplot(
    data = feature_data,
    groups = your_groups_vector,  # Replace with your actual groups vector
    feature_name = feature_name
  )
  
  # Print results
  cat("\nResults for", feature_name, ":\n")
  cat("Samples per group:", 
      paste(names(results$sample_counts), "=", results$sample_counts, collapse = ", "), 
      "\n")
  if (results$test_performed) {
    cat("Test type:", results$test_type, "\n")
    cat("P-value:", format(results$p_value, digits = 3), "\n")
    cat("Significant:", results$significant, "\n")
  } else {
    cat("No statistical test performed - insufficient samples\n")
  }
}

# For multiple features
test_multiple <- function(data_matrix, group_vector) {
  results_list <- list()
  
  for (feature in rownames(data_matrix)) {
    results_list[[feature]] <- create_stat_boxplot(
      data = as.numeric(data_matrix[feature, ]),
      groups = group_vector,
      feature_name = feature
    )
  }
  
  # Create summary of significant results
  significant_features <- names(results_list)[
    sapply(results_list, function(x) x$significant)
  ]
  
  cat("\nSignificant features (p <", alpha, "):\n")
  for (feat in significant_features) {
    cat(sprintf("%s: %s p = %.3g\n", 
                feat, 
                results_list[[feat]]$test_type,
                results_list[[feat]]$p_value))
  }
  
  return(results_list)
}

# Now use the function with named parameters:
# For protein glycan composition
glycan_results <- test_multiple(
  data_matrix = protein_gly_comp_log,
  group_vector = sub("_.*", "", colnames(protein_gly_comp_log))  # Assuming column names start with group labels
)

# For sialic acid counts
sia_results <- test_multiple(
  data_matrix = protein_gly_sia_count_log,
  group_vector = sub("_.*", "", colnames(protein_gly_sia_count_log))
)

# For fucosylation
fuc_results <- test_multiple(
  data_matrix = protein_gly_fuc_log,
  group_vector = sub("_.*", "", colnames(protein_gly_fuc_log))
)


# Create and save histogram
png("figures/gPSM_histo.png", width = 800, height = 600)
hist(log(glycoPSMs$intensity),
     breaks = 100,
     main = "Histogram of Log-Transformed GlycoPSM Intensities",
     xlab = "Log Intensity",
     ylab = "Frequency",
     col = "steelblue")
dev.off()

# Create boxplot of glycoPSM intensities
png("figures/glycoPSM_intensities_boxplot.png", width = 1000, height = 800)
boxplot(log(glycoPSMs$intensity) ~ glycoPSMs$sample,
        main = "GlycoPSM Intensities Across Samples",
        xlab = "Sample",
        ylab = "Log Intensity",
        col = "steelblue",
        las = 2,  # Rotate x-axis labels vertically
        cex.axis = 0.8,  # Reduce axis label size
        boxwex = 0.8,    # Width of boxes
        outline = TRUE,   # Show outlier points
        notch = FALSE)   # Standard boxes without notches
dev.off()

# Calculate descriptive statistics for glycoPSM intensities
glyco_stats <- data.frame(
  Raw_Intensity = c(
    Mean = mean(glycoPSMs$intensity, na.rm = TRUE),
    Median = median(glycoPSMs$intensity, na.rm = TRUE),
    SD = sd(glycoPSMs$intensity, na.rm = TRUE),
    Min = min(glycoPSMs$intensity, na.rm = TRUE),
    Max = max(glycoPSMs$intensity, na.rm = TRUE),
    Q1 = quantile(glycoPSMs$intensity, 0.25, na.rm = TRUE),
    Q3 = quantile(glycoPSMs$intensity, 0.75, na.rm = TRUE),
    N = sum(!is.na(glycoPSMs$intensity)),
    N_missing = sum(is.na(glycoPSMs$intensity))
  ),
  Log_Intensity = c(
    Mean = mean(log(glycoPSMs$intensity), na.rm = TRUE),
    Median = median(log(glycoPSMs$intensity), na.rm = TRUE),
    SD = sd(log(glycoPSMs$intensity), na.rm = TRUE),
    Min = min(log(glycoPSMs$intensity), na.rm = TRUE),
    Max = max(log(glycoPSMs$intensity), na.rm = TRUE),
    Q1 = quantile(log(glycoPSMs$intensity), 0.25, na.rm = TRUE),
    Q3 = quantile(log(glycoPSMs$intensity), 0.75, na.rm = TRUE),
    N = sum(!is.na(log(glycoPSMs$intensity))),
    N_missing = sum(is.na(log(glycoPSMs$intensity)))
  )
)

# Print the statistics
print("Descriptive Statistics for GlycoPSM Intensities:")
print(round(glyco_stats, 3))

# Calculate statistics by sample using tapply for each metric
sample_means <- tapply(glycoPSMs$intensity, glycoPSMs$sample, mean, na.rm = TRUE)
sample_medians <- tapply(glycoPSMs$intensity, glycoPSMs$sample, median, na.rm = TRUE)
sample_sds <- tapply(glycoPSMs$intensity, glycoPSMs$sample, sd, na.rm = TRUE)
sample_n <- tapply(glycoPSMs$intensity, glycoPSMs$sample, length)
sample_missing <- tapply(glycoPSMs$intensity, glycoPSMs$sample, function(x) sum(is.na(x)))

# Combine into a data frame
sample_stats_df <- data.frame(
  Sample = names(sample_means),
  Mean = as.numeric(sample_means),
  Median = as.numeric(sample_medians),
  SD = as.numeric(sample_sds),
  N = as.numeric(sample_n),
  Missing = as.numeric(sample_missing)
)


# Save statistics to file
write.csv(glyco_stats, "output_data/glycoPSM_overall_stats.csv")
write.csv(sample_stats_df, "output_data/glycoPSM_sample_stats.csv", row.names = FALSE)

# Calculate descriptive statistics for glycoPSM pep_2d
pep_2d_stats <- data.frame(
  Statistics = c(
    Mean = mean(glycoPSMs$pep_2d, na.rm = TRUE),
    Median = median(glycoPSMs$pep_2d, na.rm = TRUE),
    SD = sd(glycoPSMs$pep_2d, na.rm = TRUE),
    Min = min(glycoPSMs$pep_2d, na.rm = TRUE),
    Max = max(glycoPSMs$pep_2d, na.rm = TRUE),
    Q1 = quantile(glycoPSMs$pep_2d, 0.25, na.rm = TRUE),
    Q3 = quantile(glycoPSMs$pep_2d, 0.75, na.rm = TRUE),
    N = sum(!is.na(glycoPSMs$pep_2d)),
    N_missing = sum(is.na(glycoPSMs$pep_2d))
  )
)

# Calculate statistics by sample
sample_pep_2d_means <- tapply(glycoPSMs$pep_2d, glycoPSMs$sample, mean, na.rm = TRUE)
sample_pep_2d_medians <- tapply(glycoPSMs$pep_2d, glycoPSMs$sample, median, na.rm = TRUE)
sample_pep_2d_sds <- tapply(glycoPSMs$pep_2d, glycoPSMs$sample, sd, na.rm = TRUE)
sample_pep_2d_n <- tapply(glycoPSMs$pep_2d, glycoPSMs$sample, length)
sample_pep_2d_missing <- tapply(glycoPSMs$pep_2d, glycoPSMs$sample, function(x) sum(is.na(x)))

# Combine into a data frame
sample_pep_2d_stats <- data.frame(
  Sample = names(sample_pep_2d_means),
  Mean = as.numeric(sample_pep_2d_means),
  Median = as.numeric(sample_pep_2d_medians),
  SD = as.numeric(sample_pep_2d_sds),
  N = as.numeric(sample_pep_2d_n),
  Missing = as.numeric(sample_pep_2d_missing)
)

# Print statistics
print("Overall Descriptive Statistics for GlycoPSM pep_2d:")
print(round(pep_2d_stats, 3))
print("\nDescriptive Statistics by Sample:")
print(round(sample_pep_2d_stats, 3))

# Save statistics to files
write.csv(pep_2d_stats, "output_data/glycoPSM_pep_2d_overall_stats.csv")
write.csv(sample_pep_2d_stats, "output_data/glycoPSM_pep_2d_sample_stats.csv", row.names = FALSE)

# Create boxplot
png("figures/glycoPSM_pep_2d_boxplot.png", width = 1000, height = 800)
par(mar = c(10, 5, 4, 2))  # Increase bottom margin for labels
boxplot(log(pep_2d) ~ sample, 
        data = glycoPSMs,
        main = "Log Peptide 2D Score Across Samples",
        xlab = "",
        ylab = "Log Peptide 2D Score",
        col = "steelblue",
        las = 2,       # Rotate x-axis labels vertically
        cex.axis = 0.8,# Reduce axis label size
        boxwex = 0.8,  # Width of boxes
        outline = TRUE,# Show outlier points
        notch = FALSE  # Standard boxes without notches
)
# Add horizontal grid lines
grid(nx = 0, ny = NULL, col = "gray", lty = "dotted")
dev.off()

# Create histogram of pep_2d distribution
png("figures/glycoPSM_pep_2d_histogram.png", width = 800, height = 600)
hist(log(glycoPSMs$pep_2d),
     breaks = 50,
     main = "Distribution of log Peptide 2D Scores",
     xlab = "Peptide 2D Score",
     ylab = "Frequency",
     col = "steelblue",
     border = "white")
# Add density curve
lines(density(log(glycoPSMs$pep_2d), na.rm = TRUE), col = "red", lwd = 2)
dev.off()



## Empirical Bayes Differential Expression Analysis
# Load required packages
if (!requireNamespace("limma", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("limma")
}
library(limma)

# Create design matrix
# Extract group information (HC vs M) from column names
groups <- factor(ifelse(grepl("^HC", colnames(protein_gly_comp_log)), "HC", "M"))
print("Group levels:")
print(levels(groups))

# Create design matrix
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)
print("Design matrix column names:")
print(colnames(design))

# Create contrast matrix for HC vs M comparison
contrast.matrix <- makeContrasts(
  M_vs_HC = M - HC,
  levels = design
)

# Fit linear model
fit <- lmFit(protein_gly_comp_log, design)

# Fit contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get results
results <- topTable(fit2, 
                   coef = "M_vs_HC", 
                   number = Inf,  # Return all results
                   adjust.method = "BH")  # Benjamini-Hochberg correction

# Add more information to results
results$Feature <- rownames(results)
results$Significant <- !is.na(results$adj.P.Val) & results$adj.P.Val < 0.05
results$Direction <- ifelse(!is.na(results$logFC) & results$logFC > 0, "Up in MECFS", "Down in MECFS")

# Print summary
cat("\nDifferential Analysis Summary:\n")
cat("Total features tested:", nrow(results), "\n")
cat("Features with NA p-values:", sum(is.na(results$adj.P.Val)), "\n")
cat("Significant features (FDR < 0.05):", sum(results$Significant, na.rm = TRUE), "\n")
cat("  Up in MECFS:", sum(results$Significant & results$logFC > 0, na.rm = TRUE), "\n")
cat("  Down in MECFS:", sum(results$Significant & results$logFC < 0, na.rm = TRUE), "\n")

# Save results
write.csv(results, "output_data/protein_glycan_composition_limma_results.csv")

# Create volcano plot
png("figures/protein_glycan_composition_volcano.png", width = 1000, height = 800)
# Remove NA values for plotting
plot_data <- results[!is.na(results$P.Value) & !is.na(results$logFC), ]
plot(plot_data$logFC, -log10(plot_data$P.Value),
     pch = 20,
     col = ifelse(plot_data$Significant, "red", "grey"),
     main = "Volcano Plot: MECFS vs HC Glycan Composition",
     xlab = "Log2 Fold Change",
     ylab = "-log10(P-value)")

# Add threshold lines
abline(h = -log10(0.05), col = "blue", lty = 2)
abline(v = c(-1, 1), col = "blue", lty = 2)

# Add labels for significant features only if there are any
significant_features <- which(plot_data$Significant)
if(length(significant_features) > 0) {
    text(plot_data$logFC[significant_features],
         -log10(plot_data$P.Value)[significant_features],
         labels = rownames(plot_data)[significant_features],
         cex = 0.7, pos = 3)
} else {
    # Add a note to the plot if no significant features
    text(0, max(-log10(plot_data$P.Value)) * 0.9, 
         "No significant features found",
         cex = 1.2)
}
dev.off()

# Only create heatmap if there are significant features
sig_count <- sum(results$Significant, na.rm = TRUE)
if(sig_count > 0) {
    significant_matrix <- protein_gly_comp_log[results$Significant & !is.na(results$Significant), ]
    
    png("figures/significant_glycan_composition_heatmap.png", 
        width = 1200, height = 800)
    
    # Scale the rows
    scaled_matrix <- t(scale(t(significant_matrix)))
    
    # Create heatmap
    heatmap(scaled_matrix,
            main = "Significant Glycan Composition Changes",
            col = colorRampPalette(c("blue", "white", "red"))(100),
            scale = "none",  # already scaled
            margins = c(8,8),
            cexRow = 0.8,
            cexCol = 0.8)
    
    dev.off()
} else {
    cat("\nNo significant features found - skipping heatmap generation\n")
}

# Create boxplots only if there are significant features
top_n <- min(10, sig_count)
if(top_n > 0) {
    top_features <- rownames(results)[order(results$adj.P.Val)][1:top_n]
    
    png("figures/top_significant_glycans_boxplots.png", 
        width = 1200, height = 200 * top_n)
    par(mfrow = c(top_n, 1), mar = c(4,4,2,1))
    
    for(feature in top_features) {
        boxplot(protein_gly_comp_log[feature,] ~ groups,
                main = feature,
                xlab = "Group",
                ylab = "Log Expression",
                col = c("blue", "red"))
        stripchart(protein_gly_comp_log[feature,] ~ groups,
                   vertical = TRUE,
                   method = "jitter",
                   add = TRUE,
                   pch = 20)
    }
    
    dev.off()
} else {
    cat("\nNo significant features found - skipping boxplot generation\n")
}

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()), "output_data/limma_analysis_session_info.txt")



########################################################
cor_matrix <- cor(protein_gly_comp_log, use = "pairwise.complete.obs", method = "spearman")
# Create a heatmap of the correlation matrix
png("figures/glycan_correlation_heatmap.png", width = 1000, height = 800)
heatmap(cor_matrix, 
        main = "Protein N-Glycan Composition Correlation of Samples Heatmap from log-transformed relative intensities",
        col = colorRampPalette(c("blue", "white", "red"))(100),
        scale = "none",
        margins = c(8,8),
        cexRow = 0.8,
        cexCol = 0.8)
dev.off()

#---------------------------------------------------
protein_gly_comp_log_transposed <- t(protein_gly_comp_log)
cor_matrix_transposed <- cor(protein_gly_comp_log_transposed, use = "pairwise.complete.obs", method = "spearman")

# Remove any rows/columns that are all NA
cor_matrix_transposed <- cor_matrix_transposed[!apply(is.na(cor_matrix_transposed), 1, all), !apply(is.na(cor_matrix_transposed), 2, all)]

# Replace any remaining NA/Inf values with 0 for visualization
cor_matrix_transposed[is.na(cor_matrix_transposed)] <- 0
cor_matrix_transposed[is.infinite(cor_matrix_transposed)] <- 0

# Create a heatmap of the correlation matrix
png("figures/glycan_correlation_heatmap_transposed.png", width = 3000, height = 1200)
  heatmap(cor_matrix_transposed, 
        main = "Correlation Between Protein N-Glycan Composition Heatmap of log-transformed relative intensities",
        col = colorRampPalette(c("blue", "white", "red"))(100),
        scale = "none",
        margins = c(8,8),
        na.rm = TRUE)
dev.off()   

# Print diagnostic information
cat("\nCorrelation Matrix Summary:\n")
cat("Dimensions:", dim(cor_matrix), "\n")
cat("NA values:", sum(is.na(cor_matrix)), "\n")
cat("Inf values:", sum(is.infinite(cor_matrix)), "\n")

# Save the correlation matrix to CSV
write.csv(cor_matrix_transposed, "output_data/glycan_correlation_matrix.csv")


## Logistic Regression

# Create disease status vector (0 for Healthy, 1 for MECFS)
disease_status <- ifelse(grepl("^HC", colnames(protein_gly_comp_log)), 0, 1)

# Print initial diagnostics
cat("Initial disease status distribution:", table(disease_status), "\n")

# First transpose the original matrix
protein_gly_comp_df <- as.data.frame(t(protein_gly_comp_log))

# Add disease status to the data frame
protein_gly_comp_df$disease_status <- as.factor(disease_status)

# Print data structure before cleaning
cat("\nData structure before cleaning:\n")
print(dim(protein_gly_comp_df))
print(table(protein_gly_comp_df$disease_status))

# Remove columns with all NA values
na_cols <- colSums(is.na(protein_gly_comp_df)) == nrow(protein_gly_comp_df)
protein_gly_comp_df <- protein_gly_comp_df[, !na_cols]

# Remove columns with any NA values
na_cols <- colSums(is.na(protein_gly_comp_df)) > 0
protein_gly_comp_df <- protein_gly_comp_df[, !na_cols]

# Print data structure after cleaning
cat("\nData structure after cleaning:\n")
print(dim(protein_gly_comp_df))
print(table(protein_gly_comp_df$disease_status))

# Check if we have enough data
if(nrow(protein_gly_comp_df) < 2) {
    stop("Not enough observations after cleaning")
}


# Try fitting the model with error handling
tryCatch({
    # Fit logistic regression
    glm.fit <- glm(disease_status ~ ., 
                   data = protein_gly_comp_df, 
                   family = binomial())
    
    # Print summary of the model
    print(summary(glm.fit))
  
}, error = function(e) {
    cat("Error fitting logistic regression model:", conditionMessage(e), "\n")
    return(NULL)
})

# Linear Discriminant Analysis

library(MASS)
# Convert protein_gly_comp_log_transposed to dataframe
protein_gly_comp_df <- as.data.frame(protein_gly_comp_log_transposed)

# Remove any existing disease_status column
protein_gly_comp_df$disease_status <- NULL  # Remove if it exists

# Add disease status
protein_gly_comp_df$disease_status <- factor(ifelse(grepl("^HC", rownames(protein_gly_comp_df)), 
                                                   "Healthy", "MECFS"))

# Print initial diagnostic information
cat("Initial data dimensions:", dim(protein_gly_comp_df), "\n")
cat("Initial group sizes:", table(protein_gly_comp_df$disease_status), "\n")

# Remove columns with all NA values
na_cols <- colSums(is.na(protein_gly_comp_df)) == nrow(protein_gly_comp_df)
protein_gly_comp_df <- protein_gly_comp_df[, !na_cols]

# Remove columns with any NA values
na_cols <- colSums(is.na(protein_gly_comp_df)) > 0
protein_gly_comp_df <- protein_gly_comp_df[, !na_cols]

# Function to check if variable is constant within groups with detailed reporting
check_constant_within_groups <- function(data, response) {
    constant_vars <- character(0)
    for(col in colnames(data)[colnames(data) != response]) {
        # Check each group
        for(group in unique(data[[response]])) {
            group_data <- data[data[[response]] == group, col]
            if(length(unique(group_data)) <= 1) {
                constant_vars <- c(constant_vars, col)
                cat("Variable", col, "is constant within group", group, "\n")
                break
            }
        }
    }
    return(constant_vars)
}

# Check and remove constant variables
constant_vars <- check_constant_within_groups(protein_gly_comp_df, "disease_status")
if(length(constant_vars) > 0) {
    cat("\nRemoving", length(constant_vars), "variables that are constant within groups\n")
    protein_gly_comp_df <- protein_gly_comp_df[, !colnames(protein_gly_comp_df) %in% constant_vars]
}

# Get predictors
predictors <- protein_gly_comp_df[, !colnames(protein_gly_comp_df) %in% "disease_status"]

# Print final diagnostic information
cat("\nAfter cleaning:\n")
cat("Final data dimensions:", dim(protein_gly_comp_df), "\n")
cat("Final group sizes:", table(protein_gly_comp_df$disease_status), "\n")
cat("Number of predictors:", ncol(predictors), "\n")

# Check if we have enough data
if(nrow(protein_gly_comp_df) > 0 && 
   all(table(protein_gly_comp_df$disease_status) > 0) && 
   ncol(predictors) > 0) {
    
    # Try fitting LDA with error handling
    tryCatch({
        # Fit LDA
        lda.fit <- lda(disease_status ~ ., data = protein_gly_comp_df)
        
        # Make predictions
        lda.pred <- predict(lda.fit, protein_gly_comp_df)
        
        # Create confusion matrix
        conf_matrix <- table(Actual = protein_gly_comp_df$disease_status, 
                           Predicted = lda.pred$class)
        
        # Print results
        cat("\nLDA Results:\n")
        print(conf_matrix)
        cat("\nAccuracy:", sum(diag(conf_matrix))/sum(conf_matrix), "\n")
        
        # Print variable importance
        var_imp <- abs(lda.fit$scaling[,1])
        var_imp <- sort(var_imp, decreasing = TRUE)
        cat("\nTop 10 most important variables:\n")
        print(head(var_imp, 10))
        
        # Save results
        write.csv(var_imp, "output_data/lda_variable_importance.csv")
        
    }, error = function(e) {
        cat("\nError in LDA:", conditionMessage(e), "\n")
        cat("Problematic variables might still be present\n")
        cat("Consider checking variance of remaining variables\n")
    })
    
} else {
    cat("\nInsufficient data for LDA after cleaning.\n")
    cat("Need non-zero observations in each group and at least one predictor.\n")
}

plot(lda.fit)

# Assuming lda.fit is your LDA model from the MASS package
library(MASS)
library(ggplot2)

# Create the plot
p <- plot(lda.fit)

# To save as PNG using base R
png("lda_plot.png", width=800, height=600)
plot(lda.fit)
dev.off()

# If there's only one LD dimension, modify your plotting approach
if(is.null(dim(lda_data$x)) || ncol(lda_data$x) == 1) {
  # Create a one-dimensional plot (points along a line)
  df <- data.frame(x = lda_data$x, class = lda_data$class)
  
  p <- ggplot(df, aes(x=x, y=0, color=class)) +
    geom_point(size=3) +
    labs(title="LDA Plot", x="LD1", y="") +
    theme_minimal() +
    theme(axis.text.y=element_blank(),  # Hide y-axis text
          axis.ticks.y=element_blank()) # Hide y-axis ticks
  
} else {
  # Original 2D plot
  df <- data.frame(x = lda_data$x[,1], y = lda_data$x[,2], 
                  class = lda_data$class)
  
  p <- ggplot(df, aes(x=x, y=y, color=class)) +
    geom_point() +
    labs(title="LDA Plot", x="LD1", y="LD2") +
    theme_minimal()
}

# Density plot of the single discriminant by class
df <- data.frame(x = lda_data$x, class = lda_data$class)

p <- ggplot(df, aes(x=x, fill=class)) +
  geom_density(alpha=0.7) +
  labs(title="LDA Plot", x="LD1", y="Density") +
  theme_minimal()

# Save the plot
ggsave("lda_plot.png", plot=p, width=8, height=6, dpi=300)

# Save the plot
ggsave("lda_plot.png", plot=p, width=8, height=6, dpi=300)

# K-Nearest Neighbors

library(class)
knn.fit <- knn.cv(train = protein_gly_comp_df[, -ncol(protein_gly_comp_df)],  
               cl = protein_gly_comp_df$disease_status, 
               k = 3)
# To see the predicted classes
print(knn.fit)

# To create a confusion matrix and evaluate performance
table(knn.fit, protein_gly_comp_df$disease_status)

# To calculate accuracy
mean(knn.fit == protein_gly_comp_df$disease_status)

library(caret)
confusionMatrix(knn.fit, protein_gly_comp_df$disease_status)

##general linear model - logit regression for binary classification
library(glmnet)

# Ensure disease_status is binary (0/1)
protein_gly_comp_df$disease_status <- as.numeric(factor(protein_gly_comp_df$disease_status)) - 1

# Check for missing values
protein_gly_comp_df <- na.omit(protein_gly_comp_df)

# Fit lasso model
glmnet.fit <- cv.glmnet(x = as.matrix(protein_gly_comp_df[, -ncol(protein_gly_comp_df)]),
                        y = protein_gly_comp_df$disease_status,
                        family = "binomial",
                        alpha = 1)

# Make predictions and convert to class labels
glmnet.pred <- predict(glmnet.fit, 
                      newx = as.matrix(protein_gly_comp_df[, -ncol(protein_gly_comp_df)]), 
                      s = "lambda.min", 
                      type = "response")

glmnet.pred <- ifelse(glmnet.pred > 0.5, 1, 0)

# Print confusion matrix
conf_matrix <- table(Actual = protein_gly_comp_df$disease_status, Predicted = glmnet.pred)
print("Confusion Matrix:")
print(conf_matrix)

# Calculate accuracy
accuracy <- sum(diag(conf_matrix))/sum(conf_matrix)
cat("\nAccuracy:", round(accuracy, 3), "\n")

# Get coefficients at minimum lambda
coef_matrix <- as.matrix(coef(glmnet.fit, s = "lambda.min"))
nonzero_coef <- coef_matrix[coef_matrix != 0, , drop = FALSE]
top_coef <- head(nonzero_coef[order(abs(nonzero_coef), decreasing = TRUE), ], 5)

# Create plot
png("figures/logit_glm_protglycomp.png", 
    width = 1000, 
    height = 800, 
    res = 300)

# Set up plotting parameters
par(mar = c(5, 5, 4, 8))

# Create the plot with enhanced features
plot(glmnet.fit,
     xvar = "lambda",
     label = FALSE,  # Don't show variable labels to avoid cluttering
     main = "Logistic Regression Coefficients vs Log(Lambda)",
     sub = "Protein Glycan Composition Analysis",
     xlab = "Log(Lambda)",
     ylab = "Standardized Coefficients")

# Add grid
grid(lty = "dotted", col = "gray80")

# Add legend only if we have top coefficients
if(length(top_coef) > 0) {
    legend("topright", 
           legend = paste(names(top_coef), 
                         sprintf(": %.3f", top_coef)),
           title = "Top 5 Features",
           cex = 0.8,
           bg = "white",
           box.col = "gray")
}

# Add model performance information
text(x = par("usr")[1], 
     y = par("usr")[4],
     labels = sprintf("CV Error (min): %.3f\nAccuracy: %.3f", 
                     min(glmnet.fit$cvm),
                     accuracy),
     pos = 4,
     cex = 0.8)

dev.off()

# Save detailed results
sink("output_data/logit_glm_protglycomp_summary.txt")
cat("Logistic Regression Model Summary\n")
cat("================================\n\n")
cat("Model Information:\n")
cat("Number of features:", nrow(coef_matrix), "\n")
cat("Number of non-zero coefficients:", sum(coef_matrix != 0), "\n")
cat("Lambda sequence:", length(glmnet.fit$lambda), "values\n")
cat("Optimal lambda:", glmnet.fit$lambda[which.min(glmnet.fit$cvm)], "\n")
cat("Minimum CV error:", min(glmnet.fit$cvm), "\n")
cat("Model accuracy:", round(accuracy, 3), "\n\n")

cat("Confusion Matrix:\n")
print(conf_matrix)
cat("\n")

cat("Top Features by Absolute Coefficient Size:\n")
print(top_coef)
sink()

# Create a data frame of all coefficients
coef_df <- data.frame(
    Feature = rownames(coef_matrix),
    Coefficient = as.vector(coef_matrix),
    Abs_Coefficient = abs(as.vector(coef_matrix))
)
coef_df <- coef_df[order(-coef_df$Abs_Coefficient), ]

# Save coefficients to CSV
write.csv(coef_df, "output_data/logit_glm_protglycomp_coefficients.csv", row.names = FALSE)


# Principal Component Regression
library(pls)

pcr.fit <- pcr(disease_status ~ ., data = protein_gly_comp_df, scale = TRUE, validation = "CV")
summary(pcr.fit)

# Make predictions on the training data
pcr.pred <- predict(pcr.fit, protein_gly_comp_df, ncomp = 3)

## PLS-DA (Partial Least Squares Discriminant Analysis)
plsda.fit <- plsda(x = protein_gly_comp_df[, -ncol(protein_gly_comp_df)], 
                  y = protein_gly_comp_df$disease_status,
                  ncomp = 3, validation = "CV")
summary(plsda.fit)

# Make predictions
plsda.pred <- predict(plsda.fit, newdata = protein_gly_comp_df[, -ncol(protein_gly_comp_df)], 
                     ncomp = 3, type = "class")

# Evaluate
confusionMatrix(plsda.pred, protein_gly_comp_df$disease_status)



# Create the loadings plot
loadings_plot <- mdaplot(m$loadings, 
                         type = "p", 
                         show.labels = TRUE, 
                         show.lines = c(0, 0), 
                         cgroup = protein_gly_comp_df$disease_status)

# Close the device to save the file
dev.off()


# Open a graphics device (optional, if saving to file)
png("scores_plot.png", width = 800, height = 600)

# Create the scores plot
mdaplot(m$res$cal$scores, type = "p", 
        show.labels = TRUE, 
        show.lines = c(0, 0), 
        cgroup = protein_gly_comp_df$disease_status)

# Call plotConfidenceEllipse with only the model (ensure 'p' is defined)
plotConfidenceEllipse(m$res$cal)

# Close the device to save the file
dev.off()

# Open a graphics device (optional, if saving to file)
png("PCA_plots.png", width = 800, height = 600)
plot(m, show.labels = TRUE)
dev.off()


library(mdatools)
m <- pca(protein_gly_comp_df, 7, info = "Plasma Protein N-glycosolation HC vs MECFS PCA Model")
m = selectCompNum(m, 5)
print(m)
m$loadings[1:4, 1:4]

# Open a graphics device (optional, if saving to file)
png("loadings_plot.png", width = 800, height = 600)

# Random Forest
library(randomForest)
set.seed(20)
rf.fit <- randomForest(disease_status ~ ., data = protein_gly_comp_df)
print(rf.fit)



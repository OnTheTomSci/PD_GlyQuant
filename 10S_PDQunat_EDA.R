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
             write.csv(protein_gly_class, file = "output_data/protein_gly_class.csv")
             protein_gly_sia <- glyco_matrix(
               gpeps_dataframe = glycoPSMs,
               top_lev_group = "protein_accessions",
               glycofeature_group = "contains_NeuAc",
               value_col = "intensity",
               sample_col = "sample",
               group_col = "disease_status",
               file_prefix = "protein_gly_sia"
             )            
             write.csv(protein_gly_sia, file = "output_data/protein_gly_sia.csv")
             protein_gly_fuc <- glyco_matrix(
               gpeps_dataframe = glycoPSMs,
               top_lev_group = "protein_accessions",
               glycofeature_group = "contains_Fuc",
               value_col = "intensity",
               sample_col = "sample",
               group_col = "disease_status",
               file_prefix = "protein_gly_fuc"
             )       
             write.csv(protein_gly_fuc, file = "output_data/protein_gly_fuc.csv")
             protein_gly_sia_count <- glyco_matrix(
               gpeps_dataframe = glycoPSMs,
               top_lev_group = "protein_accessions",
               glycofeature_group = "sia_count",
               value_col = "intensity",
               sample_col = "sample",
               group_col = "disease_status",
               file_prefix = "protein_gly_sia_count"
             )            
             write.csv(protein_gly_sia_count, file = "output_data/protein_gly_sia_count.csv")
             protein_gly_comp <- glyco_matrix(
               gpeps_dataframe = glycoPSMs,
               top_lev_group = "protein_accessions",
               glycofeature_group = "glycan_composition",
               value_col = "intensity",
               sample_col = "sample",
               group_col = "disease_status",
               file_prefix = "protein_gly_comp"
             )                         
             write.csv(protein_gly_comp, file = "output_data/protein_gly_comp.csv")
             
             glycosite_gly_class <- glyco_matrix(
               gpeps_dataframe = glycoPSMs,
               top_lev_group = "gsite_ID",
               glycofeature_group = "glycan_class",
               value_col = "intensity",
               sample_col = "sample",
               group_col = "disease_status",
               file_prefix = "glycosite_gly_class"
             )            
             write.csv(glycosite_gly_class, file = "output_data/glycosite_gly_class.csv")
             glycosite_gly_sia <- glyco_matrix(
               gpeps_dataframe = glycoPSMs,
               top_lev_group = "gsite_ID",
               glycofeature_group = "contains_NeuAc",
               value_col = "intensity",
               sample_col = "sample",
               group_col = "disease_status",
               file_prefix = "glycosite_gly_sia"
             )            
             write.csv(glycosite_gly_sia, file = "output_data/glycosite_gly_sia.csv")
             glycosite_gly_fuc <- glyco_matrix(
               gpeps_dataframe = glycoPSMs,
               top_lev_group = "gsite_ID",
               glycofeature_group = "contains_Fuc",
               value_col = "intensity",
               sample_col = "sample",
               group_col = "disease_status",
               file_prefix = "glycosite_gly_fuc"
             )       
             write.csv(glycosite_gly_fuc, file = "output_data/glycosite_gly_fuc.csv")
             glycosite_gly_sia_count <- glyco_matrix(
               gpeps_dataframe = glycoPSMs,
               top_lev_group = "gsite_ID",
               glycofeature_group = "sia_count",
               value_col = "intensity",
               sample_col = "sample",
               group_col = "disease_status",
               file_prefix = "glycosite_gly_sia_count"
             )            
             write.csv(glycosite_gly_sia_count, file = "output_data/glycosite_gly_sia_count.csv")
             glycosite_gly_comp <- glyco_matrix(
               gpeps_dataframe = glycoPSMs,
               top_lev_group = "gsite_ID",
               glycofeature_group = "glycan_composition",
               value_col = "intensity",
               sample_col = "sample",
               group_col = "disease_status",
               file_prefix = "glycosite_gly_comp"
             ) 
             write.csv(glycosite_gly_class, file = "output_data/glycosite_gly_class.csv")           
             
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
  transformed_matrix <- switch(base,
    2 = log2(transformed_matrix),
    10 = log10(transformed_matrix),
    exp(1) = log(transformed_matrix),
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

protein_fuc_log <- log_transform_matrix(protein_fuc)
write.csv(protein_fuc_log, file = "output_data/protein_fuc_log.csv")

protein_sia_log <- log_transform_matrix(protein_sia) 
write.csv(protein_sia_log, file = "output_data/protein_sia_log.csv")

protein_sia_count_log <- log_transform_matrix(protein_sia_count)
write.csv(protein_sia_count_log, file = "output_data/protein_sia_count_log.csv")

# Log transform and write glycosite matrices to CSV
glycosite_gly_comp_log <- log_transform_matrix(glycosite_gly_comp)
write.csv(glycosite_gly_comp_log, file = "output_data/glycosite_gly_comp_log.csv")

glycosite_fuc_log <- log_transform_matrix(glycosite_fuc)
write.csv(glycosite_fuc_log, file = "output_data/glycosite_fuc_log.csv")

glycosite_sia_log <- log_transform_matrix(glycosite_sia)
write.csv(glycosite_sia_log, file = "output_data/glycosite_sia_log.csv")

glycosite_sia_count_log <- log_transform_matrix(glycosite_sia_count)
write.csv(glycosite_sia_count_log, file = "output_data/glycosite_sia_count_log.csv")

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
# Peptide Groups Preprocessing Functions
# This module contains functions for loading, cleaning, and preprocessing peptide groups data

library(tidyverse)
library(here)
library(janitor)
library(stringr)
library(readr)
library(dplyr)

#' Load and preprocess peptide groups data
#' 
#' @param study_info_path Path to study information file
#' @param peptide_groups_path Path to peptide groups file
#' @param glycan_class_map_path Path to glycan class mapping file
#' @param fdr_threshold FDR threshold for filtering (default: 0.001)
#' @return List containing processed data and metadata
load_and_preprocess_peptidegroups <- function(study_info_path = "input_data/10S_MECFS_GPEPS_250125_StudyInformation.txt",
                                             peptide_groups_path = "input_data/10S_MECFS_GPEPS_250125_PeptideGroups.txt",
                                             glycan_class_map_path = "input_data/glycan_class_map.csv",
                                             fdr_threshold = 0.001) {
  
  # Load data files
  StudyInformation <- read_tsv(study_info_path)
  PeptideGroups <- read_tsv(peptide_groups_path)
  glycan_class_map <- read_csv(glycan_class_map_path)
  
  # Clean and process peptide groups
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
  
  # Filter for glycosylated peptides and apply FDR threshold
  glyco_peptide_groups <- PeptideGroups %>% 
    filter(!is.na(glycan_composition)) %>%
    filter(pep_2d_by_search_engine_a2_pmi_byonic < fdr_threshold)
  
  # Clean study information
  StudyInformation <- StudyInformation %>%
    clean_names(case = "snake")
  
  # Match glycan classes
  glyco_peptide_groups <- match_glycan_class(glyco_peptide_groups, glycan_class_map)
  
  # Add glycosite ID
  glyco_peptide_groups <- glyco_peptide_groups %>%
    mutate(gsite_ID = paste(protein_accessions, protein_glycosite, sep = "_"))
  
  # Convert to long format
  glyco_peptide_groups_long <- glyco_peptide_groups %>%
    select(-any_of(starts_with("abundances_")), -any_of(starts_with("abundance_ratio")), -any_of(starts_with("found_in_"))) %>%
    pivot_longer(
      cols = starts_with("abundance_f"),
      names_to = "sample",
      values_to = "abundance"
    ) %>%
    filter(!is.na(abundance)) %>%
    mutate(sample = str_remove(sample, "^([^_]*_){4}")) %>%
    mutate(sample = str_remove(sample, "_[^_]*$"))
  
  # Calculate basic statistics
  all_samples <- unique(glyco_peptide_groups_long$sample)
  total_samples <- length(all_samples)
  
  # Find common features across samples
  common_glycopeptides <- glyco_peptide_groups_long %>%
    group_by(peptide_groups_peptide_group_id) %>%
    summarise(
      num_samples = n_distinct(sample),
      .groups = 'drop'
    ) %>%
    filter(num_samples == total_samples) %>%
    pull(peptide_groups_peptide_group_id)
  
  common_proteins <- glyco_peptide_groups_long %>%
    group_by(protein_accessions) %>%
    summarise(
      num_samples = n_distinct(sample),
      .groups = 'drop'
    ) %>%
    filter(num_samples == total_samples) %>%
    pull(protein_accessions)
  
  common_gsites <- glyco_peptide_groups_long %>%
    group_by(gsite_ID) %>%
    summarise(
      num_samples = n_distinct(sample),
      .groups = 'drop'
    ) %>%
    filter(num_samples == total_samples) %>%
    pull(gsite_ID)
  
  common_genes <- glyco_peptide_groups_long %>%
    group_by(gene_name) %>%
    summarise(
      num_samples = n_distinct(sample),
      .groups = 'drop'
    ) %>%
    filter(num_samples == total_samples) %>%
    pull(gene_name)
  
  # Print summary statistics
  cat("=== Data Loading and Preprocessing Summary ===\n")
  cat("Number of unique glycopeptides:", length(find_unique_values(glyco_peptide_groups_long, "peptide_groups_peptide_group_id")), "\n")
  cat("Number of glycopeptides common across all samples:", length(common_glycopeptides), "\n")
  cat("Number of unique glycosites:", length(find_unique_values(glyco_peptide_groups_long, "gsite_ID")), "\n")
  cat("Number of glycosites common across all samples:", length(common_gsites), "\n")
  cat("Number of unique glycoproteins:", length(find_unique_values(glyco_peptide_groups_long, "protein_accessions")), "\n")
  cat("Number of proteins common across all samples:", length(common_proteins), "\n")
  cat("Number of unique genes:", length(find_unique_values(glyco_peptide_groups_long, "gene_name")), "\n")
  cat("Number of genes common across all samples:", length(common_genes), "\n")
  cat("Total samples:", total_samples, "\n")
  cat("=============================================\n")
  
  return(list(
    glyco_peptide_groups = glyco_peptide_groups,
    glyco_peptide_groups_long = glyco_peptide_groups_long,
    StudyInformation = StudyInformation,
    glycan_class_map = glycan_class_map,
    sample_metadata = list(
      all_samples = all_samples,
      total_samples = total_samples,
      common_glycopeptides = common_glycopeptides,
      common_proteins = common_proteins,
      common_gsites = common_gsites,
      common_genes = common_genes
    )
  ))
}

#' Match glycan compositions to glycan classes
#' 
#' @param input_df Data frame with glycan_composition column
#' @param reference_df Reference data frame with glycans and glycan_class columns
#' @return Data frame with added glycan_class column
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

#' Find unique values in a data frame column
#' 
#' @param df Data frame
#' @param column_name Name of the column to extract unique values from
#' @return List of unique values
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

#' Calculate sample coverage for features
#' 
#' @param data Long format data with sample column
#' @param feature_col Column name to group by for coverage analysis
#' @return Data frame with coverage statistics
calculate_sample_coverage <- function(data, feature_col) {
  all_samples <- unique(data$sample)
  total_samples <- length(all_samples)
  
  coverage_stats <- data %>%
    group_by(!!sym(feature_col)) %>%
    summarise(
      num_samples = n_distinct(sample),
      coverage_percent = (num_samples / total_samples) * 100,
      .groups = 'drop'
    ) %>%
    arrange(desc(num_samples))
  
  return(coverage_stats)
}

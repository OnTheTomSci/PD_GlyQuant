# Debug script to check data structure for CLR analysis

library(tidyverse)

# Source functions
source("functions/peptidegroups_preprocessing.R")

# Load data
cat("Loading data...\n")
data <- load_and_preprocess_peptidegroups(
  study_info_path = "input_data/10S_MECFS_GPEPS_250125_StudyInformation.txt",
  peptide_groups_path = "input_data/10S_MECFS_GPEPS_250125_PeptideGroups.txt",
  glycan_class_map_path = "input_data/glycan_class_map.csv",
  fdr_threshold = 0.001
)

glyco_long <- data$glyco_peptide_groups_long

cat("\n=== DATA STRUCTURE CHECK ===\n")
cat("Total rows:", nrow(glyco_long), "\n")
cat("Column names:", paste(colnames(glyco_long), collapse=", "), "\n\n")

cat("=== CHECKING REQUIRED COLUMNS ===\n")
cat("Has 'glycan_composition':", "glycan_composition" %in% colnames(glyco_long), "\n")
cat("Has 'protein_accessions':", "protein_accessions" %in% colnames(glyco_long), "\n")
cat("Has 'abundance':", "abundance" %in% colnames(glyco_long), "\n")
cat("Has 'sample':", "sample" %in% colnames(glyco_long), "\n")
cat("Has 'group':", "group" %in% colnames(glyco_long), "\n\n")

cat("=== CHECKING FOR NAs ===\n")
cat("NA in glycan_composition:", sum(is.na(glyco_long$glycan_composition)), "\n")
cat("NA in protein_accessions:", sum(is.na(glyco_long$protein_accessions)), "\n")
cat("NA in abundance:", sum(is.na(glyco_long$abundance)), "\n")
cat("NA in group:", sum(is.na(glyco_long$group)), "\n\n")

cat("=== CHECKING ABUNDANCE VALUES ===\n")
cat("Rows with abundance > 0:", sum(glyco_long$abundance > 0, na.rm = TRUE), "\n")
cat("Rows with abundance = 0:", sum(glyco_long$abundance == 0, na.rm = TRUE), "\n")
cat("Min abundance:", min(glyco_long$abundance, na.rm = TRUE), "\n")
cat("Max abundance:", max(glyco_long$abundance, na.rm = TRUE), "\n\n")

cat("=== CHECKING AFTER FILTERING (like CLR function does) ===\n")
filtered_data <- glyco_long %>%
  filter(!is.na(glycan_composition), 
         !is.na(protein_accessions),
         abundance > 0)

cat("Rows after filtering:", nrow(filtered_data), "\n\n")

if (nrow(filtered_data) > 0) {
  cat("Sample of filtered data (first 10 rows):\n")
  print(filtered_data %>% 
          select(sample, group, protein_accessions, glycan_composition, abundance) %>%
          head(10))
  
  cat("\n=== UNIQUE VALUES ===\n")
  cat("Unique proteins:", n_distinct(filtered_data$protein_accessions), "\n")
  cat("Unique glycans:", n_distinct(filtered_data$glycan_composition), "\n")
  cat("Unique samples:", n_distinct(filtered_data$sample), "\n")
  cat("Unique groups:", paste(unique(filtered_data$group), collapse=", "), "\n")
} else {
  cat("\n!!! NO DATA REMAINS AFTER FILTERING !!!\n")
  cat("Checking each filter condition separately:\n\n")
  
  cat("Rows with non-NA glycan_composition:", 
      sum(!is.na(glyco_long$glycan_composition)), "\n")
  cat("Rows with non-NA protein_accessions:", 
      sum(!is.na(glyco_long$protein_accessions)), "\n")
  cat("Rows with abundance > 0:", 
      sum(glyco_long$abundance > 0, na.rm = TRUE), "\n")
  
  cat("\nSample of protein_accessions (first 20):\n")
  print(head(glyco_long$protein_accessions, 20))
  
  cat("\nSample of glycan_composition (first 20):\n")
  print(head(glyco_long$glycan_composition, 20))
}


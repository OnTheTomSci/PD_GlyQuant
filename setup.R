# Setup script for installing and loading required packages

# List of required packages
required_packages <- c(
  "tidyverse",      # For data manipulation and ggplot2
  "patchwork",      # For combining plots
  "viridis",        # For color scales
  "dplyr",          # For data manipulation
  "ggplot2",        # For plotting
  "readr",          # For reading/writing CSV files
  "tidyr",          # For data tidying
  "stringr",        # For string manipulation
  "forcats",        # For factor manipulation
  "stats"           # For statistical tests (usually included in base R)
)

# Function to install missing packages
install_missing_packages <- function() {
  new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    cat("Installing missing packages:", paste(new_packages, collapse = ", "), "\n")
    install.packages(new_packages)
  }
}

# Function to load all required packages
load_required_packages <- function() {
  # Install missing packages
  install_missing_packages()
  
  # Load all packages
  for(package in required_packages) {
    library(package, character.only = TRUE)
  }
  
  # Print package versions for reproducibility
  cat("\nLoaded package versions:\n")
  for(package in required_packages) {
    cat(sprintf("%s: %s\n", package, packageVersion(package)))
  }
}

# Check R version
check_r_version <- function(minimum_version = "4.0.0") {
  if(getRversion() < minimum_version) {
    stop("This code requires R version ", minimum_version, " or higher")
  }
  cat("R version:", getRversion(), "\n")
}

# Initialize environment
initialize_environment <- function() {
  # Check R version
  check_r_version()
  
  # Load required packages
  load_required_packages()
  
  # Create necessary directories if they don't exist
  dirs <- c("figures", "output_data", "R/functions")
  for(dir in dirs) {
    if(!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat("Created directory:", dir, "\n")
    }
  }
}

# Run initialization
initialize_environment() 
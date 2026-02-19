# Load DepMap RPPA proteomics data for AML samples
# DepMap provides hundreds of RPPA measurements across various cancer cell lines

library(tidyverse)
library(readr)
library(keras3)
# install_keras(tensorflow = "gpu")


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("depmap")
library("depmap")

# Function to download DepMap RPPA data
download_depmap_rppa <- function() {
  # DepMap RPPA data URL (this is the publicly available dataset)
  # The data includes protein expression measurements from reverse phase protein arrays
  # rppa_url <- "https://depmap.org/downloads/rppa_20Q4.csv"
  rppa_url <- "https://depmap.org/portal/data_page/?tab=allData&releasename=Harmonized%20Public%20Proteomics%2024Q4&filename=harmonized_RPPA_CCLE.csv"
  
  # Download the data
  cat("Downloading DepMap RPPA data...\n")
  rppa_data <- read_csv(rppa_url, show_col_types = FALSE)
  rppa_data <- read.csv(url(rppa_url))
  
  return(rppa_data)
}

# Function to filter for AML samples
filter_aml_samples <- function(rppa_data, sample_info) {
  # Filter RPPA data for AML samples only
  # Note: You'll need sample metadata to identify AML vs other cancer types
  
  aml_rppa <- rppa_data %>%
    filter(cell_line %in% sample_info$aml_cell_lines)
  
  return(aml_rppa)
}

# Function to format RPPA data
format_rppa_data <- function(rppa_data) {
  # Reformat from wide to long format or vice versa depending on analysis needs
  # Also standardize protein names and handle missing values
  
  rppa_formatted <- rppa_data %>%
    # Convert to long format (protein as rows, samples as columns)
    pivot_longer(cols = -cell_line, 
                 names_to = "protein", 
                 values_to = "expression_level") %>%
    # Remove missing values
    drop_na(expression_level) %>%
    # Standardize protein names
    mutate(protein = str_trim(protein)) %>%
    # Sort for consistency
    arrange(cell_line, protein)
  
  return(rppa_formatted)
}

# Main workflow
main <- function() {
  # Step 1: Download RPPA data
  rppa_raw <- download_depmap_rppa()
  
  cat("Loaded RPPA data with dimensions:", dim(rppa_raw), "\n")
  cat("Available proteins:", ncol(rppa_raw) - 1, "\n")
  
  # Step 2: Get AML-specific samples
  # First, download sample metadata to identify AML cell lines
  sample_info_url <- "https://depmap.org/downloads/sample_info.csv"
  sample_info <- read_csv(sample_info_url, show_col_types = FALSE)
  
  # Filter for AML (leukemia) samples
  aml_info <- sample_info %>%
    filter(str_detect(disease, "AML|Leukemia|acute myeloid"))
  
  cat("Found", nrow(aml_info), "AML samples in DepMap\n")
  
  # Step 3: Filter RPPA data for AML samples
  aml_rppa <- rppa_raw %>%
    filter(cell_line %in% aml_info$DepMap_ID)
  
  cat("Filtered RPPA data to", nrow(aml_rppa), "AML samples\n")
  
  # Step 4: Format the data
  aml_rppa_formatted <- format_rppa_data(aml_rppa)
  
  # Step 5: Save formatted data
  output_path <- "g:/My Drive/Projects/Kaggle_aiMATCH/Data/depmap_aml_rppa.csv"
  write_csv(aml_rppa_formatted, output_path)
  cat("Saved formatted AML RPPA data to:", output_path, "\n")
  
  return(aml_rppa_formatted)
}

# Execute main workflow
# aml_rppa_data <- main()


train_autoe <- build_and_train_autoencoder(mat = train_rppa, encoding_dim = 50, epochs = 100, batch_size = 16)



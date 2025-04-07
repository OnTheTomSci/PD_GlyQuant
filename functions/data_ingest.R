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
             
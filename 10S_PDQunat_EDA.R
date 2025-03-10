#' ---
#' title: Proteome Discover exploratory data Analysis of 10 sample ME/CFS Plasma glycoproteomics
#'   cohort
#' output: html_notebook
#'   
#' ---
#' #
#' The mass spec files of raw data are searched and analysised usising proteome discoverer that byonic as a peptide seach node. this analysis workflow uses a LFQ style quantifcation method using peak intnsities eaven when a sample does not a MS2 PSM corisponding to the peak in staed the peak is stilll quantifed on the basis that a PSM is found that in anouth sample that matches the peak's retention times.
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
Proteins <- read_tsv("10S_MECFS_GPEPS_250125_Proteins.txt")
ProteinGroups <- read_tsv("10S_MECFS_GPEPS_250125_ProteinGroups.txt")
PeptideGroups <- read_tsv("10S_MECFS_GPEPS_250125_PeptideGroups.txt")
PSMs <- read_tsv("10S_MECFS_GPEPS_250125_PSMs.txt")
ConsensusFeatures <- read_tsv("10S_MECFS_GPEPS_250125_ConsensusFeatures.txt")
InputFiles <- read_tsv("10S_MECFS_GPEPS_250125_InputFiles.txt")
PathwayProteinGroups <- read_tsv("10S_MECFS_GPEPS_250125_PathwayProteinGroups.txt")
StudyInformation <- read_tsv("10S_MECFS_GPEPS_250125_StudyInformation.txt")
ProteinAbundances <- read.csv("10S_MECFS_GPEPS_250125_Proteins_Abundances.csv")

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

png("glycan_barplot.png", width = 1600, height = 1200, res = 150)  # High-res image
par(mar = c(5, 15, 4, 2))
barplot(Glycans$count, 
        names.arg = Glycans$`Glycan Composition`, 
        horiz = TRUE,
        las = 1,
        col = "steelblue",
        main = "Glycan Composition Counts",
        xlab = "Count",
        cex.names = 0.9)
dev.off()  # Save the file


#' the following barchart of glycan PSM counts shows the plasma glycome is predominated by a few abundant glycans, thus skew the distrubution of abundances. this what we would expected and is simaly reflected in the released glycan analysis. 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::include_graphics("glycan_barplot.png")

#' 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Pivot the dataset from wide to long format
ProteinAbundances_long <- ProteinAbundances %>%
  pivot_longer(
    cols = starts_with("Abundance_"),  # Select columns dynamically
    names_to = "Sample",
    values_to = "Abundance"
  )
ProteinAbundances_long <- ProteinAbundances_long %>%
  mutate(
    Disease_Status = str_extract(Sample, "(?<=_)(Healthy|MECFS)$"),  # Extract "Healthy" or "MECFS"
    Sample = str_extract(Sample, "(?<=_)[^_]+(?=_)")  # Extract the middle part (Sample ID)
  )
ProteinAbundances_long <- ProteinAbundances_long %>%
  mutate(Abundance = replace_na(Abundance, 0))

#' 
#' Performed EDA on top 15 proteins while using glycopeptide data as as i wanted to advoid issues for my further statisical analysis that missing values could cause. 
#' 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Clean Protein_Names by removing everything including and after "OS="
ProteinAbundances_long <- ProteinAbundances_long %>%
  mutate(Protein_Names = sub("OS=.*", "", Protein_Names))  # Remove "OS=" and everything after it

# Select top 15 proteins by mean abundance
top_proteins <- ProteinAbundances_long %>%
  group_by(Protein_Names) %>%
  summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 15) %>%
  pull(Protein_Names)  # Extract protein names

# Filter dataset for only the top 15 proteins
ProteinAbundances_top15 <- ProteinAbundances_long %>%
  filter(Protein_Names %in% top_proteins)


# Create the heatmap
heatmap_plot <- ggplot(ProteinAbundances_top15, aes(x = Sample, y = Protein_Names, fill = Abundance)) +
  geom_tile(color = "white", height = 0.8) +
  scale_fill_viridis_c() +  # Corrected function
  theme_ipsum() +
  labs(
    title = "Top 15 Protein Abundance Heatmap",
    x = "Sample Name",
    y = "Protein",
    fill = "Abundance"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    panel.spacing = unit(1, "lines"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.05)))

ggsave("heatmap.png", plot = heatmap_plot, width = 10, height = 6, dpi = 300)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Perform t-tests and extract results
Top15_ProtAbun_t_test_results <- ProteinAbundances_top15 %>%
  group_by(Protein_Names) %>%
  summarise(
    t_test = list(
      tryCatch(
        t.test(
          Abundance ~ Disease_Status, 
          data = pick(everything()),  # Use `pick()` instead of `cur_data()`
          var.equal = FALSE  # Welch's t-test
        ), 
        error = function(e) NULL  # Handle errors gracefully
      )
    ),
    .groups = "drop"
  ) %>%
  mutate(
    t_test_summary = map(t_test, ~ if (!is.null(.x)) broom::tidy(.x) else NULL)
  ) %>%
  unnest(cols = t_test_summary, keep_empty = TRUE) %>% # Expand results for each protein
  select(
    Protein_Names, 
    estimate, 
    p.value, 
    conf.low, 
    conf.high, 
    method  
  )

# Display t-test results in a table
knitr::kable(
  Top15_ProtAbun_t_test_results,
  caption = "Protein Abundance T-tests",
  digits = 4
)

#' 
#' Only one of the top 15 proteins was found to be statisically significantly diffent between the Healty and ME/CFS groups. this was Alpha 2 HS glycoprotein. please not the the i used the welch t test metho as thet this method does not assume equal varrineces between the groups. 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Define significance levels for annotation
Top15_ProtAbun_t_test_results <- Top15_ProtAbun_t_test_results %>%
  mutate(
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**", 
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ),
    # Format p-values for display
    p.value_label = case_when(
      p.value < 0.001 ~ "p < 0.001",
      p.value < 0.01 ~ sprintf("p = %.3f", p.value),
      p.value < 0.05 ~ sprintf("p = %.3f", p.value),
      TRUE ~ sprintf("p = %.3f", p.value)
    )
  )

# Merge significance results with original data
ProteinAbundances_top15 <- ProteinAbundances_top15 %>%
  left_join(Top15_ProtAbun_t_test_results, by = "Protein_Names")

# Create significance data for annotation with protein-specific y-positions
significance_data <- ProteinAbundances_top15 %>%
  group_by(Protein_Names) %>%
  summarize(
    max_abundance = max(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(Top15_ProtAbun_t_test_results, by = "Protein_Names") %>%
  mutate(
    group1 = "Healthy",
    group2 = "MECFS",
    # Add some padding above the maximum value for each protein
    y_position = max_abundance * 1.1
  )

# Create the plot
p <- ggplot(ProteinAbundances_top15, aes(x = Disease_Status, y = Abundance, fill = Disease_Status)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 1) +
  scale_fill_viridis_c() + 
  # Add individual points for better data visualization
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  facet_wrap(
    ~Protein_Names, 
    ncol = 3,                  # Use 3 columns for better spacing
    scales = "free_y"          # Allow different y-scales for each protein
  ) +
  # Add significance brackets
  geom_segment(
    data = significance_data,
    aes(
      x = 1,
      xend = 2,
      y = y_position,
      yend = y_position
    ),
    inherit.aes = FALSE
  ) +
  # Add significance labels
  geom_text(
    data = significance_data,
    aes(
      x = 1.5,
      y = y_position * 1.05,
      label = paste(significance, "\n", p.value_label)
    ),
    inherit.aes = FALSE,
    size = 4  # Fix size for readability
  ) +
  theme_classic(base_size = 14) +  # Ensures a clean white background
  labs(
    title = "Box Plots of Top 15 Protein Abundances",
    x = "Disease Status",
    y = "Abundance"
  ) +
  scale_fill_manual(values = c("Healthy" = "blue", "MECFS" = "orange")) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    panel.spacing = unit(2, "lines"), # Increase space between plots
    plot.margin = margin(20, 20, 20, 20) # Add extra margin
  )

# Save the plot as a PNG file with high resolution
ggsave("Protein_Abundance_Boxplots.png", plot = p, width = 20, height = 18, dpi = 300, bg = "white")


#' 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


significance_data <- significance_data %>%
  select(Protein_Names, significance, p.value_label, y_position)

for (protein in unique(ProteinAbundances_top15$Protein_Names)) {
  clean_protein_name <- sub("OS=.*", "", protein)
  clean_protein_name <- trimws(clean_protein_name)

  protein_data <- ProteinAbundances_top15 %>% filter(Protein_Names == protein)
  sig_info <- significance_data %>% filter(Protein_Names == protein)

  p <- ggplot(protein_data, aes(x = Disease_Status, y = Abundance, fill = Disease_Status)) +
    geom_boxplot(outlier.color = "red", outlier.shape = 1) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    scale_fill_viridis_d(option = "D", labels = function(x) gsub("_", " ", x), begin = 0.2, end = 0.8) +
    theme_classic(base_size = 14) +
    labs(
      title = paste("Protein Abundance:", clean_protein_name),
      x = "Disease Status",
      y = "Abundance",
      fill = "Disease Status"
    ) +
    theme(
      strip.text = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 18, face = "bold"),
      plot.margin = margin(20, 20, 20, 20),
      panel.background = element_rect(fill = "NA", color = "NA"),
      plot.background = element_rect(fill = "NA", color = "NA")
    )

  if (nrow(sig_info) > 0) {
    p <- p + 
      geom_segment(
        data = sig_info,
        aes(x = 1, xend = 2, y = y_position, yend = y_position),
        inherit.aes = FALSE
      ) +
      geom_text(
        data = sig_info,
        aes(x = 1.5, y = y_position * 1.05, label = paste(significance, " ", p.value_label)),
        inherit.aes = FALSE,
        size = 5
      )
  }

  filename <- paste0("Protein_Abundance_", gsub("[^A-Za-z0-9]", "_", clean_protein_name), ".png")
  ggsave(filename, plot = p, width = 8, height = 6, dpi = 300, bg = NA)

  print(paste("Saved:", filename))
}

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
glycoPSMs <- glycoPSMs %>%
  mutate(
    Sample = str_extract(`Spectrum File`, "(?<=20250116_OE_TR_10S_MECFS_GPEP_)(.*?)(?=\\.raw)"),
    Disease_Status = ifelse(str_detect(Sample, "^HC"), "Healthy", 
                            ifelse(str_detect(Sample, "M"), "MECFS", NA)),
    Protein_Names = str_extract(`Master Protein Descriptions`, "^[^O]*(?=OS=)")
  )



## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
write.csv(ProteinAbundances, file = "ProteinAbundances_wide.csv")

#' 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
write.csv(ProteinAbundances_long, file = "ProteinAbundances_long.csv")


#' the following was caluatated and ploted from the glycoprotein data. As you may notice there many immugoblulin proteins found to be differenaly aubandant based on this data. It may be better use the protein groups data instead as this will group the immugoblulin family proteins together in a way that lead to less false posivies generated by the spreaded divesity of immuloglin proteins can take. 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Perform t-tests and extract results
ProtAbun_t_test_results <- ProteinAbundances_long %>%
  group_by(Protein_Names) %>%
  summarise(
    t_test = list(
      tryCatch(
        t.test(
          Abundance ~ Disease_Status, 
          data = pick(everything()),  # Use `pick()` instead of `cur_data()`
          var.equal = TRUE  #students T-test
        ), 
        error = function(e) {message(paste("Error in", Protein_Names, ":", e$message)); NULL}  # Handle errors gracefully and log message
      )
    ),
    .groups = "drop"
  ) %>%
  mutate(
    t_test_summary = map(t_test, ~ if (!is.null(.x)) broom::tidy(.x) else NULL)
  ) %>%
  unnest(cols = t_test_summary, keep_empty = TRUE) %>% # Expand results for each protein
  mutate(
    p.adjusted = p.adjust(p.value, method = "BH"),  # Benjamini-Hochberg adjustment
      `-Log10_p.value` = -log10(p.adjusted)  # -Log10 transformation
    ) %>%
  left_join(
    ProteinAbundances_long %>%
      group_by(Protein_Names, Disease_Status) %>%
      summarise(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = 'drop') %>%
      pivot_wider(names_from = Disease_Status, values_from = mean_abundance) %>%
      mutate(Fold_Change = MECFS / Healthy, Log2_FC = log2(Fold_Change)),
    by = 'Protein_Names'
  ) %>%
  select(
    Protein_Names, 
    estimate, 
    p.value, 
    p.adjusted, 
    `-Log10_p.value`,
    Fold_Change,
    Log2_FC, 
    conf.low, 
    conf.high, 
    method  
  )

#' 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
DDA_DAPs <- ProtAbun_t_test_results %>% 
  filter(Log2_FC > 0.5 & `-Log10_p.value` < 0.05) 

write.csv(DIA_DAPs, file = "DDA_DAPs.csv")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Create a basic volcano plot
ggplot(data = ProtAbun_t_test_results, aes(x = Log2_FC, y = `-Log10_p.value`)) +
  geom_vline(xintercept = c(-1.2, 1.2), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.3), col = "gray", linetype = 'dashed') + 
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
  xlim = c(-5, 5), # Adjust the range based on your data distribution
  ylim = c(0, max(ProtAbun_t_test_results$`-Log10_p.value`, na.rm = TRUE) + 0.08), # Set y-axis limits
  pCutoff = 0.3,  # Optional: Set a cutoff for significance
  FCcutoff = 1,    # Optional: Set a cutoff for fold change
  labSize = 3.0, # Optional: Adjust label size
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
wb <- wb_load("/Users/thomasreilly/Desktop/MRes Data/10S_PD_quant/Charlie_DIA/MRCFS Proteomics Results.xlsx")
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
  ylim = c(0, 6), # Adjust this based on your data
  xlim = c(-1.2, 1.2), # Adjust the range based on your data distribution
  pCutoff = 0.05,  # Optional: Set a cutoff for significance
  FCcutoff = 0.5,    # Optional: Set a cutoff for fold change
  labSize = 3.0, # Optional: Adjust label size
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
knitr::kable(
  DIA_Prot_DiffAb,
  caption = "DIA Protein Abundances",
  digits = 2,
  na = 'NA'
)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
DIA_DAPs <- DIA_Prot_DiffAb %>% 
  filter(Log2_FC > 0.5 & pval_proteins < 0.05) 

write.csv(DIA_DAPs, file = "DIA_DAPs.csv")

knitr::kable(
  DIA_DAPs,
  caption = "DIA plasma proteomics - Differentaly Aubundant Proteins",
  digits = 2,
  na = 'NA'
)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
DIA_DAPs %>%
  summarise(n_distinct(protein))

#' 
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
    disease_status = ifelse(str_detect(sample, "HC"), "Healthy", 
                            ifelse(str_detect(sample, "M"), "MECFS", NA_character_))) 
 


#' 
#' 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
glycan_class_map <- read_csv(file = "glycan_class_map.csv")


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


count_sia <- function(input_df) {
  # Ensure required column exists
  if (!"glycan_composition" %in% colnames(input_df)) {
    stop("Column 'glycan_composition' not found in input dataframe")
  }
  
  # Extract number of NeuAc occurrences and add as new column
  input_df$sia_count <- sapply(input_df$glycan_composition, function(x) {
    match <- regmatches(x, regexpr("neuac\\((\\d+)\\)", x, perl = TRUE))
    if (length(match) > 0) {
      paste("NeuAc", sub("neuac\\((\\d+)\\)", "\\1", match))
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
    match <- regmatches(x, regexpr("fuc\\((\\d+)\\)", x, perl = TRUE))
    if (length(match) > 0) {
      paste("Fucose", sub("fuc\\((\\d+)\\)", "\\1", match))
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
gly_prot_list <- glycoPSMs %>%
  dplyr::select(gene_name) %>%  # Ensure we are using dplyr's select
  distinct() %>%
  pull(gene_name) %>%
  list()  # Convert the vector to a list


  

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
  names(data_list) <- ifelse(is.na(names(data_list)) | names(data_list) == "", 
                             paste0("Sheet_", seq_along(data_list)), 
                             names(data_list))
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
  names(data_list) <- ifelse(is.na(names(data_list)) | names(data_list) == "", 
                             paste0("Sheet_", seq_along(data_list)), 
                             names(data_list))
  
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
# Check the number of groups
length(glycoPSMs %>% group_by(protein_accessions, protein_glycosite) %>% group_split())

# Check the number of unique names
glycoPSMs %>%
  distinct(protein_accessions, protein_glycosite) %>%
  mutate(name = paste(protein_accessions, protein_glycosite, sep = "_")) %>%
  pull(name) %>%
  length()


#' 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Check if lengths match
length(glycosite_list)  # Number of groups
length(unique(glycoPSMs$gsite_ID))  # Number of unique site IDs

# Preview the first few names assigned
names(glycosite_list)[1:20]
unique(glycoPSMs$gsite_ID)[1:20]


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
check_df <-  gsitePSM_list[["P01023_55"]]


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
    left_join(input_df %>% select(sample, disease_status) %>% distinct(), by = "sample") %>%
    select(sample, disease_status, glycan_composition, relative_abundance)
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
    left_join(input_df %>% select(sample, disease_status) %>% distinct(), by = "sample") %>%
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
    left_join(input_df %>% select(sample, disease_status) %>% distinct(), by = "sample") %>%
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
    left_join(input_df %>% select(sample, disease_status) %>% distinct(), by = "sample") %>%
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
    left_join(input_df %>% select(sample, disease_status) %>% distinct(), by = "sample") %>%
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
    left_join(input_df %>% select(sample, disease_status) %>% distinct(), by = "sample") %>%
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

# Function to perform analysis on a list of data frames with a filename prefix option
perform_analysis_on_list <- function(data_list, group_col, value_col, file_prefix = "Analysis_") {
  
  output_file <- paste0(file_prefix, "analysis_results.xlsx")  # Create the output filename
  wb <- wb_workbook()  # Create a new workbook for results

  for (i in seq_along(data_list)) {
    data <- data_list[[i]]  # Extract the individual data frame
    dataset_name <- paste0("Dataset_", i)  # Naming convention for each dataset

    # Count samples in each disease status group
    disease_counts <- data %>%
      group_by(!!sym(group_col)) %>%
      summarise(count = n(), .groups = "drop")

    # Check if there are at least 3 samples in each group
    if (all(c("Healthy", "MECFS") %in% data[[group_col]]) &&
        sum(data[[group_col]] == "Healthy") >= 3 &&
        sum(data[[group_col]] == "MECFS") >= 3) {

      # Calculate descriptive statistics
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

      # Perform Welch t-test
      t_test_result <- t.test(
        !!sym(value_col) ~ !!sym(group_col), 
        data = data, 
        var.equal = FALSE
      )

      # Perform ANOVA
      anova_result <- aov(!!sym(value_col) ~ !!sym(group_col), data = data)
      anova_summary <- summary(anova_result)

     # Prepare results for Excel using a stepwise approach
descriptive_stats_name <- paste0(dataset_name, "_Descriptive Stats")
t_test_name <- paste0(dataset_name, "_T-test Result")
anova_name <- paste0(dataset_name, "_ANOVA Result")

# Create each component separately
descriptive_stats_data <- descriptive_stats
t_test_data <- data.frame(
  Statistic = t_test_result$statistic,
  P_Value = t_test_result$p.value,
  DF = t_test_result$parameter,
  Confidence_Lower = t_test_result$conf.int[1],
  Confidence_Upper = t_test_result$conf.int[2]
)
anova_data <- as.data.frame(anova_summary[[1]])

# Combine into results list
results_list <- list(
  descriptive_stats_name = descriptive_stats_data,
  t_test_name = t_test_data,
  anova_name = anova_data
)

# Rename list elements properly
names(results_list) <- c(descriptive_stats_name, t_test_name, anova_name)



    # Pivot data to wide format by sample column
wide_data <- data %>%
  pivot_wider(
    names_from = sample,
    values_from = !!sym(value_col),
    values_fill = NA   # Fill missing values with NA
  )

# Add pivoted data to results list
results_list[[paste0(dataset_name, "_Wide Data")]] <- wide_data

      # If t-test is significant (p < 0.05), create a boxplot
      if (t_test_result$p.value < 0.05) {
        boxplot <- ggplot(data, aes(x = !!sym(group_col), y = !!sym(value_col))) +
          geom_boxplot(aes(fill = !!sym(group_col))) +
          theme_minimal() +
          labs(title = paste("Boxplot for", dataset_name),
               x = "Disease Status",
               y = "Measurement Value") +
          scale_fill_manual(values = c("Healthy" = "blue", "MECFS" = "red"))
        
        # Save the boxplot with prefix
        plot_filename <- paste0(file_prefix, "boxplot_", i, ".png")
        ggsave(plot_filename, plot = boxplot, width = 6, height = 4)
        
        # Add image path to results
        results_list[[paste0(dataset_name, "_Boxplot")]] <- plot_filename
      }

      # Write all results for this dataset to the workbook
      for (name in names(results_list)) {
        wb <- wb_add_worksheet(wb, name)  # Add worksheet
        wb <- wb_add_data(wb, sheet = name, x = results_list[[name]])  # Add data to the worksheet
      }
      
    } else {
      message("Skipping ", dataset_name, ": Not enough samples in each disease category (at least 3 per group required).")
    }
  }

  # Save the workbook with prefix
  wb_save(wb, output_file)
  message("Analysis complete. Results saved to: ", output_file)
}

# Example usage:
# Assume `data_list` is a list of data frames, each containing `disease_status` and `measurement_value`
# Example call with a custom prefix:

#perform_analysis_on_list(your_data_list, "disease_status", "measurement_value", "MyResults_")


#' 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
perform_analysis_on_list(Sia_glycoprot_list, "disease_status", "relative_abundance", "Sia_gprot_")
perform_analysis_on_list(glycomp_glycoprot_list, "disease_status", "relative_abundance", "Glycomp_gprot_")
perform_analysis_on_list(glyc_class_glycoprot_list, "disease_status", "relative_abundance", "glyClass_gprot_")
perform_analysis_on_list(fuc_count_glycoprot_list, "disease_status", "relative_abundance", "Fuc_count_gprot_")
perform_analysis_on_list(sia_count_glycoprot_list, "disease_status", "relative_abundance", "Sia_count_gprot_")
perform_analysis_on_list(fuc_glycoprot_list, "disease_status", "relative_abundance", "Fuc_gprot_")

perform_analysis_on_list(Sia_glycosite_list, "disease_status", "relative_abundance", "Sia_gsite_")
perform_analysis_on_list(fuc_glycosite_list, "disease_status", "relative_abundance", "Fuc_gsite_")
perform_analysis_on_list(glycomp_glycosite_list, "disease_status", "relative_abundance", "Glycomp_gsite_")
perform_analysis_on_list(glyc_class_glycosite_list, "disease_status", "relative_abundance", "glyClass_gsite_")
perform_analysis_on_list(fuc_count_glycosite_list, "disease_status", "relative_abundance", "Fuc_count_gsite_")
perform_analysis_on_list(sia_count_glycosite_list, "disease_status", "relative_abundance", "Sia_count_gsite_")


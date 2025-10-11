library(dplyr)
library(readr)
library("openxlsx2")  # for wb_* helpers
library(janitor)
library(EnhancedVolcano)



#' Have done differenial protein aubance analysis of the Charlie's DIA protemics data, as this is more correct dataset to perform analysis on then compared to my own glycoproteomics data set. this is because the glycoproteomics datasets has been biased by glycopeptide enrichment, therefore this may lead to inaccurate quants due to possible greater llels of missing values . the glycoproteomics data set also suffer from being collected by a DDA exmerimental method as this is best for glycoproteoms as the chimeric spectra would lead to signifcant issues for glycan asigments to peptides. hthe DDA method further creates issues of missing values at random due to ticinal featues inherant in the method.
             ## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             wb <- wb_load(
               "input_data/MRCFS Proteomics Results.xlsx"
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
               title = "DIA proteomics differentially abundant proteins",
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
               widthConnectors = 0.75,
               max.overlaps = 20
             )
             
             ggsave("DIA_prot_volcano.png", width = 140, height = 100, units = "mm")
             knitr::include_graphics("DIA_prot_volcano.png")


# 1) Build the significant table using the exact volcano thresholds
sig_tbl <- DIA_Prot_DiffAb %>%
  # defensively remove impossible p-values
  mutate(
    pval_proteins = ifelse(is.finite(pval_proteins) & pval_proteins > 0, pval_proteins, NA_real_),
    neglog10_p    = -log10(pval_proteins),
    regulation    = case_when(
      !is.na(pval_proteins) & pval_proteins <= 0.05 & Log2_FC >=  0.5 ~ "Up",
      !is.na(pval_proteins) & pval_proteins <= 0.05 & Log2_FC <= -0.5 ~ "Down",
      TRUE ~ "NS"
    ),
    significant   = regulation != "NS"
  ) %>%
  filter(significant) %>%
  transmute(
    protein,
    Log2_FC,
    fc_proteins,
    pval_proteins,
    neglog10_p,
    regulation
  ) %>%
  arrange(pval_proteins, desc(abs(Log2_FC)))

# 2) Preview
print(sig_tbl, n = 30, na.print = "")

# 3) Save to CSV
dir.create("output_data", showWarnings = FALSE, recursive = TRUE)
write_csv(sig_tbl, "output_data/DIA_significant_proteins.csv")

             
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
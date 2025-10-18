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

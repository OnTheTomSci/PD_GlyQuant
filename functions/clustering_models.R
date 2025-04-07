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


# Open a graphics device (optional, if saving to file)
png("scores_plot.png", width = 800, height = 600)

# Create the scores plot
mdaplot(m$res$cal$scores, type = "p", 
        show.labels = TRUE, 
        show.lines = c(0, 0), 
        cgroup = protein_gly_comp_df$disease_status)

# Call plotConfidenceEllipse with only the model (ensure 'p' is defined)
plotConfidenceEllipse(m$res$cal)

# Close the device to save the file
dev.off()

# Open a graphics device (optional, if saving to file)
png("PCA_plots.png", width = 800, height = 600)
plot(m, show.labels = TRUE)
dev.off()


library(mdatools)
m <- pca(protein_gly_comp_df, 7, info = "Plasma Protein N-glycosolation HC vs MECFS PCA Model")
m = selectCompNum(m, 5)
print(m)
m$loadings[1:4, 1:4]

# Open a graphics device (optional, if saving to file)
png("loadings_plot.png", width = 800, height = 600)

# Random Forest
library(randomForest)
set.seed(20)
rf.fit <- randomForest(disease_status ~ ., data = protein_gly_comp_df)
print(rf.fit)

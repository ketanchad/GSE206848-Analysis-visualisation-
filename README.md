# GSE206848-Analysis-visualisation-
 Dysregulated gene expression in human osteoarthritic (OA) and rheumatoid arthritis (RA) synovium using by R studio bioinformatics analysis. 
library(GEOquery)
library(limma)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(dplyr)

# Download GEO Dataset
get <- getGEO("GSE206848", GSEMatrix = TRUE)
length(get)
exprSet <- exprs(get[[1]])
head(exprSet)
dim(exprSet)


#View Sample Information
pdata <- pData(get[[1]])
head(pdata)
colnames(pdata)

# Check Expression Data   
summary(exprSet)

boxplot(exprSet,
        las = 2,
        col = "lightblue",
        main = "Raw Expression Data")

#Log2 Transformation (if needed)

exprSet_log <- log2(exprSet + 1)

boxplot(exprSet_log,
        las = 2,
        col = "orange",
        main = "Log2 Normalized Data")

# Remove zero variance genes
expr_filtered <- exprSet_log[
  apply(exprSet_log, 1, var) != 0,
]

# Run PCA
pca <- prcomp(t(expr_filtered), scale. = TRUE)

# Create dataframe
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2]
)

# PCA plot
library(ggplot2)

ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(size = 4, color = "blue") +
  theme_minimal() +
  ggtitle("PCA Plot")


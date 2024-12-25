# ---------------------------------------------------
# Single-Cell RNA-seq Analysis Workflow using Seurat
# ---------------------------------------------------

# ----------------------------
# 1. Loading Required Libraries
# ----------------------------

# Load necessary libraries
library(dplyr)       # Data manipulation
library(Seurat)
packageVersion("Seurat")

library(patchwork)   # Combining ggplot2 plots
library(ggplot2)     # Data visualization

# If any of the packages are not installed, uncomment and run the following lines:
# install.packages("dplyr")
# install.packages("patchwork")
# install.packages("ggplot2")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("Seurat")

# ----------------------------
# 2. Loading the Seurat Object
# ----------------------------

# Define file paths
input_rds <- "/home/y/local/work/Using/raw/mac.atlas.zenodo.200524.rds"
output_rds <- "/home/y/local/work/Using/final/modified_file.rds"

# Read the Seurat object from the RDS file
pbmc <- readRDS(input_rds)

# ----------------------------
# 3. Initial Exploration
# ----------------------------

# Display a summary of the Seurat object
print(pbmc)

# Inspect metadata
head(pbmc@meta.data)

# List available assays in the object
print(Assays(pbmc))

# Count the number of cells and genes
num_cells <- ncol(pbmc)
num_genes <- nrow(pbmc)
cat("Number of cells:", num_cells, "\n")
cat("Number of genes:", num_genes, "\n")

# Check the top variable features
head(VariableFeatures(pbmc))

# List available dimensional reductions
print(Reductions(pbmc))

# ---------------------------------------------------
# Assuming pbmc is your Seurat object and the 'integrated' assay is set
# ---------------------------------------------------

# Set the default assay to 'integrated'
DefaultAssay(pbmc) <- "integrated"

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))

# Generate Elbow Plot to determine the number of PCA dimensions to use
ElbowPlot(pbmc, ndims = 50) + ggtitle("Elbow Plot for PCA")

# Select PCA dimensions (based on Elbow Plot results, e.g., 1:15)
dims_use <- 1:15

# Ensure the default assay is set to 'integrated'
DefaultAssay(pbmc) <- "integrated"

# Run FindNeighbors
pbmc <- FindNeighbors(pbmc, dims = dims_use)

# Run FindClusters
pbmc <- FindClusters(pbmc, graph.name = "integrated_snn", resolution = 0.4)

# View the distribution of clusters
table(pbmc$seurat_clusters)

# Inspect the first few rows of relevant metadata columns
head(pbmc@meta.data[, c("seurat_clusters", "cluster", "integrated_snn_res.0.75", "integrated_snn_res.0.7", "integrated_snn_res.0.65")])

# ---------------------------------------------------
# Copy the 'cluster' column to the 'seurat_clusters' column
# ---------------------------------------------------
pbmc$seurat_clusters <- pbmc$cluster

# ---------------------------------------------------
# Verify that the copy was successful
# ---------------------------------------------------
table(pbmc$seurat_clusters)

# ---------------------------------------------------
# Set 'seurat_clusters' as the active identity
# ---------------------------------------------------
Idents(pbmc) <- "seurat_clusters"

# ---------------------------------------------------
# Verify the active identity
# ---------------------------------------------------
Idents(pbmc)

# ---------------------------------------------------
# Run UMAP for dimensionality reduction
# ---------------------------------------------------
pbmc <- RunUMAP(pbmc, dims = dims_use)

# ---------------------------------------------------
# Visualize the clustering results. The result lacks clusters, and the image is gray. How to fix?
# ---------------------------------------------------
DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("UMAP - Seurat Clusters")

# ---------------------------------------------------
# Identify cluster-specific marker genes
# ---------------------------------------------------
cluster_markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# ---------------------------------------------------
# View the top five marker genes for each cluster
# ---------------------------------------------------
top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

print("Top 5 Markers per Cluster:")
print(top_markers)

# ---------------------------------------------------
# Visualize marker genes using a heatmap
# ---------------------------------------------------
DoHeatmap(pbmc, features = top_markers$gene) + 
  NoLegend() +
  ggtitle("Top Markers Heatmap")

# ---------------------------------------------------
# Save the processed Seurat object
# ---------------------------------------------------
saveRDS(pbmc, file = output_rds)
cat("Processed Seurat object saved to:", output_rds, "\n")

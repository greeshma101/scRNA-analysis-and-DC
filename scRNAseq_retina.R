#Public dataset form GSE153674 "An optimized protocol for retina single-cell RNA sequencing"

# scRNA-seq Preprocessing & Clustering Pipeline (10X h5 input)

# Load packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("rhdf5", quietly = TRUE)) BiocManager::install("rhdf5")

packages <- c("Seurat", "ggplot2", "patchwork", "dplyr", "rhdf5")
lapply(packages, library, character.only = TRUE)

#-----------------------------------------------------
# Load 10X HDF5 scRNA-seq data

h5_file <- "GSM4649096_mm_ret_M1_C57OS_filtered.h5"
scRNA <- Read10X_h5(h5_file)

# Optional: check structure
head(h5read(h5_file, "matrix/barcodes"))
head(h5read(h5_file, "matrix/features/name"))
head(h5read(h5_file, "matrix/data"))

#-----------------------------------------------------
# Quality control 
# Sparsity calculation
total.zeros <- sum(scRNA == 0)
total.counts <- nrow(scRNA) * ncol(scRNA)
sparsity <- (total.zeros / total.counts) * 100
message("Matrix sparsity: ", round(sparsity, 2), "%")  # Expected ~90–95% for scRNA-seq

# Create Seurat object & basic filtering 

# Keep genes in ≥3 cells and cells with ≥200 genes
seurat <- CreateSeuratObject(counts = scRNA, min.cells = 3, min.features = 200)


# Add mitochondrial % and visualize QC metrics

seurat$mito <- PercentageFeatureSet(seurat, pattern = "^mt-")
meta <- seurat@meta.data

# Visualize QC metrics
p1 <- ggplot(meta, aes(y = nCount_RNA)) + 
  geom_violin(fill = "grey80") + geom_jitter(width = 0.3, alpha = 0.3) +
  scale_y_log10() + labs(title = "UMI Counts per Cell") + theme_classic()

p2 <- ggplot(meta, aes(y = nFeature_RNA)) +
  geom_violin(fill = "grey80") + geom_jitter(width = 0.3, alpha = 0.3) +
  scale_y_log10() + labs(title = "Genes per Cell") + theme_classic()

p3 <- ggplot(meta, aes(y = mito)) +
  geom_violin(fill = "grey80") + geom_jitter(width = 0.3, alpha = 0.3) +
  labs(title = "% Mitochondrial Content") + theme_classic()

p1 + p2 + p3 + plot_layout(nrow = 2)


# Filter out poor-quality cells

seurat <- subset(seurat, subset = mito <= 20)

#-----------------------------------------------------
# Normalize & identify variable features

seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
top20 <- head(VariableFeatures(seurat), 20)

LabelPoints(
  plot = VariableFeaturePlot(seurat),
  points = top20, ynudge = -2
)

#-----------------------------------------------------
# Scaling and PCA
#Without scaling, high-expression genes would dominate the PCA, biasing the analysis.

seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat)

ElbowPlot(seurat, ndims = 50)
DimPlot(seurat, reduction = "pca")

#-----------------------------------------------------
# Clustering & UMAP

seurat <- FindNeighbors(seurat, dims = 1:40)
seurat <- FindClusters(seurat, resolution = 0.8)
seurat <- RunUMAP(seurat, dims = 1:40)

DimPlot(seurat, reduction = "umap", label = TRUE)

#-----------------------------------------------------
# Identify marker genes

markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Top 5 markers per cluster
top5 <- markers %>%
  group_by(cluster) %>%
  top_n(5, avg_log2FC)

View(top5)

#-----------------------------------------------------
 #Manual cell type annotation

cluster.ids <- c(
  "0" = "Rods", "1" = "Rods", "2" = "Rods",
  "3" = "RGC", "4" = "Müller Glia", "5" = "Bipolar cells",
  "6" = "Cones", "7" = "Cones", "8" = "Rods",
  "9" = "Bipolar cells", "10" = "Microglia"
)

Idents(seurat) <- "seurat_clusters"
seurat <- RenameIdents(seurat, cluster.ids)
seurat$celltype <- Idents(seurat)

DimPlot(seurat, reduction = "umap", label = TRUE) + NoLegend()

#-----------------------------------------------------
# Generate reference matrix (per cell type average)
#reference_matrix <- AverageExpression(seurat, group.by = "celltype", assays = "RNA")$RNA
#dim(reference_matrix)
#head(reference_matrix[, 1:5])

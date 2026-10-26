# Companion code for Visium - Guided Clustering Tutorial
# Review paths, QC rules, and example parameters before running.

# libraries
library(Seurat)
library(ggplot2)
library(patchwork)

set.seed(42)
packageVersion("Seurat")
packageVersion("SeuratObject")

# download-data
download_visium <- function(prefix, data.dir) {
  dir.create(data.dir, recursive = TRUE, showWarnings = FALSE)
  download.file(
    paste0(prefix, "filtered_feature_bc_matrix.h5"),
    destfile = file.path(data.dir, "filtered_feature_bc_matrix.h5"),
    mode = "wb"
  )
  spatial.tar <- file.path(data.dir, "spatial.tar.gz")
  download.file(paste0(prefix, "spatial.tar.gz"), spatial.tar, mode = "wb")
  untar(spatial.tar, exdir = data.dir)
}

v1.prefix <- paste0(
  "https://cf.10xgenomics.com/samples/spatial-exp/1.0.0/",
  "V1_Mouse_Brain_Sagittal_Anterior/V1_Mouse_Brain_Sagittal_Anterior_"
)
v2.prefix <- paste0(
  "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/",
  "CytAssist_FFPE_Mouse_Brain_Rep1/CytAssist_FFPE_Mouse_Brain_Rep1_"
)

# Skip these downloads if the files are already available locally.
download_visium(v1.prefix, "data/visium_v1")
download_visium(v2.prefix, "data/visium_v2")

# load-standard
v1 <- Load10X_Spatial(
  data.dir = "data/visium_v1",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", slice = "v1", filter.matrix = TRUE
)
v2 <- Load10X_Spatial(
  data.dir = "data/visium_v2",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", slice = "v2", filter.matrix = TRUE
)

# We begin with v1, then repeat the same workflow on v2.
visium <- v1
DefaultAssay(visium) <- "Spatial"
visium

# inspect-counts
counts <- GetAssayData(visium, assay = "Spatial", layer = "counts")
counts[1:5, 1:6]
Images(visium)
head(visium[[]])

# qc-spots
VlnPlot(
  visium,
  features = c("nCount_Spatial", "nFeature_Spatial"),
  pt.size = 0, ncol = 2
)

SpatialFeaturePlot(
  visium, features = "nCount_Spatial", images = "v1"
)

# qc-mito
mt.features <- grep("^(MT-|mt-)", rownames(visium[["Spatial"]]), value = TRUE)
if (length(mt.features) > 0L) {
  visium[["percent.mt"]] <- PercentageFeatureSet(
    visium, features = mt.features, assay = "Spatial"
  )
  VlnPlot(visium, features = "percent.mt", pt.size = 0)
  SpatialFeaturePlot(visium, features = "percent.mt", images = "v1")
}

# filter-spots
min.features <- 1L
visium$qc_keep <- (
  visium$nCount_Spatial > 0 &
  visium$nFeature_Spatial >= min.features
)

table(visium$qc_keep)
SpatialDimPlot(visium, images = "v1", group.by = "qc_keep")
visium <- subset(visium, subset = qc_keep)

# normalize-spots
visium <- NormalizeData(
  visium,
  assay = "Spatial",
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# variable-spots
visium <- FindVariableFeatures(
  visium,
  assay = "Spatial",
  selection.method = "vst",
  nfeatures = 2000
)

head(VariableFeatures(visium), 10)
VariableFeaturePlot(visium)

# scale-spots
visium <- ScaleData(
  visium,
  assay = "Spatial",
  features = VariableFeatures(visium)
)

# pca-spots
visium <- RunPCA(
  visium,
  assay = "Spatial",
  features = VariableFeatures(visium),
  npcs = 50,
  seed.use = 42
)

print(visium[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(visium, dims = 1:6, cells = 500, balanced = TRUE)

# elbow-spots
ElbowPlot(visium, ndims = 50)
dims.use <- 1:30

# cluster-spots
visium <- FindNeighbors(
  visium, reduction = "pca", dims = dims.use
)
visium <- FindClusters(
  visium, resolution = 0.5, random.seed = 42
)

table(visium$seurat_clusters)

# umap-spots
visium <- RunUMAP(
  visium,
  reduction = "pca",
  dims = dims.use,
  seed.use = 42
)

p.umap <- DimPlot(
  visium, reduction = "umap",
  group.by = "seurat_clusters", label = TRUE, repel = TRUE
)
p.tissue <- SpatialDimPlot(
  visium, images = "v1", group.by = "seurat_clusters",
  label = TRUE, label.size = 3
)
p.umap | p.tissue

# standard-function
run_spot_workflow <- function(x, dims = 1:30, resolution = 0.5) {
  DefaultAssay(x) <- "Spatial"
  x$qc_keep <- x$nCount_Spatial > 0 & x$nFeature_Spatial >= 1
  x <- subset(x, subset = qc_keep)
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x, features = VariableFeatures(x))
  x <- RunPCA(x, features = VariableFeatures(x), npcs = 50, seed.use = 42)
  x <- FindNeighbors(x, reduction = "pca", dims = dims)
  x <- FindClusters(x, resolution = resolution, random.seed = 42)
  x <- RunUMAP(x, reduction = "pca", dims = dims, seed.use = 42)
  x
}

VlnPlot(v2, features = c("nCount_Spatial", "nFeature_Spatial"), pt.size = 0)
SpatialFeaturePlot(v2, features = "nCount_Spatial", images = "v2")

v2 <- run_spot_workflow(v2)
DimPlot(v2, reduction = "umap", group.by = "seurat_clusters", label = TRUE) |
  SpatialDimPlot(v2, images = "v2", group.by = "seurat_clusters", label = TRUE)

# save-standard
# For the standard Visium object analyzed above:
saveRDS(visium, "visium_v1_processed.rds")
saveRDS(v2, "visium_v2_processed.rds")
write.csv(visium[[]], "visium_v1_metadata.csv")
write.csv(Embeddings(visium, reduction = "umap"), "visium_v1_umap.csv")

# Save the corresponding object separately after running an HD branch.
# saveRDS(hd, "visium_hd_bins_processed.rds")
# saveRDS(hd.cells, "visium_hd_cells_processed.rds")

# session-command
sessionInfo()
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

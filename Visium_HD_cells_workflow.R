# Companion code for Visium - Guided Clustering Tutorial
# Review paths, QC rules, and example parameters before running.

# libraries
library(Seurat)
library(ggplot2)
library(patchwork)

set.seed(42)
packageVersion("Seurat")
packageVersion("SeuratObject")

# cells-load
stopifnot(packageVersion("Seurat") >= "5.4.0")
stopifnot(packageVersion("SeuratObject") >= "5.3.0")
stopifnot(requireNamespace("arrow", quietly = TRUE))
stopifnot(requireNamespace("sf", quietly = TRUE))

hd.cells <- Load10X_Spatial(
  data.dir = "path/to/HD_sample/outs",
  bin.size = "polygons",
  slice = "HD",
  image.name = "tissue_hires_image.png"
)
DefaultAssay(hd.cells) <- "Spatial.Polygons"
Assays(hd.cells)
Images(hd.cells)  # HD.polygons

# cells-process
VlnPlot(
  hd.cells,
  features = c("nCount_Spatial.Polygons", "nFeature_Spatial.Polygons"),
  pt.size = 0
)

qc.cells <- hd.cells[[]]
keep.cells <- rownames(qc.cells)[qc.cells$nCount_Spatial.Polygons > 0]
hd.cells <- subset(hd.cells, cells = keep.cells)

hd.cells <- NormalizeData(hd.cells, normalization.method = "LogNormalize", scale.factor = 10000)
hd.cells <- FindVariableFeatures(hd.cells, nfeatures = 2000)
hd.cells <- ScaleData(hd.cells, features = VariableFeatures(hd.cells))
hd.cells <- RunPCA(hd.cells, features = VariableFeatures(hd.cells), npcs = 50, seed.use = 42)

ElbowPlot(hd.cells, ndims = 50)
cell.dims <- 1:30
hd.cells <- FindNeighbors(hd.cells, reduction = "pca", dims = cell.dims)
hd.cells <- FindClusters(
  hd.cells, resolution = 0.5,
  cluster.name = "segmentation_clusters", random.seed = 42
)
hd.cells <- RunUMAP(hd.cells, reduction = "pca", dims = cell.dims, seed.use = 42)

# cells-plot
p.cell.umap <- DimPlot(
  hd.cells, reduction = "umap",
  group.by = "segmentation_clusters", label = TRUE
)
p.cell.tissue <- SpatialDimPlot(
  hd.cells, images = "HD.polygons",
  group.by = "segmentation_clusters", plot_segmentations = TRUE
)
p.cell.umap | p.cell.tissue

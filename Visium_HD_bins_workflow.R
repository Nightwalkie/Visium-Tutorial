# Companion code for Visium - Guided Clustering Tutorial
# Review paths, QC rules, and example parameters before running.

# libraries
library(Seurat)
library(ggplot2)
library(patchwork)

set.seed(42)
packageVersion("Seurat")
packageVersion("SeuratObject")

# hd-load
stopifnot(packageVersion("Seurat") >= "5.1.0")
stopifnot(requireNamespace("hdf5r", quietly = TRUE))
stopifnot(requireNamespace("arrow", quietly = TRUE))

hd <- Load10X_Spatial(
  data.dir = "path/to/HD_sample/outs",
  bin.size = 8,
  slice = "HD"
)

full.assay <- "Spatial.008um"
full.image <- "HD.008um"
DefaultAssay(hd) <- full.assay

Assays(hd)
Images(hd)
stopifnot(all(Cells(hd[[full.assay]]) %in% Cells(hd[[full.image]])))

# hd-qc
VlnPlot(
  hd, features = c("nCount_Spatial.008um", "nFeature_Spatial.008um"),
  pt.size = 0
)
SpatialFeaturePlot(hd, features = "nCount_Spatial.008um", images = full.image)

qc <- hd[[]]
keep.bins <- rownames(qc)[qc$nCount_Spatial.008um > 0]
hd <- subset(hd, cells = keep.bins)

hd <- NormalizeData(hd, normalization.method = "LogNormalize", scale.factor = 10000)
hd <- FindVariableFeatures(hd, selection.method = "vst", nfeatures = 2000)

# hd-sketch
hd <- SketchData(
  object = hd,
  assay = full.assay,
  ncells = min(50000L, ncol(hd[[full.assay]])),
  method = "LeverageScore",
  sketched.assay = "sketch",
  seed = 42
)

DefaultAssay(hd) <- "sketch"
hd <- FindVariableFeatures(hd, nfeatures = 2000)
hd <- ScaleData(hd, features = VariableFeatures(hd))
hd <- RunPCA(
  hd, assay = "sketch", features = VariableFeatures(hd),
  npcs = 50, reduction.name = "pca.sketch", seed.use = 42
)

ElbowPlot(hd, reduction = "pca.sketch", ndims = 50)
hd.dims <- 1:30

# hd-cluster
hd <- FindNeighbors(
  hd, assay = "sketch", reduction = "pca.sketch", dims = hd.dims,
  graph.name = c("hd_sketch_nn", "hd_sketch_snn")
)
hd <- FindClusters(
  hd, graph.name = "hd_sketch_snn", cluster.name = "hd_cluster_sketch",
  resolution = 0.5, random.seed = 42
)
hd <- RunUMAP(
  hd, reduction = "pca.sketch", dims = hd.dims,
  reduction.name = "umap.sketch", return.model = TRUE, seed.use = 42
)

hd <- ProjectData(
  object = hd,
  assay = full.assay,
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = hd.dims,
  refdata = list(hd_cluster_full = "hd_cluster_sketch")
)

# hd-plot
DefaultAssay(hd) <- full.assay
Idents(hd) <- "hd_cluster_full"
Reductions(hd)

p.hd.umap <- DimPlot(
  hd, reduction = "full.umap.sketch",
  group.by = "hd_cluster_full", label = TRUE, repel = TRUE
)
p.hd.tissue <- SpatialDimPlot(
  hd, images = full.image,
  group.by = "hd_cluster_full", label = TRUE, label.size = 3
)
p.hd.umap | p.hd.tissue

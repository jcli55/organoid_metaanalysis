library(ArchR)
library(parallel)

# After creating the ArchR project, this script runs dimentional reduction using Iterative LSI, clustering and embedding with UMAP
# Uses Harmony to batch correct, assumes already ran archr_clustering (specifically doublet removal and iterative LSI)
# 4/2/24 - Jean Li

project <- loadArchRProject("/storage/singlecell/jeanl/organoid/data/archr/proj_all")

# DimRed using Harmony
project <- addHarmony(
  ArchRProj = project,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)

# Clustering and UMAP
project <- addClusters(
  input = project,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "ClustersHarm0.8",
  resolution = 0.8
)

project <- addClusters(
  input = project,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "ClustersHarm0.5",
  resolution = 0.5
)

project <- addClusters(
  input = project,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "ClustersHarm0.2",
  resolution = 0.2
)

project <- addUMAP(
  ArchRProj = project, 
  reducedDims = "Harmony", 
  name = "HarmUMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "Sample", embedding = "HarmUMAP")
p2 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "ClustersHarm0.8", embedding = "HarmUMAP")
p3 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "ClustersHarm0.5", embedding = "HarmUMAP")
p4 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "ClustersHarm0.2", embedding = "HarmUMAP")
plotPDF(p1,p2,p3,p4, name = "Plot-UMAP-Sample-Clusters-Harmony.pdf", ArchRProj = project, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(project)
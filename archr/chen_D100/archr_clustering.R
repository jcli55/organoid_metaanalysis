library(ArchR)
library(parallel)

# After creating the ArchR project, this script runs dimentional reduction using Iterative LSI, clustering and embedding with UMAP
# 03/07/25

project <- loadArchRProject("/dfs3b/ruic20_lab/jeancl2/data/archr/chen_D100/")

idxSample <- BiocGenerics::which(project$majorclass %ni% NA)
project <- subsetArchRProject(
  ArchRProj = project,
  cells = project$cellNames[idxSample],
  outputDirectory = "/dfs3b/ruic20_lab/jeancl2/data/archr/chen_D100/",
  dropCells = TRUE,
  force = TRUE
  )

# # Remove doublets
# project <- filterDoublets(project)

# Dim Red using Iterative LSI (not deterministic!)
project <- addIterativeLSI(
  ArchRProj = project,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force = TRUE
)

#project <- addHarmony(
#  ArchRProj = project,
#  reducedDims = "IterativeLSI",
#  name = "Harmony",
#  groupBy = "Sample"
#)

# Clustering and UMAP
project <- addClusters(
  input = project,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters0.8",
  resolution = 0.8,
  force = TRUE
)

project <- addClusters(
  input = project,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters0.5",
  resolution = 0.5,
  force = TRUE
)

project <- addClusters(
  input = project,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters0.1",
  resolution = 0.1,
  force = TRUE
)

project <- addUMAP(
  ArchRProj = project, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)

saveArchRProject(project)

p1 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "maturityclass", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "majorclass", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "Clusters0.8", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "Clusters0.5", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "Clusters0.1", embedding = "UMAP")
plotPDF(p1,p2,p3,p4,p5, name = "Plot-UMAP-Clusters-Chen-D100.pdf", ArchRProj = project, addDOC = FALSE, width = 5, height = 5)

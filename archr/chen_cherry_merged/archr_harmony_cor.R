library(ArchR)
library(parallel)

# After creating the ArchR project, this script runs dimentional reduction using Iterative LSI, clustering and embedding with UMAP
# Uses Harmony to batch correct, assumes already ran archr_clustering (specifically doublet removal and iterative LSI)
# Updated to try different correlation cutoffs
# 4/23/24 - Jean Li

project <- loadArchRProject("/storage/singlecell/jeanl/organoid/data/archr/proj_all")

# DimRed using Harmony
project <- addHarmony(
  ArchRProj = project,
  reducedDims = "IterativeLSI",
  name = "Harmony0.1",
  groupBy = "Sample",
  corCutOff = 0.1
)

project <- addHarmony(
  ArchRProj = project,
  reducedDims = "IterativeLSI",
  name = "Harmony1.0",
  groupBy = "Sample",
  corCutOff = 1.0
)

project <- addUMAP(
  ArchRProj = project, 
  reducedDims = "Harmony0.1", 
  name = "HarmUMAP0.1", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

project <- addUMAP(
  ArchRProj = project, 
  reducedDims = "Harmony1.0", 
  name = "HarmUMAP1.0", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "Sample", embedding = "HarmUMAP0.1")
p2 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "predictedGroup_Un", embedding = "HarmUMAP0.1")
p3 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "Sample", embedding = "HarmUMAP1.0")
p4 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "predictedGroup_Un", embedding = "HarmUMAP1.0")
plotPDF(p1,p2,p3,p4, name = "Plot-UMAP-Harmony-Cutoffs-2.pdf", ArchRProj = project, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(project)

library(ArchR)
library(parallel)

# Runs gene integration using RNA-seq data on the project with both chen and cherry
# 4/9/24 - Jean Li

project <- loadArchRProject("/storage/singlecell/jeanl/organoid/data/archr/proj_all/")
seRNA <- readRDS("/storage/singlecell/jeanl/organoid/data/archr/integrated_rna.rds")

project <- addGeneIntegrationMatrix(
  ArchRProj = project, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  groupRNA = "majorclass",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)

p <- plotEmbedding(
  project, 
  colorBy = "cellColData", 
  name = "predictedGroup_Un"
)

p2 <- plotEmbedding(
  project,
  embedding = "HarmUMAP",
  colorBy = "cellColData",
  name = "predictedGroup_Un"
)
  
plotPDF(p,p2, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = project, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(project)

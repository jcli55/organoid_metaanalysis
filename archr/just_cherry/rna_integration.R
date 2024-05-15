library(ArchR)
library(parallel)

# Runs gene integration using RNA-seq data on the project with both chen and cherry
# 4/9/24 - Jean Li

project <- loadArchRProject("/storage/singlecell/jeanl/organoid/data/archr/cherry_project/")
seRNA <- readRDS("/storage/singlecell/jeanl/organoid/data/archr/integrated_rna.rds")
seRNA <- subset(seRNA, age == 'wk5' | age == 'wk12' | age == 'wk20' | age == 'wk28')

project <- addGeneIntegrationMatrix(
  ArchRProj = project, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = TRUE,
  force = TRUE,
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
  embedding = "UMAP",
  colorBy = "cellColData",
  name = "predictedGroup_Un"
)
  
plotPDF(p,p2, name = "Plot-UMAP-RNA-Integration-Cherry.pdf", ArchRProj = project, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(project)

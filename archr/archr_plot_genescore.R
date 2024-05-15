library(ArchR)
library(parallel)

# Plotting major class marker genescores to see if they cluster out
# 4/3/24 - Jean Li
# Plotting RPC markers to check for NRPC
# 4/25/24

project <- loadArchRProject("/storage/singlecell/jeanl/organoid/data/archr/proj_all/")

#project <- addImputeWeights(project)

#markerGenes <- c(
#  'RHO',
#  'RCVRN',
#  'CRX',
#  'PROM1',
#  'ARR3',
#  'OTX2',
#  'VSX1',
#  'VSX2',
#  'GAD1',
#  'GAD2',
#  'TFAP2A',
#  'TFAP2B',
#  'ONECUT1',
#  'ONECUT2',
#  'NEFM',
#  'RBPMS',
#  'POU4F2',
#  'RLBP1',
#  'SLC1A3',
#  'SFRP2'
#)

markerGenes <- c('SFRP2','NFIA','ASCL1','DLX1','DLX2','ATOH7','RLBP1','CRYM','CLU')

#p <- plotEmbedding(
#  ArchRProj = project, 
#  colorBy = "GeneScoreMatrix", 
#  name = markerGenes, 
#  embedding = "UMAP",
#  imputeWeights = getImputeWeights(project)
#)

p2 <- plotEmbedding(
   ArchRProj = project, 
   colorBy = "GeneScoreMatrix", 
   name = markerGenes, 
   embedding = "HarmUMAP",
   imputeWeights = getImputeWeights(project)
)

#plotPDF(plotList = p, 
#        name = "Plot-UMAP-Marker-Genes.pdf", 
#        ArchRProj = project, 
#        addDOC = FALSE, width = 5, height = 5)

plotPDF(plotList = p2, 
         name = "Plot-RPC-Marker-Genes-Harmony.pdf", 
         ArchRProj = project, 
         addDOC = FALSE, width = 5, height = 5)

#saveArchRProject(project)

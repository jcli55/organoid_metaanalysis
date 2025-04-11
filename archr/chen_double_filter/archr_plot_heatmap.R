library(ArchR)
library(parallel)

# Plotting major class marker genescores to see if they cluster out
# 4/3/24 - Jean Li
# Updated to run on just the Chen data with labels transferred, and added marker heatmap
# 4/5/24 - Jean Li

project <- loadArchRProject("/storage/singlecell/jeanl/organoid/data/archr/chen_double_filter/")

project <- addImputeWeights(project)

markerGenes <- c(
  'RHO',
  'RCVRN',
  'CRX',
  'PROM1',
  'ARR3',
  'OTX2',
  'VSX1',
  'VSX2',
  'GAD1',
  'GAD2',
  'TFAP2A',
  'TFAP2B',
  'ONECUT1',
  'ONECUT2',
  'NEFM',
  'RBPMS',
  'POU4F2',
  'RLBP1',
  'SLC1A3',
  'SFRP2'
)

p <- plotEmbedding(
  ArchRProj = project, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(project)
)

# p2 <- plotEmbedding(
#   ArchRProj = project, 
#   colorBy = "GeneScoreMatrix", 
#   name = markerGenes, 
#   embedding = "HarmUMAP",
#   imputeWeights = getImputeWeights(project)
# )

plotPDF(plotList = p, 
        name = "Plot-UMAP-Marker-Genes.pdf", 
        ArchRProj = project, 
        addDOC = FALSE, width = 5, height = 5)

# plotPDF(plotList = p2, 
#         name = "Plot-UMAP-Marker-Genes-Harmony.pdf", 
#         ArchRProj = project, 
#         addDOC = FALSE, width = 5, height = 5)

markersGS <- getMarkerFeatures(
  ArchRProj = project, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "majorclass",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = project, addDOC = FALSE)

## How I made the mclass markers heatmap like the one I have for RNA-seq, reorder rows/columns
#subset <- markersGS[which(markersGS@elementMetaData$name %in% markerGenes),]
#
#heatmapGS <- markerHeatmap(
#  seMarker = subset,
#  cutOff = "FDR <= 1 & Log2FC >= 0",
#  labelMarkers = markerGenes,
#  transpose = TRUE
#)
#heatmapGS@row_order <- c(#reorder rows)
#heatmapGS@column_order <- c(#reorder col)
#heatmapGS@column_dend_param$reorder <- FALSE
#
#plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap-Mclass", width = 8 , height = 6, ArchRProj = project, addDOC = FALSE)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
write.csv(markerList, '/storage/singlecell/jeanl/organoid/csv/gene_score_degs.csv')

saveArchRProject(project)

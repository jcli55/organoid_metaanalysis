library(ArchR)
library(parallel)

# Motif enrichment

project <- loadArchRProject("/storage/singlecell/jeanl/organoid/data/archr/chen_double_filter/")

# First make sure peaks are called based on the desired cellColData (clusters0.1 in this case)
# Start motif enrichment
project <- addMotifAnnotations(ArchRProj = project, motifSet = "cisbp", name = "Motif", force = TRUE)

# Get the DARs
dar <- getMarkerFeatures(
  ArchRProj = project, 
  useMatrix = "PeakMatrix", 
  groupBy = "majorclass",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# Get the motif enrichment
enrichMotifs <- peakAnnoEnrichment(
  seMarker = dar,
  ArchRProj = project,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = project, addDOC = FALSE)
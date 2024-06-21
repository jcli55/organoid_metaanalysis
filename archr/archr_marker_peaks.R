library(ArchR)
library(parallel)

# Peak calling and getting marker peaks
# 4/5/24 - Jean Li

project <- loadArchRProject("/storage/singlecell/jeanl/organoid/data/archr/chen_double_filter/")

project <- addGroupCoverages(ArchRProj = project, groupBy = "majorclass")

pathToMacs2 <- findMacs2()

project <- addReproduciblePeakSet(
  ArchRProj = project, 
  groupBy = "majorclass", 
  pathToMacs2 = pathToMacs2
)

project <- addPeakMatrix(project)

saveArchRProject(project)

markersPeaks <- getMarkerFeatures(
  ArchRProj = project, 
  useMatrix = "PeakMatrix", 
  groupBy = "majorclass",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

saveArchRProject(project)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerListGR <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)

write.csv(markerList, '/storage/singlecell/jeanl/organoid/csv/dars.csv')
write.csv(markerListGR, '/storage/singlecell/jeanl/organoid/csv/dars_GR.csv')

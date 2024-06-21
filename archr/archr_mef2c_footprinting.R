library(ArchR)
library(parallel)

# Motif footprinting template

project <- loadArchRProject('/storage/singlecell/jeanl/organoid/data/archr/chen_double_filter/')

motifPositions <- getPositions(project)

motifs <- c() # Enter vector of Motifs here
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs

seFoot <- getFootprints(
  ArchRProj = project, 
  positions = motifPositions[markerMotifs], 
  groupBy = "majorclass"
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = project, 
  normMethod = "Subtract",
  plotName = "Plot-Motif_Footprints.pdf",
  addDOC = FALSE,
  smoothWindow = 5
)
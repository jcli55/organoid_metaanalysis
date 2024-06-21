library(ArchR)
library(parallel)

# Motif deviations (13.1)
# Note: if errors occur during addBgdPeaks() or addDeviationsMatrix(), then try deleting the GroupCoverages and PeakCalls directories in the saved project and rerun GroupCoverages and Peak Calling

project <- loadArchRProject('/storage/singlecell/jeanl/organoid/data/archr/chen_double_filter/')

# Make sure Motif annotation was run
if("Motif" %ni% names(project@peakAnnotation)){
  project <- addMotifAnnotations(ArchRProj = project, motifSet = "cisbp", name = "Motif")
}

# Add Motif Deviations
project <- addBgdPeaks(project)

project <- addDeviationsMatrix(
  ArchRProj = project, 
  peakAnnotation = "Motif",
  force = TRUE
)

saveArchRProject(project)

# Can plot the motif deviation scores
plotVarDev <- getVarDeviations(project, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = project, addDOC = FALSE)
library(ArchR)
library(parallel)

# Computes and optionally plots co-accessibility and peak2gene linkage

project <- loadArchRProject('/storage/singlecell/jeanl/organoid/data/archr/chen_double_filter/')

# Get marker genes
markerGenes  <- c() # Enter genes of interest if want to plot browser tracks

# Co-accessibility
project <- addCoAccessibility(
  ArchRProj = project,
  reducedDims = "IterativeLSI"
)

cA <- getCoAccessibility(
    ArchRProj = project,
    corCutOff = 0.5,
    resolution = 1,
    returnLoops = FALSE
)
write.csv(cA, '/storage/singlecell/jeanl/organoid/csv/coaccessibility.csv')

cA.GR <- getCoAccessibility(
    ArchRProj = project,
    corCutOff = 0.5,
    resolution = 1,
    returnLoops = TRUE
)
write.csv(cA.GR, '/storage/singlecell/jeanl/organoid/csv/coaccessibility_GR.csv')

# Peaks-to-Genes
project <- addPeak2GeneLinks(
  ArchRProj = project,
  reducedDims = "IterativeLSI"
)

p2g <- getPeak2GeneLinks(
    ArchRProj = project,
    corCutOff = 0.45,
    resolution = 1,
    returnLoops = FALSE
)
write.csv(p2g, '/storage/singlecell/jeanl/organoid/csv/peak2genes.csv')

p2g.GR <- getPeak2GeneLinks(
    ArchRProj = project,
    corCutOff = 0.45,
    resolution = 1,
    returnLoops = TRUE
)
write.csv(p2g.GR, '/storage/singlecell/jeanl/organoid/csv/peak2genes_GR.csv')

# If marker genes are provided, plot browser tracks
if (len(markerGenes) > 0) {
  p <- plotBrowserTrack(
    ArchRProj = project, 
    groupBy = "majorclass", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = getCoAccessibility(project, resolution = 1000)
  )
  plotPDF(plotList = p, 
          name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", 
          ArchRProj = project, 
          addDOC = FALSE, width = 5, height = 5)

  p <- plotBrowserTrack(
    ArchRProj = project, 
    groupBy = "majorclass", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = getPeak2GeneLinks(project, resolution = 1000)
  )
  plotPDF(plotList = p, 
          name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.pdf", 
          ArchRProj = project, 
          addDOC = FALSE, width = 5, height = 5)
}

#saveArchRProject(project)
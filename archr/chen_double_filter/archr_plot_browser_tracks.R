library(ArchR)
library(parallel)

proj <- loadArchRProject('/dfs3b/ruic20_lab/jeancl2/data/archr/chen_double_filter/')
genes <- read.csv('/dfs3b/ruic20_lab/jeancl2/data/enhancer_assay_candidate.csv')

# Any necessary subsetting
# idxD100 <- BiocGenerics::which(proj$age %in% 'D100')
# D100 <- subsetArchRProject(
#   ArchRProj = project,
#   cells = project$cellNames[idxD100],
#   outputDirectory = "/dfs3b/ruic20_lab/jeancl2/data/archr/chen_D100_double_filter/",
#   dropCells = TRUE,
#   force = TRUE
#   )

idxD100 <- BiocGenerics::which(proj$age %in% "D100")
cellsD100 <- proj$cellNames[idxD100]
idxD183 <- BiocGenerics::which(proj$age %in% "D183")
cellsD183 <- proj$cellNames[idxD183]

p100 <- plotBrowserTrack(
    ArchRProj = proj[cellsD100, ], 
    groupBy = "maturityclass", 
    geneSymbol = genes$Gene, 
    upstream = 50000,
    downstream = 50000,
)

plotPDF(plotList = p100, 
    name = "Plot-Browser-Tracks-D100.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)

p183 <- plotBrowserTrack(
    ArchRProj = proj[cellsD183, ], 
    groupBy = "maturityclass", 
    geneSymbol = genes$Gene, 
    upstream = 50000,
    downstream = 50000,
)

plotPDF(plotList = p183, 
    name = "Plot-Browser-Tracks-D183.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)
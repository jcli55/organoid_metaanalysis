library(ArchR)
library(parallel)

# Impute weights and peak calling
# 03/07/25

project <- loadArchRProject("/dfs3b/ruic20_lab/jeancl2/data/archr/chen_D183/")

project <- addImputeWeights(project)

project <- addGroupCoverages(ArchRProj = project, groupBy = "majorclass")

pathToMacs2 <- findMacs2()
project <- addReproduciblePeakSet(
    ArchRProj = project, 
    groupBy = "majorclass", 
    pathToMacs2 = pathToMacs2
)

project <- addPeakMatrix(project)

saveArchRProject(project)
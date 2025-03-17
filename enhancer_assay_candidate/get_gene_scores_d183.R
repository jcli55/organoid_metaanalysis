library(ArchR)
library(parallel)

project <- loadArchRProject('/dfs3b/ruic20_lab/jeancl2/data/archr/chen_D183/')
gs <- getMatrixFromProject(project, 'GeneScoreMatrix')
enhancers <- read.csv("/dfs3b/ruic20_lab/jeancl2/data/enhancer_assay_candidate/enhancer_assay_candidate.csv", header = 1)

ids <- c()
meanGS <- c()
geneName <- c()

for (i in enhancers$Gene) {
    id <- which(gs@elementMetadata$name %in% i)
    ids <- append(ids, id)
    meanGS <- append(meanGS, mean(gs@assays@data$GeneScoreMatrix[id,]))
    geneName <- append(geneName, gs@elementMetadata[id, 'name'])
}

df <- data.frame(id = ids, Gene = geneName, MeanGeneScore = meanGS)
write.csv(df, '/dfs3b/ruic20_lab/jeancl2/data/enhancer_assay_candidate/enhancer_genescore_d183.csv')

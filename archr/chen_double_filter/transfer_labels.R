library(ArchR)
library(parallel)

project <- loadArchRProject('/dfs3b/ruic20_lab/jeancl2/data/archr/chen_double_filter/')
metadata <- read.csv('/dfs3b/ruic20_lab/jeancl2/data/csv/multivelo_metadata.csv')
metadata <- data.frame(metadata, row.names = 'id')

# Reformats the id's from the RNA object metadata to ArchR's id's
# Remove the '-0' tag
rownames(metadata) <- substr(rownames(metadata), 1, nchar(rownames(metadata))-2)
# Adds a '#' between the sample name and the barcode
rownames(metadata <- paste(substr(rownames(metadata), 1, nchar(rownames(metadata))-19), substr(rownames(metadata), nchar(rownames(metadata))-17, nchar(rownames(metadata))), sep="#")

majorclass <- c()
age <- c()
maturityclass <- c()
pseudotime <- c()
latenttime <- c()

for (i in project$cellNames) {
  majorclass <- append(majorclass, metadata[i, 'majorclass'])
  age <- append(age, metadata[i, 'age'])
  maturityclass <- append(maturityclass, metadata[i, 'maturityclass'])
  pseudotime <- append(pseudotime, metadata[i, 'velo_s_norm_pseudotime'])
  latenttime <- append(latenttime, metadata[i, 'latent_time'])
}

project$majorclass <- majorclass
project$age <- age
project$maturityclass <- maturityclass
project$pseudotime <- pseudotime
project$latenttime <- latenttime

saveArchRProject(project, '/dfs3b/ruic20_lab/jeancl2/data/archr/chen_double_filter/')
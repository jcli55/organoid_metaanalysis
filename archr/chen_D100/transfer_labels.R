library(ArchR)
library(parallel)

project <- loadArchRProject('/dfs3b/ruic20_lab/jeancl2/data/archr/chen_D100/')
metadata <- read.csv('/dfs3b/ruic20_lab/jeancl2/data/csv/multivelo_metadata.csv')
metadata <- data.frame(metadata, row.names = 'id')

# Reformats the id's from the RNA object metadata to ArchR's id's
# Remove the '-0' tag
rownames(metadata) <- substr(rownames(metadata), 1, nchar(rownames(metadata))-2)
# Adds a '#' between the sample name and the barcode
names <- paste(substr(rownames(metadata), 1, nchar(rownames(metadata))-19), substr(rownames(metadata), nchar(rownames(metadata))-17, nchar(rownames(metadata))), sep="#")
rownames(metadata) <- names

majorclass <- c()
maturityclass <- c()
age <- c()

for (i in project$cellNames) {
  majorclass <- append(majorclass, metadata[i, 'class'])
  maturityclass <- append(maturityclass, metadata[i, 'maturityclass'])
  age <- append(age, metadata[i, 'age'])
}

project$majorclass <- majorclass
project$maturityclass <- maturityclass
project$age <- age

saveArchRProject(project, '/dfs3b/ruic20_lab/jeancl2/data/archr/chen_D100/')
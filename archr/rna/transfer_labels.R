library(Seurat)

data <- readRDS('/storage/singlecell/jeanl/organoid/data/archr/integrated_rna.rds')
metadata <- read.csv('/storage/singlecell/jeanl/organoid/csv/metadata_with_id.csv')
metadata <- data.frame(metadata, row.names = 'id')

intersect <- intersect(rownames(metadata), Cells(data))
data <- subset(data, cells = intersect)

majorclass <- c()
age <- c()

for (i in Cells(data)) {
  majorclass <- append(majorclass, metadata[i, 'majorclass'])
  age <- append(age, metadata[i, 'age'])
}

data$majorclass <- majorclass
data$age <- age

saveRDS(data, '/storage/singlecell/jeanl/organoid/data/archr/integrated_rna.rds')

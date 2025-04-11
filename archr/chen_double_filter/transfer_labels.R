library(ArchR)
library(parallel)

project <- loadArchRProject('/storage/singlecell/jeanl/organoid/data/archr/chen_double_filter/')
metadata <- read.csv('/storage/singlecell/jeanl/organoid/csv/metadata_with_id.csv')
metadata <- data.frame(metadata, row.names = 'id')

# Reformats the id's from the RNA object metadata to ArchR's id's
# Remove the '-0' tag
rownames(metadata) <- substr(rownames(metadata), 1, nchar(rownames(metadata))-2)
# Adds a '#' between the sample name and the barcode
names <- paste(substr(rownames(metadata), 1, nchar(rownames(metadata))-19), substr(rownames(metadata), nchar(rownames(metadata))-17, nchar(rownames(metadata))), sep="#")

majorclass <- c()
age <- c()

for (i in project$sampleid) {
  majorclass <- append(majorclass, metadata[i, 'majorclass'])
  age <- append(age, metadata[i, 'age'])
}

project$majorclass <- majorclass
project$age <- age

saveArchRProject('/storage/singlecell/jeanl/organoid/data/archr/chen_double_filter/')
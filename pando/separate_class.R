library(Seurat)
library(Signac)
library(tidyverse)

# Code that generated the merged.rds file
# counts <- readRDS('/storage/singlecell/jeanl/organoid/data/Chen/rds/rna.rds') #Generated from merge_rnaseq.R
# peaks <- readRDS('/storage/singlecell/jeanl/organoid/data/Chen/rds/atac.rds') #Generated from merge_atacseq.R
# counts[['peaks']] <- peaks[['peaks']]
# saveRDS('/storage/singlecell/jeanl/organoid/data/Chen/rds/merged.rds')

# Load in the merged seurat object and remane the cellnames to match the annotated adata cell ids
data <- readRDS('/storage/singlecell/jeanl/organoid/data/Chen/rds/merged.rds')
data <- RenameCells(data, new.names = paste(colnames(data), '-0', sep=''))

# Load in the metadata that contains the cell type annotations
metadata <- read.csv('/storage/singlecell/jeanl/organoid/csv/metadata_clean.csv')

cat("Finished reading in data!")
cat("\n")

# Set the majorclasses to separate
majorclass <- c('AC','BC','HC','RGC','Cone','Rod','PRPC','NRPC','MG')

# Subset out each majorclass from the rds and save as individual rds files
for (mclass in majorclass) {
  cat('Processing ', mclass)
  cat('\n')
  barcodes <- subset(metadata, source == 'chen' & sampletype == 'organoid' & majorclass == mclass)$X
  
  subset <- subset(data, cells = barcodes)
  
  saveRDS(subset, paste0('/storage/singlecell/jeanl/organoid/data/Chen/rds/majorclass/', mclass, '.rds', sep=''))
  cat('Saved', mclass)
  cat('\n')
}

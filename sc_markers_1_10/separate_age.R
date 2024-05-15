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

# Load in the metadata
metadata <- read.csv('/storage/singlecell/jeanl/organoid/csv/metadata.csv')

cat("Finished reading in data!")
cat("\n")

# Set the ages to separate
age <- c('D35','D47','D75','D100','D123','D133','D183','D206','D243')

# Subset out each age from the rds and save as individual rds files
for (i in age) {
  cat('Processing ', i)
  cat('\n')
  barcodes <- subset(metadata, source == 'chen' & sampletype == 'organoid' & age == i)$X
  
  subset <- subset(data, cells = barcodes)
  
  saveRDS(subset, paste0('/storage/singlecell/jeanl/organoid/data/Chen/rds/age/', i, '.rds', sep=''))
  cat('Saved', i)
  cat('\n')
}

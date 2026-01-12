library(Seurat)
library(Signac)
library(tidyverse)

# Code that generated the merged.rds file
# counts <- readRDS('/storage/singlecell/jeanl/organoid/data/Chen/rds/rna.rds') #Generated from merge_rnaseq.R
# peaks <- readRDS('/storage/singlecell/jeanl/organoid/data/Chen/rds/atac.rds') #Generated from merge_atacseq.R
# counts[['peaks']] <- peaks[['peaks']]
# saveRDS('/storage/singlecell/jeanl/organoid/data/Chen/rds/merged.rds')

# Adapted to hpc3 at UCI

# Load in the merged seurat object and remane the cellnames to match the annotated adata cell ids
data <- readRDS('/dfs3b/ruic20_lab/singlecell/jeanl/organoid/data/Chen/rds/merged.rds')
data <- RenameCells(data, new.names = paste(colnames(data), '-0', sep=''))

# Load in the metadata that contains the cell type annotations
metadata <- read.csv('/dfs3b/ruic20_lab/jeancl2/data/chen_ro_clean_metadata.csv')

cat("Finished reading in data!")
cat("\n")

# Set the majorclasses to separate
class <- c('AC Precursor','BC Precursor','HC Precursor','RGC Precursor','Cone Precursor','Rod Precursor')

# Subset out each majorclass from the rds and save as individual rds files
for (mclass in class) {
  cat('Processing ', mclass)
  cat('\n')
  barcodes <- subset(metadata, source == 'chen' & sampletype == 'organoid' & subclass == mclass)$X
  
  subset <- subset(data, cells = barcodes)
  
  saveRDS(subset, paste0('/dfs3b/ruic20_lab/jeancl2/data/pando/', sub(" ", "_", mclass), '.rds', sep=''))
  cat('Saved', mclass)
  cat('\n')
}

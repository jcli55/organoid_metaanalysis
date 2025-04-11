library(Seurat)
library(Signac)

# Load the RNA seurat object and make cell names consistent with adata
project <- readRDS("/dfs3b/ruic20_lab/singlecell/jeanl/organoid/data/Chen/rds/rna.rds")
project <- RenameCells(project, new.names = paste(colnames(project), '-0', sep=''))

# Load in the metadata and subset to just Chen organoid PRPCs
data <- read.csv("/dfs3b/ruic20_lab/jeancl2/data/csv/metadata_with_id.csv", header=1)
data <- subset(data, sampletype == 'organoid' & source == 'chen' & majorclass == 'PRPC')
project <- subset(project, cells = data$id)

# Normalize the data and save
project <- NormalizeData(project)
project <- FindVariableFeatures(project)
project <- ScaleData(project)

saveRDS(project, "/dfs3b/ruic20_lab/jeancl2/data/celltran/PRPC.rds")
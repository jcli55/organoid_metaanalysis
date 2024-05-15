library(Seurat)

chen <- readRDS('/storage/singlecell/jeanl/organoid/data/Chen/rds/rna.rds')
cherry <- readRDS('/storage/singlecell/jeanl/organoid/data/Cherry/merged_rna.rds')

merge <- merge(chen, cherry)

merge <- RenameCells(merge, new.names=paste0(Cells(merge), '-0'))

saveRDS(merge, '/storage/singlecell/jeanl/organoid/data/archr/integrated_rna.rds')

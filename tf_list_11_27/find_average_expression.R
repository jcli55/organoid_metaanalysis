library(Seurat)
library(Signac)

# This script runs AverageExpression() on my organoid data by majorclass
# Used for the TF task 11/27/23

majorclass <- c('AC','BC','HC','RGC','Rod','Cone','MG','PRPC','NRPC')
gene.list <- read.csv('/storage/singlecell/jeanl/organoid/csv/tf_list.csv')
gene.list <- as.vector(gene.list$genes)

for (mclass in majorclass) {
  data <-  readRDS(paste0('/storage/singlecell/jeanl/organoid/data/Chen/rds/majorclass/', mclass, '.rds'))
  
  DefaultAssay(data) <- 'RNA'
  data <- NormalizeData(data)
  all.genes <- rownames(data)
  data <- ScaleData(data, features = all.genes)
  
  df <- AverageExpression(data, 'RNA', gene.list)
  write.csv(df, paste0('/storage/singlecell/jeanl/organoid/csv/average_expression/', mclass, 'avg_exp.csv'))
}

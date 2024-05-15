library(Seurat)
library(Signac)

# This script runs AverageExpression() on my organoid data by majorclass
# Used for the TF task 11/27/23

gene.list <- read.csv('/storage/singlecell/jeanl/organoid/csv/tf_list.csv')
gene.list <- as.vector(gene.list$genes)

data <-  readRDS(paste0('/storage/singlecell/jeanl/organoid/data/Chen/rds/majorclass/chen_all.rds'))
DefaultAssay(data) <- 'RNA'
data <- NormalizeData(data)
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)
df <- AverageExpression(data, 'RNA', gene.list)
write.csv(df, paste0('/storage/singlecell/jeanl/organoid/csv/average_expression/chen_all_avg_exp.csv'))
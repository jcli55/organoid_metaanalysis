library(Seurat)
library(Signac)

# This script runs AverageExpression() on my organoid data by age
# Used for the TF task 11/27/23

age <- c('D35','D47','D75','D100','D123','D133','D183','D206','D243')
gene.list <- c('MITF','PMEL','TYRP1','TFPI2','SLC6A15','RAB38','GJA1','PCDH7','RELN','MECOM','CPAMD8','COL9A1','ZIC1','TGFB2')

for (i in age) {
  data <-  readRDS(paste0('/storage/singlecell/jeanl/organoid/data/Chen/rds/age/', i, '.rds'))
  
  DefaultAssay(data) <- 'RNA'
  data <- NormalizeData(data)
  all.genes <- rownames(data)
  data <- ScaleData(data, features = all.genes)
  
  df <- AverageExpression(data, 'RNA', gene.list)
  write.csv(df, paste0('/storage/singlecell/jeanl/organoid/csv/sc_markers/', i, 'avg_exp.csv'))
}

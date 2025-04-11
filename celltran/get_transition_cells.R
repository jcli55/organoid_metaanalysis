library(Seurat)
library(CellTran)
library(mclust)

# Run CellTran on the PRPCs

# Load and prepare data
data <- readRDS("/dfs3b/ruic20_lab/jeancl2/data/celltran/PRPC.rds")
data <- RunPCA(data, features = rownames(data))
data <- FindNeighbors(data, dims = 1:5)
data <- FindClusters(data, resolution = 0.5)
data <- RunUMAP(data, dims=1:5)

# Calculate transition index (default params)
data <- transition_index(data)
write.csv(data$transition_index, "/dfs3b/ruic20_lab/jeancl2/data/celltran/transition_index.csv")
#FeaturePlot(data, features='transition_index', order=T, pt.size=1)
#plot(density(data$transition_index), main="Transition index distribution")

# Filter out NA values
data <- subset(data, transition_index > 0)
write.csv(data$transition_index, "/dfs3b/ruic20_lab/jeancl2/data/celltran/transition_index_na_rm.csv")

saveRDS(data, "/dfs3b/ruic20_lab/jeancl2/data/celltran/PRPC_celltran_result.rds")

# Separating cells into transiton cells and stable cells using GMM (Gaussian Mixture Model)
fit <- Mclust(data$transition_index, G=2)
transition_cells <- colnames(data)[fit$classification==2 & fit$z[,2]>0.8]
#DimPlot(data,cells.highlight = list('transition_cells'=transition_cells))

write.csv(transition_cells, "/dfs3b/ruic20_lab/jeancl2/data/celltran/transition_cells.csv")

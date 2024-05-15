library(Seurat)
library(Signac)

# Zhen's script to downsample cells in a seurat object
# He already ran this code and the rds files are made in his directory - separates the merged object by majorclass

n_features = 5000

for (x in c("RGC", "HC", "NRPC", "PRPC", "Cone", "MG", "AC", "BC", "Rod")) {
    rna <- readRDS("/storage/chentemp/zz4/adult_dev_compare/results/merged_rna/merged_rna.rds")
    meta <- read.csv("/storage/chentemp/zz4/adult_dev_compare/results/Annotation/merged_raw_filtered_umap_10000_wadult_annotated.csv")
    rownames(meta) <- meta$X
    meta <- meta[meta$majorclass==x,]
    common_cells <- intersect(colnames(rna), rownames(meta))
    rna <- subset(rna, cells = common_cells)

    atac <- readRDS("/storage/chentemp/zz4/adult_dev_compare/results/merged_atac/atac.rds")

    common_cells <- intersect(colnames(rna), colnames(atac))

    rna <- subset(rna, cells = common_cells)
    atac <- subset(atac, cells = common_cells)

    rna <- NormalizeData(rna)
    rna <- ScaleData(rna)
    rna <- FindVariableFeatures(rna, nfeatures = n_features)
    rna <- subset(rna, features = VariableFeatures(rna))

    atac <- RunTFIDF(atac)
    atac <- FindTopFeatures(atac, min.cutoff = "q5")

    atac[["RNA"]] <- rna@assays$RNA

    
    meta <- meta[common_cells, ]
    atac@meta.data <- cbind(atac@meta.data, meta[colnames(atac), ])

    saveRDS(atac, paste("/storage/chentemp/zz4/adult_dev_compare/results/pando_merged_seurat_object/pando_merged_seurat_object_",
        x, ".rds", sep = ""))
}

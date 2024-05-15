# Script by Zhen adapted for organoid data by Jean

input_path<-"/storage/singlecell/zz4/multi_organoid/data/"
output_file<-"/storage/singlecell/jeanl/organoid/data/Chen/rds/rna.rds"
# Follow vignette at
# https://satijalab.org/seurat/articles/integration_introduction.html
# And https://github.com/satijalab/seurat/issues/1720

# Import packages needed for the following analysis
suppressMessages(library(Seurat))
set.seed(0)
samples <- c("Multi_H9_D35", 
             "Multi_H9_D47", 
             "Multi_organoid_NRL_D075",
             "Multi_H9_D100",
             "Multi_NRL_GFP_D123",
             "Multi_organoid_D133",
             "Multi_Organoid_D183",
             "Multi_organoid_NRL_D206",
             "Multi_organoid_NRL_D243")

for (sam in samples) {
  # This for loop will read count matrix into Seurat object and
  # Assign Seurat object to each one of sample names.
  cat("\n")
  cat("Processing sample", sam, "now!")
  cat("\n")
  counts <- Read10X_h5(paste0(input_path, sam, "/outs/filtered_feature_bc_matrix.h5",
                              sep = ""))
  # Create a Seurat object containing the RNA data
  if (class(counts)=="list") {
    temp <- CreateSeuratObject(counts = counts$`Gene Expression`,
                               assay = "RNA")
  }
  if (class(counts)=="dgCMatrix") {
    temp <- CreateSeuratObject(counts = counts,
                               assay = "RNA")
  }
  DefaultAssay(temp) <- "RNA"
  # Added this if statement to lowercase the o in organoid to match the adata cell ids - Jean
  if (sam == "Multi_Organoid_D183") {
    temp <- RenameCells(object = temp, add.cell.id = "Multi_organoid_D183")
  } else {
    temp <- RenameCells(object = temp, add.cell.id = sam)
  }
  assign(sam, temp)
  cat("Processing sample", sam, "finished!")
  cat("\n")
  cat("\n")
}

seurat_object <- merge(Multi_H9_D35, y = c(Multi_H9_D47, 
                                           Multi_organoid_NRL_D075, 
                                           Multi_H9_D100, 
                                           Multi_NRL_GFP_D123, 
                                           Multi_organoid_D133, 
                                           Multi_Organoid_D183, 
                                           Multi_organoid_NRL_D206, 
                                           Multi_organoid_NRL_D243))

saveRDS(seurat_object, output_file)
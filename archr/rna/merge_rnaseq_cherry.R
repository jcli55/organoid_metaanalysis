# Script by Zhen adapted for organoid data by Jean

input_path<-'/storage/chentemp/u250758/mrrdir/organoid_metaanalysis/Cherry/'
output_file<-"/storage/singlecell/jeanl/organoid/data/Cherry/merged_rna.rds"
# Follow vignette at
# https://satijalab.org/seurat/articles/integration_introduction.html
# And https://github.com/satijalab/seurat/issues/1720

# Import packages needed for the following analysis
suppressMessages(library(Seurat))
set.seed(0)
samples <- c('wt5-1',
             'wt5-2',
             'wt12-1',
             'wt12-2',
             'wt12-3',
             'wt20-1',
             'wt20-2',
             'wt20-3',
             'wt28-1',
             'wt28-2',
             'wt28-3')
objects <- gsub('-','.',samples)

for (i in 1:length(samples)) {
  sam <- samples[i]
  obj <- objects[i]

  # This for loop will read count matrix into Seurat object and
  # Assign Seurat object to each one of sample names.
  cat("\n")
  cat("Processing sample", sam, "now!")
  cat("\n")
  counts <- Read10X_h5(paste0(input_path, sam, "/cellranger/outs/filtered_feature_bc_matrix.h5",
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
  #DefaultAssay(temp) <- "RNA"
  # Added this if statement to lowercase the o in organoid to match the adata cell ids - Jean
  #if (sam == "Multi_Organoid_D183") {
  #  temp <- RenameCells(object = temp, add.cell.id = "Multi_organoid_D183")
  #} else {
  #  temp <- RenameCells(object = temp, Cells(te)
  #}
  temp <- RenameCells(object = temp, add.cell.id = sam)

  assign(obj, temp)
  cat("Processing sample", sam, "finished!")
  cat("\n")
  cat("\n")
}

seurat_object <- merge(wt5.1, y = c(wt5.2,
                                    wt12.1,
                                    wt12.2,
                                    wt12.3,
                                    wt20.1,
                                    wt20.2,
                                    wt20.3,
                                    wt28.1,
                                    wt28.2,
                                    wt28.3))

saveRDS(seurat_object, output_file)

library(tidyverse)
library(Pando)
library(Seurat)
library(Signac)
library(doParallel)
library(BSgenome.Hsapiens.UCSC.hg38)

data('phastConsElements20Mammals.UCSC.hg38')
data('motifs')
data('motif2tf')

# This script runs the basic Pando workflow
# https://quadbio.github.io/Pando/articles/getting_started.html

# For loop to run the workflow through each major cell class
#majorclass <- c("RGC", "HC", "Cone", "AC", "BC", "Rod", "MG", "PRPC", "NRPC")
#majorclass <- c("PRPC", "NRPC") #Too little gene targets in Rod so execution halted, running on just progenitors and am going to go through Rod interactively
majorclass <- c('Rod')
for (mclass in majorclass) {
  # Load in the data and normalize it (Normalizing steps from these vignettes:)
  # https://satijalab.org/seurat/articles/pbmc3k_tutorial
  # https://stuartlab.org/signac/articles/pbmc_vignette 
  data <- c()
  #data <- readRDS(paste0("/storage/chentemp/zz4/adult_dev_compare/results/pando_merged_seurat_object/pando_merged_seurat_object_", mclass, ".rds"))
  data <- readRDS(paste0('/storage/singlecell/jeanl/organoid/data/Chen/rds/old_annotation/fetal_grns/', mclass, '_fet_grn.rds'))
  
  #DefaultAssay(data) <- 'RNA'
  #data <- NormalizeData(data) 
  #data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
  #all.genes <- rownames(data) 
  #data <- ScaleData(data, features = all.genes)
  
  #DefaultAssay(data) <- 'peaks'
  #data <- RunTFIDF(data)
  #data <- FindTopFeatures(data, min.cutoff = 'q0')
  #data <- RunSVD(data)
  
  # Start the Pando workflow ---------------------------------------------------
  # Initiate GRN
  #data <- initiate_grn(
  #  data,
  #  rna_assay = 'RNA',
  #  peak_assay = 'peaks',
  #  regions = phastConsElements20Mammals.UCSC.hg38 
  #)
  
  # Scanning for motifs
  #patterning_genes <- read_tsv('patterning_genes.tsv') 
    #the vignette uses patterning genes only but I use var.features (I don't have the tsv)
  
  # pattern_tfs <- patterning_genes %>%
  #   filter(type=='Transcription factor') %>%
  #   pull(symbol)
  # motif2tf_use <- motif2tf %>%
  #   filter(tf %in% pattern_tfs)
  # motifs_use <- motifs[unique(motif2tf_use$motif)]
  
  #data <- find_motifs(
  #  data, 
  #  pfm = motifs, 
  #  motif_tfs = motif2tf,
  #  genome = BSgenome.Hsapiens.UCSC.hg38
  #)
  
  # Inferring the GRN
  registerDoParallel(4)
  data <- infer_grn(
    data,
    peak_to_gene_method = 'GREAT',
    genes = data@assays$RNA@var.features,
    parallel = T
  )

  # Module discovery
  
  if (mclass == "Rod" | mclass == "MG") {
   data <- find_modules( #these are more lenient params used by the vignette, I had to use for rod and mg to match ro
     data,
     p_thresh = 0.1,
     nvar_thresh = 2,
     min_genes_per_module = 1,
     rsq_thresh = 0.05
   )
  } else {
   data <- find_modules(data) #default params
  }
  
  # Some other plots from the vignette
  #plot_gof(data, point_size=3)
  #plot_module_metrics(data)
  
  write.csv(NetworkModules(data)@meta, paste0('/storage/singlecell/jeanl/organoid/data/Chen/rds/fetal_grns/', mclass, '_fet_grn.csv'))
  saveRDS(data, paste0('/storage/singlecell/jeanl/organoid/data/Chen/rds/fetal_grns/', mclass, '_fet_grn.rds'))

  # Visualizing the GRN (The png + dev.off doesn't work in scripts)
  data <- get_network_graph(data)
  png(filename = paste0('/storage/singlecell/jeanl/organoid/figures/grns/', mclass, '_fet_grn.png'))
  plot_network_graph(data)
  dev.off()
}

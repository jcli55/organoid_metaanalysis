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
majorclass <- c('AC', 'BC', 'HC', 'RGC', 'Cone', 'Rod', 'PRPC', 'NRPC', 'MG')
for (mclass in majorclass) {
  # Load in the data and normalize it (Normalizing steps from these vignettes:)
  # https://satijalab.org/seurat/articles/pbmc3k_tutorial
  # https://stuartlab.org/signac/articles/pbmc_vignette 
  data <- c()
  data <- readRDS(paste0('/storage/singlecell/jeanl/organoid/data/Chen/rds/majorclass/', 
                         mclass, '.rds'))
  
  DefaultAssay(data) <- 'RNA'
  data <- NormalizeData(data) 
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(data) 
  data <- ScaleData(data, features = all.genes)
  
  DefaultAssay(data) <- 'peaks'
  data <- RunTFIDF(data)
  data <- FindTopFeatures(data, min.cutoff = 'q0')
  data <- RunSVD(data)
  
  # Start the Pando workflow ---------------------------------------------------
  # Initiate GRN
  data <- initiate_grn(
    data,
    rna_assay = 'RNA',
    peak_assay = 'peaks',
    regions = phastConsElements20Mammals.UCSC.hg38 
  )
  
  # Scanning for motifs
  #patterning_genes <- read_tsv('patterning_genes.tsv') 
    #the vignette uses patterning genes only but I use var.features (I don't have the tsv)
  
  # pattern_tfs <- patterning_genes %>%
  #   filter(type=='Transcription factor') %>%
  #   pull(symbol)
  # motif2tf_use <- motif2tf %>%
  #   filter(tf %in% pattern_tfs)
  # motifs_use <- motifs[unique(motif2tf_use$motif)]
  
  data <- find_motifs(
    data, 
    pfm = motifs, 
    motif_tfs = motif2tf,
    genome = BSgenome.Hsapiens.UCSC.hg38
  )
  
  # Inferring the GRN
  registerDoParallel(4)
  data <- infer_grn(
    data,
    peak_to_gene_method = 'GREAT',
    genes = data@assays$RNA@var.features,
    parallel = T
  )

  # Module discovery
  # data <- find_modules(data) #default params
  #  data <- find_modules( #these are more lenient params used by the vignette
  #    data, 
  #    p_thresh = 0.1,
  #    nvar_thresh = 2, 
  #    min_genes_per_module = 1, 
  #    rsq_thresh = 0.05
  #  )
  #  data <- find_modules( #these are even more lenient params I had to use for Rod ro
  #    data, 
  #    p_thresh = 0.2,
  #    nvar_thresh = 2, 
  #    min_genes_per_module = 1, 
  #    rsq_thresh = 0.025
  #  )
  
  # Some other plots from the vignette
  #plot_gof(data, point_size=3)
  #plot_module_metrics(data)
  
  #write.csv(NetworkModules(data)@meta, paste0('/storage/singlecell/jeanl/organoid/csv/reannotation/grns/', mclass, '_grn.csv'))
 
  # Visualizing the GRN (The png + dev.off doesn't work in scripts)
  #data <- get_network_graph(data)
  #png(filename = paste0('/storage/singlecell/jeanl/organoid/figures/grns/', mclass, '_grn.png'))
  #plot_network_graph(data)
  #dev.off()
  
  saveRDS(data, paste0('/storage/singlecell/jeanl/organoid/data/Chen/rds/grns/', mclass, '_grn.rds'))
}

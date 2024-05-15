library(ArchR)
library(parallel)
addArchRThreads(threads = 8)
addArchRGenome('hg38')

#Script on how I set up the archr project for organoid
#Updated to add Cherry's data
#Updated to just Cherry's data

# Setup the input files, sample names, and output files
outputDir <- "/storage/singlecell/jeanl/organoid/data/archr/cherry_project"

inputFiles <- c("/storage/chentemp/u250758/mrrdir/organoid_metaanalysis/Cherry/atac/wt5-1/cellranger-atac/outs/fragments.tsv.gz",
		"/storage/chentemp/u250758/mrrdir/organoid_metaanalysis/Cherry/atac/wt5-2/cellranger-atac/outs/fragments.tsv.gz",
		"/storage/chentemp/u250758/mrrdir/organoid_metaanalysis/Cherry/atac/wt12-1/cellranger-atac/outs/fragments.tsv.gz",
		"/storage/chentemp/u250758/mrrdir/organoid_metaanalysis/Cherry/atac/wt12-2/cellranger-atac/outs/fragments.tsv.gz",
		"/storage/chentemp/u250758/mrrdir/organoid_metaanalysis/Cherry/atac/wt12-3/cellranger-atac/outs/fragments.tsv.gz",
		"/storage/chentemp/u250758/mrrdir/organoid_metaanalysis/Cherry/atac/wt20-1/cellranger-atac/outs/fragments.tsv.gz",
		"/storage/chentemp/u250758/mrrdir/organoid_metaanalysis/Cherry/atac/wt20-2/cellranger-atac/outs/fragments.tsv.gz",
		"/storage/chentemp/u250758/mrrdir/organoid_metaanalysis/Cherry/atac/wt28-1/cellranger-atac/outs/fragments.tsv.gz",
		"/storage/chentemp/u250758/mrrdir/organoid_metaanalysis/Cherry/atac/wt28-2/cellranger-atac/outs/fragments.tsv.gz"
    )

names <- c("wt5-1",
	   "wt5-2",
	   "wt12-1",
           "wt12-2",
           "wt12-3",
           "wt20-1",
           "wt20-2",
           "wt28-1",
           "wt28-2"
           )

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names,
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

#doubScores <- addDoubletScores(
#    input = ArrowFiles,
#    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
#    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
#    LSIMethod = 1
#)

project <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = outputDir,
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
  
# Save the project
saveArchRProject(ArchRProj = project, outputDirectory = outputDir)

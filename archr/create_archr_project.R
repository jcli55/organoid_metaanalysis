library(ArchR)
library(parallel)
addArchRThreads(threads = 8)
addArchRGenome('hg38')

#Script on how I set up the archr project for organoid
#Updated to add Cherry's data

# Setup the input files, sample names, and output files
outputDir <- "/storage/singlecell/jeanl/organoid/data/archr/proj_all"

inputFiles <- c("/storage/singlecell/zz4/multi_organoid/data/Multi_H9_D35/outs/atac_fragments.tsv.gz", 
                "/storage/singlecell/zz4/multi_organoid/data/Multi_H9_D47/outs/atac_fragments.tsv.gz", 
                "/storage/singlecell/zz4/multi_organoid/data/Multi_organoid_NRL_D075/outs/atac_fragments.tsv.gz",
                "/storage/singlecell/zz4/multi_organoid/data/Multi_H9_D100/outs/atac_fragments.tsv.gz",
                "/storage/singlecell/zz4/multi_organoid/data/Multi_NRL_GFP_D123/outs/atac_fragments.tsv.gz",
                "/storage/singlecell/zz4/multi_organoid/data/Multi_organoid_D133/outs/atac_fragments.tsv.gz",
                "/storage/singlecell/zz4/multi_organoid/data/Multi_Organoid_D183/outs/atac_fragments.tsv.gz",
                "/storage/singlecell/zz4/multi_organoid/data/Multi_organoid_NRL_D206/outs/atac_fragments.tsv.gz",
                "/storage/singlecell/zz4/multi_organoid/data/Multi_organoid_NRL_D243/outs/atac_fragments.tsv.gz",
		"/storage/chentemp/u250758/mrrdir/organoid_metaanalysis/Cherry/atac/wt5-1/cellranger-atac/outs/fragments.tsv.gz",
		"/storage/chentemp/u250758/mrrdir/organoid_metaanalysis/Cherry/atac/wt5-2/cellranger-atac/outs/fragments.tsv.gz",
		"/storage/chentemp/u250758/mrrdir/organoid_metaanalysis/Cherry/atac/wt12-1/cellranger-atac/outs/fragments.tsv.gz",
		"/storage/chentemp/u250758/mrrdir/organoid_metaanalysis/Cherry/atac/wt12-2/cellranger-atac/outs/fragments.tsv.gz",
		"/storage/chentemp/u250758/mrrdir/organoid_metaanalysis/Cherry/atac/wt12-3/cellranger-atac/outs/fragments.tsv.gz",
		"/storage/chentemp/u250758/mrrdir/organoid_metaanalysis/Cherry/atac/wt20-1/cellranger-atac/outs/fragments.tsv.gz",
		"/storage/chentemp/u250758/mrrdir/organoid_metaanalysis/Cherry/atac/wt20-2/cellranger-atac/outs/fragments.tsv.gz",
		"/storage/chentemp/u250758/mrrdir/organoid_metaanalysis/Cherry/atac/wt28-1/cellranger-atac/outs/fragments.tsv.gz",
		"/storage/chentemp/u250758/mrrdir/organoid_metaanalysis/Cherry/atac/wt28-2/cellranger-atac/outs/fragments.tsv.gz"
    )

names <- c("Multi_H9_D35", 
           "Multi_H9_D47", 
           "Multi_organoid_NRL_D075",
           "Multi_H9_D100",
           "Multi_NRL_GFP_D123",
           "Multi_organoid_D133",
           "Multi_organoid_D183",
           "Multi_organoid_NRL_D206",
           "Multi_organoid_NRL_D243",
	         "wt5-1",
	         "wt5-2",
	         "wt12-1",
           "wt12-2",
           "wt12-3",
           "wt20-1",
           "wt20-2",
           "wt28-1",
           "wt28-2"
           )

#ArrowFiles <- createArrowFiles(
#  inputFiles = inputFiles,
#  sampleNames = names,
#  filterTSS = 4, #Dont set this too high because you can always increase later
#  filterFrags = 1000, 
#  addTileMat = TRUE,
#  addGeneScoreMat = TRUE
#)

ArrowFiles <- c("/storage/singlecell/jeanl/organoid/scripts/archr/Multi_H9_D35.arrow",
                "/storage/singlecell/jeanl/organoid/scripts/archr/Multi_H9_D47.arrow",
                "/storage/singlecell/jeanl/organoid/scripts/archr/Multi_organoid_NRL_D075.arrow",
                "/storage/singlecell/jeanl/organoid/scripts/archr/Multi_H9_D100.arrow",
                "/storage/singlecell/jeanl/organoid/scripts/archr/Multi_NRL_GFP_D123.arrow",
                "/storage/singlecell/jeanl/organoid/scripts/archr/Multi_organoid_D133.arrow",
                "/storage/singlecell/jeanl/organoid/scripts/archr/Multi_organoid_D183.arrow",
                "/storage/singlecell/jeanl/organoid/scripts/archr/Multi_organoid_NRL_D206.arrow",
                "/storage/singlecell/jeanl/organoid/scripts/archr/Multi_organoid_NRL_D243.arrow",
                "/storage/singlecell/jeanl/organoid/scripts/archr/wt5-1.arrow",
                "/storage/singlecell/jeanl/organoid/scripts/archr/wt5-2.arrow",
                "/storage/singlecell/jeanl/organoid/scripts/archr/wt12-1.arrow",
                "/storage/singlecell/jeanl/organoid/scripts/archr/wt12-2.arrow",
                "/storage/singlecell/jeanl/organoid/scripts/archr/wt12-3.arrow",
                "/storage/singlecell/jeanl/organoid/scripts/archr/wt20-1.arrow",
                "/storage/singlecell/jeanl/organoid/scripts/archr/wt20-2.arrow",
                "/storage/singlecell/jeanl/organoid/scripts/archr/wt28-1.arrow",
                "/storage/singlecell/jeanl/organoid/scripts/archr/wt28-2.arrow"
)

doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)

project <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = outputDir,
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
  
# Save the project
saveArchRProject(ArchRProj = project, outputDirectory = outputDir)

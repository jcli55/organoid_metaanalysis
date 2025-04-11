library(ArchR)
library(parallel)
addArchRThreads(threads = 8)
addArchRGenome('hg38')

#Script on how I set up the archr project for organoid

# Setup the input files, sample names, and output files
outputDir <- "/dfs3b/ruic20_lab/jeancl2/data/archr/chen_full"

inputFiles <- c("/dfs3b/ruic20_lab/jinl14/mrrdir/wkfl/hrcaprj/HECA/Organoid_retina/Lattice/Multi_H9_D35/cellranger/outs/atac_fragments.tsv.gz", 
                "/dfs3b/ruic20_lab/jinl14/mrrdir/wkfl/hrcaprj/HECA/Organoid_retina/Lattice/Multi_H9_D47/cellranger/outs/atac_fragments.tsv.gz", 
                "/dfs3b/ruic20_lab/jinl14/mrrdir/wkfl/hrcaprj/HECA/Organoid_retina/Lattice/Multi_organoid_NRL_D075/cellranger/outs/atac_fragments.tsv.gz",
                "/dfs3b/ruic20_lab/jinl14/mrrdir/wkfl/hrcaprj/HECA/Organoid_retina/Lattice/Multi_H9_D100//cellranger/Multi_H9_D100/atac_fragments.tsv.gz",
                "/dfs3b/ruic20_lab/jinl14/mrrdir/wkfl/hrcaprj/HECA/Organoid_retina/Lattice/Multi_NRL_GFP_D123/cellranger/Multi_NRL_GFP_D123/atac_fragments.tsv.gz",
                "/dfs3b/ruic20_lab/jinl14/mrrdir/wkfl/hrcaprj/HECA/Organoid_retina/Lattice/Multi_organoid_D133/cellranger/outs/atac_fragments.tsv.gz",
                "/dfs3b/ruic20_lab/jinl14/mrrdir/wkfl/hrcaprj/HECA/Organoid_retina/Lattice/Multi_Organoid_D183/cellranger/outs/atac_fragments.tsv.gz",
                "/dfs3b/ruic20_lab/jinl14/mrrdir/wkfl/hrcaprj/HECA/Organoid_retina/Lattice/Multi_organoid_NRL_D206/cellranger/outs/atac_fragments.tsv.gz",
                "/dfs3b/ruic20_lab/jinl14/mrrdir/wkfl/hrcaprj/HECA/Organoid_retina/Lattice/Multi_organoid_NRL_D243/cellranger/outs/atac_fragments.tsv.gz"
    )

names <- c("Multi_H9_D35", 
           "Multi_H9_D47", 
           "Multi_organoid_NRL_D075",
           "Multi_H9_D100",
           "Multi_NRL_GFP_D123",
           "Multi_organoid_D133",
           "Multi_organoid_D183",
           "Multi_organoid_NRL_D206",
           "Multi_organoid_NRL_D243"
           )

# If ArrowFiles are already generated, comment this block out and load in the paths of the files
ArrowFiles <- createArrowFiles(
 inputFiles = inputFiles,
 sampleNames = names,
 filterTSS = 0, # Don't filter yet
 filterFrags = 0, 
 addTileMat = TRUE,
 addGeneScoreMat = TRUE
)

# If ArrowFiles are already generated, load in the paths of those files
# ArrowFiles <- c("/dfs3b/ruic20_lab/jeancl2/data/archr/project/ArrowFiles/Multi_H9_D35.arrow",
#                 "/dfs3b/ruic20_lab/jeancl2/data/archr/project/ArrowFiles/Multi_H9_D47.arrow",
#                 "/dfs3b/ruic20_lab/jeancl2/data/archr/project/ArrowFiles/Multi_organoid_NRL_D075.arrow",
#                 "/dfs3b/ruic20_lab/jeancl2/data/archr/project/ArrowFiles/Multi_H9_D100.arrow",
#                 "/dfs3b/ruic20_lab/jeancl2/data/archr/project/ArrowFiles/Multi_NRL_GFP_D123.arrow",
#                 "/dfs3b/ruic20_lab/jeancl2/data/archr/project/ArrowFiles/Multi_organoid_D133.arrow",
#                 "/dfs3b/ruic20_lab/jeancl2/data/archr/project/ArrowFiles/Multi_organoid_D183.arrow",
#                 "/dfs3b/ruic20_lab/jeancl2/data/archr/project/ArrowFiles/Multi_organoid_NRL_D206.arrow",
#                 "/dfs3b/ruic20_lab/jeancl2/data/archr/project/ArrowFiles/Multi_organoid_NRL_D243.arrow",
# )

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
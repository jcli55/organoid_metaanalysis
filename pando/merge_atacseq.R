# Script by Zhen adapted for the organoid project by Jean

# Take all args
output_file="/storage/singlecell/jeanl/organoid/data/Chen/rds/atac_no_filter.rds"

library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(stringr)
plan("multicore", workers = 20)
set.seed(0)

mergeGRangesData <- function(bed_list) {
  peaks = data.frame()
  for (i in 1:length(bed_list)) {
    peaks = rbind(peaks, read.table(file = bed_list[i], col.names = c("chr",
                                                                      "start", "end")))
  }
  gr_list = makeGRangesFromDataFrame(peaks)
  combined.peaks <- reduce(x = gr_list)
  # Filter out bad peaks based on length
  combined.peaks@seqinfo@seqlengths <- width(combined.peaks)
  combined.peaks <- combined.peaks[combined.peaks@seqinfo@seqlengths <
                                     10000 & combined.peaks@seqinfo@seqlengths > 20]
  combined.peaks <- combined.peaks[combined.peaks@seqnames %in% c("chr1",
                                                                  "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
                                                                  "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
                                                                  "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")]
  return(combined.peaks)
}

quantify_peaks <- function(meta_list, fragment_list, combined.peaks, sams) {
  seurat_object_list <- NULL
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  # seqlevelsStyle(annotation) <- 'UCSC'
  ucsc.levels <- str_replace(string = paste("chr", seqlevels(annotation),
                                            sep = ""), pattern = "chrMT", replacement = "chrM")
  seqlevels(annotation) <- ucsc.levels
  genome(annotation) <- "hg38"
  for (i in 1:length(meta_list)) {
    md <- read.table(file = meta_list[i], stringsAsFactors = FALSE,
                     sep = ",", header = TRUE, row.names = 1)
    md$gex_barcode <- rownames(md)
    md <- md[(md$is_cell == 1) & (md$excluded_reason == 0), ]
    frags <- CreateFragmentObject(path = fragment_list[i], cells = rownames(md))
    counts <- FeatureMatrix(fragments = frags, features = combined.peaks,
                            cells = rownames(md))
    assay <- CreateChromatinAssay(counts, fragments = frags, annotation = annotation,
                                  min.cells = 10)
    seurat_object <- CreateSeuratObject(assay, assay = "peaks")
    # compute LSI seurat_object <- FindTopFeatures(seurat_object)
    # seurat_object <- RunTFIDF(seurat_object) seurat_object <-
    # RunSVD(seurat_object)
    seurat_object <- RenameCells(seurat_object, new.names = paste(sams[i],
                                                                  colnames(seurat_object), sep = "_"))
    seurat_object_list[[i]] <- seurat_object
  }
  return(seurat_object_list)
}

bed_list <- c(
  "/storage/singlecell/zz4/multi_organoid/data/Multi_H9_D35/outs/atac_peaks.bed", 
  "/storage/singlecell/zz4/multi_organoid/data/Multi_H9_D47/outs/atac_peaks.bed", 
  "/storage/singlecell/zz4/multi_organoid/data/Multi_organoid_NRL_D075/outs/atac_peaks.bed",
  "/storage/singlecell/zz4/multi_organoid/data/Multi_H9_D100/outs/atac_peaks.bed",
  "/storage/singlecell/zz4/multi_organoid/data/Multi_NRL_GFP_D123/outs/atac_peaks.bed",
  "/storage/singlecell/zz4/multi_organoid/data/Multi_organoid_D133/outs/atac_peaks.bed",
  "/storage/singlecell/zz4/multi_organoid/data/Multi_Organoid_D183/outs/atac_peaks.bed",
  "/storage/singlecell/zz4/multi_organoid/data/Multi_organoid_NRL_D206/outs/atac_peaks.bed",
  "/storage/singlecell/zz4/multi_organoid/data/Multi_organoid_NRL_D243/outs/atac_peaks.bed")
meta_list <- c(
  "/storage/singlecell/zz4/multi_organoid/data/Multi_H9_D35/outs/per_barcode_metrics.csv", 
  "/storage/singlecell/zz4/multi_organoid/data/Multi_H9_D47/outs/per_barcode_metrics.csv", 
  "/storage/singlecell/zz4/multi_organoid/data/Multi_organoid_NRL_D075/outs/per_barcode_metrics.csv",
  "/storage/singlecell/zz4/multi_organoid/data/Multi_H9_D100/outs/per_barcode_metrics.csv",
  "/storage/singlecell/zz4/multi_organoid/data/Multi_NRL_GFP_D123/outs/per_barcode_metrics.csv",
  "/storage/singlecell/zz4/multi_organoid/data/Multi_organoid_D133/outs/per_barcode_metrics.csv",
  "/storage/singlecell/zz4/multi_organoid/data/Multi_Organoid_D183/outs/per_barcode_metrics.csv",
  "/storage/singlecell/zz4/multi_organoid/data/Multi_organoid_NRL_D206/outs/per_barcode_metrics.csv",
  "/storage/singlecell/zz4/multi_organoid/data/Multi_organoid_NRL_D243/outs/per_barcode_metrics.csv")
fragment_list <- c(
  "/storage/singlecell/zz4/multi_organoid/data/Multi_H9_D35/outs/atac_fragments.tsv.gz", 
  "/storage/singlecell/zz4/multi_organoid/data/Multi_H9_D47/outs/atac_fragments.tsv.gz", 
  "/storage/singlecell/zz4/multi_organoid/data/Multi_organoid_NRL_D075/outs/atac_fragments.tsv.gz",
  "/storage/singlecell/zz4/multi_organoid/data/Multi_H9_D100/outs/atac_fragments.tsv.gz",
  "/storage/singlecell/zz4/multi_organoid/data/Multi_NRL_GFP_D123/outs/atac_fragments.tsv.gz",
  "/storage/singlecell/zz4/multi_organoid/data/Multi_organoid_D133/outs/atac_fragments.tsv.gz",
  "/storage/singlecell/zz4/multi_organoid/data/Multi_Organoid_D183/outs/atac_fragments.tsv.gz",
  "/storage/singlecell/zz4/multi_organoid/data/Multi_organoid_NRL_D206/outs/atac_fragments.tsv.gz",
  "/storage/singlecell/zz4/multi_organoid/data/Multi_organoid_NRL_D243/outs/atac_fragments.tsv.gz")

sams <- c("Multi_H9_D35", 
          "Multi_H9_D47", 
          "Multi_organoid_NRL_D075",
          "Multi_H9_D100",
          "Multi_NRL_GFP_D123",
          "Multi_organoid_D133",
          "Multi_organoid_D183",
          "Multi_organoid_NRL_D206",
          "Multi_organoid_NRL_D243")

combined.peaks <- mergeGRangesData(bed_list)
retina <- quantify_peaks(meta_list, fragment_list, combined.peaks, sams)

for (i in 1:length(sams)) {
  retina[[i]]$dataset <- sams[i]
}

seurat_object <- merge(retina[[1]], y = c(retina[2:length(sams)]), project = "retina")
seurat_object <- RunTFIDF(seurat_object)
seurat_object <- FindTopFeatures(seurat_object, min.cutoff = "q0")
seurat_object <- RunSVD(seurat_object)

saveRDS(seurat_object, output_file)

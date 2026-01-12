library(ArchR)
library(parallel)

# Peak calling and getting marker peaks
# 4/5/24 - Jean Li

project <- loadArchRProject("/storage/singlecell/jeanl/organoid/data/archr/chen_double_filter/")

project <- addGroupCoverages(ArchRProj = project, groupBy = "majorclass")

pathToMacs2 <- findMacs2()

project <- addReproduciblePeakSet(
  ArchRProj = project, 
  groupBy = "majorclass", 
  pathToMacs2 = pathToMacs2
)

project <- addPeakMatrix(project)

saveArchRProject(project)

markersPeaks <- getMarkerFeatures(
  ArchRProj = project, 
  useMatrix = "PeakMatrix", 
  groupBy = "majorclass",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

saveArchRProject(project)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerListGR <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)

write.csv(markerList, '/storage/singlecell/jeanl/organoid/csv/dars.csv')
write.csv(markerListGR, '/storage/singlecell/jeanl/organoid/csv/dars_GR.csv')

# There is a method to get up and down enriched Motifs but need to run pairwise comparison of one group to the others
 markerTest <- getMarkerFeatures(
  ArchRProj = project,
  useMatrix = "PeakMatrix",
  groupBy = "majorclass",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "PRPC",
  bgdGroups = "NRPC"
)

motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = project,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

motifsDo <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = project,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
  )

df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

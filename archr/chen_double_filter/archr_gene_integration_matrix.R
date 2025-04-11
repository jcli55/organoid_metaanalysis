library(ArchR)
library(parallel)
library(Seurat)

# Integrating snRNA-seq data to the ArchR project
# 6/21/24 - Jean Li

project <- loadArchRProject("/storage/singlecell/jeanl/organoid/data/archr/chen_double_filter/")
RNA <- readRDS('/storage/singlecell/jeanl/organoid/data/archr/integrated_rna.rds')
names <- read.csv('/storage/singlecell/jeanl/organoid/csv/ATAC_RNA_QC_intersect.csv')

intersect <- intersect(Cells(RNA), names$x)
RNA <- subset(RNA, cells = intersect)

groupList <- SimpleList(
    AC = SimpleList(
        ATAC = project$cellNames[project$majorclass %in% 'AC'],
        RNA = Cells(RNA)[RNA$majorclass %in% 'AC']
    ),
    BC = SimpleList(
        ATAC = project$cellNames[project$majorclass %in% 'BC'],
        RNA = Cells(RNA)[RNA$majorclass %in% 'BC']
    ),
    HC = SimpleList(
        ATAC = project$cellNames[project$majorclass %in% 'HC'],
        RNA = Cells(RNA)[RNA$majorclass %in% 'HC']
    ),
    Cone = SimpleList(
        ATAC = project$cellNames[project$majorclass %in% 'Cone'],
        RNA = Cells(RNA)[RNA$majorclass %in% 'Cone']
    ),
    Rod = SimpleList(
        ATAC = project$cellNames[project$majorclass %in% 'Rod'],
        RNA = Cells(RNA)[RNA$majorclass %in% 'Rod']
    ),
    MG = SimpleList(
        ATAC = project$cellNames[project$majorclass %in% 'MG'],
        RNA = Cells(RNA)[RNA$majorclass %in% 'MG']
    ),
    RGC = SimpleList(
        ATAC = project$cellNames[project$majorclass %in% 'RGC'],
        RNA = Cells(RNA)[RNA$majorclass %in% 'RGC']
    ),
    PRPC = SimpleList(
        ATAC = project$cellNames[project$majorclass %in% 'PRPC'],
        RNA = Cells(RNA)[RNA$majorclass %in% 'PRPC']
    ),
    NRPC = SimpleList(
        ATAC = project$cellNames[project$majorclass %in% 'NRPC'],
        RNA = Cells(RNA)[RNA$majorclass %in% 'NRPC']
    )
)

project <- addGeneIntegrationMatrix(
    ArchRProj = project, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = RNA,
    addToArrow = TRUE,
    force = TRUE, 
    groupList = groupList,
    groupRNA = "majorclass",
    nameCell = "predictedCell_Co",
    nameGroup = "predictedGroup_Co",
    nameScore = "predictedScore_Co"
)

project <- addImputeWeights(project)

saveArchRProject(project)

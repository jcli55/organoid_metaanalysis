#!/usr/bin/bash

SAMPLEID='Multi_organoid_D183'
OUTPUT="/storage/chentemp/u250758/organoid_metaanalysis/Chen/Multi_Organoid_D183"
BCFILE="/storage/singlecell/zz4/multi_organoid/data/Multi_Organoid_D183/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
BAMFILE="/storage/singlecell/zz4/multi_organoid/data/Multi_Organoid_D183/outs/gex_possorted_bam.bam"
GENOME="/storage/chentemp/u250758/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf"

slurmtaco.sh -p short -m 12G -- velocyto run -o $OUTPUT  -@ 12 --samtools-memory 10240 --bcfile $BCFILE -e $SAMPLEID $BAMFILE $GENOME

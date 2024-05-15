#!/usr/bin/bash

INPUT=$1

while read -r line
do

SAMPLEID=$line
OUTPUT="/storage/chentemp/u250758/organoid_metaanalysis/Chen/"$SAMPLEID
BCFILE="/storage/chentemp/u250758/organoid_metaanalysis/Chen/"$SAMPLEID"/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
BAMFILE="/storage/chentemp/u250758/organoid_metaanalysis/Chen/"$SAMPLEID"/outs/gex_possorted_bam.bam"
GENOME="/storage/chentemp/u250758/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf"

slurmtaco.sh -p short -m 12G -- velocyto run -o $OUTPUT  -@ 12 --samtools-memory 10240 --bcfile $BCFILE -e $SAMPLEID $BAMFILE $GENOME

done < $INPUT

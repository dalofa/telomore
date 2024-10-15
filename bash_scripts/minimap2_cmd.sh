#!/bin/bash
# Maps nanopore reads to a genome
REF=$1
READS=$2
THREADS=$3
OUTPUT=$4

# map reads
# It is important to use --secondary-seq and -Y in order to always
# get soft-clips in secondary reads and have their sequence included
# in the .bam file
minimap2 -a $REF $READS -t $THREADS -o $OUTPUT --secondary-seq -Y

# sort map file and assign it the intended output name
samtools sort -@ $THREADS $OUTPUT > $OUTPUT.sort.bam
mv $OUTPUT.sort.bam $OUTPUT

# index the file
samtools index $OUTPUT
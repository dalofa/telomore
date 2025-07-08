#!/bin/bash
# Maps illumina consensus to a genome
# Seems to be more accurate thatn minimap2 for the particular task
REF=$1
CONSENSUS=$2
THREADS=$3
OUTPUT=$4

# Build references for genome
bowtie2-build -q --threads $THREADS $REF $REF.bt.index

#-a to get secondary and supplementary alignments
# Local to allow soft-clips that extends beyond the genome edge (as opposed to end-to-end)
bowtie2  -f -a -p $THREADS --local --sam-no-qname-trunc -x $REF.bt.index -U $CONSENSUS -S $OUTPUT --quiet

# sort map file and assign it the intended output name
samtools sort -@ $THREADS $OUTPUT > $OUTPUT2.sort.bam
mv $OUTPUT2.sort.bam $OUTPUT

# index the file
samtools index $OUTPUT

# Remove the index-files are mapping finished
rm $REF.bt.index.*bt2

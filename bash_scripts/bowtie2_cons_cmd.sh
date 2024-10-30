#!/bin/bash
# Maps illumina consensus to a genome
REF=$1
CONSENSUS=$6
#READ1=$2
#READ2=$3
THREADS=$4
#OUTPUT=$5
OUTPUT2=$7

# Build references for genome
# For the future: Consider renaming the index-output 
bowtie2-build -q --threads $THREADS $REF $REF.bt.index

#-a to get secondary and supplementary alignments
# Local to allow soft-clips that extends beyond the genome edge (as opposed to end-to-end)
bowtie2  -a -p $THREADS --local --sam-no-qname-trunc -x $REF.bt.index -U $CONSENSUS -S $OUTPUT2 --quiet

# sort map file and assign it the intended output name
samtools sort -@ $THREADS $OUTPUT > $OUTPUT2.sort.bam
mv $OUTPUT2.sort.bam $OUTPUT2

# index the file
samtools index $OUTPUT2

# Remove the index-files are mapping finished
rm $REF.bt.index.*bt2

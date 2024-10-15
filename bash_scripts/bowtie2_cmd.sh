#!/bin/bash
# Maps illumina reads to a genome
REF=$1
READS=$2
THREADS=$3
OUTPUT=$4

# Build references for genome
# For the future: Consider renaming the index-output 
bowtie2-build -q --threads $THREADS $REF $REF.bt.index

#-a to get secondary and supplementary alignments
# Local to allow soft-clips that extends beyond the genome edge (as opposed to end-to-end)
bowtie2 --quiet -a -p $THREADS --local -x $REF.bt.index -U $READS -S $OUTPUT

# sort map file and assign it the intended output name
samtools sort -@ $THREADS $OUTPUT > $OUTPUT.sort.bam
mv $OUTPUT.sort.bam $OUTPUT

# index the file
samtools index $OUTPUT

# Remove the index-files are mapping finished
rm $REF.bt.index.*bt2
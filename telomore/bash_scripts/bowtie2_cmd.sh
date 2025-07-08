#!/bin/bash
# Maps illumina reads to a genome
REF=$1
READ1=$2
READ2=$3
THREADS=$4
OUTPUT=$5

# Build references for genome
bowtie2-build -q --threads $THREADS $REF $REF.bt.index

#-a to get secondary and supplementary alignments
# Local to allow soft-clips that extends beyond the genome edge (as opposed to end-to-end)
bowtie2  -a -p $THREADS --local --sam-no-qname-trunc -x $REF.bt.index -1 $READ1 -2 $READ2 -S $OUTPUT --quiet --no-mixed --no-discordant

# sort map file and assign it the intended output name
samtools sort -@ $THREADS $OUTPUT > $OUTPUT.sort.bam
mv $OUTPUT.sort.bam $OUTPUT

# index the file
samtools index $OUTPUT

# Remove the index-files are mapping finished
rm $REF.bt.index.*bt2

# TELOMORE
A pipeline for fishing out telomores trimmed of a linear genome (Streptomyces) by Flye.

#Process overview
The process is as follows
0. Identify the genome:
The longest contig is extracted and assumed to be the genomee.
1. Mapping of reads and extract of soft-clipped reads
Nanopore reads are mapped to the genome using minimap2 and all reads that are soft-clipped are extracted.
2. Consensus building and attachment
Using LAST and lamassemble and consensus i constructed from the soft-clip sequenced + 500 bases. The resulting consensus
is aligned to genome and any sequence that extend beyond the ends are added, creating a temporary genome+consensus.
3. Trimming of newly attached sequence
The previously identified terminal reads are mapped to genome+consensus and the sequence is trimmed. This is done by
removing bases from the ends until at least 5 reads supported the position with at least 70% of the reads supporting the
base.
4. QC and clean-up
A number of QC-maps are generated and moved to a folder called *_QC, where * is the basename of the reference supplied.
In here three bam files can be found:
- 1: The two consensuses mapped against each other
- 2: The terminal reads mapped against the final genome+trimmed_consensus
- 3: The two consensuses mapped against the genome
A file called *_consensus.log.txt can be found supplying information about the trimming of the attached sequence.


#Running on NBC-shared
Before running the script lamassemble and last must be loaded.
> module load last lamassemble
It is important to not have minimap2 loaded already as the script loads the correct version of minimap2 needed.

The script is run as follows
>python3 path_to_telomore/main.py fastq reference threads
fastq is a single file of all nanopore reads
reference is a fasta-file of the relevant genome
threads is the number of cores to use



